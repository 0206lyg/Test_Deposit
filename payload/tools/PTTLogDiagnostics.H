/* Project-local ePTT diagnostics, GPL-3.0-or-later.
 * Guards report and stop before invalid arithmetic. They never clip a field,
 * replace a coefficient, or alter the constitutive discretisation.
 */
#ifndef depositionPTTDiagnostics_H
#define depositionPTTDiagnostics_H
#include "fvCFD.H"
#include "PstreamReduceOps.H"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>
#include <sstream>

namespace Foam
{
class depositionPTTDiagnostics
{
    const volVectorField& U_;
    word mode_;
    bool enabled_;
    scalar interval_, start_, next_;
    label call_, timeIndex_, inner_;
    bool report_;
    autoPtr<volTensorField> previousGrad_;
    autoPtr<volSymmTensorField> previousTheta_;
    scalar previousGradTime_;
    label previousGradCall_;

    static string longText(const long double x)
    { std::ostringstream s; s.precision(20); s<<x; return string(s.str()); }
    static long double logMax()
    { return std::log((long double)std::numeric_limits<scalar>::max()); }
    template<class T> static bool finite(const T& value)
    {
        for (direction c=0; c<pTraits<T>::nComponents; ++c)
            if (!std::isfinite(value.component(c))) return false;
        return true;
    }
    static bool finite(const scalar value) { return std::isfinite(value); }
    template<class T> static long double largest(const T& value)
    {
        long double v=0;
        for (direction c=0; c<pTraits<T>::nComponents; ++c)
            v=std::max(v,std::fabs((long double)value.component(c)));
        return v;
    }
    static long double traceC(const tensor& d, const tensor& v)
    {
        long double t=0;
        if(d.xy()==0 && d.xz()==0 && d.yx()==0 && d.yz()==0 && d.zx()==0 && d.zy()==0)
        {
            for(direction j=0;j<3;++j)
                for(direction i=0;i<3;++i)
                    t+=(long double)v(i,j)*v(i,j)*d(j,j);
            return t;
        }
        for (direction i=0;i<3;++i)
            for (direction j=0;j<3;++j)
                for (direction k=0;k<3;++k)
                    t+=(long double)v(i,j)*(long double)d(j,k)*(long double)v(i,k);
        return t;
    }
    template<class T> static void save
    (const GeometricField<T,fvPatchField,volMesh>& f, const fileName& instance, const word& name)
    {
        GeometricField<T,fvPatchField,volMesh> out
        (IOobject(name,instance,f.mesh(),IOobject::NO_READ,IOobject::NO_WRITE,false), f);
        out.write();
    }

public:
    depositionPTTDiagnostics(const word& mode, const volVectorField& U, const dictionary& d)
    : U_(U), mode_(mode),
      enabled_(d.lookupOrDefault<Switch>("depositionDiagnostics",true)),
      interval_(d.lookupOrDefault<scalar>("depositionDiagnosticInterval",0.0001)),
      start_(d.lookupOrDefault<scalar>("depositionDiagnosticStart",0)), next_(start_),
      call_(0), timeIndex_(-1), inner_(0), report_(false),
      previousGradTime_(-1), previousGradCall_(-1)
    {
        if (interval_<=0 || start_<0)
            FatalErrorInFunction << "Invalid ePTT diagnostic interval/start" << exit(FatalError);
    }
    void begin()
    {
        ++call_;
        if (timeIndex_!=U_.time().timeIndex()) {timeIndex_=U_.time().timeIndex();inner_=0;}
        ++inner_;
        const scalar now=U_.time().value();
        report_=enabled_ && now>=next_;
        if(report_) next_=now+interval_;
    }
    // stage=0 input, 1 relaxation f/extFun, 2 solved theta before exp(eigen(theta)),
    // 3 updated conformation before conversion to tau.
    void check
    (
        const label stage,
        const volSymmTensorField& theta,
        const volSymmTensorField& tau,
        const volTensorField& eigVals,
        const volTensorField& eigVecs,
        const volScalarField& lambda,
        const volScalarField& etaP,
        const scalar epsilon, const scalar zeta,
        const volTensorField* gradient=nullptr,
        const volTensorField* extBase=nullptr
    ) const
    {
        const fvMesh& mesh=U_.mesh();
        bool bad=false; label badCell=-1,badPatch=-1,badFace=-1;
        word reason("none");
        long double badExponent=0,badMagnitude=0;
        scalar maxExponent=-GREAT,minC=GREAT,maxC=-GREAT,maxTheta=0,maxGrad=0,maxLogStiff=-GREAT,maxLogIncrement=-GREAT;
        auto inspect = [&](const label cell,const label patch,const label face,
            const symmTensor& t,const symmTensor& stress,const tensor& d,const tensor& v,
            const scalar lam,const scalar eta,const vector& vel,const tensor* g,const tensor* base)
        {
            word why("none");long double q=0, offending=0;
            if(!finite(t)||!finite(stress)||!finite(d)||!finite(v)||!finite(vel)||!finite(lam)||!finite(eta)||(g&&!finite(*g))||(base&&!finite(*base)))
                why="nonfinite_input";
            else if(lam<=0 || eta<0 || !std::isfinite(epsilon) || !std::isfinite(zeta) || zeta==1)
                why="invalid_material_coefficient";
            else if(d.xx()<=0 || d.yy()<=0 || d.zz()<=0)
                why="nonpositive_conformation_eigenvalue";
            else
            {
                const long double trace=traceC(d,v);
                q=(long double)epsilon/(1-(long double)zeta)*(trace-3);
                const long double low=std::min((long double)d.xx(),std::min((long double)d.yy(),(long double)d.zz()));
                const long double high=std::max((long double)d.xx(),std::max((long double)d.yy(),(long double)d.zz()));
                maxExponent=max(maxExponent,scalar(std::max((long double)-GREAT,std::min((long double)GREAT,q))));
                minC=min(minC,scalar(low));maxC=max(maxC,scalar(high));
                maxTheta=max(maxTheta,scalar(largest(t)));
                if(g) maxGrad=max(maxGrad,scalar(largest(*g)));
                if(report_ || stage==1)
                {
                    const long double logStiff=std::log((long double)U_.time().deltaTValue())+q-std::log((long double)lam);
                    maxLogStiff=max(maxLogStiff,scalar(std::max((long double)-GREAT,std::min((long double)GREAT,logStiff))));
                }
                if(stage==0)
                {
                    // The upstream inverse is a full tensor inverse, whose determinant
                    // can under/overflow even when each positive diagonal is finite.
                    const long double det=(long double)d.xx()*d.yy()*d.zz();
                    if(det<std::numeric_limits<scalar>::min()) {why="inverse_determinant_underflow";offending=det;}
                    else if(det>std::numeric_limits<scalar>::max()) {why="inverse_determinant_overflow";offending=det;}
                }
                if(stage==1)
                {
                    if(!std::isfinite(q)) why="nonfinite_destruction_exponent";
                    else if(q>=logMax()) {why="destruction_exp_overflow";offending=q;}
                    else if(base && largest(*base)>0 && q+std::log(largest(*base))>=logMax())
                    {why="relaxation_product_overflow";offending=q+std::log(largest(*base));}
                    if(base && largest(*base)>0)
                    {
                        const long double logInc=std::log((long double)U_.time().deltaTValue())+q+std::log(largest(*base));
                        maxLogIncrement=max(maxLogIncrement,scalar(std::max((long double)-GREAT,std::min((long double)GREAT,logInc))));
                    }
                }
                if(stage==2 && largest(t)>2*logMax())
                {why="theta_outside_exponentiable_spectrum_range";offending=largest(t);}
                else if(stage==2)
                {
                    // Gershgorin bounds certify ordinary states without repeating
                    // calcEig's full eigensolve. Only borderline spectra need it.
                    const long double r0=std::fabs((long double)t.xy())+std::fabs((long double)t.xz());
                    const long double r1=std::fabs((long double)t.xy())+std::fabs((long double)t.yz());
                    const long double r2=std::fabs((long double)t.xz())+std::fabs((long double)t.yz());
                    const long double lower=std::min((long double)t.xx()-r0,std::min((long double)t.yy()-r1,(long double)t.zz()-r2));
                    const long double upper=std::max((long double)t.xx()+r0,std::max((long double)t.yy()+r1,(long double)t.zz()+r2));
                    if(lower<=std::log((long double)std::numeric_limits<scalar>::min()) || upper>=logMax())
                    {
                    Eigen::Matrix3d a;
                    a << t.xx(),t.xy(),t.xz(),t.xy(),t.yy(),t.yz(),t.xz(),t.yz(),t.zz();
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> e(a,Eigen::EigenvaluesOnly);
                    if(e.info()!=Eigen::Success) why="theta_eigensolver_failure";
                    else
                    {
                        for(int k=0;k<3;++k)
                        {
                            const scalar x=e.eigenvalues()(k);
                            if(!std::isfinite(x)) {why="nonfinite_theta_eigenvalue";offending=x;break;}
                            if((long double)x>=logMax()) {why="conformation_exp_overflow";offending=x;break;}
                            if((long double)x<=std::log((long double)std::numeric_limits<scalar>::min()))
                            {why="conformation_exp_subnormal_or_zero";offending=x;break;}
                        }
                    }
                    }
                }
                if(stage==3)
                {
                    const long double factor=(long double)eta/((long double)lam*(1-(long double)zeta));
                    const bool diagonal=d.xy()==0 && d.xz()==0 && d.yx()==0 && d.yz()==0 && d.zx()==0 && d.zy()==0;
                    for(direction i=0;i<3;++i) for(direction j=i;j<3;++j)
                    {
                        long double c=0;
                        if(diagonal)
                        {
                            for(direction k=0;k<3;++k)
                                c+=(long double)v(i,k)*d(k,k)*v(j,k);
                        }
                        else
                        {
                            for(direction k=0;k<3;++k) for(direction l=0;l<3;++l)
                                c+=0.5L*((long double)v(i,k)*d(k,l)*v(j,l)+(long double)v(j,k)*d(k,l)*v(i,l));
                        }
                        c-=i==j ? 1 : 0;
                        const long double result=factor*c;
                        if(!std::isfinite(result)||std::fabs(result)>std::numeric_limits<scalar>::max())
                        {why="stress_conversion_overflow";offending=result;}
                    }
                }
            }
            if(why!="none"&&!bad)
            {bad=true;badCell=cell;badPatch=patch;badFace=face;reason=why;badExponent=q;badMagnitude=offending;}
        };
        forAll(theta,c)
            inspect(c,-1,-1,theta[c],tau[c],eigVals[c],eigVecs[c],lambda[c],etaP[c],U_[c],gradient?&(*gradient)[c]:nullptr,extBase?&(*extBase)[c]:nullptr);
        // calcEig updates internal cells only. Patch theta is not exponentiated there.
        forAll(theta.boundaryField(),p)
            forAll(theta.boundaryField()[p],f)
                if(stage!=2)
                    inspect(mesh.boundary()[p].faceCells()[f],p,f,theta.boundaryField()[p][f],tau.boundaryField()[p][f],eigVals.boundaryField()[p][f],eigVecs.boundaryField()[p][f],lambda.boundaryField()[p][f],etaP.boundaryField()[p][f],U_.boundaryField()[p][f],gradient?&gradient->boundaryField()[p][f]:nullptr,extBase?&extBase->boundaryField()[p][f]:nullptr);
        bool anyBad=bad;reduce(anyBad,orOp<bool>());
        if(report_||anyBad)
        {
            reduce(maxExponent,maxOp<scalar>());reduce(minC,minOp<scalar>());reduce(maxC,maxOp<scalar>());
            reduce(maxTheta,maxOp<scalar>());reduce(maxGrad,maxOp<scalar>());reduce(maxLogStiff,maxOp<scalar>());reduce(maxLogIncrement,maxOp<scalar>());
            Info<<"EPTT_DIAGNOSTIC model="<<mode_<<" time="<<U_.time().value()<<" dt="<<U_.time().deltaTValue()
                <<" inner="<<inner_<<" call="<<call_<<" stage="<<stage<<" maxExponent="<<maxExponent
                <<" minEigC="<<minC<<" maxEigC="<<maxC<<" maxAbsTheta="<<maxTheta
                <<" maxAbsGradU="<<maxGrad<<" maxLogDtFOverLambda="<<maxLogStiff
                <<" maxLogDtRelaxationComponent="<<maxLogIncrement<<endl;
        }
        if(anyBad)
        {
            // Lazy geometry creation can communicate across processor patches.
            // Every rank must materialize it before the rank-local report branch.
            const volVectorField& cellCentres=mesh.C();
            const surfaceVectorField& faceCentres=mesh.Cf();
            const fileName instance=fileName("postProcessing")/"epttDiagnostics"/mode_/(U_.time().timeName()+"_call"+Foam::name(call_)+"_stage"+Foam::name(stage));
            mkDir(mesh.time().path()/instance);
            OFstream record(mesh.time().path()/instance/"failure.txt");
            record.precision(17);
            record<<"EPTT_DIAGNOSTIC_FAILURE\nmode "<<mode_<<"\nstage "<<stage<<"\ntime "<<U_.time().value()<<"\ndeltaT "<<U_.time().deltaTValue()<<"\ninnerIteration "<<inner_<<"\ncall "<<call_<<"\nrank "<<Pstream::myProcNo()<<"\nlocalFailure "<<bad<<"\nreason "<<reason<<"\ncell "<<badCell<<"\npatch "<<badPatch<<"\nface "<<badFace<<"\n";
            if(bad)
            {
                const vector xyz=badPatch<0 ? cellCentres[badCell] : faceCentres.boundaryField()[badPatch][badFace];
                record<<"xyz "<<xyz<<"\nexponent "<<longText(badExponent)<<"\nlogScalarMax "<<longText(logMax())<<"\noffendingValue "<<longText(badMagnitude)<<"\n";
                if(badPatch>=0) record<<"patchName "<<mesh.boundary()[badPatch].name()<<"\n";
                const symmTensor& bt=badPatch<0 ? theta[badCell] : theta.boundaryField()[badPatch][badFace];
                const symmTensor& bs=badPatch<0 ? tau[badCell] : tau.boundaryField()[badPatch][badFace];
                const tensor& bd=badPatch<0 ? eigVals[badCell] : eigVals.boundaryField()[badPatch][badFace];
                const tensor& bv=badPatch<0 ? eigVecs[badCell] : eigVecs.boundaryField()[badPatch][badFace];
                const scalar bl=badPatch<0 ? lambda[badCell] : lambda.boundaryField()[badPatch][badFace];
                const scalar be=badPatch<0 ? etaP[badCell] : etaP.boundaryField()[badPatch][badFace];
                const vector& bu=badPatch<0 ? U_[badCell] : U_.boundaryField()[badPatch][badFace];
                record<<"theta "<<bt<<"\ntau "<<bs<<"\neigVals "<<bd<<"\neigVecs "<<bv
                    <<"\nU "<<bu<<"\nlambda "<<bl<<"\netaP "<<be<<"\nepsilon "<<epsilon<<"\nzeta "<<zeta<<"\n";
                if(finite(bd)&&finite(bv)) record<<"traceC "<<longText(traceC(bd,bv))<<"\n";
                if(finite(bl)&&bl>0)
                    record<<"logDtFOverLambda "<<longText(std::log((long double)U_.time().deltaTValue())+badExponent-std::log((long double)bl))<<"\n";
                if(gradient) record<<"gradU "<<(badPatch<0 ? (*gradient)[badCell] : gradient->boundaryField()[badPatch][badFace])<<"\n";
                if(extBase) record<<"relaxationBase "<<(badPatch<0 ? (*extBase)[badCell] : extBase->boundaryField()[badPatch][badFace])<<"\n";
                if(previousGrad_.valid()) record<<"gradU_previousSuccessfulCorrection "<<(badPatch<0 ? previousGrad_()[badCell] : previousGrad_().boundaryField()[badPatch][badFace])<<"\n";
                if(previousTheta_.valid()) record<<"theta_previousSuccessfulCorrection "<<(badPatch<0 ? previousTheta_()[badCell] : previousTheta_().boundaryField()[badPatch][badFace])<<"\n";
                if(theta.nOldTimes()) record<<"theta_previousTime "<<(badPatch<0 ? theta.oldTime()[badCell] : theta.oldTime().boundaryField()[badPatch][badFace])<<"\n";

                Pout<<"EPTT_DIAGNOSTIC_FAILURE model="<<mode_<<" time="<<U_.time().value()<<" inner="<<inner_<<" stage="<<stage<<" reason="<<reason<<" cell="<<badCell<<" patch="<<badPatch<<" face="<<badFace<<" xyz="<<xyz<<" exponent="<<longText(badExponent)<<" snapshot="<<mesh.time().path()/instance<<endl;
            }
            record<<"previousSuccessfulGradTime "<<previousGradTime_<<"\npreviousSuccessfulGradCall "<<previousGradCall_<<"\nSnapshot is NOT a restart checkpoint. theta_previousTime is from the previous physical time step; current fields are at the indicated constitutive stage.\n";
            record.flush();
            save(theta,instance,"theta_current");save(tau,instance,"tau_current");save(eigVals,instance,"eigVals_current");save(eigVecs,instance,"eigVecs_current");save(U_,instance,"U_current");save(lambda,instance,"lambda_current");save(etaP,instance,"etaP_current");
            if(theta.nOldTimes()) save(theta.oldTime(),instance,"theta_previousTime");
            if(gradient) save(*gradient,instance,"gradU_current");
            if(extBase) save(*extBase,instance,"relaxationBase_current");
            if(previousGrad_.valid()) save(previousGrad_(),instance,"gradU_previousSuccessfulCorrection");
            if(previousTheta_.valid()) save(previousTheta_(),instance,"theta_previousSuccessfulCorrection");
            // All ranks enter this branch and finish writing before collective exit.
            label sync=1;reduce(sync,sumOp<label>());
            FatalErrorInFunction<<"ePTT diagnostic guard stopped before invalid arithmetic; inspect processor*/"<<instance<<"/failure.txt. No stress was clipped."<<exit(FatalError);
        }
    }
    void success(const volTensorField& grad, const volSymmTensorField& theta)
    {
        if(enabled_)
        {
            if(!previousGrad_.valid()) previousGrad_.reset(new volTensorField(IOobject("depPreviousGrad"+mode_,U_.time().timeName(),U_.mesh(),IOobject::NO_READ,IOobject::NO_WRITE,false),grad));
            else previousGrad_()=grad;
            if(!previousTheta_.valid()) previousTheta_.reset(new volSymmTensorField(IOobject("depPreviousTheta"+mode_,U_.time().timeName(),U_.mesh(),IOobject::NO_READ,IOobject::NO_WRITE,false),theta));
            else previousTheta_()=theta;
            previousGradTime_=U_.time().value();previousGradCall_=call_;
        }
    }
};
}
#endif
