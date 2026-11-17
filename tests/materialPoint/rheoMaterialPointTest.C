// GPL-3.0-or-later: a material-point regression driver for the project add-on.
// Uses the actual OpenFOAM-compiled constitutive library with a prescribed
// homogeneous shear gradient and zero transport, starting at theta = 0.
#include "fvCFD.H"
#include "constitutiveModel.H"
#include "PTTLogDiagnostics.H"
#include <Eigen/Dense>

namespace
{
Eigen::Matrix3d matrix(const Foam::tensor& t)
{
    Eigen::Matrix3d a;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) a(i,j)=t(i,j);
    return a;
}
Foam::tensor tensorValue(const Eigen::Matrix3d& a)
{
    return Foam::tensor(a(0,0),a(0,1),a(0,2),a(1,0),a(1,1),a(1,2),a(2,0),a(2,1),a(2,2));
}
void requireClose(const double value,const double expected,const char* what)
{
    if (!std::isfinite(value) || std::abs(value-expected)>2e-11*(1+std::abs(expected)))
        FatalErrorInFunction<<what<<": got "<<value<<", expected "<<expected
            <<Foam::exit(Foam::FatalError);
}
void tensorAlgebraChecks(Foam::fvMesh& mesh)
{
    using namespace Foam;
    // Actual uploaded failure eigenbasis; all three directions are mixed.
    const tensor r(-.43355922891605519,-.30476922132258921,.8480224742046274,
        -.87462354084456284,.3688711874100255,-.31459133633848968,
        .21693330039953435,.87809439631120545,.42648584307644366);
    const tensor spectra[3]={tensor::I,tensor(.2,0,0,0,1.7,0,0,0,12),
        tensor(3.7230316007289579e-6,0,0,0,.2615797065245114,0,0,0,618196.71650896512)};
    const scalar lam=.023488735159794109,eta=.15192835835914401,zeta=.17;
    volTensorField R(IOobject("auditR",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedTensor("R",dimless,r));
    volTensorField D(IOobject("auditD",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedTensor("D",dimless,tensor::I));
    volScalarField lambda(IOobject("auditLambda",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedScalar("lambda",dimTime,lam));
    volScalarField etaP(IOobject("auditEta",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedScalar("eta",dimPressure*dimTime,eta));
    volSymmTensorField tau(IOobject("auditTau",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedSymmTensor("tau",dimPressure,symmTensor::zero));
    const symmTensor prescribed(7,1,2,8,3,9);
    volSymmTensorField fixedTau(IOobject("auditFixedTau",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedSymmTensor("tau",dimPressure,prescribed),"fixedValue");
    auto check=[&](const tensor& rotation,const tensor& diagonal,const scalar l,const scalar e,
                   const tensor& base,const scalar trace,const symmTensor& stress)
    {
        const Eigen::Matrix3d v=matrix(rotation),d=matrix(diagonal),identity=Eigen::Matrix3d::Identity();
        const Eigen::Matrix3d c=v*d*v.transpose();
        const Eigen::Matrix3d b=(v*(d.inverse()-identity)*v.transpose())/l;
        const Eigen::Matrix3d t=e/(l*(1-zeta))*(c-identity);
        requireClose(trace,c.trace(),"conformation trace");
        const tensor fullStress(stress);
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        {
            requireClose(base(i,j),b(i,j),"relaxation base");
            requireClose(fullStress(i,j),t(i,j),"polymer stress");
        }
        const tensor primitiveC=depositionPTTAlgebra::rotateDiagonal(rotation,diagonal);
        const tensor primitiveB=depositionPTTAlgebra::relaxationBase(rotation,diagonal,l);
        const tensor primitiveTau(depositionPTTAlgebra::stress(rotation,diagonal,e,l,zeta));
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        {
            requireClose(primitiveC(i,j),c(i,j),"primitive conformation");
            requireClose(primitiveB(i,j),b(i,j),"primitive relaxation");
            requireClose(primitiveTau(i,j),t(i,j),"primitive stress");
        }
    };
    // Reorder first-cell spectra; patch inputs intentionally differ from cells.
    for (label order=0;order<3;++order)
    {
        forAll(D,i)
        {
            D[i]=spectra[(i+order)%3];lambda[i]=lam*(i+1);etaP[i]=eta/(i+1);
        }
        forAll(D.boundaryField(),p) forAll(D.boundaryField()[p],f)
        {
            D.boundaryFieldRef()[p][f]=spectra[(f+p+order+1)%3];
            R.boundaryFieldRef()[p][f]=r.T();
            lambda.boundaryFieldRef()[p][f]=lam*(p+2);
            etaP.boundaryFieldRef()[p][f]=eta/(f+2);
        }
        tmp<volTensorField> base=depositionPTTAlgebra::relaxationBase(R,D,lambda);
        tmp<volScalarField> trace=depositionPTTAlgebra::traceConformation(R,D);
        depositionPTTAlgebra::updateStress(tau,R,D,etaP,lambda,zeta);
        depositionPTTAlgebra::updateStress(fixedTau,R,D,etaP,lambda,zeta);
        forAll(fixedTau.boundaryField(),p) forAll(fixedTau.boundaryField()[p],f)
            for (direction c=0;c<6;++c)
                requireClose(fixedTau.boundaryField()[p][f].component(c),prescribed.component(c),"prescribed boundary stress");
        forAll(D,i) check(R[i],D[i],lambda[i],etaP[i],base()[i],trace()[i],tau[i]);
        forAll(D.boundaryField(),p) forAll(D.boundaryField()[p],f)
            check(R.boundaryField()[p][f],D.boundaryField()[p][f],lambda.boundaryField()[p][f],
                etaP.boundaryField()[p][f],base().boundaryField()[p][f],trace().boundaryField()[p][f],tau.boundaryField()[p][f]);
    }
    Info<<"AUDIT_ALGEBRA_PASS failure-spectrum, 3D rotations, first-cell ordering, boundary values"<<endl;
}
}
int main(int argc,char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    tensorAlgebraChecks(mesh);
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("zero",dimVolume/dimTime,0.));
    const scalar g=runTime.controlDict().lookupOrDefault<scalar>("shearRate",1000.);
    const label nc=runTime.controlDict().lookupOrDefault<label>("innerCorrections",4);
    const tensor rotation=runTime.controlDict().lookupOrDefault<tensor>("rotation",tensor::I);
    const Eigen::Matrix3d q=matrix(rotation);
    const Eigen::Matrix3d gradient=q*matrix(tensor(0,0,0,g,0,0,0,0,0))*q.transpose();
    volTensorField L(IOobject("L",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedTensor("L",dimless/dimTime,tensorValue(gradient)));
    constitutiveModel model(U,phi);
    while(runTime.loop())
    {
        for(label i=0;i<nc;++i) model.correct(nullptr,&L);
        if(runTime.timeIndex()%1000==0)
        {
            const symmTensor t=model.tau()()[0];
            Info<<"AUDIT "<<runTime.value()<<' '<<t.xy()<<' '<<t.xx()-t.yy()<<endl;
        }
    }
    const symmTensor t=model.tau()()[0];
    Info<<"AUDIT_FINAL "<<runTime.value()<<' '<<t.xy()<<' '<<t.xx()-t.yy()<<endl;
    Info<<"AUDIT_TENSOR "<<t.xx()<<' '<<t.xy()<<' '<<t.xz()<<' '<<t.yy()<<' '<<t.yz()<<' '<<t.zz()<<endl;
    runTime.writeNow();
}
