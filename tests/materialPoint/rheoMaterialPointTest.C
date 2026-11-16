// GPL-3.0-or-later: a material-point regression driver for the project add-on.
// Uses the actual OpenFOAM-compiled constitutive library with a prescribed
// homogeneous shear gradient and zero transport, starting at theta = 0.
#include "fvCFD.H"
#include "constitutiveModel.H"
int main(int argc,char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("zero",dimVolume/dimTime,0.));
    const scalar g=runTime.controlDict().lookupOrDefault<scalar>("shearRate",1000.);
    const label nc=runTime.controlDict().lookupOrDefault<label>("innerCorrections",4);
    volTensorField L(IOobject("L",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedTensor("L",dimless/dimTime,tensor(0,0,0,g,0,0,0,0,0)));
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
    runTime.writeNow();
}
