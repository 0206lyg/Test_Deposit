volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA) + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);
adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);
if (simple.dict().lookupOrDefault<Switch>("consistent", true))
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}
tUEqn.clear();
constrainPressure(p, U, phiHbyA, rAtU());
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p, "laplacian(p|(ap-H1))") == fvc::div(phiHbyA)
    );
    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();
    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}
#include "continuityErrs.H"
p.relax();
U = HbyA - rAtU()*fvc::grad(p);
U.correctBoundaryConditions();
