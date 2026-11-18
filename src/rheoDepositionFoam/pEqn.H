volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
adjustPhi(phiHbyA, U, p);
tUEqn.clear();
constrainPressure(p, U, phiHbyA, rAU);
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();
    if (simple.finalNonOrthogonalIter()) phi = phiHbyA - pEqn.flux();
}
p.relax();
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
