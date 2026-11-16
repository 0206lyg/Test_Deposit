Info<< "Reading p (kinematic pressure) and U\n" << endl;
volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
if (p.dimensions() != dimVelocity*dimVelocity || U.dimensions() != dimVelocity)
{
    FatalErrorInFunction << "p must be kinematic pressure [0 2 -2 0 0 0 0] and U velocity."
        << exit(FatalError);
}
#include "createPhi.H"
label pRefCell = 0;
scalar pRefValue = 0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

const dictionary* coupled = mesh.solutionDict().subDictPtr("coupledSolvers");
if (coupled)
{
    const dictionary* uptau = coupled->subDictPtr("Uptau");
    if (uptau &&
        (uptau->lookupOrDefault<Switch>("solveCoupledUp", false) ||
         uptau->lookupOrDefault<Switch>("solveCoupledTau", false)))
    {
        FatalErrorInFunction << "This executable supports native segregated solvers only."
            << exit(FatalError);
    }
}
constitutiveModel constEq(U, phi);
volSymmTensorField tauTotal
(
    IOobject("tauTotal", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    constEq.tauTotal()
);
