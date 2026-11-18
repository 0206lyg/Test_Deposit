// Particle-free planar Stokes flow with rheoTool constitutive stresses.
// Based on rheoTool rheoFoam and OpenFOAM Foundation SIMPLE, GPLv3 or later.
#include "fvCFD.H"
#include "simpleControl.H"
#include "constitutiveModel.H"
#include "planarNormalVelocity.H"
#include "meshSearch.H"
#include "fixedGradientFvPatchFields.H"
#include "Residuals.H"
#include <fstream>
#include <cmath>
#include "planarProfiles.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    simpleControl simple(mesh);
    #include "createFields.H"
    #include "createDepositionFields.H"

    if (!mesh.steady())
        FatalErrorInFunction << "A steadyState time discretisation is required" << exit(FatalError);
    Info<< "Steady creeping flow: each time index is one SIMPLE iteration." << nl;
    while (simple.loop(runTime))
    {
        Info<< "Iteration = " << runTime.timeName() << nl;
        #include "updateOutlet.H"
        #include "UEqn.H"
        #include "pEqn.H"
        constEq.correct();
        #include "writeDepositionFields.H"
        runTime.write();
    }
    Info<< "End" << endl;
    return 0;
}
