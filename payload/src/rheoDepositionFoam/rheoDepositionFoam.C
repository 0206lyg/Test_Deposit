/*---------------------------------------------------------------------------*\
  rheoDepositionFoam: native-only, static-mesh transient SIMPLEC solver.

  Derived from rheoTool rheoFoam (fppimenta/rheoTool, GPL-3.0-or-later)
  and OpenFOAM Foundation's SIMPLE/PIMPLE pressure correction algorithms.
  The solver is a project add-on; it does not modify the OpenFOAM installation.

  Copyright (C) OpenFOAM Foundation and the original rheoTool authors.
  This derivative is distributed under the GNU General Public License,
  version 3 or later, WITHOUT ANY WARRANTY. See the bundled COPYING.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "constitutiveModel.H"
#include "interpolationCellPoint.H"
#include "meshSearch.H"
#include "wallFvPatch.H"
#include "DepositionConvergence.H"
#include <fstream>

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);
    #include "createFields.H"
    #include "createDepositionFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nTransient native SIMPLEC iteration toward a steady solution.\n"
        << "No particle transport equation is solved.\n" << endl;

    while (simple.run(runTime))
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const label nInIter = simple.dict().lookupOrDefault<label>("nInIter", 3);
        if (nInIter < 1)
        {
            FatalErrorInFunction << "SIMPLE.nInIter must be at least 1" << exit(FatalError);
        }
        for (label inner = 0; inner < nInIter; ++inner)
        {
            Info<< "Inner iteration " << inner + 1 << '/' << nInIter << nl;
            #include "UEqn.H"
            #include "pEqn.H"
            constEq.correct();
        }

        tauTotal = constEq.tauTotal();
        // A short diagnostic run can end before its first scheduled write.
        // Save the final completed state without changing write-time alignment.
        const bool depFinalStep = runTime.value() >=
            runTime.endTime().value() - 0.5*runTime.deltaTValue();
        #include "writeDepositionFields.H"
        if (depFinalStep) runTime.writeNow();
        else runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }
    Info<< "End\n" << endl;
    return 0;
}
