if (depFinalStep || runTime.outputTime() || runTime.value()-depLastMonitor >= depMonitorInterval)
{
    // The Time registry may already have re-read this dictionary; therefore
    // compare the resolved values, not only readIfModified()'s return value.
    scalarField depOldControls(13, 0);
    depOldControls[0]=depMinSteadyTime; depOldControls[1]=depSteadyWindow;
    depOldControls[2]=depTol; depOldControls[3]=depMassTol;
    depOldControls[4]=depMonitorInterval; depOldControls[5]=depRequiredChecks;
    depOldControls[6]=depAbsU; depOldControls[7]=depAbsTau;
    depOldControls[8]=depAbsWall; depOldControls[9]=depAbsApproach;
    depOldControls[10]=depAbsDP; depOldControls[11]=depAbsMass;
    depOldControls[12]=depAutoStop;
    depDict.readIfModified();
    #include "readDepositionControls.H"
    const bool depControlsChanged =
        depOldControls[0]!=depMinSteadyTime || depOldControls[1]!=depSteadyWindow ||
        depOldControls[2]!=depTol || depOldControls[3]!=depMassTol ||
        depOldControls[4]!=depMonitorInterval || depOldControls[5]!=depRequiredChecks ||
        depOldControls[6]!=depAbsU || depOldControls[7]!=depAbsTau ||
        depOldControls[8]!=depAbsWall || depOldControls[9]!=depAbsApproach ||
        depOldControls[10]!=depAbsDP || depOldControls[11]!=depAbsMass ||
        depOldControls[12]!=scalar(depAutoStop);
    if (depControlsChanged)
    {
        depConvergence.reset();
        depHavePrevious = false;
        Info<< "Deposition convergence controls reloaded; starting fresh comparisons." << nl;
    }
    // Output writes between scheduled monitors never advance the stop counter
    // or shorten the baseline used to measure field changes.
    const bool depMonitorDue = !depHavePrevious ||
        runTime.value()-depLastMonitor >= depMonitorInterval;
    const tmp<volSymmTensorField> depPolymer = constEq.tau();
    N1 = depPolymer().component(symmTensor::XX)-depPolymer().component(symmTensor::YY);
    pPa = p*constEq.rho();
    wallShearStressTotal = dimensionedVector("zero", dimPressure, vector::zero);
    scalar depWallIntegral = 0, depWallArea = 0, depWallMax = 0;
    scalar depDownIntegral = 0, depDownArea = 0;
    scalar depQin = 0, depQout = 0, depPIn = 0, depPOut = 0;
    scalar depInArea = 0, depOutArea = 0;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const scalarField& area = patch.magSf();
        if (patch.name() == "inlet")
        {
            depQin -= sum(phi.boundaryField()[patchi]);
            depPIn += sum(pPa.boundaryField()[patchi]*area);
            depInArea += sum(area);
        }
        if (patch.name() == "outlet")
        {
            depQout += sum(phi.boundaryField()[patchi]);
            depPOut += sum(pPa.boundaryField()[patchi]*area);
            depOutArea += sum(area);
        }
        if (!isA<wallFvPatch>(patch)) continue;
        const vectorField normals(patch.nf());
        const vectorField traction(tauTotal.boundaryField()[patchi] & normals);
        vectorField tangential(traction - normals*(normals & traction));
        wallShearStressTotal.boundaryFieldRef()[patchi] == tangential;
        const vectorField& centres = patch.Cf();
        forAll(area, facei)
        {
            const vector& centre = centres[facei];
            const scalar absY = mag(centre.y());
            const bool awayFromHeightEdges = mag(centre.z()) < depHalfHeight-depEdge;
            const scalar stress = mag(tangential[facei]);
            if (patch.name() == "wallsContraction" && awayFromHeightEdges &&
                absY > depHalfDown+depEdge && absY < depHalfUp-depEdge)
            {
                depWallIntegral += stress*area[facei];
                depWallArea += area[facei];
                depWallMax = max(depWallMax, stress);
            }
            if (patch.name() == "wallsDownstream" && awayFromHeightEdges &&
                centre.x() > depEdge && centre.x() < depROILength)
            {
                depDownIntegral += stress*area[facei];
                depDownArea += area[facei];
            }
        }
    }
    reduce(depQin, sumOp<scalar>()); reduce(depQout, sumOp<scalar>());
    reduce(depPIn, sumOp<scalar>()); reduce(depPOut, sumOp<scalar>());
    reduce(depInArea, sumOp<scalar>()); reduce(depOutArea, sumOp<scalar>());
    reduce(depWallIntegral, sumOp<scalar>()); reduce(depWallArea, sumOp<scalar>());
    reduce(depDownIntegral, sumOp<scalar>()); reduce(depDownArea, sumOp<scalar>());
    reduce(depWallMax, maxOp<scalar>());
    const scalar depQScale = max(mag(depQin), scalar(1e-30));
    const scalar depMass = mag(depQout-depQin)/depQScale;
    const scalar depDP = depPIn/max(depInArea,VSMALL)-depPOut/max(depOutArea,VSMALL);
    const scalar depWallMean = depWallIntegral/max(depWallArea,VSMALL);
    const scalar depDownMean = depDownIntegral/max(depDownArea,VSMALL);
    interpolationCellPoint<vector> depInterpolate(U);
    List<vector> depSampleU(depSamplePoints.size(), vector::zero);
    forAll(depSamplePoints, pointi)
        if (depSampleOwner[pointi] == Pstream::myProcNo())
            depSampleU[pointi] = depInterpolate.interpolate(depSamplePoints[pointi], depSampleCell[pointi]);
    Pstream::listCombineGather(depSampleU, plusEqOp<vector>());
    Pstream::listCombineScatter(depSampleU);
    scalar depApproach = 0, depNetApproach = 0;
    forAll(depSampleU, pointi)
    {
        depApproach += max(depSampleU[pointi].x(), scalar(0))*depDY*depDZ;
        depNetApproach += depSampleU[pointi].x()*depDY*depDZ;
    }
    // Mixed absolute/relative tolerances retain meaning near zero. Field
    // changes use ROI RMS values, so their absolute tolerances have SI units.
    const scalar depURMS = Foam::sqrt(gSum(depROIVolume*magSqr(U.primitiveField()))/depROITotalVolume);
    const scalar depTRMS = Foam::sqrt(gSum(depROIVolume*magSqr(tauTotal.primitiveField()))/depROITotalVolume);
    const scalar depPrevURMS = Foam::sqrt(gSum(depROIVolume*magSqr(depPrevU))/depROITotalVolume);
    const scalar depPrevTRMS = Foam::sqrt(gSum(depROIVolume*magSqr(depPrevTau))/depROITotalVolume);
    const scalar depChangeU = Foam::sqrt(gSum(depROIVolume*magSqr(U.primitiveField()-depPrevU))/depROITotalVolume);
    const scalar depChangeTau = Foam::sqrt(gSum(depROIVolume*magSqr(tauTotal.primitiveField()-depPrevTau))/depROITotalVolume);
    const scalar depChangeWall = mag(depWallMean-depPrevWall);
    const scalar depChangeApproach = mag(depApproach-depPrevApproach);
    const scalar depChangeDP = mag(depDP-depPrevDP);
    const scalar depUScale = max(depURMS, depPrevURMS);
    const scalar depTauScale = max(depTRMS, depPrevTRMS);
    const scalar depWallScale = max(mag(depWallMean),mag(depPrevWall));
    const scalar depApproachScale = max(mag(depApproach),mag(depPrevApproach));
    const scalar depDPScale = max(mag(depDP),mag(depPrevDP));
    const scalar depRelU = depHavePrevious ? depChangeU/max(depUScale,scalar(1e-30)) : 1;
    const scalar depRelTau = depHavePrevious ? depChangeTau/max(depTauScale,scalar(1e-30)) : 1;
    const scalar depRelWall = depHavePrevious ? depChangeWall/max(depWallScale,scalar(1e-30)) : 1;
    const scalar depRelApproach = depHavePrevious ? depChangeApproach/max(depApproachScale,scalar(1e-30)) : 1;
    const scalar depRelDP = depHavePrevious ? depChangeDP/max(depDPScale,scalar(1e-30)) : 1;
    const bool depStable = depHavePrevious &&
        depositionConvergence::withinTolerance(depChangeU,depUScale,depTol,depAbsU) &&
        depositionConvergence::withinTolerance(depChangeTau,depTauScale,depTol,depAbsTau) &&
        depositionConvergence::withinTolerance(depChangeWall,depWallScale,depTol,depAbsWall) &&
        depositionConvergence::withinTolerance(depChangeApproach,depApproachScale,depTol,depAbsApproach) &&
        depositionConvergence::withinTolerance(depChangeDP,depDPScale,depTol,depAbsDP) &&
        depositionConvergence::withinTolerance(mag(depQout-depQin),max(mag(depQin),mag(depQout)),depMassTol,depAbsMass);
    const bool depStationary = depMonitorDue && depConvergence.observe
    (
        depStable, runTime.value(), depRequiredChecks, depMinSteadyTime, depSteadyWindow
    );
    const scalar depDuration = depConvergence.duration(runTime.value());
    const scalar depUMax = gMax(mag(U.primitiveField()));
    const scalar depNMax = gMax(N1.primitiveField());
    if (Pstream::master())
    {
        std::ofstream out((depOutDir/"metrics.csv").c_str(), std::ios::app);
        out.precision(16);
        out << runTime.value() << ',' << depQin << ',' << depQout << ',' << depMass << ',' << depDP << ','
            << depWallMean << ',' << depWallMax << ',' << depDownMean << ','
            << depApproach << ',' << depNetApproach << ',' << depOffset << ','
            << depRelU << ',' << depRelTau << ',' << depRelWall << ',' << depRelApproach << ',' << depRelDP << ','
            << depUMax << ',' << depNMax << ',' << depDuration << '\n';
        if (!out) FatalErrorInFunction << "Could not write metrics.csv" << exit(FatalError);
        Info<< "Deposition diagnostics: dP=" << depDP << " Pa, step wall traction=" << depWallMean
            << " Pa, approach=" << depApproach << " m3/s, mass imbalance=" << depMass
            << ", relative U/tau=" << depRelU << '/' << depRelTau
            << ", independent check=" << depMonitorDue
            << ", stable checks=" << depConvergence.count() << '/' << depRequiredChecks
            << ", stable duration=" << depDuration << " s" << nl;
        if (depStationary)
        {
            std::ofstream marker((depOutDir/"stationary.json").c_str());
            marker.precision(16);
            marker << "{\"status\":\"stationary\",\"time\":" << runTime.value()
                << ",\"relativeU\":" << depRelU << ",\"relativeTau\":" << depRelTau
                << ",\"relativeWall\":" << depRelWall << ",\"relativeApproach\":" << depRelApproach
                << ",\"relativePressure\":" << depRelDP << ",\"massImbalance\":" << depMass
                << ",\"stableDuration_s\":" << depDuration
                << ",\"consecutiveChecks\":" << depConvergence.count()
                << ",\"requiredConsecutiveChecks\":" << depRequiredChecks
                << ",\"monitorInterval_s\":" << depMonitorInterval
                << ",\"minimumSteadyTime_s\":" << depMinSteadyTime
                << ",\"steadyWindow_s\":" << depSteadyWindow
                << ",\"stopCriterion\":\"independent_field_and_observable_checks\"}\n";
        }
        // The grid and signed Ux values make the approach integral auditable.
        std::ofstream samples((depOutDir/"approach_latest.csv").c_str());
        samples.precision(16);
        samples << "x_m,y_m,z_m,Ux_m_s,Uy_m_s,Uz_m_s,area_m2\n";
        forAll(depSamplePoints, pointi)
        {
            const vector& pt = depSamplePoints[pointi];
            const vector& vel = depSampleU[pointi];
            samples << pt.x() << ',' << pt.y() << ',' << pt.z() << ',' << vel.x() << ','
                << vel.y() << ',' << vel.z() << ',' << depDY*depDZ << '\n';
        }
    }
    if (depMonitorDue)
    {
        depPrevU = U.primitiveField(); depPrevTau = tauTotal.primitiveField();
        depPrevWall = depWallMean; depPrevApproach = depApproach; depPrevDP = depDP;
        depHavePrevious = true;
        depLastMonitor = runTime.value();
    }
    if (depStationary && depAutoStop)
    {
        Info<< "Stationarity confirmed by " << depConvergence.count()
            << " independent field/observable checks (interval " << depMonitorInterval
            << " s, sustained " << depDuration << " s); writing final fields and stopping." << nl;
        runTime.writeAndEnd();
    }
}
