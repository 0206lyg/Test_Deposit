if (runTime.outputTime() || runTime.value()-depLastMonitor >= depMonitorInterval)
{
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
    // Ratios are changes over one output interval, with physical floors near zero.
    const scalar depUNorm = gSum(depROIVolume*magSqr(U.primitiveField()));
    const scalar depTNorm = gSum(depROIVolume*magSqr(tauTotal.primitiveField()));
    const scalar depRelU = depHavePrevious ? Foam::sqrt
    (
        gSum(depROIVolume*magSqr(U.primitiveField()-depPrevU))/max(depUNorm, scalar(1e-40))
    ) : 1;
    const scalar depRelTau = depHavePrevious ? Foam::sqrt
    (
        gSum(depROIVolume*magSqr(tauTotal.primitiveField()-depPrevTau))/max(depTNorm, scalar(1e-30))
    ) : 1;
    const scalar depRelWall = depHavePrevious ? mag(depWallMean-depPrevWall)/max(mag(depWallMean),scalar(1e-6)) : 1;
    const scalar depRelApproach = depHavePrevious ? mag(depApproach-depPrevApproach)/max(mag(depApproach),depQScale*1e-6) : 1;
    const scalar depRelDP = depHavePrevious ? mag(depDP-depPrevDP)/max(mag(depDP),scalar(1e-6)) : 1;
    const bool depStable = depHavePrevious && runTime.value() >= depMinSteadyTime &&
        depRelU < depTol && depRelTau < depTol && depRelWall < depTol &&
        depRelApproach < depTol && depRelDP < depTol && depMass < depMassTol;
    if (depStable)
    {
        if (depStableSince < 0) depStableSince = runTime.value();
    }
    else depStableSince = -1;
    const scalar depDuration = depStableSince < 0 ? 0 : runTime.value()-depStableSince;
    const bool depStationary = depStable && depDuration+SMALL >= depSteadyWindow;
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
            << ", relative U/tau=" << depRelU << '/' << depRelTau << nl;
        if (depStationary)
        {
            std::ofstream marker((depOutDir/"stationary.json").c_str());
            marker.precision(16);
            marker << "{\"status\":\"stationary\",\"time\":" << runTime.value()
                << ",\"relativeU\":" << depRelU << ",\"relativeTau\":" << depRelTau
                << ",\"relativeWall\":" << depRelWall << ",\"relativeApproach\":" << depRelApproach
                << ",\"relativePressure\":" << depRelDP << ",\"massImbalance\":" << depMass
                << ",\"stableDuration_s\":" << depDuration << "}\n";
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
    depPrevU = U.primitiveField(); depPrevTau = tauTotal.primitiveField();
    depPrevWall = depWallMean; depPrevApproach = depApproach; depPrevDP = depDP;
    depHavePrevious = true;
    depLastMonitor = runTime.value();
    if (depStationary && depAutoStop)
    {
        Info<< "Stationarity criteria sustained for " << depDuration << " s; writing final fields and stopping." << nl;
        runTime.writeAndEnd();
    }
}
