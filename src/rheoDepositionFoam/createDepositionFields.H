// Diagnostics are project-local. All stresses below are in Pa.
IOdictionary depDict
(
    IOobject("depositionProperties", runTime.constant(), mesh,
             IOobject::MUST_READ, IOobject::NO_WRITE)
);
const scalar depHalfUp = depDict.lookupOrDefault<scalar>("upstreamHalfWidth", 90e-6);
const scalar depHalfDown = depDict.lookupOrDefault<scalar>("downstreamHalfWidth", 22.5e-6);
const scalar depHalfHeight = depDict.lookupOrDefault<scalar>("halfHeight", 50e-6);
const scalar depOffset = depDict.lookupOrDefault<scalar>("approachOffset", 2e-6);
const scalar depEdge = depDict.lookupOrDefault<scalar>("edgeExclusion", 2e-6);
const scalar depROILength = depDict.lookupOrDefault<scalar>("roiHalfLength", 180e-6);
const label depNy = depDict.lookupOrDefault<label>("sampleNy", 24);
const label depNz = depDict.lookupOrDefault<label>("sampleNz", 32);
const scalar depMinSteadyTime = depDict.lookupOrDefault<scalar>("minimumSteadyTime", 0.4);
const scalar depSteadyWindow = depDict.lookupOrDefault<scalar>("steadyWindow", 0.2);
const scalar depTol = depDict.lookupOrDefault<scalar>("relativeTolerance", 1e-4);
const scalar depMassTol = depDict.lookupOrDefault<scalar>("massTolerance", 1e-5);
const scalar depMonitorInterval = depDict.lookupOrDefault<scalar>("monitorInterval", 0.01);
const Switch depAutoStop = depDict.lookupOrDefault<Switch>("stopWhenStationary", true);
if (depNy < 2 || depNz < 2 || depOffset <= 0 || depEdge < 0 ||
    depHalfUp-depHalfDown-2*depEdge <= 0 || depHalfHeight-depEdge <= 0)
{
    FatalErrorInFunction << "Invalid deposition sampling geometry" << exit(FatalError);
}
volVectorField wallShearStressTotal
(
    IOobject("wallShearStressTotal", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimensionedVector("zero", dimPressure, vector::zero)
);
volScalarField N1
(
    IOobject("N1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimensionedScalar("zero", dimPressure, 0)
);
volScalarField pPa
(
    IOobject("p_Pa", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    p*constEq.rho()
);
// State change is checked on the contraction ROI, not diluted by long ducts.
vectorField depPrevU(U.primitiveField());
symmTensorField depPrevTau(tauTotal.primitiveField());
scalar depPrevWall = 0, depPrevApproach = 0, depPrevDP = 0;
scalar depStableSince = -1;
bool depHavePrevious = false;
scalar depLastMonitor = -GREAT;
scalarField depROIVolume(mesh.nCells(), 0.0);
forAll(depROIVolume, celli)
{
    if (mag(mesh.C()[celli].x()) <= depROILength)
        depROIVolume[celli] = mesh.V()[celli];
}
pointField depSamplePoints(2*depNy*depNz);
const scalar depDY = (depHalfUp-depHalfDown-2*depEdge)/depNy;
const scalar depDZ = 2*(depHalfHeight-depEdge)/depNz;
label depPointIndex = 0;
for (label side = -1; side <= 1; side += 2)
    for (label jy = 0; jy < depNy; ++jy)
        for (label jz = 0; jz < depNz; ++jz)
            depSamplePoints[depPointIndex++] = point
            (
                -depOffset,
                side*(depHalfDown+depEdge+(jy+0.5)*depDY),
                -depHalfHeight+depEdge+(jz+0.5)*depDZ
            );
meshSearch depSearch(mesh);
labelList depSampleCell(depSamplePoints.size(), -1);
labelList depSampleOwner(depSamplePoints.size(), labelMax);
forAll(depSamplePoints, pointi)
{
    depSampleCell[pointi] = depSearch.findCell(depSamplePoints[pointi]);
    if (depSampleCell[pointi] >= 0) depSampleOwner[pointi] = Pstream::myProcNo();
}
Pstream::listCombineGather(depSampleOwner, minEqOp<label>());
Pstream::listCombineScatter(depSampleOwner);
forAll(depSampleOwner, pointi)
{
    if (depSampleOwner[pointi] == labelMax)
        FatalErrorInFunction << "Approach sample outside mesh: " << depSamplePoints[pointi]
            << exit(FatalError);
}
const fileName depGlobalCase = Pstream::parRun() ? runTime.path()/".." : runTime.path();
const fileName depOutDir = depGlobalCase/"postProcessing"/"deposition";
if (Pstream::master())
{
    mkDir(depOutDir);
    const fileName csvPath = depOutDir/"metrics.csv";
    if (!isFile(csvPath))
    {
        std::ofstream out(csvPath.c_str());
        out << "time,Qin_m3_s,Qout_m3_s,massImbalance,deltaP_Pa,"
            << "wallShearMean_Pa,wallShearMax_Pa,downstreamWallShearMean_Pa,"
            << "approachFlux_m3_s,approachNetFlux_m3_s,approachOffset_m,"
            << "relativeU,relativeTau,relativeWall,relativeApproach,relativePressure,"
            << "maxVelocity_m_s,maxN1_Pa,stableDuration_s\n";
    }
}
