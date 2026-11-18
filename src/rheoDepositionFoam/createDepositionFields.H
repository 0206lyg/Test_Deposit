IOdictionary depDict
(
    IOobject("depositionProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);
const word depModel(depDict.lookup("model"));
const scalar halfUp = readScalar(depDict.lookup("upstreamHalfWidth"));
const scalar offset = readScalar(depDict.lookup("approachOffset"));
const scalar depth = readScalar(depDict.lookup("meshDepth"));
const scalar qHalf = readScalar(depDict.lookup("qHalfPerDepth"));
const scalar depEta = readScalar(depDict.lookup("outletSolventViscosity"));
const scalar depTol = depDict.lookupOrDefault<scalar>("relativeTolerance", 1e-3);
const scalar resTol = depDict.lookupOrDefault<scalar>("residualTolerance", 1e-5);
const scalar massTol = depDict.lookupOrDefault<scalar>("massTolerance", 1e-5);
const label monitorEvery = depDict.lookupOrDefault<label>("monitorEvery", 10);
const label stableChecks = depDict.lookupOrDefault<label>("stableChecks", 3);
const scalar outletRelax = depDict.lookupOrDefault<scalar>("outletRelaxation", 0.3);
const bool elastic = depModel == "peo";
const label inletID = mesh.boundaryMesh().findPatchID("inlet");
const label outletID = mesh.boundaryMesh().findPatchID("outlet");
if (depth <= 0 || qHalf <= 0 || depEta <= 0 || offset <= 0 || monitorEvery < 1 ||
    stableChecks < 1 || depTol <= 0 || resTol <= 0 || massTol <= 0 ||
    outletRelax <= 0 || outletRelax > 1 || mesh.nGeometricD() != 2 ||
    inletID < 0 || outletID < 0)
    FatalErrorInFunction << "Invalid 2D mesh or steady calculation parameters" << exit(FatalError);
if (!p.boundaryField()[outletID].fixesValue())
    FatalErrorInFunction << "Outlet p must use fixedValue" << exit(FatalError);

volScalarField pPa
(
    IOobject("p_Pa", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    p*constEq.rho()
);
volScalarField N1
(
    IOobject("N1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimensionedScalar("zero", dimPressure, 0)
);
DynamicList<WallRow> localWalls;
label wallPatchID[3];
for (label w=0; w<3; ++w)
{
    const label k = mesh.boundaryMesh().findPatchID(wallPatches[w]);
    if (k<0) FatalErrorInFunction << "Missing " << wallPatches[w] << exit(FatalError);
    wallPatchID[w]=k;
    const fvPatch& patch = mesh.boundary()[k];
    forAll(patch.Cf(), i)
    {
        const point& pFace = patch.Cf()[i];
        WallRow row;
        row[0]=w; row[1]=w==1 ? halfUp-pFace.y() : pFace.x();
        row[2]=pFace.x(); row[3]=pFace.y(); row[4]=patch.magSf()[i]/depth;
        localWalls.append(row);
    }
}
List<List<WallRow>> byProc = collectRows<WallRow>(localWalls);
Pstream::scatterList(byProc);
label wallStart=0;
DynamicList<WallRow> wallRows;
forAll(byProc, proc)
{
    if (proc==Pstream::myProcNo()) wallStart=wallRows.size();
    wallRows.append(byProc[proc]);
}
const label nWall=wallRows.size();
pointField samplePoints(nWall);
meshSearch search(mesh);
labelList sampleCell(nWall,-1), owner(nWall,labelMax);
forAll(wallRows,i)
{
    const WallRow& r=wallRows[i];
    const bool step=label(r[0])==1;
    samplePoints[i]=point(r[2]-(step?offset:0),r[3]-(step?0:offset),0);
    sampleCell[i]=search.findCell(samplePoints[i]);
    if (sampleCell[i]>=0) owner[i]=Pstream::myProcNo();
}
Pstream::listCombineGather(owner,minEqOp<label>());
Pstream::listCombineScatter(owner);
forAll(owner,i) if (owner[i]==labelMax)
    FatalErrorInFunction << "Wall-offset sample outside mesh at " << samplePoints[i] << exit(FatalError);
const fileName casePath = Pstream::parRun() ? runTime.path()/".." : runTime.path();
const fileName depOutDir = casePath/"postProcessing"/"deposition";
if (Pstream::master())
{
    mkDir(depOutDir);
    const fileName history = depOutDir/"metrics.csv";
    const bool newHistory = !isFile(history) || fileSize(history)==0;
    std::ofstream metrics(history.c_str(),std::ios::app);
    if (newHistory) metrics << "iteration,Qin_m2_s,Qout_m2_s,massImbalance,residualU,residualP,"
            << "residualTheta1,residualTheta2,relativeWall,relativeApproach,relativeApproachProfile,"
            << "stepApproach_m2_s,stepApproachFraction,outletTractionMax_Pa,stableChecks\n";
}
scalarField prevWall(nWall,0), prevNormal(nWall,0);
FixedList<scalar,3> prevQ(scalar(0));
bool havePrevious=false;
label consecutive=0;
