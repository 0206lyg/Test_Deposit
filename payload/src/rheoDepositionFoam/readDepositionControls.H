// Included at creation and again only at monitor/write events.
// Sampling geometry is fixed for a run; these numerical stop controls are mutable.
depMinSteadyTime = depDict.lookupOrDefault<scalar>("minimumSteadyTime", 0);
depSteadyWindow = depDict.lookupOrDefault<scalar>("steadyWindow", 0);
depTol = depDict.lookupOrDefault<scalar>("relativeTolerance", 1e-4);
depMassTol = depDict.lookupOrDefault<scalar>("massTolerance", 1e-5);
depMonitorInterval = depDict.lookupOrDefault<scalar>("monitorInterval", 0.01);
depRequiredChecks = depDict.lookupOrDefault<label>("requiredConsecutiveChecks", 3);
depAutoStop = depDict.lookupOrDefault<Switch>("stopWhenStationary", true);
depAbsU = depDict.lookupOrDefault<scalar>("absoluteVelocityTolerance", 1e-9);
depAbsTau = depDict.lookupOrDefault<scalar>("absoluteStressTolerance", 1e-6);
depAbsWall = depDict.lookupOrDefault<scalar>("absoluteWallTolerance", 1e-6);
depAbsApproach = depDict.lookupOrDefault<scalar>("absoluteApproachTolerance", 1e-20);
depAbsDP = depDict.lookupOrDefault<scalar>("absolutePressureTolerance", 1e-3);
depAbsMass = depDict.lookupOrDefault<scalar>("absoluteMassTolerance", 1e-20);
if (depMinSteadyTime < 0 || depSteadyWindow < 0 || depTol <= 0 ||
    depMassTol <= 0 || depMonitorInterval <= 0 || depRequiredChecks < 2 ||
    depAbsU < 0 || depAbsTau < 0 || depAbsWall < 0 || depAbsApproach < 0 ||
    depAbsDP < 0 || depAbsMass < 0 ||
    !std::isfinite(depMinSteadyTime+depSteadyWindow+depTol+depMassTol+
                   depMonitorInterval+depAbsU+depAbsTau+depAbsWall+
                   depAbsApproach+depAbsDP+depAbsMass))
{
    FatalErrorInFunction << "Invalid convergence controls in depositionProperties. "
        << "Use nonnegative time windows/absolute tolerances, positive relative "
        << "tolerances and monitorInterval, and requiredConsecutiveChecks >= 2."
        << exit(FatalError);
}
