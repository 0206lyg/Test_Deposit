// Stokes balance: div(extra stress)/rho = grad(p). No convective inertia.
tmp<fvVectorMatrix> tUEqn(-constEq.divTau(U));
fvVectorMatrix& UEqn = tUEqn.ref();
UEqn.relax();
solve(UEqn == -fvc::grad(p));
