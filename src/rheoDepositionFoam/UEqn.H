// Original rheoTool divTau supplies solvent stress, polymer stress divergence,
// and the selected BSD/coupling stress--velocity stabilization.
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U) == constEq.divTau(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();
UEqn.relax();
solve(UEqn == -fvc::grad(p));
