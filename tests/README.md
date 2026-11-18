# Constitutive regression

After the project library has been built and the project's OpenFOAM-8 environment is loaded:

```bash
python3 tests/run_material_point.py
```

This compiles a small optional test executable and creates homogeneous-shear fixtures. It applies a prescribed homogeneous shear gradient at 1000 s^-1 to the **actual compiled PTTLog/multiMode library**, with zero advection and initially zero log-conformation. It compares polymer shear stress and first normal-stress difference after 0.3 s against the exact steady exponential-PTT shear solution. The selected PEO case is read from `config/experiment.json`. A second coefficient set tests a strongly nonlinear relaxation term; it is a numerical regression fixture, not an additional experimental fluid.

The steadyState fixture uses 5,000 iterations starting from the identity conformation. It tests the steady constitutive diagonal at the same prescribed shear rate with zero convective flux.

The transient fixtures use dt = 1e-5 s and four nonlinear corrections per time step. Default acceptance is 0.1% relative error for both stresses. Fixtures, compact progress logs and `report.json` are generated in `validation/material_point`. It is a constitutive-equation regression, not evidence that a contraction flow is spatially converged or physically steady.

The generator's repeated-eigenvalue limit is essential for the A = I initial state: the exact log-conformation strain source is L + L^T, so initially pure shear must generate a nonzero shear stress. The upstream denominator regularization by itself omits this degenerate symmetric block. The patch restores that mathematical limit without changing the ePTT constitutive equation.

The OpenFOAM-dependent driver is GPL-3.0-or-later, consistent with the parent solver and upstream library.
