# Planar steady workflow verification

Base revision: `df22a9f6bd6607754c1fa10ade7272133e6f71b6`. The checks below use Foundation OpenFOAM 8 on a Linux CPU host. Local wall times depend on the concurrent workload. No Slurm allocation at Gautschi was used.

| Check | Observed result |
| --- | --- |
| Four one-layer meshes | blockMesh/checkMesh pass with two solution directions and zero non-orthogonality. Final meshes contain 807, 6,552, 14,796 and 33,510 cells. |
| Inlet | Mesh-face averages conserve the prescribed half flow, including 6.94444444×10⁻⁷ m²/s at 0.5 mL/h. The two inlet stress, log-conformation and spectral states agree with the prescribed steady shear state. |
| Steady constitutive equation | Prescribed shear at 100, 1000 and 4000 s⁻¹ from rest reproduces the analytical two-mode shear and first normal stress difference. At 1000 s⁻¹ the polymer shear is 29.44964815 Pa and N1 is 875.15057528 Pa. |
| Two-dimensional sampling | An affine divergence-free field on the Cartesian mesh is recovered from oriented face fluxes by the offset sampler to 2.5×10⁻¹⁴ m/s in serial and two-rank MPI tests. Integrals use line lengths and units of m²/s. |
| Resume | A deliberately capped Cross run reported unconverged. The workflow restarted from its last complete common checkpoint, converged and reused that checkpoint on a second invocation. Wrong-iteration and truncated diagnostic files failed validation. |
| Slit reference | For 0.5 mL/h the 6,552-cell result at x = +500 µm gives Cross/Gly75 = 0.6120. The fully developed planar reference is 0.6113. |

A clean checkout of the base commit with the overlay was built through `build_dep.sbatch`. The compiled material-point fixtures for standard shear, the strong-nonlinearity regression, rotated 3D shear, and the new steadyState case passed. A normal two-rank `run_dep.sbatch peo 0.5 smoke` completed all five stages and produced reconstructed fields. With identical initialization of the final stage, the serial and MPI step approach fractions were 0.590846% and 0.590844% of the imposed half flow, respectively. This small mesh serves as a workflow check.

All three materials converged on 6,552 cells at 0.3 and 0.5 mL/h and on 14,796 cells at 0.5 mL/h. At 0.5 mL/h, the dominant ePTT band spans the last 1.79 and 1.82 µm of the step on the two grids. Using the same lower-resolution initial field, step approach fractions were 0.304% and 0.276% of the imposed half flow. The traction at the face with the highest approach speed changed from 56 to 79 Pa. Its absolute magnitude remains grid-sensitive. The distributed launcher completed all five stages at 0.3 and 0.5 mL/h on one local CPU rank in approximately 26 minutes per complete ePTT workflow. At 0.5 mL/h the stages used 1,430, 2,280, 4,890, 6,860 and 6,910 iterations. Same-grid continuation gave 0.398% and 0.325% at 0.3 and 0.5 mL/h, versus 0.378% and 0.304% after same-flow interpolation onto that mesh. The dominant band reached the lip in every comparison.

The comparison data accompanying the overlay identify the flow, grid, iteration, residual and inlet/outlet balance for each published solution. Grid sensitivity of near-lip flow and stress is material-specific. The stated stop criteria concern the discrete steady equations on a given mesh.

## Earlier 3D work (historical)

This is software/build and short-run verification. No production steady contraction solution or mesh-independent deposition-mechanism result is claimed.

## Environment actually tested

- Foundation OpenFOAM-8, official Ubuntu focal package `openfoam8_20210316_amd64.deb`, build `8-1c9b5879390b`, extracted into a separate scratch tree.
- Debian 13 host, GCC 14.2.0, OpenMPI 5.0.7, Eigen 3.4.0.
- Project-local source generation, library and solver compilation through the delivered `build_dep.sbatch` script.
- Local MPI tests used OpenMPI's shared-memory transport (`btl=self,sm`, `pml=ob1`). These test-only settings are not forced onto Gautschi's existing OpenMPI installation.
- Gautschi itself, its GCC 11.4.1/OpenMPI 4.1.6 modules, Slurm allocation and network were not accessed.

## Passed checks

| Check | Evidence / outcome |
|---|---|
| Full build script | Completed source verification, compilation, material regression and build manifest |
| Selected two-mode ePTT at gamma=1000/s | Polymer shear 29.4496479716 Pa, N1=875.150573678 Pa; maximum relative error 6.07e-9 versus the exact steady-shear solution |
| Strong-nonlinearity numerical fixture | Maximum relative error 1.32e-12; explicitly not an alternative proposed fluid |
| Initial log conformation A=I | Actual compiled pure-shear startup generates stress and reaches the correct analytical state |
| Mesh smoke/coarse | Actual blockMesh and checkMesh: Mesh OK, 9440/268992 cells, non-orthogonality 0 |
| Pressure solver | On the smoke PEO case, PCG with nested GAMG converged in all 30 tested inner pressure solves; final flux imbalance 2.68e-6 |
| Gly75 serial launcher | New case and actual restart from 1e-4 to 2e-4 s completed; final flux imbalance 9.31e-10 |
| PEO two-rank MPI | Fresh launcher, field writes, reconstruction and a separate actual restart from 2e-5 to 4e-5 s completed; restart final flux imbalance 1.35e-7 |
| Optional Cross control | Actual generated case and short native-solver execution completed |
| Equilibrium auto-stop | A separate exact zero-flow fixture stopped at 3e-6 s and wrote consistent final fields plus marker; this fixture is not an experimental case |
| Failure classification | Short finite-flow tests correctly report unconverged, never steady |
| Restart safeguards | Complete field footer check, tau/theta/eigenstate requirement, common MPI checkpoint, rank mismatch and stale marker checks |
| Base/build safeguards | Missing/stale compiler, MPI, source or add-on hashes block execution; existing case creation is refused |
| Comparison plots | Actual CFD CSV histories rendered; step-face and downstream-sidewall regions labeled separately |
| ParaView pipeline | ParaView 5.13.2 loaded actual U, N1 and total wall traction on seven wall patches and wrote an editable state; differing final times preserved |

Compact run-status/metrics evidence and the material-point report are under `validation/evidence/`. The short-run stress and pressure values in those files are transient plumbing-test results, not scientific predictions of the experimental steady state.

## Not yet established

- Long-time convergence or sustained physical unsteadiness at each of the three experimental flow rates.
- Mesh and timestep convergence of the corner flow, mean wall traction and finite-offset approach integral.
- Performance or full 16-rank execution on the user's Gautschi allocation.
- Actual three-dimensional PNG rendering by ParaView: the tested installation required an unavailable display. State creation and the complete data/filter pipeline passed; CSV PNG rendering passed.
- Quantitative agreement with experimental vortices or deposition, which is a scientific validation step after the runs.

The time cap does not assert steady state. Every completed run is classified by saved fields, finite metrics and a matching stationarity marker. Changing linear residual tolerances, the diagnostic sampling interval or the stationarity thresholds changes the numerical acceptance test and should be recorded.
