# Validation — 2026-11-16 update

This update changes Newtonian termination and adds PEO diagnosis. It does not claim to cure or to identify the initiating cause of the Gautschi PEO instability.

Environment: isolated official Foundation OpenFOAM-8 (20210316), GCC 14.2.0, OpenMPI 5.0.7, Eigen 3.4.0. The user's Gautschi installation was not accessed or modified.

## Executed checks

- Compiled the actual updated OpenFOAM-8 solver and constitutive library.
- Compiled the pure C++ convergence policy and exercised first baseline, independent sample count, failed-check reset, duplicate-time rejection, near-zero mixed tolerances, nonfinite rejection and the retained PEO time-window policy.
- Actual Newtonian equilibrium fixture: zero flow (test only, not experimental fluid data), 1e-6 s monitor interval, stopped after the first baseline plus three independent successful comparisons at t=4e-6 s. Saved fields and stationarity marker agreed; run checker returned stationary.
- Actual homogeneous-shear constitutive regression: selected two-mode model polymer shear 29.4496479716 Pa, N1 875.150573678 Pa; maximum relative error 6.07e-9 versus analytical steady shear. A separate strong-nonlinearity fixture passed at 1.32e-12. This tests preservation of the constitutive response, not contraction-flow stability.
- Artificial destruction-exponential overflow: exponent 830.892 exceeds double-precision log(max)=709.783. Guard identified mode 2, stage 1, coordinates and input fields before the exponential executed.
- Asymmetric four-rank failure fixture: only rank 2 was given the invalid state. All four ranks saved failure records and snapshots, then exited with a controlled nonzero status. No FPE or MPI error occurred in the final test. This test exposed and verified the correction of lazy mesh-geometry initialization in the diagnostic branch.
- Actual existing-case convergence dictionary update preserved a backup and applied minTime=0, window=0 and three independent checks.
- Shell syntax, Python compilation, case-generation checks, wrapper argument/path dispatch, resume dispatch, preservation of diagnostic failure exit code and evidence collection.
- Installer tests: preserved customized existing launch/configuration scripts; created source backups; repeated install was a no-op; unexpected source modifications were rejected before any target bytes changed.

## Final output at an end-time cap

The short diagnostic launcher retains the original 0.1 s write interval to avoid changing time-step alignment. A diagnostic cap of 0.003 s can therefore precede its first scheduled write. The updated solver explicitly computes diagnostic fields and writes the final completed state at an end-time cap. This fixes a bug discovered in the actual MPI launcher test: normal solver completion previously left no evolved checkpoint when endTime < writeInterval. It does not count that extra write as an independent convergence check.

The final native MPI launcher and material-test results are recorded in the accompanying evidence summary.

## Limits

A bounded four-rank coarse PEO run exercised the real 268,992-cell geometry and was stopped cleanly before the user's failure time. It did not reproduce or identify the initiating Gautschi failure. The actual local startup pressure solve also showed slower convergence than the supplied Gautschi log. Do not use that brief probe as a performance estimate or stability certification.

The real diagnostic run on Gautschi is still required. Its compact diagnostic_bundle.tar.gz should establish the first captured failing stage/cell, surrounding state and preceding nonlinear-relaxation history. Mesh convergence, physical validation and a production-ready remedy for the PEO instability are outside the completed diagnosis update.
