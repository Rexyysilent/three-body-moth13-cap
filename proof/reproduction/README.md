# Reproducing the interval certificate

## Fast archived-log verification

From this directory under WSL Ubuntu, run:

    bash scripts/verify_existing.sh

This recompiles the point preconditioner and high-precision Krawczyk checker, consumes the archived center and outer interval logs, and must finish with:

    PROOF_OK 1
    certificate verified

## Full regeneration

Run:

    JOBS=8 bash scripts/rerun_all.sh

The full route compiles every verifier, generates 244 center jobs and 244 outer-box jobs, validates each segment tube and first variational map, assembles the logs, and executes the final Krawczyk test. Set JOBS to fit available memory and CPU.

## Inputs and outputs

- data/mesh244.certmesh: exact-decimal mesh, schedule, center, and radii.
- data/mpfi_center.log: 80-decimal center enclosures.
- data/outer_box.log: full-box state, derivative, and collision-tube enclosures.
- logs/final_certificate.txt: publication certificate.
- source/: C++ interval integrators, refinement utilities, and Krawczyk verifiers.
- scripts/: job generation, assembly, quick verification, and full rerun.

The scripts prefer a local dependency prefix at ../../work/deps when present and otherwise use system GMP/MPFR/MPFI packages. Do not enable fast-math.
