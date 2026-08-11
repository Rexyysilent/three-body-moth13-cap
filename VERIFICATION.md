# Independent verification protocol

Please report the release tag, archive SHA-256, operating system, compiler version, library versions, commands used, exit codes, and complete stdout/stderr. Do not report a successful verification from a trajectory plot or floating-point return residual alone.

## Level 0: payload integrity

```bash
cd proof
sha256sum -c checksums.sha256
```

This proves only that the frozen files match the released payload.

## Level 1: archived-enclosure certificate check

```bash
cd proof
bash reproduction/scripts/verify_existing.sh
```

Success requires `PROOF_OK 1` and `certificate verified`. This recompiles the final point-preconditioner and MPFI Krawczyk stages but reuses the archived interval flow enclosures.

## Level 2: complete enclosure regeneration

```bash
cd proof
JOBS=8 bash reproduction/scripts/rerun_all.sh
```

This regenerates the 244 center integrations and 244 full-box state/C1 tubes before rebuilding the certificate. A Level 2 report is the minimum meaningful independent rerun of the supplied implementation.

## Level 3: independent implementation

Reimplement the reduced equations, reversing boundary map, multiple-shooting system, Krawczyk inclusion, and full-time collision tubes with an independent validated ODE package. A Level 3 agreement is the strongest verification because it does not share the supplied interval integrator.

## Required checks

- Exact masses and initial symmetry are read as rationals/intervals, not binary literals asserted as exact reals.
- Directed rounding is active; fast-math is disabled.
- Every local flow and first variational map encloses its entire input box and time slab.
- `PRECONDITIONER_NONSINGULAR 1`, `KRAWCZYK_STRICT_INCLUSION 1`, and `PROOF_OK 1` are present.
- The full-tube squared-separation lower bound is positive.
- Minimal period is not reported as rigorous unless all earlier returns are independently interval-excluded.

Use the GitHub issue template **Independent verification report** so positive and negative reruns remain auditable.
