# Contributing

The most valuable contributions are independent verification reports, portability fixes that preserve directed rounding, and independently implemented interval checks.

1. Verify the frozen payload checksums.
2. State whether you ran Level 1, Level 2, or Level 3 from `VERIFICATION.md`.
3. Preserve raw logs and environment details.
4. Open an issue before changing mathematical conventions, coordinate order, or certificate inequalities.

Do not weaken compiler rounding flags, replace interval bounds with ordinary floating-point tolerances, or describe numerical primitivity as a theorem.

By contributing, you agree that your contribution is made available under the
license assigned to the target file class in `LICENSES.md`: Apache-2.0 for code
and automation, or CC-BY-4.0 for documentation and research data.
