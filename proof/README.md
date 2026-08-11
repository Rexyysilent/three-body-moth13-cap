# Validated three-body orbit bundle

This bundle certifies an apparently uncatalogued collisionless periodic orbit for masses **(1, 1, 1001/1000)** in the planar Newtonian three-body problem.

- [Technical proof report](proof_report.md)
- [Machine-readable candidate](candidate.json)
- [Orbit geometry](orbit_geometry.svg)
- [Final interval certificate](reproduction/logs/final_certificate.txt)
- [Reproduction guide](reproduction/README.md)
- [SHA-256 checksums](checksums.sha256)

A 1,955-dimensional outward-rounded Krawczyk calculation proves strict inclusion with maximum normalized radius **0.2296323**. Validated full-time tubes prove every squared pair separation exceeds **0.003**. The period interval is recorded in candidate.json.

Claim boundaries:

- rigorous: existence, collisionlessness, labeled periodicity, and local uniqueness in the box;
- numerical: Moth-I¹³ topology and minimal-period primitivity;
- catalog-conditional: apparent novelty.

Start with proof_report.md. Run reproduction/scripts/verify_existing.sh for the quick archived-log recheck.
