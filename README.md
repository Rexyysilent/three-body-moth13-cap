# Three-body Moth-I^13 computer-assisted proof

Release candidate `v1.0.0-rc1` for a rigorously validated collisionless periodic orbit in the planar Newtonian three-body problem with masses `(1, 1, 1001/1000)`.

[![Verify release](https://github.com/Rexyysilent/three-body-moth13-cap/actions/workflows/verify.yml/badge.svg)](https://github.com/Rexyysilent/three-body-moth13-cap/actions/workflows/verify.yml)

## Result

The bundled interval certificate establishes:

- existence of a labeled inertial periodic orbit;
- period in the interval recorded in `proof/candidate.json`;
- collisionlessness, with every squared pair separation greater than `0.003`;
- local uniqueness inside the 1,955-dimensional shooting box;
- a local, presently unquantified, C1 continuation in the central mass.

The Moth-I^13 topology and minimal-period primitivity are numerical characterizations. The novelty statement is deliberately limited to **apparently uncatalogued** after the documented finite catalog audit.

## Verify the release candidate

On Ubuntu 24.04 or WSL Ubuntu 24.04:

```bash
bash scripts/check_release.sh
```

That command verifies the immutable proof payload, checks the scripts and metadata, recompiles the Krawczyk checkers, and requires both:

```text
PROOF_OK 1
certificate verified
```

For a complete regeneration of all 244 interval flow and derivative enclosures:

```bash
cd proof
JOBS=8 bash reproduction/scripts/rerun_all.sh
```

The full regeneration is substantially more expensive and is the preferred independent verification target.

## Container check

```bash
docker build -t three-body-moth13-cap:rc1 .
docker run --rm three-body-moth13-cap:rc1
```

## Manuscript

The `manuscript/` directory contains the reviewed LaTeX source, generated BBL,
vector figure, and compiled PDF for the research-facing paper. The prepared
arXiv bundle consists of:

```text
main.tex
references.bib
main.bbl
orbit_geometry.pdf
```

The manuscript proves existence, collisionlessness, and local uniqueness, and
states the unquantified local C1 mass-family consequence. It keeps the syzygy
word, minimal-period primitivity, and catalog novelty explicitly outside the
interval theorem.

## Contents

- `proof/`: frozen scientific payload and its original SHA-256 manifest.
- `proof/proof_report.md`: technical proof report.
- `proof/reproduction/`: C++ sources, interval data, logs, and rerun scripts.
- `manuscript/`: arXiv-ready TeX source, bibliography, figure, PDF, and metadata.
- `AUTHORSHIP.md`: accountable-authorship requirement and AI contribution record.
- `DISCLOSURE.md`: draft computational/tool-use disclosure.
- `CITATION.cff`: citation metadata for GitHub and Zenodo.
- `LICENSES.md`: Apache-2.0/CC-BY-4.0 file-class license map.
- `VERIFICATION.md`: verification levels and reporting format.
- `RELEASE_CHECKLIST.md`: completed checks and publication blockers.
- `.github/`: CI and independent-verification issue template.

## Publication status

OpenAI Codex / GPT-5.6 Sol (`gpt-5.6-sol`) is recorded as an AI computational
contributor, not as the scholarly or legal author. `Rexyysilent` is the approved
GitHub account and public project identity; **Souparna Majumder** is the
accountable scholarly author and release creator. This release candidate is
public for independent
verification; it is not yet the final `v1.0.0`/Zenodo record.

## Licensing

- Code, build files, and automation: Apache-2.0.
- Documentation, figures, certificates, logs, and research data: CC-BY-4.0.

See `LICENSES.md` for the authoritative path map and canonical license texts.

## Claim boundaries

| Status | Claims |
|---|---|
| Rigorous | Existence, labeled periodicity, collisionlessness, local uniqueness in the certified box |
| Analytic consequence | Unquantified local C1 mass family by the implicit-function theorem |
| Numerical | Energy, Moth-I^13 syzygy class, sampled clearance, minimal-period evidence |
| Catalog-conditional | Apparent novelty |

See `proof/proof_report.md` for the mathematical argument and `VERIFICATION.md` before reporting a reproduction.
