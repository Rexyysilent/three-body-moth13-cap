# Manuscript package

This directory contains the arXiv-ready source for:

> A computer-assisted proof of an unequal-mass Moth-I^13 periodic orbit in the planar Newtonian three-body problem

Author: **Rexyysilent**

## Build

The release-preparation build used Tectonic 0.17.0:

```text
tectonic main.tex --keep-logs --keep-intermediates
```

The arXiv upload bundle should contain `main.tex`, `references.bib`, and
`orbit_geometry.pdf`. A generated `main.bbl` is also included in the prepared
source archive for portability.

## Scope

The manuscript proves existence, collisionlessness, and local uniqueness in
the certified shooting box. It states an unquantified local C1 mass-family
consequence. It does not claim an interval proof of minimal period, syzygy
itinerary, or exhaustive catalog novelty.

## Licensing

The manuscript source, figure, and associated documentation are licensed under
CC BY 4.0. Code in the parent repository is Apache-2.0.
