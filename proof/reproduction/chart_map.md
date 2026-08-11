# Chart and visual QA notes

## Orbit geometry

- Question: what paths do the bodies trace over one reconstructed period?
- Takeaway: the certified center follows three intricate bounded paths with repeated close passages.
- Form: static two-dimensional path plot, x-position versus y-position, one labeled series per body.
- Data: 245 half-orbit mesh centers plus the exact reversing reconstruction.
- Palette: blue and orange-red for the equal outer masses; olive-green for the central mass.
- Trust boundary: center visualization only; it proves neither closure nor collision clearance.
- QA: self-contained SVG, white background, 1000×760 view box, explicit legend, no external assets, finite in-frame coordinates.

## Certificate ratios

- Question: how far are the normalized validation bounds from the strict threshold one?
- Takeaway: the largest displayed ratio is the Krawczyk inclusion ratio, approximately 0.2296.
- Form: native vertical bar chart in the report artifact.
- Fields: certificate component, normalized ratio, strict threshold, interpretation, and source field.
- Palette: single blue root; no redundant legend.
- Trust boundary: values derive directly from the archived final certificate.
