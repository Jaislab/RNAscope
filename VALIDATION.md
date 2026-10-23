# Version 0.2 validation and comparison

The workflow was revised using the newer uploaded `rnascope_pnoc_cre(1).py`.
The user selected nucleus-plus-2 µm assignment. The ≥3 detected-peak criterion
per gene remains unchanged. No overlap target enters the algorithm.

## Current example results

Both profiles were executed on the original `20X_L_CZI.czi`, using identical
StarDist/CLAHE nuclear segmentation and intensity-peak detection. Only the
assignment radius differs.

| Measurement | Nucleus only | Nucleus + 2 µm |
|---|---:|---:|
| Accepted nuclei | 1,554 | 1,554 |
| Rejected candidate masks | 216 | 216 |
| Pnoc-positive | 291 | 371 |
| Cre-positive | 89 | 121 |
| Double positive | 78 | 109 |
| Cre-positive / Pnoc-negative | 11 | 12 |
| Cre-positive also Pnoc-positive | 87.64% | **90.08%** |
| Pnoc-positive also Cre-positive | 26.80% | 29.38% |

These are provisional **detected-peak** results. Every Cre-only object is
included in the corresponding review CSV/PDF/PNG pages, with separate raw
channel panels, outlines and detections. Nuclear-only measurements are retained
in the expanded run's cell table as well.

The new DAPI masks follow broader nuclear outlines and reduce the fragmented
masks visible in the previous result. The representative comparison image
`results_diagnostics/segmentation_comparison.png` uses the same raw DAPI crop
and only accepted objects in both panels. It was not selected by overlap score.
This visual check is not a manually annotated accuracy benchmark.

## What changed in the newer uploaded script

- Local DAPI contrast enhancement before StarDist and a nucleus/background
  contrast filter, with 20–300 µm² nuclear area limits.
- A different peak detector, using morphological background correction and
  intensity maxima, rather than multiscale LoG.
- Its default assignment includes a 2 µm surrounding region. This changes the
  biological measurement and was explicitly selected by the user.
- It can substitute saturated area / 3.3 µm² for peak counts. This is an
  unvalidated area-based estimate, not an observation of individual puncta.

Separately executing the newer reference's relevant functions gave:

| Uploaded reference variant | Cre+ | Both+ | Both / Cre+ |
|---|---:|---:|---:|
| Nucleus-only detected peaks | 89 | 78 | 87.64% |
| +2 µm detected peaks | 122 | 110 | 90.16% |
| +2 µm including saturated-area estimates | 122 | 111 | 90.98% |

The revised production script retains detected peaks and saturation flags,
without area-derived count substitution. It also blocks assignment inside
rejected DAPI objects. That safeguard explains the slight difference from the
reference's expanded result. The script uses metadata/type-derived saturation
limits, rather than assuming an image's observed maximum must be saturation.

The earlier uploaded script was superseded during work. Its trial results are
not the basis of the current profile. The preceding v0.1 result (75%) and its
settings/results remain available separately for reproducibility; they used a
different segmentation and detector and should not be compared as though only
the assignment radius changed.

## Software and runtime evidence

- All **19 tests pass** in the installed Python 3.12 environment. Tests cover
  both detectors, background noise, DAPI contrast filtering, explicit backend
  selection, nuclear/expanded assignment, excluded-object blocking, physical
  distances, independent gene thresholds, denominator calculations and exports.
- StarDist inference was actually run on CPU with the official
  `2D_versatile_fluo` model. The registered archive checksum and TLS verification
  were preserved. Model file hashes and package versions are saved per scene.
- Both saved profiles were exercised end to end with successful exit status.
  Their nuclear label arrays are identical, confirming a controlled comparison.
- Every exported spot's nuclear and assignment IDs were checked against the
  saved masks. Every cell's raw nuclear and assigned counts were recomputed
  from the spot CSV. The expanded masks were checked for the 2 µm distance
  bound and exclusion of rejected-nucleus pixels.
- The manifests' script hashes match the delivered source. The dependency/model
  setup script was rerun successfully and executed the 19-test suite.

There are no manual ground-truth labels or negative-control images in this
session. Saturated clusters, background peaks, remaining segmentation errors,
and RNA in neighboring processes can affect positivity. The current >90%
result does not establish biological sensitivity/specificity. Review the masks
and the 12 Cre-only cells before interpreting these results across animals.
