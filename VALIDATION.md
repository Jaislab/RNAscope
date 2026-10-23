# Methods validation and reference results

## Endpoint and scope

The intended primary endpoint is the fraction of Cre-positive nuclei with
detectable Pnoc signal in the same nuclear mask. The motivating images show
nucleus-enriched Cre fluorescence and predominantly cytoplasmic Pnoc signal.
Nuclear detection can therefore miss Pnoc expression outside the nucleus.
Neither total Pnoc-expressing cells nor the fraction of all Pnoc cells expressing
Cre can be inferred from these nuclear measurements alone.

The optional nucleus-plus-2 µm region includes additional surrounding RNA.
This is approximate nucleus-associated signal, not a validated whole-cell
boundary or a purely nuclear-overlap measurement. The region applies to both
RNA channels, so the Cre-positive denominator can also change.

## Reference image and profiles

The reference `20X_L_CZI.czi` is a single 1024 × 1024 section, four 16-bit
channels, with X/Y sampling 0.690535481770835 µm/pixel. The source CZI and runtime
outputs are not distributed with this repository. The numerical examples are
reference checks, not a public regression dataset or manual accuracy benchmark.

Both current profiles use identical StarDist/CLAHE nuclear segmentation and
intensity-peak detection, with an independent ≥3 peak criterion for each gene.
Only the assignment radius differs. The supplied parameter files reproduce the
following reference results:

| Measurement | Nuclear mask | Nuclear mask + 2 µm |
|---|---:|---:|
| Accepted nuclei | 1,554 | 1,554 |
| Rejected candidate masks | 216 | 216 |
| Pnoc-positive objects | 291 | 371 |
| Cre-positive objects | 89 | 121 |
| Double-positive objects | 78 | 109 |
| Cre-positive / Pnoc-below-cutoff objects | 11 | 12 |
| Cre-positive with detectable Pnoc | **87.64%** | **90.08%** |
| Detected Pnoc-positive objects also Cre-positive | 26.80% | 29.38% |

These are provisional **detected-peak** results. A Pnoc-below-cutoff nucleus is
not necessarily a Pnoc-negative cell. The reverse percentage refers only to
the observed region-positive Pnoc subset. No expected overlap percentage is an
analysis parameter; the expanded result exceeding 90% does not validate accuracy.

## Reference implementation comparison

The v0.2.0 method was compared with a separate reference implementation using
CLAHE/StarDist, local DAPI-contrast filtering and morphological-background
intensity peaks. Results from that comparison were:

| Reference variant | Cre-positive | Double positive | Both / Cre-positive |
|---|---:|---:|---:|
| Nuclear-mask detected peaks | 89 | 78 | 87.64% |
| +2 µm detected peaks | 122 | 110 | 90.16% |
| +2 µm with saturated-area estimates | 122 | 111 | 90.98% |

This repository preserves detected peaks instead of replacing counts with
saturated area / 3.3 µm². It also blocks assignment inside rejected DAPI objects.
These documented choices account for the small differences in the expanded
comparison. The separate reference is a software comparison, not independent
biological ground truth. The earlier watershed/LoG result is documented in
[the historical v0.1 note](docs/archive/validation-v0.1.md).

## Software evidence

All **26 tests pass** in the tested Python 3.12/Linux environment. Tests cover
known peaks and background noise, explicit profiles, every shipped JSON, default
profile export, missing distribution metadata, early environment failure,
production-profile exports with mocked nucleus inference, thresholds, region
assignment, excluded-object blocking, percentages, failures and file outputs.

Actual CPU inference with official `2D_versatile_fluo` weights was executed on
the reference CZI. Registered archive checksum and TLS verification were
preserved; per-scene metadata record model file hashes and package versions.
Both profiles completed end to end with the reference counts above. Their
nuclear label arrays match, supporting a controlled region comparison.

During the v0.2.1 audit fix, the nucleus-only run was also repeated with the
`tensorflow` distribution-version lookup deliberately raising
`PackageNotFoundError`. Inference and exports still completed using the imported
module's version. This simulates the provenance failure; it does not claim an
actual installation test of every alternative TensorFlow package or platform.

Exported spot IDs/counts were checked against saved nuclear and assignment
labels. The expanded masks respect the 2 µm distance limit and rejected-object
blocking. Script hashes in the run manifests match the source used. Active
methods/radius are printed and recorded; missing `--config` now fails before
analysis, instead of selecting historical parameters.

## Experimental work still needed

No expert-annotated ground truth or matched negative-control images were used
for these reference runs. Review DAPI masks, negative controls, near-cutoff
calls, and saturated RNA clusters. The dedicated review outputs include every
accepted Cre-positive/Pnoc-below-cutoff object. QC flags remain broad and are
not automatic biological exclusions or calibrated class-error probabilities.

Further work includes matched negative controls, sample/animal grouping,
cytoplasmic cell-boundary information when whole-cell claims are intended,
and separate validation of any future Z-stack or chance-overlap analysis.
An independent software audit and matching reference counts do not establish
biological sensitivity, specificity or a complete Pnoc-cell census.
