# Methods validation and reference results

## Endpoint and scope

The intended primary endpoint is the fraction of Cre-positive nuclei with
detectable Pnoc signal in the same nuclear mask. Detected Cre peaks are fewer
than Pnoc peaks; low-Cre cells can miss the ≥3 criterion and inflate the apparent
Pnoc-only class. This detection-sensitivity limitation motivates the chosen
direction of analysis. Conditioning on detected Cre is not proof of unbiased
population sensitivity, nor does it recover undetected Cre-expressing cells.
Neither gene's complete expressing-cell population is established by these counts.

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

## Detected-peak distribution

The exported v0.2.1 peak tables were recomputed during round-two audit review:

| Measurement | Cre | Pnoc |
|---|---:|---:|
| All detected peaks | 1,209 | 3,531 |
| In any candidate nuclear mask, including rejected masks | 806 (66.7%) | 2,311 (65.4%) |
| In accepted nuclear masks | 778 (64.4%) | 2,195 (62.2%) |
| In accepted nucleus-plus-2 µm regions | 1,013 (83.8%) | 2,869 (81.3%) |

The two probes have similar spatial fractions among **detected** peaks.
These values do not establish a strong Cre-nuclear versus Pnoc-cytoplasmic
localization contrast. They also do not prove identical biological localization:
peak thresholds, optical resolution, background and the mask definition affect
the observed distribution. The earlier strong localization wording was removed.

Two researchers in the lab performed manual coexpression counts of approximately
90%; the auditor relayed those results and did not perform the counts.
That is useful context for the Cre-conditioned endpoint. Individual annotations,
matching object IDs, cell denominators, and numbers of sections/animals were
not supplied, so no cell-by-cell sensitivity/specificity or
formal uncertainty comparison is claimed. The reverse automated percentage is
particularly vulnerable to under-detection of low-Cre cells; its full discrepancy
from reported manual counting has not been assigned to one proven mechanism.

In the expanded run's 371 detected Pnoc-positive objects, the median assigned
Cre peak count is 1: 267/371 (72.0%) have ≥1 Cre peak, 179/371 (48.2%) have ≥2,
and 109/371 (29.4%) meet the ≥3 criterion. These observed counts illustrate
cutoff sensitivity; they do not establish the true Cre-expression distribution
or the cause of low counts. RNA coexpression is not a functional recombination
or reporter-labeling assay.

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

All **40 tests pass** in the tested Python 3.12/Linux environment. Tests cover
known peaks and background noise, explicit profiles, every shipped JSON, default
profile export, missing distribution metadata, early environment failure,
production-profile exports with mocked nucleus inference, thresholds, region
assignment, excluded-object blocking, percentages, failures and file outputs.
New tests separate priority and informational flags, exercise component-bit
saturation, verify ImageJ calibration and exact IDs at dtype boundaries, and
check that the annotated overlay retains at least native image resolution.
Cre-endpoint triage is tested for both genes' cutoffs, custom thresholds,
saturation below cutoff, excluded objects, ordering, summary counts and export.

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

No object-level expert annotations or matched negative-control images were
provided for these reference runs; the lab researchers' approximate manual totals
are described above. Review DAPI masks, negative controls, near-cutoff
calls, and saturated RNA clusters. The dedicated review outputs include every
accepted Cre-positive/Pnoc-below-cutoff object. Informational QC flags remain
available independently of the priority queue and are not automatic biological
exclusions or calibrated class-error probabilities.

Further work includes matched negative controls, sample/animal grouping,
cytoplasmic cell-boundary information when whole-cell claims are intended,
and separate validation of any future Z-stack or chance-overlap analysis.
An independent software audit and matching reference counts do not establish
biological sensitivity, specificity or a complete Pnoc-cell census.

## Round-two export and review changes (v0.2.2)

The reference image's segmentation, detected peak positions, assigned counts
and expression classes are unchanged. The prioritization rule changes the
expanded run's `nuclei_needing_review` from 955 to 413 of 1,554 accepted objects;
this is a queue-size change, not an improvement in measured biological accuracy.
All original reasons remain available in per-reason flags. DAPI saturation alone
no longer causes priority review. Every one of the 12 expanded-region Cre-only
objects remains included in the dedicated review pages.

The single-image annotated overlay retains at least native resolution and IDs
matching the cell table. ImageJ companion masks use direct µm calibration and
losslessly preserve IDs; native uint32 masks are retained. The reference CZI
reports 16 component bits, so metadata-based saturation resolves to the previous
65535 ceiling and does not change its counts.

## Round-three attribution and review tiers (v0.2.3)

The lab researchers' manual counts are correctly attributed above. Segmentation,
peak detection, assignments and expression classes remain unchanged.

The expanded `review_queue.csv` contains the existing 413 queued objects,
ordered with 77 currently Cre-positive objects (`review_priority=cre_endpoint`)
before 336 others (`general`). The first tier contains 67 double-positive and
10 Cre-only objects. The round-three report described all 77 as near the Pnoc
cutoff; the CSV shows only 22 meet that condition (17 double-positive, 5 Cre-only).
The broader 77 includes objects queued for near-Cre counts or saturated Pnoc
below cutoff. Near-Cre counts affect the observed denominator, and the general
tier can also contain missed Cre-positive objects. This triage does not claim
that only the first tier can change the endpoint.

The nucleus-only profile queues 362 objects, of which 63 are in the Cre-endpoint
tier; 20 of those are near the Pnoc cutoff. All Cre-only objects still receive
dedicated review pages. No queue/tier flag changes a cell's automated class.

## Parameter sensitivity

The subsequent [324-setting sweep](docs/SENSITIVITY.md) reuses these exported
masks and peaks while changing independent per-gene cutoffs (≥1–6) and counting
radius (0–4 µm). It reproduces baseline per-object and summary counts before
comparisons, preserves the production settings and reports numerator/denominator
counts plus a fixed baseline Cre group. Four additional tests cover sweep
behavior and multi-scene aggregation. Higher fractions are not taken as proof
of more accurate segmentation or positivity thresholds.

## Cellpose second opinion

The optional [Cellpose report](docs/CELLPOSE_COMPARISON.md) records five actual
CPU model/input variants, their method-agreement statistics and cell-region
coverage. All use the original ≥3 cutoffs and frozen peak detections. Baseline
StarDist counts reproduce. A blinded eight-tile review package is prepared;
expert object-level labels have not been supplied and no accuracy winner is
claimed. Production masks/settings remain unchanged. Six additional tests
cover comparison logic and blinded-package behavior.
