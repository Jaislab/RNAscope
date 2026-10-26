# Methods and interpretation

## Intended measurement

The primary endpoint is Pnoc signal associated with Cre-positive nuclei.
Cre has fewer detected peaks than Pnoc in the reference image. Under-detection
of weak Cre signal can move Cre-expressing cells into the Pnoc-only class under
the ≥3-peak rule, affecting the reverse conditional percentage. Conditioning on
detected Cre-positive objects addresses the chosen endpoint but does not remove
all selection bias or recover cells missed by the detector.

The similar detected-peak fractions inside candidate DAPI masks (66.7% Cre,
65.4% Pnoc) do not establish stronger nuclear localization of one probe. These
are thresholded-detection and mask statistics, not a biochemical localization
assay. A DAPI mask defines a nuclear region, not an entire cell. Neither nuclear
nor expanded-region counts establish the total population expressing a gene.

The two channel thresholds are applied independently (≥3 peaks in each region).
Double positivity requires both genes in the same region; it does not require
individual peaks to overlap. Counts are detections, not calibrated transcript
abundances. Different gains/laser settings prevent direct abundance comparison
across channels without additional calibration.

The conditional Cre-positive → Pnoc-positive fraction can inform driver-line
characterization as detected RNA coexpression. Cre mRNA is not a direct assay of
functional recombination or of labeling efficiency across the Pnoc population.
A Cre-dependent reporter is a distinct readout that requires separate validation
and probe-specific names in outputs; current Cre outputs must not be relabeled
as reporter measurements merely by changing a channel index.

## Image reading and calibration

The reader reconstructs individual CZI scenes with C/Y/X axes. Non-singleton
Z, T or other extra axes are rejected; no max projection or first-plane
selection occurs. Input images remain read-only. Multiple scenes are processed
separately, and failed items do not prevent processing later files.

| Role | Metadata match | Fallback 0-based index |
|---|---|---:|
| Cre | Cy5-T1 / Cy5 | 0 |
| Pnoc | At520-T3 / ATTO 520 | 2 |
| DAPI | DAPI-T4 / DAPI | 3 |

Metadata take precedence. Fallback indices are used only when channel names are
missing, with a warning. Ambiguous/unrecognized names require an explicit
`channels.<role>.index` after confirming probe identity. Dye metadata alone do
not establish which biological probe was used. Cy3 is only a QC reference.

X/Y sampling comes from CZI scaling metadata. Missing calibration requires
`pixel_size_um: [y_um, x_um]`; no assumption of square pixels is made silently.
Coordinates are local to the scene's top-left (x=column, y=row), not stage
coordinates. Physical masks and assignment distances account for anisotropy.

## Analysis steps

### 1. Nuclear segmentation

RNA channels never determine the masks. The method is selected explicitly by
the profile and is never replaced silently; a missing model or environment
aborts the run during preflight.

- **Cellpose `nuclei` (main profile).** Raw DAPI is handed, as a numeric NPY
  array, to an isolated subprocess running pinned Cellpose 3.1.1.3 on CPU with
  the legacy `nuclei` weights (SHA-256 verified before loading). The model
  receives an exact 17-pixel diameter, its own 1st/99th-percentile
  normalization, flow threshold 0.4, cell-probability threshold 0, minimum mask
  size 15 pixels, resampling and no augmentation. The worker echoes the keyword
  set it ran; masks whose echo differs from the request are refused. Square XY
  sampling is required.
- **StarDist `2D_versatile_fluo` (historical profiles).** Smooth DAPI at 1.4 µm,
  normalize to the 1st/99.8th percentiles and apply CLAHE with a 25 µm window
  and clip 0.01. StarDist predicts outlines at a target sampling of
  0.69 µm/pixel with its default probability threshold unless overridden.
- **Watershed (legacy profile).** Li or Otsu threshold of background-corrected
  DAPI, enclosed-hole filling, distance-transform seeds at a 3 µm minimum
  separation and watershed splitting. An optional local-watershed backend uses
  local thresholding and solidity filtering.

### 2. Acceptance

Accept candidates of 20–300 µm² (15–250 µm² in the legacy profile), excluding
image-border objects and, where `min_dapi_contrast` is set, objects whose DAPI
mean, smoothed at `nuclei.contrast_smooth_um` (1.4 µm, independent of the
StarDist preprocessing scale), is below 1.15 times the mean in an available
1.5–3.5 µm surrounding ring. An unavailable ring is unknown, not grounds for
rejection. Rejected candidates remain in the CSV and label image with their
reasons and never enter summary totals.

### 3. Counting region

Use nuclear masks directly, or expand accepted masks by the configured radius
with nearest-label assignment and no overlap. Pixels belonging to rejected DAPI
objects are blocked. The 2 µm expansion is an approximate nucleus-associated
region and can include another cell's cytoplasm or processes.

### 4. RNA peak detection

The shipped Cellpose and StarDist profiles independently remove each RNA
channel's morphological background at an 8 µm radius, smooth at about 0.55 µm,
and detect intensity maxima above the median plus 10 × robust noise
(`1.4826 × MAD`), with a 1 µm minimum separation. The legacy profile instead
uses multiscale, scale-normalized LoG responses with per-scale signed noise
estimates; each exported peak reports `sigma_um`/`response` from the scale with
the maximal normalized response (a size estimate) and
`detection_snr`/`detection_sigma_um` from the scale whose SNR triggered the
detection, and the two can differ. Peaks within three filter widths of the
image edge are not reported by the LoG detector.

### 5. Classification

Count peak centers strictly inside each selected region and apply the
per-gene ≥3 criterion. Nuclear counts are preserved alongside expanded counts.
All profile parameters are saved with each run.

## Controls and quality limits

Calibrate masks and spot thresholds with manual annotations, representative
positive fields, and negative-control slides such as dapB under matched imaging
conditions. A near-zero negative-control positive rate supports a threshold but
is not, by itself, a complete accuracy assessment. Fix parameters across
comparable samples; never tune the desired overlap percentage.

There is no automated cell-boundary segmentation, spectral unmixing, chromatic
registration or Cy3 subtraction. Subpixel/merged puncta and saturated clusters
cannot be resolved reliably. Counts include detected peaks, which can include
artifact or cluster maxima. No saturated-area-to-transcript substitution is used.

Saturation first uses `qc.saturation_value` when set; otherwise per-channel
CZI `ComponentBitCount`, falling back to image-level component bits and finally
the integer storage ceiling. Valid 12/14-bit metadata therefore give 4095/16383
for 16-bit containers. Missing/invalid bit metadata fall back to storage range;
floating-point images require an explicit physical ceiling. Metadata and the
resolved threshold are recorded. No intensity values or peak thresholds change.

All QC reasons remain in `qc_flags` and separate boolean `flag_*` columns.
`needs_review` is restricted to accepted objects with a count within ±1 of a
gene cutoff, or saturated RNA for a gene whose count is below cutoff.
`review_reasons` lists the actionable conditions. DAPI saturation, irregular
shape, or a saturated RNA-positive object above cutoff remain informational
unless a priority condition also applies. DAPI saturation is evaluated inside
the nuclear mask; RNA saturation and cluster flags inside the counting region. This is a triage rule, not a claim
that other masks/classes cannot be wrong. Review all Cre-positive/Pnoc-below-
cutoff objects in the dedicated pages, plus representative positives/negatives.

Queued objects have `review_priority=cre_endpoint` when currently Cre-positive,
or `general` otherwise; unqueued/excluded objects have `none`. `review_queue.csv`
lists all queued objects in that order, with descending combined peak count
and then nucleus ID for ties. Summary field `nuclei_cre_endpoint_review` counts
the first tier. This rule covers near-Cre denominator changes as well as Pnoc
classification uncertainty. The general tier also matters because below-cutoff
Cre objects can cross into the denominator after review. Class-example crops
are retained and ordered by tier; every Cre-only object still gets dedicated
review pages. Triage flags never alter assignments, counts or expression calls.

Native TIFF labels retain unsigned 32-bit IDs and pixels/cm calibration.
ImageJ companions (`*_labels_imagej.tif`) use direct µm calibration: uint16
when IDs fit, otherwise float32 for IDs up to 2^24, whose integer values are
exactly representable. Beyond that bound the optional companion is omitted
with a logged diagnostic; native uint32 masks still preserve IDs. No renumbering
or rounding is performed. These are label masks, not intensity images. Spot CSV `nucleus_id` is the raw
containing nucleus, `assignment_id` is the region owner, and `counted` indicates
contribution to an accepted object. Excluded objects never enter summary totals.

## Runtime and reproducibility

The main workflow uses two pinned environments on Linux CPU with Python 3.12:
the analysis environment (`requirements-lock.txt`, NumPy ≥2.1 for the CZI
reader) and the Cellpose environment (`requirements-cellpose.txt`, NumPy <2.1),
connected through numeric NPY files and a JSON parameter/metadata exchange.
The unit suite also passes on Python 3.13 in the analysis environment. The
historical StarDist environment was tested on Linux CPU and Python 3.12.
Windows and Apple-silicon macOS have not been validated here. The pinned
TensorFlow 2.21 dependency has no Intel-macOS wheel, so that environment is not
an installation target for Intel Macs.
Compatible TensorFlow installations may use different distribution names;
version recording falls back to the imported module's version and then
`unknown`, rather than aborting analysis. Missing capabilities/model failures
are reported once during preflight, before reading image pixels.

Model downloads use the publisher's registered archive checksum and normal
TLS verification. `KERAS_HOME` selects the cache; otherwise `.cache/keras/`
beside the script is used. Models, input hashes, source hash, resolved profiles,
dependency versions and numerical detector diagnostics are recorded. Decode
and TensorFlow thread defaults are bounded; model tiling remains configurable
through `nuclei.stardist_tiles` and is not selected automatically. The Cellpose
worker inherits thread limits from the runner's environment and is bounded by
`nuclei.cellpose_timeout_s`; CZI scaling entries without a usable value are
treated as missing calibration.

QC rendering works from 8-bit display copies of each channel and hands
matplotlib only crop windows or block-reduced overview panels, so its memory
follows the native pixel count modestly (about 0.2 GB peak for a 2048 x 2048
scene) rather than dominating the analysis. Rendering never feeds back into
segmentation, detection or counts.

Corrupt/unsupported images produce failed summary rows and a nonzero batch
exit, while remaining images continue. Verify completeness before using totals.
Batch percentages use summed counts; zero denominators are NA/null. Overlapping
fields are not deduplicated, and biological grouping must come from the study.
