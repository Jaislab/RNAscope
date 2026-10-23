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

## StarDist profiles

1. Smooth DAPI at 1.4 µm, normalize and apply CLAHE with a 25 µm window and
   clip 0.01. StarDist `2D_versatile_fluo` predicts nucleus outlines at target
   sampling 0.69 µm/pixel. Its default probability threshold is used unless
   explicitly overridden. RNA channels do not determine the masks.
2. Accept nuclei of 20–300 µm², excluding image-border objects and objects whose
   smoothed DAPI mean is below 1.15 times the mean in an available 1.5–3.5 µm
   surrounding ring. An unavailable ring is unknown, not grounds for rejection.
   Rejected candidates remain in the CSV and label image with their reasons.
3. Use nuclear masks directly, or expand accepted masks by the configured
   radius with nearest-label assignment and no overlap. Pixels belonging to
   rejected DAPI objects are blocked. The 2 µm expansion is an approximate
   nucleus-associated region and can include another cell's cytoplasm/processes.
4. Independently remove each RNA channel's morphological background at an 8 µm
   radius, smooth at about 0.55 µm, and detect intensity maxima above median plus
   10 × robust noise (`1.4826 × MAD`). Minimum separation uses a 1 µm footprint.
5. Count peak centers strictly inside each selected region and apply the
   per-gene ≥3 criterion. Preserve nuclear counts alongside expanded counts.

StarDist is explicitly selected and never silently replaced by thresholding.
The optional local-watershed backend uses local DAPI thresholding, watershed
and solidity filtering. The historical profile uses Li/Otsu segmentation and
multiscale, scale-normalized LoG responses with per-scale signed noise estimates.
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
unless a priority condition also applies. This is a triage rule, not a claim
that other masks/classes cannot be wrong. Review all Cre-positive/Pnoc-below-
cutoff objects in the dedicated pages, plus representative positives/negatives.

Native TIFF labels retain unsigned 32-bit IDs and pixels/cm calibration.
ImageJ companions (`*_labels_imagej.tif`) use direct µm calibration: uint16
when IDs fit, otherwise float32 for IDs up to 2^24, whose integer values are
exactly representable. Beyond that bound the optional companion is omitted
with a logged diagnostic; native uint32 masks still preserve IDs. No renumbering
or rounding is performed. These are label masks, not intensity images. Spot CSV `nucleus_id` is the raw
containing nucleus, `assignment_id` is the region owner, and `counted` indicates
contribution to an accepted object. Excluded objects never enter summary totals.

## Runtime and reproducibility

The pinned StarDist environment was tested on Linux CPU and Python 3.12.
Windows and Apple-silicon macOS have not been validated here. The pinned
TensorFlow 2.21 dependency has no Intel-macOS wheel, so this environment is not
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
through `nuclei.stardist_tiles` and is not selected automatically.

Corrupt/unsupported images produce failed summary rows and a nonzero batch
exit, while remaining images continue. Verify completeness before using totals.
Batch percentages use summed counts; zero denominators are NA/null. Overlapping
fields are not deduplicated, and biological grouping must come from the study.
