# RNAscope: Pnoc association in Cre-positive nuclei

Batch analysis of single-section Zeiss CZI images. The script segments DAPI
nuclei, detects Pnoc (ATTO 520) and Cre (Cy5) RNA peaks, counts peaks in each
nucleus-associated region, and reports **the fraction of Cre-positive objects
with detectable Pnoc RNA signal**:

```text
pct_cre_expressing_pnoc = 100 × double_positive / cre_positive
```

Cre is the **sparser detected signal** in the reference image: 1,209 Cre peaks
versus 3,531 Pnoc peaks. Cells expressing Cre weakly can fall below the ≥3-peak
Cre cutoff, increasing the apparent Pnoc-only group and depressing the reverse
percentage. The Cre-conditioned endpoint focuses on the detected Cre-positive
subset; it does not recover missed Cre cells or guarantee an unbiased census.

Detected peak positions do **not** establish a strong localization contrast:
66.7% of Cre peaks and 65.4% of Pnoc peaks lie inside candidate DAPI nuclear
masks. DAPI supports a nuclear-region measurement; assigning whole cytoplasms
to individual cells requires additional boundary information. See the
[measured fractions and interpretation](VALIDATION.md#detected-peak-distribution).

## Status

The main profile, `settings_cellpose.json`, segments raw DAPI with the pinned
legacy Cellpose 3.1.1.3 `nuclei` model and counts peaks in a non-overlapping
2 µm expansion of each accepted nucleus. On the reference image it reproduces
409 accepted nuclei, 89 Pnoc-positive, 38 Cre-positive and 36 double-positive
objects (94.74% of Cre-positive objects).

**The main profile's object count is unvalidated against expert labels.** On
the same image the StarDist profile accepts 1,554 nuclei; only 124 of those
match a Cellpose object at IoU ≥ 0.5, and 285 of the 409 Cellpose objects match
no StarDist nucleus. The 94.74% value rests on 38 Cre-positive objects. No
expert annotation has established which segmentation is closer to the true
nuclei, and the
[blinded review package](docs/CELLPOSE_COMPARISON.md#blinded-review-package)
has not been scored. A higher coexpression percentage does not establish better
accuracy. Treat the main profile's counts and percentages as provisional, and
complete expert review on representative images before interpreting scientific
results.

## Installation

Linux CPU with Python 3.12 is the supported setup; the unit suite also passes
on Python 3.13 in the analysis environment. Cellpose 3.1.1.3 requires NumPy
<2.1 and the CZI reader requires NumPy ≥2.1, so the main workflow uses two
environments: the batch runner exchanges numeric NPY arrays with an isolated
Cellpose subprocess. Do not install both dependency sets together.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
python3 -m venv .venv-cellpose
.venv-cellpose/bin/python -m pip install -r requirements-cellpose.txt
```

| Variable | Purpose | Default |
|---|---|---|
| `CELLPOSE_PYTHON` | Absolute path of the Cellpose environment's interpreter | `.venv-cellpose/bin/python` beside `rnascope.py` |
| `CELLPOSE_LOCAL_MODELS_PATH` | Cellpose weight cache | `.cache/cellpose` beside the scripts |
| `KERAS_HOME` | StarDist/Keras model cache (historical profiles only) | `.cache/keras` beside the scripts |

The `nuclei` weights are downloaded once over verified HTTPS and must match the
SHA-256 fingerprint recorded in [the comparison report](docs/CELLPOSE_COMPARISON.md);
that fingerprint is a reproducibility check, not a publisher-issued checksum.
The historical StarDist profiles need `requirements-stardist.txt` in a third,
separate environment. Missing dependencies, weights or worker failures abort
the run; no method is substituted silently.

## Running an analysis

```bash
.venv/bin/python rnascope.py --input /path/to/czi_folder --output results_cellpose_main --config settings_cellpose.json
```

- `--config` is required; no profile is selected implicitly.
  `--write-default-config my_settings.json` writes a copy of the main profile
  to edit.
- `--output` must be a new or empty folder. Input CZI files are opened
  read-only. Add `--recursive` to include subfolders; their relative paths are
  kept in the output tree.
- Each scene is processed separately. A corrupt or unsupported file or scene
  produces a failed summary row and a nonzero exit status while the remaining
  images continue, so check `batch_summary.json` for `complete: true` before
  using totals.
- Only single optical sections are accepted; Z, T or other extra axes of size
  >1 are rejected rather than projected. The Cellpose profile additionally
  requires square XY sampling (within 1%).
- X/Y calibration comes from the CZI scaling metadata; when it is missing or
  unusable, set `pixel_size_um: [y_um, x_um]` in the profile.

| Profile | Nuclei | Counting region | Role |
|---|---|---|---|
| `settings_cellpose.json` | Cellpose `nuclei`, raw DAPI | Nucleus + 2 µm | Main workflow |
| `settings.json` | StarDist `2D_versatile_fluo`, CLAHE | Nucleus + 2 µm | Historical reference (1,554 nuclei) |
| `settings_nucleus_only.json` | StarDist `2D_versatile_fluo`, CLAHE | Nuclear mask only | Historical nuclear-overlap measurement |
| `settings_legacy.json` | Li/Otsu watershed | Nuclear mask only | Historical v0.1 method with LoG peaks |

All profiles use the same channel mapping, the independent ≥3-peak criterion
per gene, the same QC and the same exports. The active method and expansion are
printed at startup and recorded in every summary and manifest.

## Cellpose profile settings

The `nuclei` section of `settings_cellpose.json` controls segmentation and
acceptance; the RNA detector and assignment settings are shared with the
StarDist +2 µm profile.

| Key | Value | Meaning |
|---|---:|---|
| `cellpose_diameter_px` | 17 | Nuclear diameter passed to Cellpose in pixels, never via a rounded µm conversion. Values outside 3–500 are rejected as probable unit mistakes. The µm equivalent (about 11.74 µm on the reference image) is recorded in scene metadata. |
| `cellpose_threads` | 2 | CPU threads; also exported to the worker's environment as OMP, OpenBLAS, MKL and Numba limits unless those variables are already set. |
| `cellpose_timeout_s` | 600 | Limit for one worker call, including the first weight download. Raise it for large scenes. Partial downloads left by a killed worker are removed on the next start. |
| `min_area_um2` / `max_area_um2` | 20 / 300 | Accepted nuclear area. |
| `min_dapi_contrast` | 1.15 | Minimum ratio of the nucleus's smoothed DAPI mean to the mean of a 1.5–3.5 µm surrounding ring; an empty ring is unknown, not a rejection. |
| `contrast_smooth_um` | 1.4 | Smoothing scale for that contrast ratio, independent of StarDist preprocessing. |
| `exclude_border` | true | Reject objects touching the image edge. |

Fixed Cellpose inference settings (1st/99th-percentile normalization, flow
threshold 0.4, cell-probability threshold 0, minimum mask size 15 pixels, no
augmentation) are defined once in `cellpose_backend.py`. The batch runner sends
the complete keyword set to the worker, the worker echoes exactly what it ran,
and masks whose echo differs from the request are refused. Scene metadata
therefore records the model, loaded-weight fingerprint, package versions,
diameter and the inference parameters that produced the masks; run manifests
record the hashes of both scripts. Different image sampling may require a
different pixel diameter; these are starting settings, not validated universal
parameters.

## Endpoint and biological limits

In each counting region (the nucleus, or the nucleus plus 2 µm in the main
profile), an object is Cre-positive with **≥3 detected Cy5 peaks** and
Pnoc-positive with **≥3 detected ATTO 520 peaks**. Both genes must meet their
thresholds in the same region; their individual puncta do not need to coincide.
This measures nucleus-associated RNA signal, not molecular colocalization or
complete whole-cell expression.

A below-cutoff object can belong to a cell expressing that gene weakly or
outside the region. The script therefore **does not estimate the total
Pnoc-positive or Cre-positive cell population**. The reverse percentage,
`pct_pnoc_expressing_cre`, is particularly sensitive to missed low-Cre cells
and should not be used as the fraction of all Pnoc cells expressing Cre.
Neither conditional percentage is automatically an unbiased population estimate.

The 2 µm expansion applies to **both genes**. It can include cytoplasm or
neighboring processes, and it can change which objects meet the Cre cutoff.
Its result must not be labeled nuclear overlap; nuclear-only counts are kept
in the cell table for comparison. Rejected nuclear pixels are never assigned
to a neighboring accepted nucleus.

The Cre-conditioned fraction describes **detected RNA coexpression**. It can
inform Cre-driver characterization, but does not directly measure functional
recombination, reporter specificity, or the efficiency of labeling all Pnoc cells.

## Outputs

Batch level:

- `summary.csv`: one row per scene with counts, both percentages, the active
  profile, assignment radius, warnings and errors; rewritten after every file.
- `batch_summary.csv` / `.json`: summed counts, count-weighted percentages,
  numbers of successful and failed items and `complete`.
- `config_used.json`, `run_manifest.json`, `run.log`: the merged profile, script
  and worker hashes, package versions, input file hashes, timings and per-scene
  diagnostics.

Per scene, under `images/<CZI path>/scene_<id>/`:

- `nuclei.csv`: every candidate mask with centroid, area, shape, DAPI contrast,
  acceptance and exclusion reason, nuclear and assigned-region peak counts,
  classification, review tier and one `flag_*` column per QC reason.
- `spots.csv`: every detected peak with position, scale (`sigma_um`), raw
  response, detection SNR, the scale whose SNR triggered the detection
  (`detection_sigma_um`), contrast, containing nucleus, assigned object and
  whether it was counted.
- `review_queue.csv`: queued objects, `cre_endpoint` tier first.
- `nuclei_labels.tif`, `assignment_labels.tif`: uint32 label masks with
  pixels-per-centimeter calibration; `*_labels_imagej.tif`: Fiji/ImageJ
  companions calibrated directly in µm with unchanged IDs.
- `annotated_overlay.png` / `.pdf`: one full-image overlay at native resolution
  or better, with accepted nucleus IDs matching the CSV. The PNG draws IDs and
  markers into the pixels; the PDF embeds the native-resolution image with
  vector, zoomable labels. `annotated_overlay_legend.json` explains the colors.
- `qc_overview.png`, `review_crops.png`, `review_crop_ids.json`: channel panels,
  boundaries, detections and class examples. Overview panels of scenes larger
  than about 1,200 pixels are block-averaged to the panel size; crops and
  overlays use full-resolution pixels.
- `cre_only_review.csv` / `.pdf` / PNG pages: every accepted Cre-positive object
  below the Pnoc cutoff, for inspection of missed signal, background and saturation.
- `metadata.json`: pixel size, channel metadata and resolved indices,
  segmentation and detection diagnostics, assignment and review policies.

All QC images use 8-bit display scaling; analysis uses the original
intensities. Zero denominators are `NA` in CSV and `null` in JSON. Images and
nuclei are not automatically independent biological replicates; group results
by animal or sample according to the experimental design.

## Review workflow

`needs_review` marks accepted objects within ±1 peak of either cutoff, or
below a gene's cutoff with saturated RNA signal for that gene. `review_reasons`
lists those conditions. Informational flags (DAPI saturation, irregular mask
shape, bright clusters, saturation above cutoff) are kept in the `flag_*`
columns but do not queue an object on their own. DAPI flags are evaluated
inside the nuclear mask; RNA saturation and cluster flags inside the counting
region. This is a triage heuristic, not an error probability or exclusion
rule, and it does not capture every possible segmentation problem.

`review_priority` is `cre_endpoint` for queued Cre-positive objects, `general`
for other queued objects, and `none` otherwise. On the StarDist +2 µm reference
run, 77 of 413 queued objects are in the Cre-endpoint tier and 22 of those are
near the Pnoc cutoff. Near-Cre-cutoff calls affect the denominator, and the
general tier matters because below-cutoff Cre objects can become positive after
review. `nuclei_cre_endpoint_review` in the summaries counts the first tier.
Review crops keep one example per class and then follow the queue order. Every
Cre-only object gets a dedicated review page whether or not it is queued.

Saturation uses `qc.saturation_value` when set, otherwise the CZI
`ComponentBitCount`, otherwise the integer storage limit. Saturated clusters
are flagged, never converted to estimated counts.

## Historical StarDist reference results

The reference single-section `20X_L_CZI.czi` (1024 × 1024 pixels, 0.69 µm per
pixel) produced, with the StarDist profiles:

| Measurement | Nuclear mask | Nuclear mask + 2 µm |
|---|---:|---:|
| Accepted nuclei | 1,554 | 1,554 |
| Cre-positive objects | 89 | 121 |
| Double-positive objects | 78 | 109 |
| Cre-positive with detectable Pnoc | **87.64%** | **90.08%** |

These are provisional detected-peak results. The two percentages use different
Cre-positive denominators; the expanded value is not a more accurate
nuclear-overlap result merely because it exceeds 90%. No target overlap
percentage enters the algorithm. The StarDist profiles and ≥3 criteria are
unchanged since v0.2.0. See [methods](docs/METHODS.md) and
[validation](VALIDATION.md) for controls, assumptions and reproducibility.

## Comparison and sensitivity tools

**Cellpose comparison.** [docs/CELLPOSE_COMPARISON.md](docs/CELLPOSE_COMPARISON.md)
runs the legacy `nuclei` and `cyto3` models on the reference image with frozen
peak detections and the ≥3 criteria, including DAPI preprocessing controls,
Pnoc/DAPI cell-outline estimates, per-object disagreements and a blinded
annotation package. `prepare_cellpose_input.py` and `compare_cellpose.py` run
in the analysis environment and accept scene folders from `--recursive` runs;
`cellpose_infer.py` runs in the Cellpose environment. Pnoc-guided cell
estimates link only 20 of the 121 StarDist Cre-positive objects; they are not
validated whole-cell boundaries and the main profile does not use them.

**Parameter sweep.** `sweep_parameters.py` reuses a complete run's masks and
peaks to vary the per-gene cutoffs and counting radius without rerunning
segmentation:

```bash
.venv/bin/python sweep_parameters.py --results results_expanded --output results_sweep
```

The [324-setting sensitivity analysis](docs/SENSITIVITY.md) shows why a higher
percentage can reflect a smaller Cre-positive denominator. Pnoc ≥2 / Cre ≥3 at
2 µm gives 114/121 (94.21%); that is a candidate for manual validation, not an
accuracy optimum or a change to any shipped profile.

## Tests and continuous integration

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python tests/smoke_cellpose.py
```

The unit suite covers configuration validation and every shipped profile,
channel resolution, axis handling, calibration metadata, nuclear segmentation
and acceptance filters, both peak detectors, assignment and classification,
review tiers and flags, saturation handling, exports including ImageJ
companions and the annotated overlay, a complete batch through the reader
interface with nested input folders followed by a sweep of its exports, the
Cellpose worker contract (parameter echo, weight fingerprint, timeout and
thread limits) with mocked inference, the comparison tools and the sweep. It
runs on every push through GitHub Actions
(`.github/workflows/tests.yml`, Linux, Python 3.12).

The smoke check needs the Cellpose environment. It downloads the verified
weights if necessary, runs real inference on four synthetic nuclei and checks
classifications, the 2 µm expansion, provenance and CSV/TIFF/PDF exports. It
does not validate real-image segmentation accuracy or CZI decoding. Passing
tests are software evidence, not a biological sensitivity or specificity
benchmark.

## Further documentation

- [docs/METHODS.md](docs/METHODS.md): measurement definition, image reading,
  segmentation and detection steps, controls, QC rules and reproducibility.
- [VALIDATION.md](VALIDATION.md): reference results, detected-peak
  distributions, software evidence and the experimental work still needed.
- [docs/SENSITIVITY.md](docs/SENSITIVITY.md): cutoff and radius sweep.
- [docs/CELLPOSE_COMPARISON.md](docs/CELLPOSE_COMPARISON.md): Cellpose variants,
  method agreement and the blinded review package.
- [docs/AUDIT_REVIEW.md](docs/AUDIT_REVIEW.md): audit findings and dispositions.
- [CHANGELOG.md](CHANGELOG.md).
