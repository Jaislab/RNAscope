# RNAscope: DAPI-defined Pnoc and Cre cells

Batch analysis of all `.czi` / `.CZI` images in one folder. The supplied
`settings.json` uses **StarDist nuclei plus a non-overlapping 2 µm surrounding
region**. A cell is Pnoc-positive with **at least 3 detected ATTO
520 peaks** and Cre-positive with **at least 3 detected Cy5 peaks** in its region.
Both criteria must be satisfied for a double positive; the two genes' peaks do
not need to overlap each other.

These are **nucleus-defined regions**, not measured whole-cell boundaries.
The region can include cytoplasm and neighboring processes. Review this
assignment in crowded tissue. Nuclear-only counts are retained alongside the
region counts in each cell's CSV row. Cy3 is a visual reference only.

## Example result

On `20X_L_CZI.czi`, using the supplied profile:

| Measurement | Result |
|---|---:|
| Accepted DAPI nuclei | 1,554 |
| Pnoc-positive cells | 371 |
| Cre-positive cells | 121 |
| Double-positive cells | 109 |
| Cre-positive cells also expressing Pnoc | **90.08%** |
| Pnoc-positive cells also expressing Cre | 29.38% |

With the same masks and detector, nucleus-only assignment gives 87.64% of
Cre-positive nuclei also positive for Pnoc. These are **provisional automated
results**, requiring review of segmentation, background and saturated clusters.
See [validation and comparison](VALIDATION.md) for the full evidence and limits.

## Install and run

Use **Python 3.12** (the tested version). From this folder on Linux/macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-stardist.txt
python rnascope.py --input /path/to/czi_folder --output /path/to/new_results --config settings.json
```

On Windows, create the environment with `py -3.12 -m venv .venv` and activate
with `.venv\Scripts\Activate.ps1` in PowerShell. Quote paths containing spaces.
The pinned StarDist/TensorFlow stack was executed on Linux CPU; the script does
not require a GPU. First model use downloads the official
`2D_versatile_fluo` weights from the StarDist GitHub release. The publisher's
registered archive checksum and TLS verification are preserved. Model files
are cached under `.cache/keras/` beside the script unless `KERAS_HOME` is set.
The code records hashes of the model files used. TensorFlow installation is
substantially larger than the base image-analysis dependencies.

Use `--recursive` for subfolders. Input images are read-only. Multiple scenes
are analyzed separately. Non-singleton Z, T or other extra dimensions are
rejected, never projected or silently discarded. Image decoding uses at most
two threads; TensorFlow defaults to two intra/inter-op threads. No web server,
database, API credentials or GPU are needed.

Use a new or empty output folder for each run. A bad file/scene is logged and
remaining files continue. Exit 0 means all selected scenes were processed;
exit 1 means one or more failed; exit 2 means invalid arguments/settings.
Processing success is not biological validation.

## Available settings

- **`settings.json`**: newer uploaded script's approach to DAPI preprocessing
  (smoothing, local contrast enhancement, StarDist, brightness relative to a
  surrounding ring), intensity peak detection, and 2 µm assignment. No saturated
  area-to-punctum conversion. This is the current workflow.
- **`settings_nucleus_only.json`**: the same segmentation and detection with
  assignment expansion set to zero, for a direct comparison.
- **`settings_legacy.json`**: the preceding Li/watershed and multiscale LoG
  nucleus-only analysis, preserved for reproducibility.

For the legacy workflow, `requirements-lock.txt` supplies the smaller dependency
set. Omitting `--config`, or using `--write-default-config my_settings.json`,
selects the lightweight built-in watershed/LoG, nuclear-only defaults, **not**
the supplied StarDist/2 µm profile. Specify `--config settings.json` for the
current workflow. Partial JSON overrides are supported; unknown keys fail.
StarDist is selected explicitly and never silently replaced by a fallback.

`nuclei.method` can be `watershed`, `local_watershed`, or `stardist`.
`spots.<gene>.method` can independently be `log` or `intensity`.
`assignment.expand_um` sets the assignment radius; zero uses nuclear masks
only. Per-gene `positive_min_spots` values remain 3 in all supplied profiles.
There is no target overlap percentage in the analysis.

## Channels and pixel size

| Role | Metadata name / fluorophore | Fallback 0-based index |
|---|---|---:|
| Cre | Cy5-T1 / Cy5 | 0 |
| Pnoc | At520-T3 / ATTO 520 | 2 |
| DAPI | DAPI-T4 / DAPI | 3 |

Channel names/fluorophores take precedence over fallback indices. Fallbacks are
used only when names are missing, with a warning. Unrecognized or ambiguous
named channels fail; after checking probe identity, set `channels.<role>.index`
explicitly if needed. Metadata identifying a dye do not prove the biological
probe identity.

X/Y pixel sizes come from CZI metadata. Missing calibration requires an
explicit `pixel_size_um: [y_um, x_um]` override. Overrides are recorded.
Coordinates are local to the scene's top-left, not microscope stage coordinates;
x is column and y is row. Assignment distances account for anisotropic pixels.

## Current algorithm

1. Smooth DAPI by 1.4 µm, normalize, and enhance local contrast with CLAHE
   (25 µm window, clip 0.01). StarDist predicts individual nucleus outlines at
   a target sampling of 0.69 µm/pixel. Its default confidence cutoff is used;
   a configured override is recorded. Neither RNA channel influences this step.
2. Accept nuclei of 20–300 µm² that do not touch the image border and whose
   smoothed DAPI intensity is at least 1.15 times the local background in an
   available 1.5–3.5 µm ring. A missing ring is unknown and does not reject the
   object. Rejected candidates remain in the CSV/label image with reasons.
3. Expand accepted nuclei by 2 µm using nearest-label assignment, without
   overlap. Pixels inside rejected DAPI objects are blocked: a rejected nucleus
   is not treated as cytoplasm of an accepted neighbor. This safeguard can make
   results differ slightly from the uploaded script.
4. For each RNA channel, remove morphological background at an 8 µm radius,
   smooth at approximately 0.55 µm, then find intensity maxima above median plus
   10 times a robust noise estimate (`1.4826 × MAD`). Minimum peak separation
   uses a 1 µm physical footprint. At the example sampling these reproduce the
   uploaded script's smoothing and 3×3 local-maximum neighborhood.
5. Count detected peak centers in each assignment region. Apply the independent
   ≥3 criteria. Also count the same peaks strictly inside each nucleus.

The legacy option uses background-corrected Li/Otsu segmentation, distance
watershed, and scale-normalized multiscale LoG peak detection, with separate
signed-filter noise estimates. Its parameters remain unchanged. The optional
local-watershed method follows the uploaded threshold fallback and excludes
low-solidity candidates. The new StarDist preprocessing was adopted from the
**newer** uploaded script, not its superseded version.

## Saturation, background and interpretation

All counts are provisional until compared with manual annotations and assay
controls. Three peaks is an expression rule, not background rejection. Check
DAPI masks, noise-only fields, negative controls, crowded fields, bright blobs,
and near-cutoff cells before processing an experimental cohort. Fix one set of
parameters for comparable acquisitions; do not tune an overlap percentage to
match biological expectations.

The example has approximately 0.69 µm pixels and saturated signal in both RNA
channels. Individual subpixel or merged puncta cannot be resolved reliably.
The script flags saturated or clustered regions and counts detected peaks;
these can still include artifact/cluster maxima and are not guaranteed
individual transcripts. It does **not** infer hidden puncta from saturated area.
The newer uploaded script's area/3.3 µm² estimate was tested separately for
comparison and is labeled an estimate in the validation report. Observed image
maximum is not automatically considered detector saturation: the detector limit
is taken from the integer data type, or `qc.saturation_value` explicitly.

Flags include irregular nuclear shape, saturation, large clusters and counts
within one peak of a cutoff. Flagged objects stay in provisional summaries.
No flag means only that these automated checks found none. There is no automatic
spectral unmixing, channel registration, or Cy3 subtraction. Different gains and
laser settings prevent direct comparison of absolute transcript abundance
across RNA channels without additional calibration.

## Outputs

At the results root:

- `summary.csv`: scene counts, both overlap percentages, assignment radius,
  warnings and failures. Failed counts are `NA`.
- `batch_summary.csv` / `.json`: sums across successful scenes, completeness
  and failures. Partial batches must not be presented as complete datasets.
- `config_used.json`, `run_manifest.json`, `run.log`: settings, input hashes,
  script hash/version, dependency versions, run times and diagnostics.

Under `images/<relative-CZI-filename>/scene_<id>/`:

- `nuclei.csv`: ID, centroid, nuclear area, shape, DAPI contrast ratio,
  acceptance/exclusion reason, `pnoc_nuclear_spots` / `cre_nuclear_spots`,
  assigned-region `pnoc_spots` / `cre_spots`, positivity, class and QC flags.
  Excluded candidates retain measurements but never enter summary counts.
- `spots.csv`: all detections, including unassigned ones. `nucleus_id` is the
  raw nuclear mask containing the peak (0 = none); `assignment_id` is the owner
  after expansion. `counted` says whether the peak contributed to an accepted
  region. `response` and `detection_snr` depend on the selected detector.
  `sigma_um` is available for LoG only, not an intensity-peak size measurement.
- `nuclei_labels.tif`: 32-bit IDs including rejected candidates; background 0.
  `assignment_labels.tif` is also saved when expansion is nonzero. TIFF
  calibration uses pixels/cm, with physical µm sizes in metadata.
- `qc_overview.png`, `review_crops.png`: nuclear and assignment boundaries,
  detected spots, classes, and Cy3 reference. White boundaries are nuclei and
  cyan boundaries are expanded regions. Cre is red, Pnoc green and DAPI blue.
  Display contrast does not change analysis intensities. Review crops are a
  sample of the dataset, with IDs in `review_crop_ids.json`.
- `cre_only_review.csv`, `cre_only_review.pdf`, `cre_only_review_XX.png`:
  **every accepted Cre+/Pnoc− object**, with separate grayscale DAPI/Pnoc/Cre
  panels, nuclear/assignment outlines and counted peaks. PDF/PNG files are
  created only when such objects exist. Their outlines use yellow for nuclei
  and cyan for assignment regions.
- `metadata.json`, `summary.json`: resolved channels, model hashes, calibration,
  numerical noise/threshold diagnostics, assignment policy and scene totals.

Both positive populations include double positives:

```text
Pnoc expressing Cre (%) = 100 × double_positive / pnoc_positive
Cre expressing Pnoc (%) = 100 × double_positive / cre_positive
```

Zero denominators are `NA` in CSV / `null` in JSON. Batch percentages use summed
counts, not the average of image percentages. Images are not automatically
biological replicates; group by animal/sample using your study design.
Overlapping fields are not deduplicated.

## Verification

```bash
python -m unittest discover -s tests -v
```

The 19 tests cover known peaks and noise, contrast filtering, explicit backend
selection, touching/border nuclei, channel identity, stack rejection, thresholds,
nuclear versus expanded assignments, rejected-object blocking, anisotropic
spacing, percentages, failures and exports. StarDist routing is mocked in the
unit tests; actual model inference was separately executed on the uploaded CZI.
See `VALIDATION.md` for current results and the method comparison. Software
checks do not establish biological sensitivity or specificity.
