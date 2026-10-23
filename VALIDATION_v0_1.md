# Validation of this implementation

The Python 3.12 environment installs successfully from `requirements-lock.txt`.
Rerunning the saved installation script succeeds without changing the source
or settings. All **15 software tests pass**, including a synthetic image with
known nuclei, known spot locations, and known overlap percentages.

The uploaded `20X_L_CZI.czi` was processed with the supplied `settings.json`.
Its metadata identify a single 1024 × 1024 section, four 16-bit channels, and
0.690535481770835 µm pixels in X and Y. Every exported spot-to-nucleus assignment
was checked against the saved label image. Per-nucleus spot counts, aggregate
counts, TIFF calibration, configuration identity, and the script hash in the
run manifest were also checked.

These are **provisional automated nuclear-associated results**, before manual
annotation or negative-control calibration:

| Measurement | Value |
|---|---:|
| Accepted nuclear masks | 1,700 |
| Excluded watershed candidates | 2,069 |
| Pnoc-positive (≥3 detected spots) | 279 |
| Cre-positive (≥3 detected spots) | 60 |
| Double positive | 45 |
| Pnoc only | 234 |
| Cre only | 15 |
| Neither | 1,406 |
| Pnoc-positive expressing Cre | 16.13% |
| Cre-positive expressing Pnoc | 75.00% |
| Accepted nuclei with one or more review flags | 1,062 |

Excluded candidates are predominantly small fragments/debris; they are not an
estimate of additional biological cells. Nuclear segmentation and true punctum
counts remain to be checked by an expert. Both RNA channels contain saturated
pixels. Review flags include near-cutoff counts as well as image-quality issues;
flagged nuclei remain included in the provisional totals.

Visual inspection compared Li and Otsu nuclear thresholding on the example.
Li followed the broader DAPI signal better, while Otsu missed dim nuclear
regions. This supports the provisional default; it is not an expert-annotated
accuracy benchmark. No biological sensitivity, specificity, or absolute
transcript-count accuracy is claimed. Cy3 is a reference panel only, with no
automatic artifact exclusion or spectral correction.

A separate recursive batch smoke check used two copies of the same CZI, one
with uppercase `.CZI`, plus an intentionally corrupt file. Both valid images
completed with identical counts, the invalid image was reported, aggregate
counts equaled twice the single-image counts, completeness was false, and the
process returned exit 1. These duplicate-image outputs were used only to check
batch behavior and are not included in the example biological results.

Reproduce the example on this cloud machine using a new output directory:

```bash
cd /workspace/RNAscope
MPLCONFIGDIR=/tmp/rnascope-matplotlib .venv/bin/python rnascope.py \
  --input /workspace/attachments/b4183988-0b21-4f53-9763-7c1e1d17dae7 \
  --output results_example_repeat \
  --config settings.json
```

The bundle includes code, settings, dependencies, tests and example results;
the original CZI and installed virtual environment are not included.
