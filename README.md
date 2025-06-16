# Preclinical Radiomics Pipeline

A fully reproducible, open‑source workflow for small‑animal MRI/CT **radiomics**, **dose mapping**, and **statistical analysis**. The code accompanies the manuscript **“Automated pre‑clinical radiomics pipeline for hippocampal dose‑response studies”** and lets anyone clone the repo, install one environment, and reproduce every table/figure in the paper.

---

## 📁 Repository Layout

```
preclinical-radiomics-pipeline/
├── src/rad_pipeline/           # importable Python package (pip‑installable)
│   ├── io/
│   │   ├── dose_conversion.py  # LET→RT‑Dose DICOM
│   │   └── dicom_slices.py     # split enhanced CT volume
│   ├── preprocessing/
│   ├── features/radiomics.py   # feature extraction pipeline
│   ├── segmentation/mr_dose.py # MRI + dose segmentation & stats
│   └── analysis/pca.py         # PCA scatter/loadings figure
├── scripts/                    # thin CLI wrappers (installed as console‑scripts)
├── tests/                      # unit + regression tests (pytest)
├── environment.yml             # mamba/conda lock‑file
├── Dockerfile                  # one‑command reproducible container
└── README.md                   # you are here
```

---

## 🚀 Quick Start

```bash
# 1 — clone & install (Conda)
git clone https://github.com/<your_user>/preclinical-radiomics-pipeline.git
cd preclinical-radiomics-pipeline
mamba env create -f environment.yml
conda activate preclinical-radiomics
pip install -e .   # editable install for scripts/rad_pipeline

# 2 — run unit tests (optional)
pytest -q

# 3 — slice an enhanced CT DICOM volume
slice-volume  input_CT_enhanced.dcm  --orientations original ray

# 4 — convert LET map to RT‑Dose DICOM
dose-convert   template_RD.dcm dose.npy let.npy  out_let_map.dcm

# 5 — extract radiomic features from two mice
extract-features  \
  --mice-dir data/Manually_preprocessed_and_verified \
  --mice Mouse27_Verified Mouse30_Verified           \
  --out-dir results/features

# 6 — segment dose and compute DVH stats
segment-dose  \
  --mr Mouse33_MR_BSpline_to_CT.nii                \
  --dose Mouse33_Dose_Mousehead.nii                \
  --label Mouse33_Atlas_Registered_to_MR.nii        \
  --hierarchy-csv DSURQE_mapping.csv               \
  --export-nifti                                   \
  --out results/dose_stats.csv

# 7 — recreate PCA figure from manuscript data
pca-plot --csv results/features/combined_features_all_mice.csv \
         --out results/pca_scatter_and_loadings.png
```

---

## 🛠️ Installation Options

### 1. Conda/Mamba (recommended)

```bash
mamba env create -f environment.yml
conda activate preclinical-radiomics
pip install -e .
```

### 2. Docker

```bash
docker build -t preclinical-radiomics .
docker run --rm -it -v $PWD:/workspace preclinical-radiomics
```

---

## 📄 License

Released under the **MIT License**. See `LICENSE` for details.

---

## ✏️ Citation

If you use this pipeline in academic work, please cite both the paper and the repository:

```bibtex
@article{your2025radiomics,
  title   = {Automated pre‐clinical radiomics pipeline for hippocampal dose‐response studies},
  author  = {Your Name and Co‑authors},
  journal = {To be announced},
  year    = {2025},
  note    = {Code available at https://github.com/<your_user>/preclinical-radiomics-pipeline}
}
```

---

## 📧 Contact

Questions or suggestions? Open an issue or drop us an email at **[your.email@university.edu](mailto:your.email@university.edu)**.

---

