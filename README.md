# Preclinical Radiomics Pipeline

A fully reproducible, open‑source workflow for small‑animal MRI/CT **radiomics**, **dose mapping**, and **statistical analysis**. The code accompanies the manuscript **“An Open-Source Irradiation and Data-Handling Framework for Pre-Clinical Ion-Beam Research”**.

---

## Repository Layout

```
preclinical-radiomics-pipeline/
├── src/rad_pipeline/           
│   ├── io/
│   │   ├── dose_conversion.py  # LET array →RT‑Dose DICOM
│   │   └── dicom_slices.py     # split CT volume to slices for the TPS
│   ├── features/radiomics.py   # feature extraction pipeline
│   ├── segmentation/mr_dose.py # MRI + dose segmentation & statistics
│   ├── segmentation/atlas_segment_export.py # Region export from Digimouse or DSURQE atlases
│   └── analysis/pca.py         # PCA scatter/loadings figure
├── environment.yml             
├── Dockerfile                  
└── README.md                   
```

---

## Quick Start

```bash
# 1. clone & install
git clone https://github.com/<your_user>/preclinical-radiomics-pipeline.git
cd preclinical-radiomics-pipeline

# 2. create & activate a clean virtual environment
python3 -m venv .venv
source .venv/bin/activate          # Windows PowerShell:  .venv\Scripts\Activate

# 3. upgrade pip & install the package (editable mode recommended)
python -m pip install --upgrade pip
pip install -e . 

# 4. slice an enhanced CT DICOM volume
slice-volume  input_CT_enhanced.dcm  --orientations original ray
```python .\src\io\dicom_slices.py .\data\slice-volume\20241031094809_CT_ISRATV_0.dcm```


# 5. convert LET map to RT‑Dose DICOM
dose-convert   template_RD.dcm dose.npy let.npy  out_let_map.dcm
```python .\src\io\dose_conversion.py .\data\dose-convert\template_dose.dcm  .\data\dose-convert\dose_npy.npy  .\data\dose-convert\let_npy.npy  .\data\dose-convert\out_let_map.dcm```

# 6. extract radiomic features from two mice
extract-features  \
  --mice-dir data/Manually_preprocessed_and_verified \
  --mice Mouse27_Verified Mouse30_Verified           \
  --out-dir results/features

# 7. segment dose and compute DVH stats
segment-dose  \
  --mr Mouse33_MR_BSpline_to_CT.nii                \
  --dose Mouse33_Dose_Mousehead.nii                \
  --label Mouse33_Atlas_Registered_to_MR.nii       \
  --hierarchy-csv DSURQE_mapping.csv               \
  --export-nifti                                   \
  --out results/dose_stats.csv

# 8. extract Hippocampus from DSURQE atlas (left+right merged)
atlas-segment \
  --atlas dsurqe \
  --label Mouse33_Atlas_Registered_to_MR.nii \
  --hierarchy-csv DSURQE_mapping.csv \
  --region Hippocampus \
  --merge-sides \
  --out-dir results/segments

# 9. extract Cerebellum from Digimouse atlas
atlas-segment \
  --atlas digimouse \
  --label atlas_380x992x208.img \
  --digimouse-map atlas_380x992x208.txt \
  --region cerebellum \
  --out-dir results/segments

# 10. recreate PCA figure from manuscript data
pca-plot --csv results/features/combined_features_all_mice.csv \
         --out results/pca_scatter_and_loadings.png
```

---

## Atlas Resources

This pipeline supports the following preclinical mouse atlases:

### DSURQE Mouse Brain Atlas

- **Source**: [Rodare Repository – Record 915](https://rodare.hzdr.de/record/915)
- **Citation**:  
  Dorr AE, Lerch JP, Spring S, Kabani N, Henkelman RM.  
  *High resolution three-dimensional brain atlas using an average magnetic resonance image of 40 adult C57Bl/6J mice.*  
  **Neuroimage. 2008; 42(1):60–69.**  
  DOI: [10.1016/j.neuroimage.2008.03.037](https://doi.org/10.1016/j.neuroimage.2008.03.037)
  
  Müller, Johannes, Suckert, Theresa, Beyreuther, Elke, Schneider, Moritz, Boucsein, Marc,
  Bodenstein, Elisabeth, … Dietrich, Antje. (2021). Slice2Volume: Fusion of multimodal medical imaging 
  and light microscopy data of irradiation-injured brain tissue in 3D. 
  (Version 0.3.1) [Data set]. Rodare. http://doi.org/10.14278/rodare.915
  

### Digimouse Whole-Body Atlas

- **Source**: [USC Neuroimage – Digimouse Project](https://neuroimage.usc.edu/neuro/Digimouse)
- **Citation**:  
  Dogdas B, Stout D, Chatziioannou A, Leahy RM.  
  *Digimouse: A 3D Whole Body Mouse Atlas from CT and Cryosection Data.*  
  **Phys Med Biol. 2007; 52:577–587.**  
  DOI: [10.1088/0031-9155/52/3/003](http://dx.doi.org/10.1088%2F0031-9155%2F52%2F3%2F003)

  D. Stout, P. Chow, R. Silverman, R. M. Leahy, X. Lewis, S. Gambhir, 
  A. Chatziioannou, Creating a whole body digital mouse atlas with PET, CT and cryosection images,
  Molecular Imaging and Biology.2002; 4(4): S27

---

## Installation Options

### 1.Venv/pip

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows PowerShell:  .venv\Scripts\Activate
python -m pip install --upgrade pip
pip install -e . 
```

### 2.Docker

```bash
docker build -t preclinical-radiomics:latest .
docker run --rm -it -v %cd%:/workspace preclinical-radiomics:latest
```

---

## License

Released under the **MIT License**. See `LICENSE` for details.

---

## Citation

If you use this pipeline in academic work, please cite both the paper and the repository:

```bibtex
@article{your2025radiomics,
  title   = {An Open-Source Irradiation and Data-Handling Framework for Pre-Clinical Ion-Beam Research},
  author  = {Isselmou Abdarahmane, Lorenz Wolf, Peter Kuess, Gerd Heilemann, Silvia Stocchiero, Barbara Knäusl, Ingo Feinerer, Markus Zeilinger, Dietmar Georg},
  journal = {***},
  year    = {2025},
  note    = {***}
}
```


