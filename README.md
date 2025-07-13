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
│   ├── features/radiomics_pipeline.py   # feature extraction pipeline
│   ├── segmentation/mr_dose.py # MRI + dose segmentation & statistics
│   ├── segmentation/atlas_segment_export.py # Region export from Digimouse or DSURQE atlases
│   └── analysis/pca.py         # PCA scatter/loadings figure
├── environment.yml             
├── Dockerfile                  
└── README.md                   
```

---

## Quick Start


#### 1. clone & install
```bash
git clone https://github.com/<your_user>/preclinical-radiomics-pipeline.git
cd preclinical-radiomics-pipeline
```
### 2. create & activate a clean virtual environment
```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows PowerShell:  .venv\Scripts\Activate
```
### 3. upgrade pip & install the package
```bash
python -m pip install --upgrade pip
pip install -e . 
```
### 4. slice an enhanced CT DICOM volume
```bash
python .\src\io\dicom_slices.py .\data\slice-volume\20241031094809_CT_ISRATV_0.dcm  CT_Volume.dcm  --orientations original ray
```
### 5. convert LET map to RT‑Dose DICOM
```bash
python .\src\io\dose_conversion.py .\data\dose-convert\template_dose.dcm  .\data\dose-convert\dose_npy.npy  .\data\dose-convert\let_npy.npy  .\data\dose-convert\out_let_map.dcm
```
### 6. extract radiomic features from two mice
```bash
python .\src\features\radiomics_pipeline.py --mice-dir data/MR --mice Mouse_01 --out-dir data/Radiomics_Features
```
### 7. segment dose and compute DVH stats
```bash
python .\src\segmentation\mr_dose.py --mr .\data\segment-dose\mr\mr_volume.nii --dose .\data\segment-dose\dose\dose_volume.nii --atlas .\data\segment-dose\atlas\registered_atlas.nii --hierarchy-csv .\data\segment-dose\atlas\registered_atlas_labels.csv
```
### 8. extract Hippocampus from DSURQE atlas (left+right merged)
```bash
python .\src\segmentation\atlas_regions_nifti.py --atlas dsurqe --label .\data\segment-dose\atlas\registered_atlas.nii --hierarchy-csv .\data\segment-dose\atlas\registered_atlas_labels.csv --region "Hippocampal region" --side both --merge-sides --smooth-radius 0 --out-dir .\data\segment-atlas
```
### 9. extract Cerebellum from Digimouse atlas
```bash
python .\src\segmentation\atlas_regions_nifti.py --atlas digimouse --label .\data\segment-dose\atlas\atlas_380x992x208.img --hierarchy-csv .\data\segment-dose\atlas\atlas_380x992x208.txt --region "cerebellum" -- --smooth-radius 0 --out-dir .\data\segment-atlas
```
### 10. Create PCA Analysis
```bash
python .\src\analysis\pca.py --csv data/Radiomics_Features/combined_features_all_mice.csv
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


