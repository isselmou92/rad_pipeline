#!/usr/bin/env python3
"""
mr_dose.py
──────────
• Loads MRI, dose, and atlas‐label NIfTI volumes
• Builds left / right ROIs from a CSV mapping (label numbers ↔ structures)
• Combines ROIs into hierarchy-level masks, smooths them, and aligns them to
  the MRI or dose grid on demand
• Segments both MR and dose volumes, then computes dose statistics
  (Dmin / Dmax / Dmean / D2 / D10 / D50 / D90 / D95 / D98, HI, etc.)
• Saves one CSV with all hierarchy-side statistics and (optionally) exports the
  segmented NIfTI volumes.

Identical maths, identical statistics — just parameterised paths and a small
`argparse` CLI wrapper for reproducibility.

Example
-------
python -m rad_pipeline.segmentation.mr_dose \
    --mr Mouse33_MR_BSpline_to_CT.nii \
    --label Mouse33_Atlas_Registered_to_MR.nii \
    --dose Mouse33_Dose_Mousehead.nii \
    --hierarchy-csv DSURQE_mapping.csv \
    --out dose_stats.csv \
    --export-nifti
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import SimpleITK as sitk


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers (unchanged from the original script)
# ─────────────────────────────────────────────────────────────────────────────
def load_nifti_volume(path: str | Path) -> sitk.Image:
    """Read a floating-point NIfTI volume."""
    return sitk.ReadImage(str(path), sitk.sitkFloat32)


def load_labelmap(path: str | Path) -> sitk.Image:
    """Read the atlas labelmap (integer image)."""
    return sitk.ReadImage(str(path), sitk.sitkUInt16)


# ─────────────────────────────────────────────────────────────────────────────
# Core processing helpers (verbatim logic)
# ─────────────────────────────────────────────────────────────────────────────
def segment_volume(volume: sitk.Image, roi_mask: sitk.Image) -> sitk.Image:
    """Mask *volume* by *roi_mask* (labels > 0 retained)."""
    return volume * sitk.Cast(roi_mask, volume.GetPixelID())


def compute_dose_statistics(segmented_dose: sitk.Image) -> Dict[str, float] | None:
    """Return common DVH-style statistics for a masked dose image."""
    arr = sitk.GetArrayFromImage(segmented_dose)
    arr = arr[arr > 0]  # drop background

    if arr.size == 0:
        return None

    return {
        "voxels": int(arr.size),
        "dmin": float(arr.min()),
        "dmean": float(arr.mean()),
        "dmax": float(arr.max()),
        "d2_gy": float(np.percentile(arr, 98)),
        "d10_gy": float(np.percentile(arr, 90)),
        "d50_gy": float(np.percentile(arr, 50)),
        "d90_gy": float(np.percentile(arr, 10)),
        "d95_gy": float(np.percentile(arr, 5)),
        "d98_gy": float(np.percentile(arr, 2)),
        "HI": float(
            (np.percentile(arr, 90) - np.percentile(arr, 10))
            / np.percentile(arr, 50)
        ),
    }


def join_rois_to_single_image(rois: List[sitk.Image]) -> sitk.Image | None:
    """Logical OR of a list of ROI masks, preserving geometry."""
    if not rois:
        return None
    base = sitk.Image(rois[0].GetSize(), rois[0].GetPixelID())
    base.CopyInformation(rois[0])
    for roi in rois:
        base = sitk.Or(base, roi)
    return base


def smooth_segmentation_closing(seg: sitk.Image, radius: int = 2) -> sitk.Image:
    """Binary closing (fill small holes) on a mask."""
    filt = sitk.BinaryMorphologicalClosingImageFilter()
    filt.SetKernelRadius(radius)
    return filt.Execute(seg)


def match_physical_space(
    reference: sitk.Image, moving: sitk.Image
) -> sitk.Image:
    """Resample *moving* into *reference*’s grid using nearest-neighbour."""
    res = sitk.ResampleImageFilter()
    res.SetReferenceImage(reference)
    res.SetInterpolator(sitk.sitkNearestNeighbor)
    res.SetTransform(sitk.Transform())
    return res.Execute(moving)


# ─────────────────────────────────────────────────────────────────────────────
# ROI extraction
# ─────────────────────────────────────────────────────────────────────────────
def extract_bilateral_rois(
    labelmap: sitk.Image, csv_df: pd.DataFrame
) -> Dict[str, sitk.Image]:
    """
    Build one binary mask per {Structure}-(left/right) using the CSV that
    contains columns:  Structure | hierarchy | left label | right label
    """
    rois: Dict[str, sitk.Image] = {}
    for _, row in csv_df.iterrows():
        for side in ("left", "right"):
            lbl = row[f"{side} label"]
            if pd.notna(lbl):
                mask = sitk.BinaryThreshold(
                    labelmap,
                    lowerThreshold=int(lbl),
                    upperThreshold=int(lbl),
                    insideValue=1,
                    outsideValue=0,
                )
                rois[f"{row['Structure']} ({side})"] = mask
    return rois


def group_rois_by_hierarchy(
    csv_df: pd.DataFrame, rois: Dict[str, sitk.Image]
) -> Dict[str, Dict[str, List[sitk.Image]]]:
    """
    Return a nested dict:  {hierarchy: {'left': [roi1, …], 'right': [roi1, …]}}
    """
    grouped: Dict[str, Dict[str, List[sitk.Image]]] = {}
    for hierarchy, rows in csv_df.groupby("hierarchy"):
        for _, row in rows.iterrows():
            for side in ("left", "right"):
                roi_key = f"{row['Structure']} ({side})"
                if roi_key in rois:
                    grouped.setdefault(hierarchy, {"left": [], "right": []})
                    grouped[hierarchy][side].append(rois[roi_key])
    return grouped


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────
def process(
    mr_path: Path,
    dose_path: Path,
    label_path: Path,
    csv_path: Path,
    export_nifti: bool,
    out_dir: Path,
) -> pd.DataFrame:
    """Run the full segmentation/statistics workflow."""
    out_dir.mkdir(parents=True, exist_ok=True)

    print("• loading volumes …")
    mr = load_nifti_volume(mr_path)
    dose = load_nifti_volume(dose_path)
    atlas = load_labelmap(label_path)
    mapping = pd.read_csv(csv_path)

    print("• building ROIs …")
    bilateral_rois = extract_bilateral_rois(atlas, mapping)
    grouped = group_rois_by_hierarchy(mapping, bilateral_rois)

    rows: List[Tuple[str, str, Dict[str, float]]] = []

    print("• processing hierarchies …")
    for hierarchy, sides in grouped.items():
        print(f"  – {hierarchy}")
        # merge, smooth, align
        merged = {
            side: smooth_segmentation_closing(join_rois_to_single_image(rois))
            for side, rois in sides.items()
        }

        for side, roi in merged.items():
            if roi is None:
                continue

            aligned_dose_roi = match_physical_space(dose, roi)
            aligned_mr_roi = match_physical_space(mr, roi)

            seg_dose = segment_volume(dose, aligned_dose_roi)
            stats = compute_dose_statistics(seg_dose)
            if stats:
                rows.append((hierarchy, side, stats))

            if export_nifti:
                seg_dir = out_dir / "segments"
                seg_dir.mkdir(exist_ok=True)
                sitk.WriteImage(seg_dose, seg_dir / f"dose_{hierarchy}_{side}.nii")
                mr_seg = segment_volume(mr, aligned_mr_roi)
                sitk.WriteImage(mr_seg, seg_dir / f"mr_{hierarchy}_{side}.nii")

    # flatten into DataFrame
    df = (
        pd.DataFrame(rows, columns=["hierarchy", "side", "stats"])
        .join(pd.json_normalize(rows, record_path="stats", errors="ignore"))
        .drop(columns="stats")
    )
    df.to_csv(out_dir / "dose_stats.csv", index=False)
    print("✓ written", out_dir / "dose_stats.csv")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Segment dose & MR volumes by atlas hierarchy and compute dose statistics."
    )
    p.add_argument("--mr", required=True, help="MRI NIfTI volume")
    p.add_argument("--dose", required=True, help="Dose NIfTI volume")
    p.add_argument("--label", required=True, help="Atlas labelmap NIfTI")
    p.add_argument(
        "--hierarchy-csv",
        required=True,
        help="CSV with columns Structure / hierarchy / left label / right label",
    )
    p.add_argument(
        "--out",
        default="dose_stats.csv",
        help="Output CSV (will also create segments/ folder if --export-nifti)",
    )
    p.add_argument(
        "--export-nifti",
        action="store_true",
        help="Write segmented MR & dose volumes to NIfTI files.",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out_path = Path(args.out).resolve()
    process(
        mr_path=Path(args.mr),
        dose_path=Path(args.dose),
        label_path=Path(args.label),
        csv_path=Path(args.hierarchy_csv),
        export_nifti=args.export_nifti,
        out_dir=out_path.parent,
    )


if __name__ == "__main__":
    main()
