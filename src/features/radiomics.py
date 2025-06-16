#!/usr/bin/env python3
"""
radiomics.py
────────────
End-to-end feature-extraction pipeline:

1. Resample      – SimpleITK BSpline to 0.2 mm³
2. Normalise     – Z-score intensity normalisation
3. Filter        – clamp outliers ±3 σ
4. Discretise    – 32-level fixed-bin width
5. Pyradiomics   – extract all default features per ROI
6. CSV output    – one file per ROI + combined table

The code is almost line-for-line from the original
*radiomics_feature_extraction_with_image_normalization.py*
so existing results remain reproducible.

CLI example
------------
python -m rad_pipeline.features.radiomics \
    --mice-dir data/Manually_preprocessed_and_verified \
    --mice Mouse27_Verified Mouse30_Verified \
    --out-dir data/radiomics_output
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor


# -------------------------------------------------------------------------
# Augmentation helpers
# -------------------------------------------------------------------------

# def elastic_deformation(volume, alpha, sigma, random_seed=42):
#     np.random.seed(random_seed)  # Set seed for reproducibility
#     shape = volume.shape
#     dx = ndi.gaussian_filter((np.random.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha
#     dy = ndi.gaussian_filter((np.random.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha
#     dz = ndi.gaussian_filter((np.random.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha
#
#     z, y, x = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing='ij')
#     indices = np.reshape(z + dz, (-1, 1)), np.reshape(y + dy, (-1, 1)), np.reshape(x + dx, (-1, 1))
#
#     deformed_volume = ndi.map_coordinates(volume, indices, order=1, mode='reflect').reshape(shape)
#     return deformed_volume
#
#
# def adjust_brightness(volume, factor):
#     normalized_volume = volume / np.max(volume)  # Normalize to avoid exceeding intensity range
#     adjusted_volume = normalized_volume * factor
#     return np.clip(adjusted_volume, 0, 1)  # Ensure values remain in valid range
#
#
# def add_gaussian_noise(volume, mean=0, std_factor=0.05, random_seed=42):
#     np.random.seed(random_seed)  # Set seed for reproducibility
#     std = std_factor * (np.max(volume) - np.min(volume))
#     noise = np.random.normal(mean, std, volume.shape)
#     noisy_volume = volume + noise
#     return np.clip(noisy_volume, 0, 1)  # Clip values to valid range


# -------------------------------------------------------------------------
# Pre-processing helpers
# -------------------------------------------------------------------------


def resample_image(image: sitk.Image, new_spacing: List[float] = [0.2, 0.2, 0.2]) -> sitk.Image:
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    new_size = [
        int(round(original_size[i] * (original_spacing[i] / new_spacing[i])))
        for i in range(3)
    ]

    resampler = sitk.ResampleImageFilter()
    resampler.SetInterpolator(sitk.sitkBSpline)
    resampler.SetOutputSpacing(new_spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputDirection(image.GetDirection())
    resampler.SetOutputOrigin(image.GetOrigin())
    return resampler.Execute(image)


def normalize_intensity(image: sitk.Image) -> sitk.Image:
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    return sitk.ShiftScale(image, shift=-stats.GetMean(), scale=1 / stats.GetSigma())


def filter_outliers(image: sitk.Image, num_stddev: int = 3) -> sitk.Image:
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    lower = stats.GetMean() - num_stddev * stats.GetSigma()
    upper = stats.GetMean() + num_stddev * stats.GetSigma()
    return sitk.Clamp(image, sitk.sitkFloat32, lower, upper)


def discretize_image(image: sitk.Image, num_bins: int = 32) -> sitk.Image:
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    bin_width = (stats.GetMaximum() - stats.GetMinimum()) / num_bins
    discretised = sitk.IntensityWindowing(
        image,
        windowMinimum=stats.GetMinimum(),
        windowMaximum=stats.GetMaximum(),
        outputMinimum=0,
        outputMaximum=num_bins - 1,
    )
    return sitk.Cast(discretised, sitk.sitkUInt8)


def preprocess_image(image: sitk.Image) -> sitk.Image:
    return discretize_image(
        filter_outliers(
            normalize_intensity(
                resample_image(image)
            )
        )
    )


# -------------------------------------------------------------------------
# I/O helpers (unchanged)
# -------------------------------------------------------------------------


def load_nifti_volume(path: str | Path) -> sitk.Image:
    return preprocess_image(sitk.ReadImage(str(path), sitk.sitkFloat32))


def load_labelmap(path: str | Path) -> sitk.Image:
    return sitk.ReadImage(str(path), sitk.sitkFloat32)


def resample_labelmap_to_image(label: sitk.Image, image: sitk.Image) -> sitk.Image:
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(image)
    resample.SetInterpolator(sitk.sitkBSpline)
    resample.SetOutputDirection(image.GetDirection())
    resample.SetOutputSpacing(image.GetSpacing())
    resample.SetOutputOrigin(image.GetOrigin())
    resample.SetSize(image.GetSize())
    return resample.Execute(label)


# -------------------------------------------------------------------------
# Feature normalisation & CSV helpers (unchanged)
# -------------------------------------------------------------------------


def normalize_features_z_score(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    df = df.dropna(subset=["Value"])
    df.loc[:, "Value"] = (df["Value"] - df["Value"].mean()) / df["Value"].std()
    return df


def save_features_to_csv(df: pd.DataFrame, out_dir: Path, image_name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{image_name}_features.csv", index=False)


# -------------------------------------------------------------------------
# Core extraction workflow (unchanged)
# -------------------------------------------------------------------------


def extract_features_from_nifti(
    nifti_file: Path,
    label: sitk.Image,
    out_dir: Path,
    image_type: str,
    image_name: str,
    combined_df: pd.DataFrame,
    extractor: featureextractor.RadiomicsFeatureExtractor,
) -> pd.DataFrame:
    image = load_nifti_volume(nifti_file)
    resampled_label = resample_labelmap_to_image(label, image)
    result = extractor.execute(image, resampled_label)

    df = pd.DataFrame(
        [
            (k, np.ndarray.item(v) if isinstance(v, np.ndarray) else v)
            for k, v in result.items()
        ],
        columns=["Feature", "Value"],
    )
    df["Image"] = image_type
    df["FileName"] = image_name
    df = normalize_features_z_score(df)

    combined_df = pd.concat([combined_df, df], ignore_index=True)
    save_features_to_csv(df, out_dir, image_name)
    print(f"✓ features saved for {image_name}")
    return combined_df


def process_mouse(
    mouse_dir: Path,
    combined_df: pd.DataFrame,
    out_dir: Path,
    extractor: featureextractor.RadiomicsFeatureExtractor,
) -> pd.DataFrame:
    data_file = mouse_dir / "MR_BSpline_Registered_to_CT.nii"
    rois_dir = mouse_dir / "rois"

    regions = {
        "Hippocampal Region Left": "Hippocampal_Region_Left.nii",
        "Hippocampal Region Right": "Hippocampal_Region_Right.nii",
    }

    for region_name, labelmap_file in regions.items():
        label_path = rois_dir / labelmap_file
        if label_path.exists():
            labelmap = load_labelmap(label_path)
            combined_df = extract_features_from_nifti(
                data_file,
                labelmap,
                out_dir / region_name.replace(" ", "_"),
                region_name,
                mouse_dir.name,
                combined_df,
                extractor,
            )
        else:
            print(f"⚠ {labelmap_file} not found in {rois_dir}; skipping.")
    return combined_df


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run radiomics feature extraction for a list of mice."
    )
    p.add_argument(
        "--mice-dir",
        required=True,
        help="Root directory that contains one sub-folder per mouse.",
    )
    p.add_argument(
        "--mice",
        nargs="+",
        required=True,
        help="Space-separated list of mouse sub-folder names to process.",
    )
    p.add_argument(
        "--out-dir",
        required=True,
        help="Where to write per-ROI CSVs and the combined CSV.",
    )
    p.add_argument(
        "--bin-width",
        type=int,
        default=32,
        help="Number of intensity bins (default: 32).",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    # configure Pyradiomics extractor – only geometry tolerance kept from original
    settings = {"geometryTolerance": 2e-1}
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)

    mice_dir = Path(args.mice_dir)
    out_dir = Path(args.out_dir)
    mice_list = args.mice

    combined_df = pd.DataFrame()

    for mouse in mice_list:
        mouse_path = mice_dir / mouse
        if mouse_path.is_dir():
            combined_df = process_mouse(mouse_path, combined_df, out_dir, extractor)
        else:
            print(f"⚠ mouse folder {mouse_path} missing; skipping.")

    combined_df.to_csv(out_dir / "combined_features_all_mice.csv", index=False)
    print("✓ feature extraction finished for all mice.")


if __name__ == "__main__":
    main()
