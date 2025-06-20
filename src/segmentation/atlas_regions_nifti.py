#!/usr/bin/env python3
"""
atlas_segment_export.py
───────────────────────
Export binary region masks from either **DSURQE** or **Digimouse** mouse brain
atlases.

Features
~~~~~~~~
* Choose atlas type with `--atlas {dsurqe|digimouse}`
* For **DSURQE**
  • needs a CSV that lists *Structure / hierarchy / left label / right label*
  • supports left / right export or merged sides
* For **Digimouse**
  • needs a plain-text table (`label_id[+label_id…] --> region_name`)
  • regions are not lateralised; the `--side` flag is ignored
* `--list-regions` ( `-l` ) prints all available regions and exits
* Masks can be smoothed with binary closing before saving

Quick examples
~~~~~~~~~~~~~~
```bash
# DSURQE Hippocampus, merged sides
python -m atlas_segment_export \
  --atlas dsurqe \
  --label Atlas_DSURQE.nii \
  --hierarchy-csv DSURQE_mapping.csv \
  --region Hippocampus \
  --merge-sides

# Digimouse Cerebellum (from supplied txt table)
python -m atlas_segment_export \
  --atlas digimouse \
  --label digimouse_atlas.nii \
  --digimouse-map atlas_380x992x208.txt \
  --region cerebellum

# List all Digimouse regions
python -m atlas_segment_export --atlas digimouse --digimouse-map atlas_380x992x208.txt -l
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd
import SimpleITK as sitk

# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────
def load_labelmap(path: str | Path) -> sitk.Image:
    """Read a label map as UInt16 NIfTI/Analyze image."""
    return sitk.ReadImage(str(path), sitk.sitkUInt16)


# ─────────────────────────────────────────────────────────────────────────────
# DSURQE helpers (same logic as previous version)
# ─────────────────────────────────────────────────────────────────────────────
def extract_bilateral_rois(
    labelmap: sitk.Image, csv_df: pd.DataFrame
) -> Dict[str, sitk.Image]:
    rois: Dict[str, sitk.Image] = {}
    for _, row in csv_df.iterrows():
        for side in ("left", "right"):
            lbl = row.get(f"{side} label")
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
    grouped: Dict[str, Dict[str, List[sitk.Image]]] = {}
    for hierarchy, rows in csv_df.groupby("hierarchy"):
        for _, row in rows.iterrows():
            for side in ("left", "right"):
                key = f"{row['Structure']} ({side})"
                if key in rois:
                    grouped.setdefault(hierarchy, {"left": [], "right": []})
                    grouped[hierarchy][side].append(rois[key])
    return grouped


# ─────────────────────────────────────────────────────────────────────────────
# Digimouse helpers
# ─────────────────────────────────────────────────────────────────────────────
def parse_digimouse_table(txt_path: Path) -> Dict[str, List[int]]:
    """Parse lines like `4+5+6 --> whole brain` into {region: [4,5,6]}"""
    mapping: Dict[str, List[int]] = {}
    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "-->" not in line:
                continue
            left, right = [part.strip() for part in line.split("-->")]
            ids = [int(t.strip()) for t in left.split("+") if t.strip().isdigit()]
            if ids:
                mapping[right.lower()] = ids  # store lowercase for ci lookup
    return mapping


def mask_from_label_ids(labelmap: sitk.Image, ids: List[int]) -> sitk.Image:
    """Return union‐mask of all *ids* in the labelmap."""
    base = sitk.Image(labelmap.GetSize(), sitk.sitkUInt8)
    base.CopyInformation(labelmap)
    for lbl in ids:
        part = sitk.BinaryThreshold(
            labelmap, lowerThreshold=lbl, upperThreshold=lbl, insideValue=1, outsideValue=0
        )
        base = sitk.Or(base, part)
    return base


# ─────────────────────────────────────────────────────────────────────────────
# Common helpers
# ─────────────────────────────────────────────────────────────────────────────
def join_masks(rois: List[sitk.Image]) -> sitk.Image | None:
    if not rois:
        return None
    base = sitk.Image(rois[0].GetSize(), rois[0].GetPixelID())
    base.CopyInformation(rois[0])
    for roi in rois:
        base = sitk.Or(base, roi)
    return base


def smooth_mask(mask: sitk.Image, radius: int) -> sitk.Image:
    if radius <= 0:
        return mask
    flt = sitk.BinaryMorphologicalClosingImageFilter()
    flt.SetKernelRadius(radius)
    return flt.Execute(mask)


# ─────────────────────────────────────────────────────────────────────────────
# Core exporter
# ─────────────────────────────────────────────────────────────────────────────
def export_dsurqe(
    labelmap: sitk.Image,
    mapping_csv: Path,
    region: str,
    side: str | None,
    merge_sides: bool,
    smooth_radius: int,
    out_dir: Path,
) -> None:
    mapping_df = pd.read_csv(mapping_csv)
    bilateral_rois = extract_bilateral_rois(labelmap, mapping_df)
    grouped = group_rois_by_hierarchy(mapping_df, bilateral_rois)

    match = next((h for h in grouped if h.lower() == region.lower()), None)
    if match is None:
        print(f"Error: hierarchy '{region}' not found (try --list-regions)", file=sys.stderr)
        sys.exit(1)

    sides = grouped[match]
    out_dir.mkdir(parents=True, exist_ok=True)

    def _write(mask: sitk.Image, name: str) -> None:
        if mask is None:
            print(f"Warning: no voxels for '{name}' – skipping.")
            return
        sitk.WriteImage(mask, out_dir / name)
        print("✓ written", out_dir / name)

    if side in ("left", "right"):
        mask = smooth_mask(join_masks(sides[side]), smooth_radius)
        _write(mask, f"{region}_{side}.nii")
    else:
        if merge_sides:
            mask = smooth_mask(join_masks(sides["left"] + sides["right"]), smooth_radius)
            _write(mask, f"{region}_merged.nii")
        else:
            for s in ("left", "right"):
                mask = smooth_mask(join_masks(sides[s]), smooth_radius)
                _write(mask, f"{region}_{s}.nii")


def export_digimouse(
    labelmap: sitk.Image,
    table_path: Path,
    region: str,
    smooth_radius: int,
    out_dir: Path,
) -> None:
    mapping = parse_digimouse_table(table_path)
    region_key = region.lower()
    if region_key not in mapping:
        print(f"Error: region '{region}' not found in Digimouse table (try --list-regions)", file=sys.stderr)
        sys.exit(1)

    mask = smooth_mask(mask_from_label_ids(labelmap, mapping[region_key]), smooth_radius)
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{region.replace(' ', '_')}.nii"
    sitk.WriteImage(mask, out_dir / fname)
    print("✓ written", out_dir / fname)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Export binary masks from DSURQE or Digimouse mouse brain atlases.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Atlas choice
    p.add_argument("--atlas", choices=["dsurqe", "digimouse"], default="dsurqe", help="Atlas type")

    # Common inputs
    p.add_argument("--label", required=True, help="Labelmap image (NIfTI or Analyze)")

    # DSURQE-specific
    p.add_argument("--hierarchy-csv", help="DSURQE mapping CSV (required if --atlas=dsurqe)")

    # Digimouse-specific
    p.add_argument("--digimouse-map", help="Digimouse label table TXT (required if --atlas=digimouse)")

    # Region selection / listing
    p.add_argument("--region", help="Region name / hierarchy to extract")
    p.add_argument("-l", "--list-regions", action="store_true", help="List available regions and exit")

    # Export options (DSURQE only)
    p.add_argument("--side", choices=["left", "right", "both"], default="both", help="(DSURQE) Which side(s) to export")
    p.add_argument("--merge-sides", action="store_true", help="Merge left+right (DSURQE, when --side=both)")

    # Common options
    p.add_argument("--smooth-radius", type=int, default=2, help="Binary closing radius (voxels)")
    p.add_argument("--out-dir", default="./segments", help="Output directory")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Entrypoint
# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:
    args = _parse_args()

    labelmap = load_labelmap(args.label)

    # Region listing shortcut
    if args.list_regions:
        if args.atlas == "dsurqe":
            if not args.hierarchy_csv:
                print("Error: --hierarchy-csv is required to list DSURQE regions", file=sys.stderr)
                sys.exit(1)
            mapping_df = pd.read_csv(args.hierarchy_csv)
            regions = sorted(set(mapping_df["hierarchy"].dropna()))
        else:  # digimouse
            if not args.digimouse_map:
                print("Error: --digimouse-map is required to list Digimouse regions", file=sys.stderr)
                sys.exit(1)
            regions = sorted(parse_digimouse_table(Path(args.digimouse_map)).keys())
        print("\nAvailable regions:\n" + "\n".join(regions))
        sys.exit(0)

    # Region name required for export
    if not args.region:
        print("Error: --region is required (unless --list-regions)", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.out_dir)

    if args.atlas == "dsurqe":
        if not args.hierarchy_csv:
            print("Error: --hierarchy-csv is required for DSURQE export", file=sys.stderr)
            sys.exit(1)
        export_dsurqe(
            labelmap=labelmap,
            mapping_csv=Path(args.hierarchy_csv),
            region=args.region,
            side=(None if args.side == "both" else args.side),
            merge_sides=args.merge_sides,
            smooth_radius=max(args.smooth_radius, 0),
            out_dir=out_dir,
        )
    else:  # digimouse
        if not args.digimouse_map:
            print("Error: --digimouse-map is required for Digimouse export", file=sys.stderr)
            sys.exit(1)
        export_digimouse(
            labelmap=labelmap,
            table_path=Path(args.digimouse_map),
            region=args.region,
            smooth_radius=max(args.smooth_radius, 0),
            out_dir=out_dir,
        )


if __name__ == "__main__":
    main()
