"""
Utility functions to split an Enhanced CT DICOM volume into a folder of
single‑slice CT images.

This version ports **all DICOM‑handling logic** from the PyInstaller GUI script
(`volume_to_slice_gui_rayCommandSupport.py`) into a lean, CLI‑only helper.  All
GUI, Nextcloud and threading code has been stripped, but every tag, UID and
naming rule needed for RayStation/LUT mapping is preserved.
"""
from __future__ import annotations

import logging
import os
import random
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import UID

__all__ = [
    "volume_to_slices",
    "create_dicom_slices_originalOrientation",
    "create_dicom_slices_rayCommandOrientation",
]

# ---------------------------------------------------------------------------
# Helper class – harvests and normalises metadata once per input volume
# ---------------------------------------------------------------------------


class metadataHelper:
    """Gather all information we later replicate to the 2‑D slices."""

    #: UID root strings we occasionally use when generating new identifiers
    STUDY_UID_ROOT = "1.2.826.0.1.3680043.10.543."
    DICOM_OBJECTS_UID = UID("1.2.826.0.1.3680043.1.2.100.6.40.0.76")

    def __init__(self, file_path: str | os.PathLike):
        self.filepath = str(file_path)
        self.originalDICOM: pydicom.FileDataset = pydicom.dcmread(self.filepath)
        self._extract_metadata()

    # ------------------------------------------------------------------
    # Metadata extraction (largely lifted 1‑to‑1 from the GUI code)
    # ------------------------------------------------------------------

    def _extract_metadata(self):  # noqa: C901  (yes, it is long – mirrors original logic)
        ds = self.originalDICOM  # shorthand

        # ---------- dates & times --------------------------------------------------
        self.StudyDate = ds.get("StudyDate", "")
        self.StudyTime = ds.get("StudyTime", "")
        self.SeriesDate = ds.get("SeriesDate", "")
        self.SeriesTime = ds.get("SeriesTime", "")
        # human‑readable YYYY‑MM‑DD used in naming later
        self.SeriesDateFormatted = (
            f"{self.SeriesDate[:4]}-{self.SeriesDate[4:6]}-{self.SeriesDate[6:]}"
            if self.SeriesDate
            else ""
        )

        # keep ContentDate/Time from the original file if present – fall back to *now*
        self.ContentDate = ds.get("ContentDate") or datetime.now().strftime("%Y%m%d")
        self.ContentTime = ds.get("ContentTime") or datetime.now().strftime("%H%M%S")

        # ---------- patient --------------------------------------------------------
        self.PatientName = ds.get("PatientName", "")
        self.PatientID = ds.get("PatientID", "")
        self.PatientBirthDate = ds.get("PatientBirthDate", "") or "20000101"
        self.PatientSex = "O"  # explicitly set to Other (RayStation requirement in original script)

        # ---------- manufacturer ---------------------------------------------------
        self.Manufacturer = ds.get("Manufacturer", "")
        self.ManufacturerModelName = ds.get("ManufacturerModelName", "")
        self.ReferringPhysicianName = ds.get("ReferringPhysicianName", "")

        # ---------- UIDs & SOP class ----------------------------------------------
        self.modality = ds.get("Modality", "")
        if self.modality == "CT":
            self.SOPClassUID = pydicom.uid.CTImageStorage
        else:
            # fall back to whatever the file declares (could be an MR‑Enhanced, etc.)
            self.SOPClassUID = ds.SOPClassUID
        self.StudyInstanceUID = ds.StudyInstanceUID

        # ---------- comments & descriptions ---------------------------------------
        self.ImageComments = ds.get("ImageComments", "")
        self.StudyDescription = ds.get("StudyDescription", "")
        self.StudyID = ds.get("StudyID", "")

        # ---------- acquisition parameters ----------------------------------------
        fg = ds.SharedFunctionalGroupsSequence[0]
        # (0018,0060)  KVP
        self.KVP = fg[(0x0018, 0x9325)][0][(0x0018, 0x0060)].value
        # (0018,9328)  exposure time [ms]
        self.ExposureTimeInms = fg[(0x0018, 0x9321)][0][(0x0018, 0x9328)].value
        # (0018,9330)  tube current [mA]
        self.XrayTubeCurrentMa = fg[(0x0018, 0x9321)][0][(0x0018, 0x9330)].value
        self.XRayTubeCurrentInuA = self.XrayTubeCurrentMa * 1000  # μA
        # (0018,9332)  mAs
        self.ExposureInmAs = fg[(0x0018, 0x9321)][0][(0x0018, 0x9332)].value

        # ---------- reconstruction parameters -------------------------------------
        parts = os.path.splitext(os.path.basename(self.filepath))[0].split("_")
        self.ReconstructionAlgorithm = parts[2] if len(parts) > 2 else ""
        self.ReconstructionNumber = parts[3] if len(parts) > 3 else ""

        # ---------- protocol name (RayStation LUT mapping) ------------------------
        in_plane_um = int(fg[(0x0028, 0x9110)][0][(0x0028, 0x0030)].value[0] * 1000)
        series_descr_clean = (
            str(ds.SeriesDescription).split("/")[1].strip().replace(" ", "").replace("-", "")
            if ds.SeriesDescription and "/" in ds.SeriesDescription
            else str(ds.SeriesDescription).replace(" ", "").replace("-", "")
        )
        self.ProtocolName = (
            f"{series_descr_clean}-{self.KVP}kVp"
            f"{int(round(self.XRayTubeCurrentInuA, 0))}uA-{in_plane_um}um-"
            f"{self.ReconstructionAlgorithm}"
        )

        # ---------- bit & pixel info ---------------------------------------------
        self.BitsAllocated = ds.BitsAllocated
        self.SamplesPerPixel = ds.SamplesPerPixel
        self.BitsStored = ds.BitsStored
        self.HighBit = ds.HighBit

        self.PixelSpacing = fg[(0x0028, 0x9110)][0][(0x0028, 0x0030)].value
        self.SliceThickness = fg[(0x0028, 0x9110)][0][(0x0018, 0x0050)].value

        # orientation & position
        self.ImageOrientationPatient = fg[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value
        self.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        self.PhotometricInterpretation = "MONOCHROME2"

        # ---------- rescale --------------------------------------------------------
        ip_fg = ds.SharedFunctionalGroupsSequence[0]
        self.RescaleIntercept_init = ip_fg[(0x0028, 0x9145)][0][(0x0028, 0x1052)].value
        self.RescaleSlope_init = ip_fg[(0x0028, 0x9145)][0][(0x0028, 0x1053)].value

        self.arr = ds.pixel_array  # shape [z, y, x]
        self.arr_rescaled = (
            np.round(
                self.arr.astype(np.float64) * float(self.RescaleSlope_init) + float(self.RescaleIntercept_init),
                0,
            ).astype(np.int16)
        )
        self.NumberOfFrames = self.arr_rescaled.shape[0]
        self.PixelRepresentation = 1  # signed after rescaling
        self.RescaleSlope = 1
        self.RescaleIntercept = 0

        # ---------- window preset --------------------------------------------------
        self.WindowCenter = -300
        self.WindowWidth = 1300

        # ---------- text for SeriesDescription ------------------------------------
        try:
            pat_part = str(ds.PatientName).split("/")[2].strip()
        except Exception:
            pat_part = str(ds.PatientName).strip()
        try:
            ser_part = str(ds.SeriesDescription).split("/")[1].strip()
        except Exception:
            ser_part = str(ds.SeriesDescription).strip()
        # Clean series descriptor without back‑slashes inside an f‑string expression
        ser_part_clean = re.sub(r"[ \-]", "", ser_part)
        self.SeriesDescription = f"{pat_part}_{ser_part_clean}"

        # ---------- StationName (critical for LUT mapping) ------------------------
        # The original GUI script copies ManufacturerModelName when StationName
        # is absent – we mimic that to avoid AttributeError.
        self.StationName = getattr(ds, "StationName", self.ManufacturerModelName)

        # ---------- misc UID roots for later use ----------------------------------
        self.SeriesInstanceUID = ds.SeriesInstanceUID
        self.FrameOfReferenceUID = ds.FrameOfReferenceUID

        self.SeriesNumber = str(random.randint(0, 9999))
        self.SOPSeriesRoot = f"1.2.276.0.7230010.3.{random.randint(1000, 9999)}"
        self.ImageSeriesUIDRoot = f"1.2.276.0.7230010.3.{random.randint(1000, 9999)}"

        # ---------- file‑meta copy -------------------------------------------------
        self.file_meta = FileMetaDataset()
        src_fm = ds.file_meta
        self.file_meta.FileMetaInformationVersion = src_fm.FileMetaInformationVersion
        self.file_meta.MediaStorageSOPClassUID = src_fm.MediaStorageSOPClassUID
        self.file_meta.FileMetaInformationGroupLength = src_fm.FileMetaInformationGroupLength
        self.file_meta.ImplementationVersionName = "DicomObjects.NET"
        self.file_meta.ImplementationClassUID = self.DICOM_OBJECTS_UID
        self.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian

        self.is_little_endian = True
        self.is_implicit_VR = True

    # ------------------------------------------------------------------
    # Convenience wrappers --------------------------------------------------------

    def new_filename(self, out_dir: str | os.PathLike, idx: int) -> str:
        return os.path.join(out_dir, f"CT{idx}.dcm")


# ---------------------------------------------------------------------------
# Study‑level metadata helper (unchanged except for new tags)
# ---------------------------------------------------------------------------

def set_metadata(meta: metadataHelper, filename: str) -> FileDataset:
    """Populate tags that stay constant across all slices."""

    ds = FileDataset(filename, {}, file_meta=meta.file_meta, preamble=b"\0" * 128)
    ds.file_meta = meta.file_meta

    ds.is_little_endian = meta.is_little_endian
    ds.is_implicit_VR = meta.is_implicit_VR

    # bits / samples
    ds.BitsAllocated = meta.BitsAllocated
    ds.SamplesPerPixel = meta.SamplesPerPixel
    ds.PixelRepresentation = meta.PixelRepresentation
    ds.BitsStored = meta.BitsStored
    ds.HighBit = meta.HighBit

    # IDs & descriptions
    ds.SOPClassUID = meta.SOPClassUID
    ds.Modality = meta.modality
    ds.StudyInstanceUID = meta.StudyInstanceUID
    ds.StudyDescription = meta.StudyDescription
    ds.StudyID = meta.StudyID

    # patient
    ds.PatientName = meta.PatientName
    ds.PatientID = meta.PatientID
    ds.PatientBirthDate = meta.PatientBirthDate
    ds.PatientSex = meta.PatientSex
    ds.ImageComments = meta.ImageComments

    # dates / times
    ds.StudyDate = meta.StudyDate
    ds.StudyTime = meta.StudyTime
    ds.SeriesDate = meta.SeriesDate
    ds.SeriesTime = meta.SeriesTime
    ds.ContentDate = meta.ContentDate
    ds.ContentTime = meta.ContentTime

    # manufacturer & acquisition
    ds.Manufacturer = meta.Manufacturer
    ds.ManufacturerModelName = meta.ManufacturerModelName
    ds.ReferringPhysicianName = meta.ReferringPhysicianName

    ds.KVP = meta.KVP
    ds.ExposureTimeInms = meta.ExposureTimeInms
    ds.ExposureInmAs = meta.ExposureInmAs
    ds.XRayTubeCurrentInuA = meta.XRayTubeCurrentInuA
    ds.ReconstructionAlgorithm = meta.ReconstructionAlgorithm
    ds.ProtocolName = meta.ProtocolName
    ds.StationName = meta.StationName

    # photometric & rescale
    ds.PhotometricInterpretation = meta.PhotometricInterpretation
    ds.RescaleIntercept = meta.RescaleIntercept
    ds.RescaleSlope = meta.RescaleSlope

    # windowing
    ds.WindowCenter = meta.WindowCenter
    ds.WindowWidth = meta.WindowWidth

    return ds


# ---------------------------------------------------------------------------
# Slice creation helpers
# ---------------------------------------------------------------------------

def _write_slice(
    ds: FileDataset,
    CTslice: np.ndarray,
    meta: metadataHelper,
    idx: int,
    orientation: str,
    series_uid: str,
    frame_uid: str,
    series_number: str,
    image_orientation_patient,
    patient_position,
    filename: str,
):
    """Populate series‑ & instance‑level tags and write to *filename*."""

    ds.SeriesInstanceUID = series_uid
    ds.FrameOfReferenceUID = frame_uid
    ds.SeriesNumber = series_number
    ds.SeriesDescription = (
        meta.SeriesDescription if orientation == "original" else meta.SeriesDescription + "_rayOrientation"
    )

    ds.Rows, ds.Columns = CTslice.shape
    ds.PixelSpacing = meta.PixelSpacing
    ds.SliceThickness = meta.SliceThickness
    ds.NumberOfFrames = 1

    ds.PatientPosition = patient_position
    ds.ImageOrientationPatient = image_orientation_patient

    # instance‑specific
    sop_instance_uid = pydicom.uid.generate_uid()
    ds.file_meta.MediaStorageSOPInstanceUID = sop_instance_uid
    ds.SOPInstanceUID = sop_instance_uid
    ds.InstanceNumber = str(idx)
    ds.ImagePositionPatient = [
        meta.originalDICOM.PerFrameFunctionalGroupsSequence[0][(0x0020, 0x9113)][0][(0x0020, 0x0032)].value[0],
        meta.originalDICOM.PerFrameFunctionalGroupsSequence[0][(0x0020, 0x9113)][0][(0x0020, 0x0032)].value[1],
        meta.originalDICOM.PerFrameFunctionalGroupsSequence[0][(0x0020, 0x9113)][0][(0x0020, 0x0032)].value[2]
        - (idx - 1) * meta.SliceThickness,
    ]

    ds.PixelData = CTslice.tobytes()
    ds.save_as(filename)


# ---------------------- orientation helpers ---------------------------------

def create_dicom_slices_originalOrientation(meta: metadataHelper, cb=None):
    out_dir = os.path.join(
        os.path.dirname(meta.filepath),
        os.path.splitext(os.path.basename(meta.filepath))[0] + "_originalOrientation_slices",
    )
    os.makedirs(out_dir, exist_ok=True)

    series_uid = pydicom.uid.generate_uid()
    frame_uid = pydicom.uid.generate_uid()
    series_number = str(random.randint(0, 9999))

    for idx, slice_arr in enumerate(meta.arr_rescaled, start=1):
        filename = meta.new_filename(out_dir, idx)
        ds = set_metadata(meta, filename)
        _write_slice(
            ds,
            slice_arr,
            meta,
            idx,
            "original",
            series_uid,
            frame_uid,
            series_number,
            meta.ImageOrientationPatient,
            "HFP",
            filename,
        )
        if cb:
            cb(idx / meta.NumberOfFrames * 100)


def create_dicom_slices_rayCommandOrientation(meta: metadataHelper, cb=None):
    """RayStation‑compatible orientation (HFS, [1,0,0,0,1,0])."""
    out_dir = os.path.join(
        os.path.dirname(meta.filepath),
        os.path.splitext(os.path.basename(meta.filepath))[0] + "_rayOrientation_slices",
    )
    os.makedirs(out_dir, exist_ok=True)

    series_uid = pydicom.uid.generate_uid()
    frame_uid = pydicom.uid.generate_uid()
    series_number = str(random.randint(0, 9999))

    # flip AP axis + rotate 180° around IS
    arr_hfs = np.rot90(np.flip(meta.arr_rescaled, axis=1), k=2, axes=(1, 2))

    for idx, slice_arr in enumerate(arr_hfs, start=1):
        filename = meta.new_filename(out_dir, idx)
        ds = set_metadata(meta, filename)
        _write_slice(
            ds,
            slice_arr,
            meta,
            idx,
            "ray",
            series_uid,
            frame_uid,
            series_number,
            [1, 0, 0, 0, 1, 0],
            "HFS",
            filename,
        )
        if cb:
            cb(idx / meta.NumberOfFrames * 100)


# ---------------------------------------------------------------------------
# Public convenience wrapper
# ---------------------------------------------------------------------------

def volume_to_slices(
    volume_path: str | Path,
    orientations: tuple[str, ...] = ("original", "ray"),
    progress: bool = True,
):
    volume_path = Path(volume_path)
    meta = metadataHelper(volume_path)

    def _print(pct: float):
        if progress:
            print(f"\rProgress: {pct:6.2f}%", end="", flush=True)

    if "original" in orientations:
        create_dicom_slices_originalOrientation(meta, _print if progress else None)
    if "ray" in orientations:
        create_dicom_slices_rayCommandOrientation(meta, _print if progress else None)

    if progress:
        print("\nDone.")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Split an Enhanced CT DICOM volume into single‑slice DICOM images.",
    )
    parser.add_argument("dicom", help="Path to the input Enhanced CT DICOM file")
    parser.add_argument(
        "--orientations",
        nargs="+",
        default=["original", "ray"],
        choices=["original", "ray"],
        help="Which orientations to export (default: both)",
    )
    parser.add_argument(
        "--no‑progress",
        dest="no_progress",
        action="store_true",
        help="Disable textual progress output",
    )

    args = parser.parse_args()
    volume_to_slices(args.dicom, tuple(args.orientations), progress=not args.no_progress)
