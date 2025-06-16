"""
dicom_slices.py
---------------
Utility functions to split an Enhanced CT DICOM volume into a folder of
single-slice CT images.  All graphical-user-interface (tkinter) and
Nextcloud-upload functionality present in the original script have been removed,
keeping the module focused purely on DICOM I/O.  Core metadata handling and
slice-export logic are preserved verbatim so the output remains identical to
the original workflow.
"""

from __future__ import annotations

import os
import random
from datetime import datetime
from pathlib import Path
import logging

import numpy as np
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import UID

# ---------------------------------------------------------------------------
# Helper class
# ---------------------------------------------------------------------------


class metadataHelper:
    def __init__(self, file_path):
        self.originalDICOM = pydicom.dcmread(file_path)
        self.filepath = file_path
        self.getAndAdjustMetadata()

    def getAndAdjustMetadata(self):

        dicom_volume = self.originalDICOM

        # Series and Study  (description and IDs set end of function)
        self.StudyDate = dicom_volume.StudyDate
        self.StudyTime = dicom_volume.StudyTime
        self.SeriesDate = dicom_volume.SeriesDate
        self.SeriesDateFormatted = f"{self.SeriesDate[:4]}-{self.SeriesDate[4:6]}-{self.SeriesDate[6:]}"
        self.SeriesTime = dicom_volume.SeriesTime

        # ContentDate & –Time (= modification date/time)
        self.ContentDate = datetime.now().strftime("%Y%m%d")
        self.ContentTime = datetime.now().strftime("%H%M%S")

        # Patient
        self.PatientName = dicom_volume.PatientName
        self.PatientID = dicom_volume.PatientID
        self.PatientBirthDate = dicom_volume.PatientBirthDate
        self.PatientSex = "O"

        # Manufacturer details
        self.Manufacturer = dicom_volume.Manufacturer
        self.ManufacturerModelName = dicom_volume.ManufacturerModelName
        self.ReferringPhysicianName = dicom_volume.ReferringPhysicianName

        # SOP class, Study, Series, Frame UIDs
        self.SOPClassUID = dicom_volume.SOPClassUID
        self.StudyInstanceUID = dicom_volume.StudyInstanceUID

        # Comments
        self.ImageComments = dicom_volume.ImageComments

        # Study description / ID
        self.StudyDescription = dicom_volume.StudyDescription
        self.StudyID = dicom_volume.StudyID

        # File-meta information (copied 1-to-1)
        self.file_meta = FileMetaDataset()
        self.file_meta.FileMetaInformationVersion = dicom_volume.file_meta.FileMetaInformationVersion
        self.file_meta.MediaStorageSOPClassUID = dicom_volume.file_meta.MediaStorageSOPClassUID
        self.file_meta.FileMetaInformationGroupLength = dicom_volume.file_meta.FileMetaInformationGroupLength
        self.file_meta.FileMetaInformationVersion = dicom_volume.file_meta.FileMetaInformationVersion
        self.file_meta.ImplementationVersionName = "DicomObjects.NET"
        self.file_meta.ImplementationClassUID = UID("1.2.826.0.1.3680043.1.2.100.6.40.0.76")
        self.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian

        self.is_little_endian = True
        self.is_implicit_VR = True

        # X-ray acquisition parameters
        self.KVP = dicom_volume.SharedFunctionalGroupsSequence[0][0x0018, 0x9325][0][0x0018, 0x0060].value
        self.ExposureTimeInms = dicom_volume.SharedFunctionalGroupsSequence[0][0x0018, 0x9321][0][0x0018, 0x9328].value
        self.XrayTubeCurrentMa = dicom_volume.SharedFunctionalGroupsSequence[0][0x0018, 0x9321][0][
            0x0018, 0x9330
        ].value
        self.XRayTubeCurrentInuA = self.XrayTubeCurrentMa * 1000
        self.ExposureInmAs = dicom_volume.SharedFunctionalGroupsSequence[0][0x0018, 0x9321][0][0x0018, 0x9332].value

        # Reconstruction parameters
        self.ReconstructionAlgorithm = os.path.splitext(os.path.split(self.filepath)[1])[0].split("_")[2]
        self.ReconstructionNumber = os.path.splitext(os.path.split(self.filepath)[1])[0].split("_")[3]

        # ProtocolName (important for CT-LUT mapping in RayStation)
        self.ProtocolName = (
            f"{str(dicom_volume.SeriesDescription).split('/')[1].strip()}_"
            f"{int(dicom_volume.SharedFunctionalGroupsSequence[0][0x0028,0x9110][0][0x0028,0x0030].value[0]*1000)}um"
            f"-{self.ReconstructionAlgorithm}"
        )

        # Bit information
        self.BitsAllocated = dicom_volume.BitsAllocated
        self.SamplesPerPixel = dicom_volume.SamplesPerPixel
        self.BitsStored = dicom_volume.BitsStored
        self.HighBit = dicom_volume.HighBit

        # Pixel dimensions
        self.PixelSpacing = dicom_volume.SharedFunctionalGroupsSequence[0][0x0028, 0x9110][0][0x0028, 0x0030].value
        self.SliceThickness = dicom_volume.SharedFunctionalGroupsSequence[0][0x0028, 0x9110][0][0x0018, 0x0050].value

        # Image orientation and patient position
        self.ImageOrientationPatient = dicom_volume.SharedFunctionalGroupsSequence[0][0x0020, 0x9116][0][
            0x0020, 0x0037
        ].value

        self.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]  # RayStation requirement
        self.PhotometricInterpretation = "MONOCHROME2"

        # Rescale pixel data
        self.RescaleIntercept_init = dicom_volume.SharedFunctionalGroupsSequence[0][0x0028, 0x9145][0][
            0x0028, 0x1052
        ].value
        self.RescaleSlope_init = dicom_volume.SharedFunctionalGroupsSequence[0][0x0028, 0x9145][0][
            0x0028, 0x1053
        ].value

        self.arr = dicom_volume.pixel_array  # [z=IS, y=AP, x=RL]
        self.arr_rescaled = np.round(
            self.arr.astype(np.float64) * float(self.RescaleSlope_init) + float(self.RescaleIntercept_init), 0
        ).astype(np.int16)
        self.NumberOfFrames = self.arr_rescaled.shape[0]

        self.PixelRepresentation = 1  # signed int now (initially unsigned)
        self.RescaleSlope = 1
        self.RescaleIntercept = 0

        # Window preset
        self.WindowCenter = -300
        self.WindowWidth = 1300

        # Series description / ID
        self.SeriesDescription = (
            str(dicom_volume.PatientName).split("/")[2].strip()
            + "_"
            + str(dicom_volume.SeriesDescription).split("/")[1].strip().replace(" ", "").replace("-", "")
        )

        # Station name (CT-LUT mapping in RayStation)
        self.StationName = dicom_volume.StationName

        self.SeriesInstanceUID = dicom_volume.SeriesInstanceUID
        self.FrameOfReferenceUID = dicom_volume.FrameOfReferenceUID

        # Study UID and Root
        self.StudyUIDRoot = "1.2.826.0.1.3680043.10.543."
        self.StudyRoot = dicom_volume.StudyInstanceUID.split(".")[-1]

        # Conversion for DICOMDIR creation (SeriesNumber etc.)
        self.SeriesNumber = str(random.randint(0, 10000))

        # SOP series root (for individual slices)
        self.SOPSeriesRoot = f"1.2.276.0.7230010.3.{random.randint(1000,9999)}"

        # UID root used for instances
        self.ImageSeriesUIDRoot = f"1.2.276.0.7230010.3.{random.randint(1000,9999)}"

        # SOP instances will be created later for every slice


# ---------------------------------------------------------------------------
# Study-level metadata helper
# ---------------------------------------------------------------------------


def set_metadata(metadataHelper, filename):
    """
    Set tags that are constant across the study.
    Series- and instance-specific tags are added outside this helper.
    """
    ds = FileDataset(filename, {}, file_meta=metadataHelper.file_meta, preamble=b"\0" * 128)
    ds.file_meta = metadataHelper.file_meta

    ds.is_little_endian = metadataHelper.is_little_endian
    ds.is_implicit_VR = metadataHelper.is_implicit_VR

    # Bits / samples
    ds.BitsAllocated = metadataHelper.BitsAllocated
    ds.SamplesPerPixel = metadataHelper.SamplesPerPixel
    ds.PixelRepresentation = metadataHelper.PixelRepresentation
    ds.BitsStored = metadataHelper.BitsStored
    ds.HighBit = metadataHelper.HighBit

    # Identifiers
    ds.SOPClassUID = metadataHelper.SOPClassUID
    ds.StudyInstanceUID = metadataHelper.StudyInstanceUID
    ds.StudyDescription = metadataHelper.StudyDescription
    ds.StudyID = metadataHelper.StudyID

    # Patient
    ds.ImageComments = metadataHelper.ImageComments
    ds.PatientName = metadataHelper.PatientName
    ds.PatientID = metadataHelper.PatientID
    ds.PatientBirthDate = metadataHelper.PatientBirthDate
    ds.PatientSex = metadataHelper.PatientSex

    # Dates / times
    ds.StudyDate = metadataHelper.StudyDate
    ds.StudyTime = metadataHelper.StudyTime
    ds.SeriesDate = metadataHelper.SeriesDate
    ds.SeriesTime = metadataHelper.SeriesTime
    ds.ContentDate = metadataHelper.ContentDate
    ds.ContentTime = metadataHelper.ContentTime

    # Manufacturer
    ds.Manufacturer = metadataHelper.Manufacturer
    ds.ManufacturerModelName = metadataHelper.ManufacturerModelName
    ds.ReferringPhysicianName = metadataHelper.ReferringPhysicianName

    # X-ray acquisition parameters
    ds.KVP = metadataHelper.KVP
    ds.ExposureTimeInms = metadataHelper.ExposureTimeInms
    ds.ExposureInmAs = metadataHelper.ExposureInmAs
    ds.XRayTubeCurrentInuA = metadataHelper.XRayTubeCurrentInuA
    ds.ReconstructionAlgorithm = metadataHelper.ReconstructionAlgorithm
    ds.ProtocolName = metadataHelper.ProtocolName
    ds.StationName = metadataHelper.StationName

    # Photometric interpretation & rescale
    ds.PhotometricInterpretation = metadataHelper.PhotometricInterpretation
    ds.RescaleIntercept = metadataHelper.RescaleIntercept
    ds.RescaleSlope = metadataHelper.RescaleSlope

    # Window preset
    ds.WindowWidth = metadataHelper.WindowWidth
    ds.WindowCenter = metadataHelper.WindowCenter

    return ds


# ---------------------------------------------------------------------------
# Slice-creation helpers (progress_callback made optional)
# ---------------------------------------------------------------------------


def create_dicom_slices_originalOrientation(metadataHelper, progress_callback=None):
    # Output directory
    dirName = (
        os.path.split(metadataHelper.filepath)[0]
        + "/"
        + os.path.splitext(os.path.split(metadataHelper.filepath)[1])[0]
        + "_originalOrientation_slices"
    )

    if not os.path.exists(dirName):
        os.makedirs(dirName)

    SeriesInstanceUID = pydicom.uid.generate_uid()
    FrameOfReferenceUID = pydicom.uid.generate_uid()
    SeriesNumber = str(random.randint(0, 10000))

    for index, CTslice in enumerate(metadataHelper.arr_rescaled):
        filename = dirName + f"/CT{str(index + 1)}.dcm"

        ds = set_metadata(metadataHelper, filename)

        # Series-level tags
        ds.SeriesInstanceUID = SeriesInstanceUID
        ds.FrameOfReferenceUID = FrameOfReferenceUID
        ds.SeriesDescription = metadataHelper.SeriesDescription
        ds.SeriesNumber = SeriesNumber
        ds.Rows, ds.Columns = CTslice.shape
        ds.PixelSpacing = metadataHelper.PixelSpacing
        ds.SliceThickness = metadataHelper.SliceThickness
        ds.NumberOfFrames = metadataHelper.NumberOfFrames
        ds.PatientPosition = "HFP"
        ds.ImageOrientationPatient = metadataHelper.ImageOrientationPatient

        # Instance-level tags
        SOPInstanceUID = pydicom.uid.generate_uid()
        ds.file_meta.MediaStorageSOPInstanceUID = SOPInstanceUID
        ds.SOPInstanceUID = SOPInstanceUID
        ds.PixelData = CTslice.tobytes()
        ds.InstanceNumber = str(index + 1)
        ds.NumberOfFrames = 1
        ds.ImagePositionPatient = [
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[0],
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[1],
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[2]
            - index * metadataHelper.SliceThickness,
        ]

        ds.save_as(filename)

        progress = (index + 1) / metadataHelper.NumberOfFrames * 100
        if progress_callback is not None:
            progress_callback(progress)


def create_dicom_slices_rayCommandOrientation(metadataHelper, progress_callback=None):
    """
    Create a RayStation-compatible orientation:
    ImageOrientationPatient = [1,0,0,0,1,0] and PatientPosition = 'HFS'.
    """
    dirName = (
        os.path.split(metadataHelper.filepath)[0]
        + "/"
        + os.path.splitext(os.path.split(metadataHelper.filepath)[1])[0]
        + "_rayOrientation_slices"
    )

    if not os.path.exists(dirName):
        os.makedirs(dirName)

    SeriesInstanceUID = pydicom.uid.generate_uid()
    FrameOfReferenceUID = pydicom.uid.generate_uid()
    SeriesNumber = str(random.randint(0, 10000))

    # Flip AP axis & rotate 180° around IS to match orientation
    arr_HFS = np.flip(metadataHelper.arr_rescaled, axis=1)

    for index, CTslice in enumerate(arr_HFS):
        filename = dirName + f"/CT{str(index + 1)}.dcm"

        ds = set_metadata(metadataHelper, filename)

        # Series-level tags
        ds.SeriesInstanceUID = SeriesInstanceUID
        ds.FrameOfReferenceUID = FrameOfReferenceUID
        ds.SeriesDescription = metadataHelper.SeriesDescription + "_rayOrientation"
        ds.SeriesNumber = SeriesNumber
        ds.Rows, ds.Columns = CTslice.shape
        ds.PixelSpacing = metadataHelper.PixelSpacing
        ds.SliceThickness = metadataHelper.SliceThickness
        ds.NumberOfFrames = metadataHelper.NumberOfFrames
        ds.PatientPosition = "HFS"
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]

        # Instance-level tags
        SOPInstanceUID = pydicom.uid.generate_uid()
        ds.file_meta.MediaStorageSOPInstanceUID = SOPInstanceUID
        ds.SOPInstanceUID = SOPInstanceUID
        ds.PixelData = CTslice.tobytes()
        ds.InstanceNumber = str(index + 1)
        ds.NumberOfFrames = 1
        ds.ImagePositionPatient = [
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[0],
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[1],
            metadataHelper.originalDICOM.PerFrameFunctionalGroupsSequence[0][0x0020, 0x9113][0][
                0x0020, 0x0032
            ].value[2]
            - index * metadataHelper.SliceThickness,
        ]

        ds.save_as(filename)

        progress = (index + 1) / metadataHelper.NumberOfFrames * 100
        if progress_callback is not None:
            progress_callback(progress)


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------


def volume_to_slices(
    volume_path: str | Path,
    orientations: tuple[str, ...] = ("original", "ray"),
    progress: bool = True,
) -> None:
    """
    Convert an Enhanced CT DICOM volume into 2-D slices.

    Parameters
    ----------
    volume_path
        Path to the Enhanced CT DICOM file.
    orientations
        Tuple containing "original" and/or "ray".
    progress
        If ``True`` a simple textual progress bar is printed to stdout.
    """
    volume_path = Path(volume_path)
    helper = metadataHelper(str(volume_path))

    def _print_progress(pct: float):
        if progress:
            print(f"\rProgress: {pct:6.2f}%", end="", flush=True)

    if "original" in orientations:
        create_dicom_slices_originalOrientation(helper, _print_progress if progress else None)
    if "ray" in orientations:
        create_dicom_slices_rayCommandOrientation(helper, _print_progress if progress else None)
    if progress:
        print("\nDone.")


# ---------------------------------------------------------------------------
# Simple CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Split an Enhanced CT DICOM volume into single-slice DICOM images."
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
        "--no-progress",
        action="store_true",
        help="Disable textual progress output",
    )

    args = parser.parse_args()
    volume_to_slices(
        args.dicom,
        tuple(args.orientations),
        progress=not args.no_progress,
    )
