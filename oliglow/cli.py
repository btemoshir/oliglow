from pathlib import Path
from typing import Literal
from oliglow import deisotope, utils, averageine
from tap import Tap
import polars as pl
from clr_loader import get_mono


class Settings(Tap):
    rawfile: Path  # path to the .raw file to be deisotoped

    outfile: Path  # path to the deisotoped output file .tsv file

    mono: Path = "/opt/homebrew/Cellar/mono/6.12.0.206/lib/libmono-2.0.1.dylib"
    # Path to the mono installation.
    # On arm OSX (Mac) use the installation in the homebrew directory and not the default
    # (the default is only for the intel acrchitecture)

    # Deistoping options:
    avgine: str = "RNA_with_backbone"
    # Averagine to use for deisotoping.
    # Options are "RNA", "DNA" (to be implemented), "RNA_with_backbone", "RNA_with_thiophosphate_backbone"
    min_score: float = 150.0
    minimum_intensity: float = 0.0
    mass_error_tol: float = 0.02
    truncate_after: float = 0.8
    scale: Literal["sum"] = "sum"
    max_missed_peaks: int = 0
    error_tol: float = 2e-5
    scale_method: Literal["sum", "max"] = "sum"
    min_num_peaks_per_ms2_scan: int = 0


def main():
    settings = Settings().parse_args()
    # print(settings)

    # Check if the input file exists
    if not settings.rawfile.exists():
        raise FileNotFoundError(
            f"The specified raw file does not exist: {settings.rawfile}"
        )

    # Check if the output directory exists, create it if it doesn't
    if not settings.outfile.parent.exists():
        settings.outfile.parent.mkdir(parents=True, exist_ok=True)

    # Check if the mono library exists
    if not Path(settings.mono).exists():
        raise FileNotFoundError(
            f"The specified mono library does not exist: {settings.mono}"
        )

    # Load the mono library
    get_mono(libmono=settings.mono)

    if settings.avgine == "RNA":
        avgine = averageine.average_nucleoside_rna
    elif settings.avgine == "RNA_with_backbone":
        avgine = averageine.averagine_rna_with_backbone
    elif settings.avgine == "RNA_with_thiophosphate_backbone":
        avgine = averageine.averagine_rna_with_thiophosphate_backbone
    else:
        raise ValueError(
            f"Invalid averagine option: {settings.avgine}. Choose from 'RNA', 'RNA_with_backbone', or 'RNA_with_thiophosphate_backbone'."
        )

    # Create a parameter dictionary for the deisotope function:
    params_dict = {
        "avgine": avgine,
        "min_score": settings.min_score,
        "mass_error_tol": settings.mass_error_tol,
        "truncate_after": settings.truncate_after,
        "scale": settings.scale,
        "max_missed_peaks": settings.max_missed_peaks,
        "error_tol": settings.error_tol,
        "scale_method": settings.scale_method,
        "minimum_intensity": settings.minimum_intensity,
    }

    # Call the deisotope function with the parameters
    df = (
        deisotope.deisotope_all_ms2_scans(
            raw_file_path=settings.rawfile,
            deisotoping_parameters=params_dict,
            aggregate_masses=None,
            min_num_peaks_per_ms2_scan=settings.min_num_peaks_per_ms2_scan,
        )
        .sort("precursor_neutral_mass")
        .filter(pl.col("ms1_deisotope_success"))
    )

    df.write_csv(settings.outfile, separator="\t")
