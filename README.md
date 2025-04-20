<!-- # Oliglow

A Python package for deisotoping and deconvoluting nucleic acid sequences.

## Installation
```bash
pip install oliglow -->

<p align="center">
  <img src="https://img.shields.io/badge/Oliglow-%F0%9F%94%8D%20Deisotope%20Nucleic%20Acid%20Oligos-brightgreen?style=for-the-badge" alt="Oliglow Badge">
</p>

<h1 align="center">🧬 Oliglow</h1>

<p align="center"><i>Glow your oligos from the dark spectra</i></p>

<!-- <p align="center"><i>Illuminate your oligos. Strip away the noise.</i></p> -->

---

## 🌟 What is Oliglow?

**Oliglow** is a Python package for **deisotoping** and **deconvoluting** *mass spectrometry* of nucleic acid sequences and oligonucleotide spectra.

Whether you're cleaning mass spectrometry data, analysing complex oligo mixtures, or developing new bioinformatics tools — Oliglow helps you extract meaningful sequence signals with minimal fuss. 

It is designed to extract fragment masses (MS1) from tandem mass spectrometry data and map them to the corresponding oligonucleotide (MS1) masses for each fragment. This process serves as a critical pre-processing step in utilizing tandem mass spectrometry for oligonucleotide sequencing. By leveraging the `lionelmssq` library, **Oliglow** ensures accurate and efficient extraction of these masses, enabling researchers to streamline their workflows and gain deeper insights into nucleic acid sequences. This functionality is particularly valuable for applications in bioinformatics, molecular biology, and analytical chemistry, where precise mass spectrometry data is essential for downstream analysis.

This package leverages the powerful [`ms_deisotope`](https://github.com/mobiusklein/ms_deisotope) library as its backend for performing deisotoping and deconvolution tasks. However, **Oliglow** extends its functionality by generating a comprehensive table that aggregates fragments from each raw file. 

The resulting table provides detailed insights, allowing users to trace each fragment back to its original peaks and the corresponding scan from which it was derived. This feature ensures transparency and reproducibility in the analysis, making it easier to validate results and integrate them into downstream workflows.

---

## 🔧 Features

- ✅ Deisotoping of noisy nucleic acid signal data  
- ✅ Deconvolution of overlapping mass/charge peaks  
- ✅ Easy-to-use Python API  
- ✅ Backtracking support for each output fragment  
- ✅ Command line tools
- ✅ Easy-to-use Python API  
- ✅ Lightweight and dependency-aware  

---

## 📦 Installation

First, make sure you have a Python 3.7+ environment (Micromamba, Conda, or virtualenv).

### Option 1: Pip (Editable Mode)

```bash
# Inside your active environment:
pip install -e .
```

### Option 2: Micromamba

```bash
micromamba activate your-env
cd /path/to/oliglow
pip install -e .
```

> 🔁 This installs the package in **editable mode**, so your code changes reflect immediately.

---

## 🚀 Usage

### Command Line Interface (CLI)

The `oliglow` CLI provides tools for deisotoping and deconvolution of tandem mass spectrometry (LCMS) data for nucleic acid oligos.

#### Example Usage

```bash
oliglow --rawfile /path/to/input.raw \
        --outfile /path/to/output.tsv \
        --avgine "RNA_with_backbone" \
        --min-score 150 \
        --mass-error-tol 0.02 \
        --truncate-after 0.8 \
        --scale "sum" \
        --max-missed-peaks 0 \
        --error-tol 2e-5 \
        --scale-method "sum" \
        --minimum-intensity 0.0 \
        --min-num-peaks-per-ms2-scan 0
```
```
Options
--rawfile: Path to the input raw file (required).
--outfile: Path to the output file (required).
--avgine: Averagine model to use (RNA, RNA_with_backbone, or RNA_with_thiophosphate_backbone).
--min-score: Minimum score for deisotoping.
--mass-error-tol: Mass error tolerance in ppm.
--truncate-after: How much of the deisoptoping envelope to consider. For MS2 scans, typical value should be 0.8.
--scale: Scaling factor for intensities.
--max-missed-peaks: Maximum number of missed peaks allowed.
--error-tol: Error tolerance for deisotoping (in ppm).
--scale-method: Method for scaling intensities ("sum", "max",...).
--minimum-intensity: Minimum intensity threshold for any peak to be considered.
--min-num-peaks-per-ms2-scan: Minimum number of peaks required per MS2 scan for the scan to be considered.
```

### Python

```python
from oliglow import deisotope, averageine

# Example usage with your data
params_dict = {
    "avgine": averageine.averagine_rna_with_backbone,
    "min_score": 150.0,
    "mass_error_tol": 0.02,
    "truncate_after": 0.8,
    "scale": "sum",
    "max_missed_peaks": 0,
    "error_tol": 2e-5,
    "scale_method": "sum",
    "minimum_intensity": 0.0,
}

df = deisotope.deisotope_all_ms2_scans("/path/to/input.raw",
    deisotoping_parameters=params_dict,
    aggregate_masses=None,
    min_num_peaks_per_ms2_scan=0,
)

#Export this to a tsv file!
df.write_csv("/path/to/input.tsv",separator="\t")

```

You can also build your own analysis pipeline using these functions as components.

---

## 📁 Project Structure

```
Oliglow/
├── oliglow/            # Package code
│   ├── __init__.py
│   └── averagine.py
    └── cli.py
    └── deisotope.py
    └── utils.py
├── tests/              # Unit tests
├── README.md
├── LICENSE
├── pyproject.toml
```

---

## 📜 License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**.  
See the [LICENSE](./LICENSE) file for full details.

---

## 📚 Citation

If you use **Oliglow** in your research or software, please cite it as follows:

```
@software{oliglow,
  author       = {Moshir Harsh},
  title        = {Oliglow: A Python package for deisotoping and deconvoluting nucleic acid mass spectrometry spectra},
  year         = {2025},
  url          = {https://github.com/yourusername/oliglow},
  version      = {0.1.0},
  license      = {GPL-3.0}
}
```

<!-- > 💡 Tip: You can also create a DOI by archiving the package on [Zenodo](https://zenodo.org/) for easier citation. -->

---

## 🤝 Contributions

Pull requests and issues are welcome!  
Feel free to fork the repo, open an issue, or start a discussion if you have ideas or run into problems.

---

<p align="center">
  <b>🔬 Oliglow – Because your oligos deserve clarity.</b>
</p>
"""