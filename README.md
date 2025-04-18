<!-- # Oliglow

A Python package for deisotoping and deconvoluting nucleic acid sequences.

## Installation
```bash
pip install oliglow -->

<p align="center">
  <img src="https://img.shields.io/badge/Oliglow-%F0%9F%94%8D%20Clean%20Nucleic%20Data%20Tool-brightgreen?style=for-the-badge" alt="Oliglow Badge">
</p>

<h1 align="center">🧬 Oliglow</h1>

<p align="center"><i>Glow your oligos from the dark spectra</i></p>

<!-- <p align="center"><i>Illuminate your oligos. Strip away the noise.</i></p> -->

---

## 🌟 What is Oliglow?

**Oliglow** is a Python package for **deisotoping** and **deconvoluting** *mass spectrometry* of nucleic acid sequences and oligonucleotide spectra.

Whether you're cleaning mass spectrometry data, analysing complex oligo mixtures, or developing new bioinformatics tools — Oliglow helps you extract meaningful sequence signals with minimal fuss. 

It is meant to extract the fragment masses (MS1) from tandem mass spectrometry and the oligo (MS1) mass corresponding to each fragment. This is first pre-processing step in using tandem mass spectrometry for oligo sequencing using `lionelmssq`.

This package uses [`ms_deisotope`](https://github.com/mobiusklein/ms_deisotope) at its backend but creates a table aggregating fragments from each raw file where its possible to backtrack each fragment to its origin peaks and the origin scan.

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
        --avgine RNA \
        --min-score 50 \
        --mass-error-tol 10 \
        --truncate-after 1000 \
        --scale 1.0 \
        --max-missed-peaks 2 \
        --error-tol 0.01 \
        --scale-method intensity \
        --minimum-intensity 100 \
        --min-num-peaks-per-ms2-scan 5
```
```
Options
--rawfile: Path to the input raw file (required).
--outfile: Path to the output file (required).
--avgine: Averagine model to use (RNA, RNA_with_backbone, or RNA_with_thiophosphate_backbone).
--min-score: Minimum score for deisotoping.
--mass-error-tol: Mass error tolerance in ppm.
--truncate-after: Maximum number of peaks to consider.
--scale: Scaling factor for intensities.
--max-missed-peaks: Maximum number of missed peaks allowed.
--error-tol: Error tolerance for deisotoping.
--scale-method: Method for scaling intensities (intensity, etc.).
--minimum-intensity: Minimum intensity threshold for peaks.
--min-num-peaks-per-ms2-scan: Minimum number of peaks required per MS2 scan.
```

### Python

```python
from oliglow import deisotope

# Example usage with your data
#TODO

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