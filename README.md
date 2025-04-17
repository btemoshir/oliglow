<!-- # Oliglow

A Python package for deisotoping and deconvoluting nucleic acid sequences.

## Installation
```bash
pip install oliglow -->

<p align="center">
  <img src="https://img.shields.io/badge/Oliglow-%F0%9F%94%8D%20Clean%20Nucleic%20Data%20Tool-brightgreen?style=for-the-badge" alt="Oliglow Badge">
</p>

<h1 align="center">🧬 Oliglow</h1>

<p align="center"><i>Illuminate your oligos. Strip away the noise.</i></p>

---

## 🌟 What is Oliglow?

**Oliglow** is a Python package for **deisotoping** and **deconvoluting** nucleic acid sequences and oligonucleotide spectra.

Whether you're cleaning mass spectrometry data, analysing complex oligo mixtures, or developing new bioinformatics tools — Oliglow helps you extract meaningful sequence signals with minimal fuss.

---

## 🔧 Features

- ✅ Deisotoping of noisy nucleic acid signal data  
- ✅ Deconvolution of overlapping mass/charge peaks  
- ✅ Easy-to-use Python API  
- ✅ Lightweight and dependency-aware  
- ✅ Ideal for bioinformaticians and analytical labs

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

```python
from oliglow import deisotope, deconvolute

# Example usage with your data
cleaned_data = deisotope(raw_signal)
final_output = deconvolute(cleaned_data)

print(final_output)
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
  title        = {Oliglow: A Python package for deisotoping and deconvoluting nucleic acid spectra},
  year         = {2025},
  url          = {https://github.com/yourusername/oliglow},
  version      = {0.1.0},
  license      = {GPL-3.0}
}
```

> 💡 Tip: You can also create a DOI by archiving the package on [Zenodo](https://zenodo.org/) for easier citation.

---

## 🤝 Contributions

Pull requests and issues are welcome!  
Feel free to fork the repo, open an issue, or start a discussion if you have ideas or run into problems.

---

<p align="center">
  <b>🔬 Oliglow – Because your oligos deserve clarity.</b>
</p>
"""