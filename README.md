# SpotSweeper_py_paper

**Data analysis and figure generation for the *SpotSweeper-py* manuscript**

This repository contains the **data analysis Jupyter notebooks** used in the SpotSweeper-py manuscript.  
It is intended for **reproducibility of the manuscript results** and **does not contain the SpotSweeper-py software package itself**.

---

## Manuscript

**Title:** *SpotSweeper-py: spatially-aware quality control metrics for spatial omics data in the Python ecosystem*  
**Authors:** Xingyi Chen, Michael Totty, Stephanie C. Hicks  
**Venue:** bioRxiv (2025)  
**DOI:** https://doi.org/10.64898/2025.12.06.692760

---

## SpotSweeper-py Package

The SpotSweeper-py software package is maintained in a **separate repository**:  
https://github.com/danielchen05/spotsweeper_py

It is also available on PyPI:  
https://pypi.org/project/spotsweeper/

---

## Analysis Notebooks

All analysis is organized as Jupyter notebooks.

| Notebook | Description |
|--------|------------|
| `Visium.ipynb` | Analysis of 10x Visium datasets; comparison of global vs local QC metrics |
| `VisiumHD.ipynb` | Analysis of 10x Visium HD data at higher spatial resolution |

Each notebook is intended to be run **top-to-bottom** after downloading the required data and setting appropriate paths.

---

## Data Availability & Download Instructions

Raw datasets used in the manuscript are **not stored directly in this repository** due to size constraints.

### Data Sources

The datasets can be downloaded from the 10x Genomics website:

- **Visium:** https://www.10xgenomics.com/datasets/gene-and-protein-expression-library-of-human-breast-cancer-cytassist-ffpe-2-standard  
- **Visium HD:** https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ffpe-if  

### How to Download

1. Download datasets from the sources listed above.
2. Extract them locally if needed.
3. Update file paths in the corresponding notebooks if necessary.


> Notes on preprocessing (e.g., filtering, binning, normalization) are documented **within each notebook**.  
> Downstream analyses expect an **AnnData** object; feel free to use your preferred method to load the data. (especially for the VisiumHD data)

---

## Questions?

For questions related to the analysis or manuscript figures, please contact:  
Xingyi (Daniel) Chen — xchen274@jh.edu

For questions related to the SpotSweeper-py software package, see:  
https://github.com/danielchen05/spotsweeper_py
