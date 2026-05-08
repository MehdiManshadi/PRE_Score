# PRE Score: Predicted Regulatory Effects of GWAS Variants

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![RegulomeDB](https://img.shields.io/badge/Data-RegulomeDB-orange)](https://regulomedb.org/)

A Python pipeline to compute **Predicted Regulatory Effects (PRE)** of GWAS-associated SNPs and their linkage disequilibrium (LD) proxies, by querying the [RegulomeDB](https://regulomedb.org/) API and aggregating cell-type-specific regulatory evidence.

This tool implements the scoring framework described in:

> *Guo et al. (2019). A systems biology approach uncovers cell-specific gene regulatory effects of genetic associations in multiple sclerosis.* Nature Communications. [DOI: 10.1038/s41467-019-09773-y](https://doi.org/10.1038/s41467-019-09773-y)

---

## Table of Contents

- [Background](#background)
- [Workflow Overview](#workflow-overview)
- [Installation](#installation)
- [Input Format](#input-format)
- [Usage](#usage)
- [Output](#output)
- [RegulomeDB API Response Structure](#regulomedb-api-response-structure)
- [Citation](#citation)

---

## Background

GWAS studies identify SNPs associated with disease, but these variants are often in non-coding regions. Interpreting their functional impact requires linking them to regulatory elements and target genes in a cell-type-specific manner.

This pipeline:
1. Takes a list of GWAS lead SNPs
2. Expands them with LD proxy SNPs (r² ≥ 0.5)
3. Queries each SNP against RegulomeDB to retrieve overlapping regulatory features (ChIP-seq peaks, histone marks, eQTLs, DNase sites, PWMs, etc.)
4. Computes a **Sum of Weighted Weights (SWW)** per gene per cell type — the **PRE score**

---

## Workflow Overview
