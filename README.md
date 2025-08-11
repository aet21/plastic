---
title: "Introduction to PLASTIC"
author:
- name: "Andrew E. Teschendorff"
  affiliation: 
  - Shanghai Institute for Nutrition and Health, SINH, CAS.
date: "2025-08-11"
package: plastic
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

PLASTIC (Predictive Landscape Analysis of Single-cell Transcriptomes for Identifying Cancer) is an R-package containing functions for estimating dedifferentiation (epigenetic reprogramming) and cellular plasticity at single-cell resolution from single-cell or single-nucleus RNA-Seq data. It is aimed mainly at multi-stage scRNA-Seq or snRNA-Seq studies of cancer development that have profiled cells from normal and precancerous lesions, in addition to cancer itself. PLASTIC can be used to help explore the cellular states of preneoplastic lesions and to identify those cells at higher cancer-risk.

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/plastic")
```

# References

Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of differentiation potency from a cell's transcriptome. Nat Commun. 2017 Jun 1;8:15599. doi: 10.1038/ncomms15599.

Teschendorff AE, Wang N. Improved detection of tumor suppressor events in single-cell RNA-Seq data. NPJ Genom Med. 2020 Oct 7;5:43. doi: 10.1038/s41525-020-00151-y .

Teschendorff AE, Maity AK, Hu X, Weiyan C, Lechner M. Ultra-fast scalable estimation of single-cell differentiation potency from scRNA-Seq data. Bioinformatics. 2021 Jul 12;37(11):1528-1534. doi: 10.1093/bioinformatics/btaa987 .

Liu T, Zhao X, Lin Y, Luo Q, Zhang S, Xi Y, Chen Y, Lin L, Fan W, Yang J, Ma Y, Maity AK, Huang Y, Wang J, Chang J, Lin D, Teschendorff AE, Wu C. Computational Identification of Preneoplastic Cells Displaying High Stemness and Risk of Cancer Progression. Cancer Res. 2022 Jul 18;82(14):2520-2537. doi: 10.1158/0008-5472.CAN-22-0668 .

Teschendorff AE, Feinberg AP. Statistical mechanics meets single-cell biology. Nat Rev Genet. 2021 Jul;22(7):459-476. doi: 10.1038/s41576-021-00341-z .

Chang J, Lu J, Liu Q, Xiang T, Zhang S, Yi Y, Li D, Liu T, Liu Z, Chen X, Dong Z, Li C, Yi H, Yu S, Huang L, Qu F, Wang M, Wang D, Dong H, Cheng G, Zhu L, Li J, Li C, Wu P, Xie X, Teschendorff AE, Lin D, Wang X, Wu C. Single-cell multi-stage spatial evolutional map of esophageal carcinogenesis. Cancer Cell. 2025 Mar 10;43(3):380-397.e7. doi: 10.1016/j.ccell.2025.02.00