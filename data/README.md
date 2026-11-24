# Input data

The two expression matrices are **not redistributed in this repository**. One of them
(`Leukemia_GSE28497.csv`, 106 MB) exceeds GitHub's 100 MB per-file limit, and both are already
publicly available from their original sources.

## How to obtain them

Download both curated matrices from the Curated Microarray Database (CuMiDa):

- <https://sbcb.inf.ufrgs.br/cumida>
- Datasets: **Leukemia_GSE28497** and **Leukemia_GSE71935**

Place the two CSV files in this directory, or edit the `D` variable at the top of each script in
`code/` to point at wherever you have saved them.

## What the files contain

| File | Samples | Classes | Probe sets | Platform |
|---|---|---|---|---|
| `Leukemia_GSE28497.csv` | 281 | 7 | 22,283 | GPL96 (Affymetrix U133A) |
| `Leukemia_GSE71935.csv` | 46 | 2 | 54,675 | GPL570 (Affymetrix U133 Plus 2.0) |

Column 1 is the sample identifier, column 2 (`type`) is the class label, and the remaining
columns are log2 expression values by probe set.

## Primary sources

The underlying raw data are deposited in NCBI GEO and are in the public domain:

- **GSE28497** — Coustan-Smith E, *et al.* New markers for minimal residual disease detection in
  acute lymphoblastic leukemia. *Blood.* 2011;117(23):6267–76. doi:10.1182/blood-2010-12-324004
- **GSE71935** — Helsmoortel HH, *et al.* LIN28B overexpression defines a novel fetal-like
  subgroup of juvenile myelomonocytic leukemia. *Blood.* 2016;127(9):1163–72.
  doi:10.1182/blood-2015-09-667808

## Note on the curated labels

CuMiDa provides no sample-level mapping back to GSM accessions. GEO lists 288 samples for
GSE28497 and 47 for GSE71935; the CuMiDa curations used here contain 281 and 46. The arrays
excluded during curation cannot be identified from the curated files. CuMiDa also labels 46
GSE28497 samples `B-CELL_ALL_T-ALL`; this is read as a prefix artefact and the class is treated
as T-ALL throughout. Both points are discussed in Section 2.1 of the paper.
