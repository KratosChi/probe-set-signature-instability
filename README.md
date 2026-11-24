# Probe-set signature instability across five feature selection algorithms

Analysis code and complete numerical results for:

> Singh C, Majji SSR, Singh S, Chaudhary P. **Quantifying probe-set signature instability
> across five feature selection algorithms in curated leukemia microarray cohorts:
> a comparative benchmark.** *Gene Expression* (under review).

Five standard feature selection algorithms were applied to two curated leukemia microarray
cohorts to test whether they are reproducible **with themselves** under resampling, and whether
they agree **with each other** on the probe sets they select. They are, and they largely do not:
at a signature size of 50, the median pairwise Jaccard index between algorithms was 0.212, and
59.2% of the pooled signature was selected by a single algorithm and no other.

## Repository layout

```
code/      analysis and figure scripts, numbered in execution order
results/   complete numerical results as JSON (every method, signature size and cohort)
data/      instructions for obtaining the two input matrices (not redistributed here)
```

## Scripts

| Script | What it does |
|---|---|
| `code/01_stability_analysis.py` | Within-method stability (Kuncheva, Jaccard), between-method agreement and selection multiplicity for both cohorts. Writes `results_A.json`, `results_B.json`. |
| `code/02_subtype_and_classification.py` | Per-subtype one-versus-rest stability, class separability (silhouette), and cross-validated classification performance. Writes `results_perclass.json`, `results_perf.json`. |
| `code/03_lasso_revision_analysis.py` | Added at revision: counts non-zero LASSO coefficients per fit and re-tunes the L1 penalty by inner five-fold cross-validation. Writes `lasso_rev_A.json`, `lasso_rev_B.json`. |
| `code/04_review_statistics.py` | Added at revision: reproduces every statistic reported in the response to reviewers, including the Kuncheva recomputation with P = 2,000, bootstrap confidence intervals and the Nadeau–Bengio corrected t-test. Reads the JSON files only — no input data needed. |
| `code/fig1_workflow.py` | Figure 1 (study design). |
| `code/fig23_stability_agreement.py` | Figures 2 and 3. |
| `code/fig45_subtype_accuracy.py` | Figures 4 and 5, plus the statistics in Results 3.5. |
| `code/graphical_abstract.py` | Graphical abstract. |

## Reproducing the analysis

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# obtain the two input matrices first (see data/README.md), then:
python code/01_stability_analysis.py
python code/02_subtype_and_classification.py
python code/03_lasso_revision_analysis.py
python code/04_review_statistics.py        # runs from results/ alone
python code/fig23_stability_agreement.py
python code/fig45_subtype_accuracy.py
```

Each script reads the input matrices from the path set in its `D` variable near the top; edit
that to point at your local copy. All random seeds are fixed (seed 42), so re-running reproduces
the published values exactly.

`04_review_statistics.py` needs no input data and is the quickest way to verify the reported
numbers.

## Results files

All results are plain JSON.

| File | Contents |
|---|---|
| `results_A.json` | GSE28497: within-method Kuncheva and Jaccard at k = 10/25/50/100, pairwise between-method Jaccard, consensus signatures, selection multiplicity, probe selection frequencies. |
| `results_B.json` | The same for GSE71935. |
| `results_perclass.json` | Per-subtype stability, silhouette widths and Spearman correlations. |
| `results_perf.json` | Balanced accuracy and macro-F1 per algorithm, including all 15 individual fold scores. |
| `results_stats.json` | Correlations and tests reported in Results 3.5. |
| `lasso_rev_A.json`, `lasso_rev_B.json` | Non-zero coefficient counts, selected penalties and stability under fixed versus tuned C. |

## Software

Python 3.11. Exact versions in `requirements.txt`. The analysis runs on a laptop; the full
pipeline takes roughly 90 minutes, dominated by ReliefF and mRMR.

## Licence

Code is released under the MIT Licence (see `LICENSE`). The input data are redistributed by
CuMiDa and NCBI GEO under their own terms; see `data/README.md`.

## Citation

If you use this code, please cite the paper above and the CuMiDa database:

> Feltes BC, Chandelier EB, Grisci BI, Dorn M. CuMiDa: an extensively curated microarray database
> for benchmarking and testing of machine learning approaches in cancer research.
> *J Comput Biol.* 2019;26(4):376–86. doi:10.1089/cmb.2018.0238
