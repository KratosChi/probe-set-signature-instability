"""
04_review_statistics.py
Recomputations performed in response to peer review.

Reproduces, from the stored results files:
  * within-method Jaccard alongside Kuncheva                      (Reviewer 2, point 7)
  * Kuncheva recomputed with P = 2,000 instead of P = 22,215      (Reviewer 2, point 6)
  * conversion of Jaccard indices to shared probe-set counts      (Reviewer 2, point 9)
  * bootstrap 95% CIs for median Jaccard and balanced accuracy    (Editor, point 2)
  * Nadeau-Bengio corrected resampled t-test                      (Reviewer 2, point 8)
  * exact (unrounded) stability and accuracy ranges               (Editor, point 4)

Run from the repository root:  python code/04_review_statistics.py
"""
import json, itertools, numpy as np
from pathlib import Path
from scipy import stats

RES = Path(__file__).resolve().parent.parent / "results"
A  = json.loads((RES / "results_A.json").read_text())["A"]
B  = json.loads((RES / "results_B.json").read_text())["B"]
PF = json.loads((RES / "results_perf.json").read_text())
M  = ["ANOVA", "mRMR", "ReliefF", "LASSO", "RF"]
KS = ["10", "25", "50", "100"]
P_UNIVERSE, P_PREFILTER = 22215, 2000


def kuncheva_from_r(r, k, P):
    return (r - k * k / P) / (k - k * k / P)


def r_from_kuncheva(I, k, P):
    """Invert the index to recover the mean intersection size."""
    return I * (k - k * k / P) + k * k / P


print("=== Within-method Jaccard and Kuncheva (Supplementary Table S1) ===")
for name, R in [("GSE28497", A), ("GSE71935", B)]:
    print(name)
    for m in M:
        cells = " ".join(
            f"k={k}: J={R['within'][m][k]['jaccard']:.3f} I={R['within'][m][k]['kuncheva']:.3f}"
            for k in KS)
        print(f"  {m:8s} {cells}")

print("\n=== Kuncheva recomputed with P = 2,000 ===")
worst = 0.0
for name, R in [("GSE28497", A), ("GSE71935", B)]:
    for m in M:
        for k in KS:
            kk = int(k)
            I = R["within"][m][k]["kuncheva"]
            I2 = kuncheva_from_r(r_from_kuncheva(I, kk, P_UNIVERSE), kk, P_PREFILTER)
            worst = max(worst, abs(I2 - I))
print(f"  maximum absolute change across all 40 cells: {worst:.4f}")

print("\n=== Jaccard -> shared probe sets, k = 50 (GSE28497) ===")
cons = {m: set(A["consensus"]["50"][m]) for m in M}
pairs = []
for a, b in itertools.combinations(M, 2):
    inter = len(cons[a] & cons[b])
    J = inter / len(cons[a] | cons[b])
    pairs.append((f"{a}|{b}", J, inter))
for lbl, J, inter in sorted(pairs, key=lambda x: x[1]):
    print(f"  {lbl:16s} J={J:.3f}  shared={inter:2d}/50 ({100*inter/50:.0f}% of each signature)")
Js = [p[1] for p in pairs]
print(f"  median J = {np.median(Js):.3f}")

print("\n=== Bootstrap 95% confidence intervals (20,000 resamples) ===")
rng = np.random.default_rng(42)
def boot_ci(v, fn):
    bs = np.array([fn(rng.choice(v, len(v), replace=True)) for _ in range(20000)])
    return np.percentile(bs, 2.5), np.percentile(bs, 97.5)

lo, hi = boot_ci(np.array(Js), np.median)
print(f"  median pairwise Jaccard GSE28497 = {np.median(Js):.3f} (95% CI {lo:.3f}-{hi:.3f})")
consB = {m: set(B["consensus"]["50"][m]) for m in M}
JsB = [len(consB[a] & consB[b]) / len(consB[a] | consB[b]) for a, b in itertools.combinations(M, 2)]
lo, hi = boot_ci(np.array(JsB), np.median)
print(f"  median pairwise Jaccard GSE71935 = {np.median(JsB):.3f} (95% CI {lo:.3f}-{hi:.3f})")
for m in M:
    v = np.array(PF[m]["bacc_all"])
    lo, hi = boot_ci(v, np.mean)
    print(f"  balanced accuracy {m:8s} {v.mean():.3f} (95% CI {lo:.3f}-{hi:.3f})")

print("\n=== Nadeau-Bengio corrected resampled t-test (5-fold x 3 repeats) ===")
def nadeau_bengio(x, y, test_over_train=0.25):
    d = np.asarray(x) - np.asarray(y)
    n = len(d)
    se = np.sqrt((1 / n + test_over_train) * d.var(ddof=1))
    t = d.mean() / se
    return d.mean(), t, 2 * stats.t.sf(abs(t), n - 1)

for a, b, label in [("RF", "mRMR", "best vs worst accuracy"),
                    ("ReliefF", "RF", "most vs least stable")]:
    diff, t, p = nadeau_bengio(PF[a]["bacc_all"], PF[b]["bacc_all"])
    w = stats.wilcoxon(PF[a]["bacc_all"], PF[b]["bacc_all"]).pvalue
    print(f"  {label:24s} {a} vs {b}: diff={diff:+.4f} corrected P={p:.3f} (uncorrected Wilcoxon P={w:.3f})")

print("\n=== Exact ranges (Editor, point 4) ===")
st = [A["within"][m]["50"]["kuncheva"] for m in M]
ba = [PF[m]["bacc_mean"] for m in M]
print(f"  stability        {min(st):.4f} to {max(st):.4f}  range {max(st)-min(st):.4f}")
print(f"  balanced accuracy {min(ba):.4f} to {max(ba):.4f}  range {max(ba)-min(ba):.4f}")
print(f"  ratio {(max(st)-min(st))/(max(ba)-min(ba)):.1f}x")
