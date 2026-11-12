#!/usr/bin/env python3
"""Recompute the COVID-19 common-effect signed study-terminal currents b_{ab,k}
from the arm counts in figures/data/covid_AE_study_polygon_currents.csv and
compare with the displayed A:E values; also report A:B and B:E currents and the
all-target fitted log odds ratios (corrected logits, arm variances
1/(r+1/2)+1/(n-r+1/2))."""
import csv, sys, numpy as np, itertools
p = sys.argv[1]
rows = list(csv.DictReader(open(p)))
studies = {}
for r in rows:
    studies.setdefault(r["study"], []).append((r["t"], int(r["r"]), int(r["n"]), float(r["signed_current"])))
trt = sorted({r["t"] for r in rows}); n = len(trt); ix = {t: i for i, t in enumerate(trt)}
L = np.zeros((n, n)); rhs = np.zeros(n); blocks = {}
for s, arms in studies.items():
    m = len(arms); T = np.zeros((m, n)); v = np.zeros(m); mu = np.zeros(m)
    for i, (t, rr, nn, _) in enumerate(arms):
        T[i, ix[t]] = 1; v[i] = 1/(rr+.5) + 1/(nn-rr+.5); mu[i] = np.log((rr+.5)/(nn-rr+.5))
    g = 1/v; K = np.diag(g) - np.outer(g, g)/g.sum()
    blocks[s] = (T, K, mu, [a[0] for a in arms], [a[3] for a in arms])
    L += T.T @ K @ T; rhs += T.T @ K @ mu
Lp = np.linalg.pinv(L); theta = Lp @ rhs
def currents(a, b):
    q = np.zeros(n); q[ix[a]] = 1; q[ix[b]] = -1; phi = Lp @ q
    return {s: (blk[1] @ (blk[0] @ phi)) for s, blk in blocks.items()}
cur = currents("A", "E"); maxerr = 0
for s, blk in blocks.items():
    for t, disp, got in zip(blk[3], blk[4], cur[s]):
        maxerr = max(maxerr, abs(disp - got))
print("A:E max |displayed - recomputed| =", maxerr)
print("fitted log OR (second minus first):")
for a, b in itertools.combinations(trt, 2):
    print(f"  {a}:{b} {theta[ix[b]]-theta[ix[a]]: .6f}  var {((np.eye(n)[ix[a]]-np.eye(n)[ix[b]]) @ Lp @ (np.eye(n)[ix[a]]-np.eye(n)[ix[b]])):.6f}")
for tgt in [("A","B"),("A","E"),("B","E")]:
    c = currents(*tgt); print(tgt, {s: np.round(c[s], 6).tolist() for s in c})
