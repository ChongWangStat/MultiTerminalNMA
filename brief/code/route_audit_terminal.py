#!/usr/bin/env python3
"""Route audit on the terminal-current graph (study level) and on the aggregate
hat-matrix flow graph (treatment level) for every treatment pair of the COVID-19
network, common-effect model.

Study level: one vertex per treatment and per study; an edge treatment -> study
where a terminal current enters the study and study -> treatment where it leaves.
Treatment level: one edge per treatment pair with positive aggregate hat-matrix
flow w_jl (phi_j - phi_l).
For each graph: number of vertices |V|, edges |E|, complete routes M, rank of the
edge-by-route incidence matrix A, dim Pi(f) = M - rank(A) (= M - |E| + |V| - 2),
and whether route weights are unique.

Usage: python3 code/route_audit_terminal.py figures/data/covid_AE_study_polygon_currents.csv audit_output
"""
import csv, sys, itertools, os, numpy as np
from collections import defaultdict
src, outdir = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src)))
studies = defaultdict(list)
for r in rows: studies[r["study"]].append((r["t"], int(r["r"]), int(r["n"])))
trt = sorted({r["t"] for r in rows}); n = len(trt); ix = {t: i for i, t in enumerate(trt)}
L = np.zeros((n, n)); W = np.zeros((n, n)); blocks = {}
for s, arms in studies.items():
    v = np.array([1/(r+.5)+1/(nn-r+.5) for _, r, nn in arms]); g = 1/v
    K = np.diag(g)-np.outer(g, g)/g.sum(); idx = [ix[t] for t, _, _ in arms]
    L[np.ix_(idx, idx)] += K; blocks[s] = (idx, K)
    for i1, i2 in itertools.combinations(range(len(arms)), 2):
        w = g[i1]*g[i2]/g.sum(); W[idx[i1], idx[i2]] += w; W[idx[i2], idx[i1]] += w
Lp = np.linalg.pinv(L, rcond=1e-12); TOL = 1e-12

def audit(edges, src_v, snk_v):
    """edges: list of (u, v, f). Returns |V|, |E|, M, rank, dim, unique, formula check."""
    out = defaultdict(list)
    for i, (u, v, f) in enumerate(edges): out[u].append((i, v))
    routes = []
    def dfs(u, path):
        if u == snk_v: routes.append(list(path)); return
        for i, v in out[u]: path.append(i); dfs(v, path); path.pop()
    dfs(src_v, [])
    M = len(routes); A = np.zeros((len(edges), M))
    for r, p in enumerate(routes):
        for i in p: A[i, r] = 1
    rank = int(np.linalg.matrix_rank(A)); V = len({u for u, _, _ in edges} | {v for _, v, _ in edges})
    assert rank == len(edges) - V + 2, "rank formula"
    return V, len(edges), M, rank, M - rank, M == rank

lines = []; n_u_study = n_u_trt = 0
for a, b in itertools.combinations(trt, 2):
    q = np.zeros(n); q[ix[a]] = 1; q[ix[b]] = -1; phi = Lp @ q
    # study level: terminal currents
    e_s = []
    for s, (idx, K) in blocks.items():
        cur = K @ phi[idx]
        for j, c in zip(idx, cur):
            if c > TOL: e_s.append(("T"+trt[j], "S"+s, c))
            elif c < -TOL: e_s.append(("S"+s, "T"+trt[j], -c))
    Vs, Es, Ms, rs, ds, us = audit(e_s, "T"+a, "T"+b)
    # treatment level: aggregate hat-matrix flow
    e_t = [("T"+trt[j], "T"+trt[l], W[j, l]*(phi[j]-phi[l])) for j in range(n) for l in range(n)
           if W[j, l] > 0 and phi[j]-phi[l] > TOL]
    Vt, Et, Mt, rt, dt, ut = audit(e_t, "T"+a, "T"+b)
    n_u_study += us; n_u_trt += ut
    lines.append((f"{a}:{b}", Vs, Es, Ms, ds, us, Vt, Et, Mt, dt, ut))
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "route_audit_terminal_covid.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["target", "study_V", "study_E", "study_routes", "study_dim", "study_unique",
                "trt_V", "trt_E", "trt_routes", "trt_dim", "trt_unique"])
    for l in lines: w.writerow(l)
with open(os.path.join(outdir, "route_audit_terminal_rows.tex"), "w") as fh:
    for t, Vs, Es, Ms, ds, us, Vt, Et, Mt, dt, ut in lines:
        fh.write(f"{t} & {Vs} & {Es} & {Ms} & {ds} & {'yes' if us else 'no'} & {Et} & {Mt} & {dt} & {'yes' if ut else 'no'} \\\\\n")
for l in lines: print(*l)
print("unique at study level:", n_u_study, "of", len(lines), "; unique at treatment level:", n_u_trt, "of", len(lines))
