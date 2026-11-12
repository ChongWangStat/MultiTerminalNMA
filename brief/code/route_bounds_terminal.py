#!/usr/bin/env python3
"""Recompute full B:E route allocation intervals on the terminal-current graph.

Inputs are the fixed five-study illustrative arm counts from CSP (not a
source-faithful clinical NMA). Distinguishes the intact-study bipartite
graph from the CSP edge graph; no within-study pairing is imposed.

Usage: python3 code/route_bounds_terminal.py \
    figures/data/covid_AE_study_polygon_currents.csv audit_output

Writes a CSV and TeX rows with the flow-proportional Markov selection and
sharp LP bounds for every complete B:E route, plus an independent
Kirchhoff/flow-margins residual check.
"""
import csv
import itertools
import os
import sys
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog

source = sys.argv[1] if len(sys.argv) > 1 else 'figures/data/covid_AE_study_polygon_currents.csv'
outdir = sys.argv[2] if len(sys.argv) > 2 else 'audit_output'
rows = list(csv.DictReader(open(source, newline='', encoding='utf-8')))
studies = defaultdict(list)
for row in rows:
    studies[row['study']].append((row['t'], int(row['r']), int(row['n'])))
trt = sorted({t for arms in studies.values() for t, _, _ in arms})
idx = {t: i for i, t in enumerate(trt)}
L = np.zeros((len(trt), len(trt)))
blocks = {}
for study, arms in studies.items():
    s = np.array([1 / (r + .5) + 1 / (n-r+.5) for _, r, n in arms])
    g = 1/s
    K = np.diag(g) - np.outer(g, g)/g.sum()
    ix = [idx[t] for t, _, _ in arms]
    L[np.ix_(ix, ix)] += K
    blocks[study] = (ix, K)
q = np.zeros(len(trt))
q[idx['B']] = 1.
q[idx['E']] = -1.
phi = np.linalg.pinv(L, rcond=1e-12) @ q
edges = []
for study, (ix, K) in blocks.items():
    current = K @ phi[ix]
    for j, c in zip(ix, current):
        if c > 1e-12:
            edges.append((trt[j], study, float(c)))
        elif c < -1e-12:
            edges.append((study, trt[j], float(-c)))
source_node = 'B'
sink_node = 'E'
next_edges = defaultdict(list)
for e, (u, v, flow) in enumerate(edges):
    next_edges[u].append(e)
paths = []
def dfs(node, route, nodes):
    if node == sink_node:
        paths.append((route[:], nodes[:]))
        return
    for e in next_edges[node]:
        u, v, _ = edges[e]
        assert v not in nodes, 'The positive support should be acyclic'
        dfs(v, route+[e], nodes+[v])
dfs(source_node, [], [source_node])
assert len(paths) == 8, 'Expected eight B:E study-terminal routes'
A = np.zeros((len(edges), len(paths)))
for r, (path, nodes) in enumerate(paths):
    for e in path:
        A[e, r] = 1
f = np.array([float(flow) for _, _, flow in edges])
outflow = defaultdict(float)
for u, v, fe in edges:
    outflow[u] += fe
markov = np.array([
    np.prod([edges[e][2] / outflow[edges[e][0]] for e in path])
    for path, _ in paths
])
residual = float(np.max(np.abs(A @ markov-f)))
assert residual < 1e-10, f'Inconsistent Markov allocation: {residual}'
assert abs(sum(markov)-1) < 1e-10
mins, maxs = [], []
for j in range(len(paths)):
    c = np.zeros(len(paths)); c[j] = 1.
    lo = linprog(c, A_eq=A, b_eq=f, bounds=(0, None), method='highs')
    hi = linprog(-c, A_eq=A, b_eq=f, bounds=(0, None), method='highs')
    assert lo.success and hi.success, (lo.message, hi.message)
    mins.append(max(float(lo.fun), 0))
    maxs.append(float(-hi.fun))
    assert mins[-1]-1e-10 <= markov[j] <= maxs[-1]+1e-10
assert np.linalg.matrix_rank(A) == 6
os.makedirs(outdir, exist_ok=True)
short = {'RECOVERY':'Rec','SOLIDARITY':'Sol','ACTT-1':'Act','CoDEX':'Cdx','REMAP-CAP':'Rmp'}
records = []
for j, (path, nodes) in enumerate(paths):
    lbl = '-'.join(short.get(node, node) for node in nodes)
    records.append({'route':lbl, 'markov':markov[j], 'lower':mins[j], 'upper':maxs[j]})
records = sorted(records, key=lambda row: row['route'])
with open(os.path.join(outdir, 'route_bounds_terminal_BE.csv'), 'w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['route', 'markov', 'lower', 'upper'])
    writer.writeheader()
    writer.writerows(records)
with open(os.path.join(outdir, 'route_bounds_terminal_BE_rows.tex'), 'w', encoding='utf-8') as fh:
    for r in records:
        fh.write(f"\\texttt{{{r['route']}}} & {r['markov']:.6f} & "
                 f"{r['lower']:.6f} & {r['upper']:.6f} \\\\\n")
with open(os.path.join(outdir, 'route_bounds_terminal_BE_rows.tex'), 'a', encoding='utf-8') as fh:
    fh.write('\\bottomrule\n')
print('Study-terminal B:E routes:', len(paths), '; dim Pi(f):', len(paths)-np.linalg.matrix_rank(A))
print('Max margin reconstruction residual:', residual)
for r in records:
    print(r['route'], f"Markov={r['markov']:.6f}",
          f"bounds=[{r['lower']:.6f},{r['upper']:.6f}]")
