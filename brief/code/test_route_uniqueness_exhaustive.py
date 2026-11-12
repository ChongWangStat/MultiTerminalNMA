#!/usr/bin/env python3
"""Exhaustively check the private-edge/full-rank equivalence on small DAGs."""
from itertools import combinations
import numpy as np

def paths(n, edges):
    adj = {v: [] for v in range(n)}
    for u, v in edges:
        adj[u].append(v)
    out = []
    def visit(v, es):
        if v == n-1:
            out.append(tuple(es))
            return
        for w in adj[v]:
            visit(w, es+[(v,w)])
    visit(0, [])
    return out

checked = 0
for n in (4, 5, 6):
    possible = list(combinations(range(n), 2))
    for mask in range(1 << len(possible)):
        edges = [e for i, e in enumerate(possible) if (mask >> i) & 1]
        pp = paths(n, edges)
        if len(pp) <= 1:
            continue
        used = sorted(set(e for p in pp for e in p))
        if len(pp) > len(used):
            full_rank = False
        else:
            a = np.array([[int(e in p) for p in pp] for e in used], dtype=float)
            full_rank = np.linalg.matrix_rank(a) == len(pp)
        private = all(any(sum(e in q for q in pp)==1 for e in p) for p in pp)
        assert full_rank == private, (n, edges, pp, full_rank, private)
        checked += 1
    print('checked n=', n)
print('Total enumerated DAG flow supports:', checked)
