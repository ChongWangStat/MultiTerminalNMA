#!/usr/bin/env python3
"""Recompute the depression-network reconstruction diagnostic from the public
Davies et al. (2022) repository files (AnnieDavies/NMA_and_RW).

Inputs (downloaded verbatim from the repository, main branch):
  Rucker2014.txt   contrast-level data: TE, seTE, treat1, treat2, studlab
  Hflow13.txt      aggregate hat-matrix row (evidence flow) for comparison 1 vs 3
  Short-Prop.txt   shortest-path contribution proportions for 1 vs 3
  RW_prop.txt      random-walk contribution proportions for 1 vs 3

Method. Arm variances are recovered within each trial from the pairwise
variances (v_i + v_j = seTE_ij^2). The common-effect NMA is the GLS fit with
study information matrices K_k = S_k^{-1} - g_k g_k^T / d_k (positive arm
variance class), which is also the netmeta multi-arm adjusted analysis. Pooled
direct effects use the Kron-reduced within-trial weights w_{k,jl} =
g_kj g_kl / d_k. The aggregate hat row is h_{jl} = w_{jl}(phi_j - phi_l) with
phi = L^+ (e_1 - e_3). The substitution diagnostic is
  theta_sub = sum_j s_j p_j ybar_j,  s_j = sign(h_j),  p_j = published share.
"""
import sys, json, itertools
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else "."
rows = []
with open(f"{d}/Rucker2014.txt") as f:
    hdr = f.readline().split()
    for line in f:
        p = line.split()
        if len(p) < 6:
            continue
        rows.append((float(p[1]), float(p[2]), int(p[3]), int(p[4]), p[5]))
TE = np.array([r[0] for r in rows]); SE = np.array([r[1] for r in rows])
T1 = np.array([r[2] for r in rows]); T2 = np.array([r[3] for r in rows])
ST = np.array([r[4] for r in rows])
n = 11
studies = sorted(set(ST), key=lambda s: int(s))
print("contrast records", len(rows), "studies", len(studies))

# --- arm variances per study
L = np.zeros((n, n))
Kblocks = {}; arms_of = {}; svec = {}
rhs = np.zeros(n)
w_direct = {}   # (j,l) -> list of (weight, effect l vs j i.e. TE with treat1=j, treat2=l)
for s in studies:
    idx = np.where(ST == s)[0]
    arms = sorted(set(T1[idx]) | set(T2[idx]))
    m = len(arms); pos = {a: i for i, a in enumerate(arms)}
    # solve v_i + v_j = se^2
    A = np.zeros((len(idx), m)); b = np.zeros(len(idx))
    for r, i in enumerate(idx):
        A[r, pos[T1[i]]] = 1; A[r, pos[T2[i]]] = 1; b[r] = SE[i] ** 2
    if m == 2:
        v = np.array([b[0] / 2, b[0] / 2])  # split arbitrarily; only the sum matters
    else:
        v, res, rk, _ = np.linalg.lstsq(A, b, rcond=None)
        assert np.allclose(A @ v, b, atol=1e-10), (s, A @ v - b)
    g = 1 / v; dk = g.sum()
    K = np.diag(g) - np.outer(g, g) / dk
    Kblocks[s] = K; arms_of[s] = arms; svec[s] = v
    Tk = np.zeros((m, n))
    for a in arms: Tk[pos[a], a - 1] = 1
    L += Tk.T @ K @ Tk
    # contrast-space contribution to the normal equations: sum_k X_k^T V_k^-1 y_k
    # baseline contrasts relative to first arm; y = TE for (first, other)
    base = arms[0]; others = arms[1:]
    y = np.zeros(m - 1); C = np.zeros((m - 1, m))
    for q, o in enumerate(others):
        C[q, pos[base]] = -1; C[q, pos[o]] = 1
        # find TE for (base, o): TE is effect of treat2 relative to treat1? take sign accordingly
        hit = [i for i in idx if (T1[i] == base and T2[i] == o)]
        if hit:
            y[q] = TE[hit[0]]
        else:
            hit = [i for i in idx if (T1[i] == o and T2[i] == base)]
            y[q] = -TE[hit[0]]
    V = C @ np.diag(v) @ C.T
    X = C @ Tk
    rhs += X.T @ np.linalg.solve(V, y)
    # Kron-reduced direct weights
    for i in idx:
        j, l = T1[i], T2[i]
        w = g[pos[j]] * g[pos[l]] / dk
        w_direct.setdefault((j, l), []).append((w, TE[i]))

Lp = np.linalg.pinv(L)
theta = Lp @ rhs   # arm-level treatment effects (identified up to constant)
q = np.zeros(n); q[0] = 1; q[2] = -1
est_1minus3 = q @ theta
print("theta_1 - theta_3 =", est_1minus3, "   theta_3 - theta_1 =", -est_1minus3)

# pooled direct effects (TE orientation: as recorded, treat1 -> treat2)
pairs_all = list(itertools.combinations(range(1, n + 1), 2))
assert len(pairs_all) == 55
ybar = {}; wsum = {}
for (j, l), lst in w_direct.items():
    W = sum(w for w, _ in lst); ybar[(j, l)] = sum(w * e for w, e in lst) / W; wsum[(j, l)] = W
print("distinct direct comparisons", len(ybar))

phi = Lp @ q
hagg = {}
for (j, l) in ybar:
    hagg[(j, l)] = wsum[(j, l)] * (phi[j - 1] - phi[l - 1])
# exact aggregate identity: sum_j h_j * (effect of l relative to j) ... orientation check
recon_a = sum(hagg[p] * ybar[p] for p in ybar)
print("sum h*ybar =", recon_a, " (compare with +/- estimate)")

# published vectors
def readvec(fn):
    return np.array([float(x) for x in open(f"{d}/{fn}").read().split()])
hflow = readvec("Hflow13.txt"); short = readvec("Short-Prop.txt"); rw = readvec("RW_prop.txt")
if len(hflow) == 54: hflow = np.append(hflow, 0.0)
print("sums of shares:", short.sum(), rw.sum())
pub = {pairs_all[i]: (hflow[i], short[i], rw[i]) for i in range(55)}
nz = {p for p, v in pub.items() if abs(v[0]) > 0 or v[1] > 0 or v[2] > 0}
print("published nonzero pairs:", len(nz), " observed:", len(ybar), " match:", nz <= set(ybar))

# compare our hat row with published flow (up to a global sign)
ours = np.array([hagg.get(p, 0.0) for p in pairs_all])
pubf = np.array([pub[p][0] for p in pairs_all])
sgn = np.sign(ours @ pubf)
print("max |ours - sgn*published flow| =", np.max(np.abs(ours - sgn * pubf)), " sign", sgn)

# substitution diagnostic with the published orientation s_j = sign(published flow)
def subst(shares):
    tot = 0.0
    for p in ybar:
        f, sh, rwv = pub[p]
        s = np.sign(f) if f != 0 else 0.0
        tot += s * shares[p] * ybar[p]
    return tot
theta_short = subst({p: pub[p][1] for p in pairs_all})
theta_rw = subst({p: pub[p][2] for p in pairs_all})
fitted = sgn * est_1minus3  # in the orientation of the published flow row
# the published flow row reconstructs the fitted estimate exactly:
recon_pub = sum(pub[p][0] * ybar[p] for p in ybar)
out = {
    "fitted_NMA_common_effect_1v3_published_orientation": float(recon_pub),
    "fitted_from_GLS_theta": float(fitted),
    "theta_sub_shortest_path": float(theta_short),
    "theta_sub_random_walk": float(theta_rw),
    "discrepancy_shortest": float(theta_short - recon_pub),
    "discrepancy_random_walk": float(theta_rw - recon_pub),
    "sum_shares_shortest": float(short.sum()),
    "sum_shares_random_walk": float(rw.sum()),
    "max_abs_diff_hat_row_vs_published_flow": float(np.max(np.abs(ours - sgn * pubf))),
    "n_trials": len(studies), "n_contrast_records": len(rows), "n_direct_comparisons": len(ybar),
}
print(json.dumps(out, indent=1))
json.dump(out, open(f"{d}/depression_diagnostic_output.json", "w"), indent=1)

# --- per-comparison table (LaTeX rows) -------------------------------------
with open(f"{d}/depression_percomparison_rows.tex", "w") as f:
    for p in pairs_all:
        if p not in ybar:
            continue
        fl, sh, rwv = pub[p]
        f.write("%d:%d & %.6f & %.6f & %.6f & %.6f \\\\\n" % (p[0], p[1], ybar[p], fl, sh, rwv))
print("wrote per-comparison rows")
