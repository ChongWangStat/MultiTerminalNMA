#!/usr/bin/env python3
"""Multi-terminal study resistance for trials with correlated arms.

Recomputes Supplement Section S3 from the contrast standard errors published in
Rucker G, Schwarzer G. Trials and triangles: Network meta-analysis of multi-arm
trials with correlated arms. Res Synth Methods 2025;16(6):961-974
(Table 3: real example; Table D1: fictitious four-trial network).
Only standard errors are used; the terminal conductance matrix K_k, the network
variances and the terminal currents do not depend on the effect estimates.

Usage: python3 code/correlated_arm_examples.py > audit_output/correlated_arm_examples.txt
"""
import itertools, numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=140)

def K_three_arm(v12, v13, v23):
    """Terminal conductance of a three-arm trial from its three contrast variances."""
    cov = (v12 + v13 - v23) / 2                      # cov(y2 - y1, y3 - y1)
    S = np.array([[v12, cov], [cov, v13]]); C = np.array([[-1., 1, 0], [-1., 0, 1]])
    return C.T @ np.linalg.inv(S) @ C, S

def implied_arm_variances(v12, v13, v23):
    return np.array([(v12 + v13 - v23) / 2, (v12 + v23 - v13) / 2, (v13 + v23 - v12) / 2])

def assemble(trts, blocks):
    n = len(trts); ix = {t: i for i, t in enumerate(trts)}; L = np.zeros((n, n))
    for arms, K in blocks.values():
        idx = [ix[t] for t in arms]; L[np.ix_(idx, idx)] += K
    return ix, np.linalg.pinv(L, rcond=1e-12)

def report(trts, blocks, targets):
    ix, Lp = assemble(trts, blocks); n = len(trts)
    print("network variances:", {f"{a}:{b}": round(float((np.eye(n)[ix[a]] - np.eye(n)[ix[b]]) @ Lp @ (np.eye(n)[ix[a]] - np.eye(n)[ix[b]])), 4)
                                 for a, b in itertools.combinations(trts, 2)})
    for a, b in targets:
        q = np.zeros(n); q[ix[a]] = 1; q[ix[b]] = -1; phi = Lp @ q
        print(f"comparison {a}:{b}  potentials (sink = 0):", {t: round(float(phi[ix[t]] - phi[ix[b]]), 5) for t in trts})
        tot = np.zeros(n)
        for s, (arms, K) in blocks.items():
            idx = [ix[t] for t in arms]; cur = K @ phi[idx]; tot[idx] += cur
            print(f"    {s:8s} terminal currents", {t: round(float(c), 4) for t, c in zip(arms, cur)})
        assert np.allclose(tot, q, atol=1e-9), "Kirchhoff balance"

print("=== (a) Real example (Table 3): Batch7 A-Control, Batch9 B-Control, Batch10 (Control, B, C) correlated ===")
se = dict(BCtl=0.1700709, CCtl=0.2739463, BC=0.2132226)
v = {k: x * x for k, x in se.items()}
print("contrast variances of Batch10:", {k: round(x, 6) for k, x in v.items()},
      " triangle check V(C-Ctl) <= V(B-Ctl) + V(B-C):", round(v["CCtl"], 6), "<=", round(v["BCtl"] + v["BC"], 6), v["CCtl"] <= v["BCtl"] + v["BC"])
K10, S10 = K_three_arm(v["BCtl"], v["CCtl"], v["BC"])       # terminals (Control, B, C)
print("implied arm variances (Control, B, C):", implied_arm_variances(v["BCtl"], v["CCtl"], v["BC"]))
print("det Sigma =", np.linalg.det(S10), " K (Control, B, C) =\n", K10)
k2 = lambda s: np.array([[1., -1], [-1, 1]]) / s ** 2
blocks = {"Batch7": (["Control", "A"], k2(0.3316958)), "Batch9": (["Control", "B"], k2(0.2240514)), "Batch10": (["Control", "B", "C"], K10)}
report(["Control", "A", "B", "C"], blocks, [("B", "C")])

print("\n=== (b) Fictitious network (Table D1): four three-arm trials with correlated arms ===")
se = {"ABC": (1.0, 1.7, 1.0), "ABD": (1.0, 6.5, 6.9), "ACD": (1.7, 6.5, 6.5), "BCD": (1.0, 6.9, 6.5)}   # (12, 13, 23)
blocks = {}
for s, (a, b, c) in se.items():
    K, _ = K_three_arm(a * a, b * b, c * c); blocks[s] = (list(s), K)
    print(s, "implied arm variances", np.round(implied_arm_variances(a * a, b * b, c * c), 3), " K =", np.round(K, 4).tolist())
print("published NMA variances (Table D1): A:B 0.364, A:C 0.754, B:C 0.364, A:D 13.502, B:D 14.019, C:D 13.502")
report(list("ABCD"), blocks, [("B", "D"), ("A", "D")])

print("\n=== (c) Four terminals: every triangle inequality holds, yet one off-diagonal entry of K is positive ===")
n = 4; K = np.zeros((n, n))
for i, j in itertools.combinations(range(n), 2):
    c = -0.05 if (i, j) == (0, 1) else 1.0
    K[i, i] += c; K[j, j] += c; K[i, j] -= c; K[j, i] -= c
Kp = np.linalg.pinv(K); R = lambda i, j: Kp[i, i] + Kp[j, j] - 2 * Kp[i, j]
print("eigenvalues of K:", np.round(np.linalg.eigvalsh(K), 4), " K[1,2] =", K[0, 1])
print("contrast variances (effective resistances):", {f"{i+1}{j+1}": round(float(R(i, j)), 4) for i, j in itertools.combinations(range(n), 2)})
viol = [(i, j, k) for i, j in itertools.combinations(range(n), 2) for k in range(n) if k not in (i, j) and R(i, j) > R(i, k) + R(k, j) + 1e-12]
print("triangle violations:", viol)

print("\n=== (d) Redundant contrasts: all three pairwise contrasts with the pseudoinverse give the same K ===")
D = np.array([[-1., 1, 0], [-1., 0, 1], [0, -1., 1]]); C = D[:2]; H = D @ np.linalg.pinv(C)
Kfull = D.T @ np.linalg.pinv(H @ S10 @ H.T, rcond=1e-12) @ D
print("max |K_full - K| for Batch10:", float(np.abs(Kfull - K10).max()))
