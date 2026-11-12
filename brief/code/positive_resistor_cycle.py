#!/usr/bin/env python3
"""Verify a directed study-terminal cycle even when all resistor weights are positive."""
import numpy as np

K = np.zeros((4, 4))
for i in range(4):
    for j in range(i + 1, 4):
        g = 10.0 if {i, j} in ({0, 1}, {2, 3}) else 0.3
        assert g > 0
        K[i, i] += g
        K[j, j] += g
        K[i, j] -= g
        K[j, i] -= g
H = np.zeros((4, 4))
H[1, 1] = H[2, 2] = 1.0
H[1, 2] = H[2, 1] = -1.0
q = np.array([1., 0., 0., -1.])
potentials = np.linalg.pinv(K + H) @ q
terminal = K @ potentials
two_arm = H @ potentials
assert np.max(np.abs(terminal + two_arm - q)) < 1e-12
assert np.all(K[np.triu_indices(4, 1)] < 0)
assert np.allclose(terminal, [1., -0.417025, 0.417025, -1.], atol=1e-5)
assert terminal[1] < 0 < two_arm[1] and terminal[2] > 0 > two_arm[2]
print("Four-arm terminal currents:", np.array2string(terminal, precision=6))
print("Two-arm terminal currents:", np.array2string(two_arm, precision=6))
print("Cycle: 2 -> two-arm study -> 3 -> four-arm study -> 2")
print("All pairwise resistor conductances are positive.")
