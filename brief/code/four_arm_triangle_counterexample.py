#!/usr/bin/env python3
"""Four-arm example: all variance triangles hold but positive resistors fail."""
from itertools import permutations
import numpy as np

k = np.zeros((4, 4))
for i in range(4):
    for j in range(i + 1, 4):
        conductance = -0.05 if (i, j) == (0, 1) else 1.0
        k[i, i] += conductance
        k[j, j] += conductance
        k[i, j] -= conductance
        k[j, i] -= conductance

c = np.column_stack((np.ones(3) * -1, np.eye(3)))
pinv = np.linalg.pinv(k)
sigma = c @ pinv @ c.T
svals = np.linalg.eigvalsh(sigma)
r = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        d = np.eye(4)[i] - np.eye(4)[j]
        r[i, j] = d @ pinv @ d
triangles = [r[i, j] <= r[i, m] + r[m, j] + 1e-12
             for i, j, m in permutations(range(4), 3)]

assert np.all(np.linalg.eigvalsh(k)[1:] > 0)
assert np.all(svals > 0)
assert np.all(triangles)
assert k[0, 1] > 0
np.testing.assert_allclose(c.T @ np.linalg.inv(sigma) @ c, k, atol=1e-12)
print('K =')
print(np.array2string(k, precision=6))
print('Sigma eigenvalues:', np.array2string(svals, precision=6))
print('Within-study contrast variances =')
print(np.array2string(r, precision=6))
print('All 24 oriented triangle checks:', all(triangles))
print('Positive off-diagonal K[1,2] =', k[0, 1])
print('Ordinary positive resistors cannot realize this study.')
