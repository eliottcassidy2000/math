# Direct finite-n test of THM-4755 statement 1 (the local-limit derivation), independent of Wormald's Hessian:
#   E[H'(T(T*, n))] / g'(dn) = P_walk(balanced) / P_iid(balanced),
# P_walk = Markov chain with P = (I + 2A)/d from the uniform start, 'balanced' = n visits to every state.
# (The identity itself: a path's block word w has weight prod(1/2 stay, 1 forward, 0 backward) * (n!)^d.)
# Computed exactly in floating point for growing n, then extrapolated polynomially in 1/n.
import sys
import numpy as np
from math import lgamma, log, exp
from fractions import Fraction
def ratio(A, n):
    d = len(A)
    P = (np.eye(d) + 2 * np.array(A, dtype=float)) / d
    shape = (n + 1,) * d + (d,)
    f = np.zeros(shape)
    for j in range(d):
        idx = [0] * d; idx[j] = 1
        f[tuple(idx) + (j,)] = 1.0 / d
    for _ in range(d * n - 1):
        g = np.zeros(shape)
        fp = f @ P  # fp[..., j] = sum_i f[..., i] P[i, j]
        for j in range(d):
            src = [slice(None)] * d; dst = [slice(None)] * d
            src[j] = slice(0, n); dst[j] = slice(1, n + 1)
            g[tuple(dst) + (j,)] += fp[tuple(src) + (j,)]
        f = g
    pw = f[(n,) * d].sum()
    log_piid = lgamma(d * n + 1) - d * lgamma(n + 1) - d * n * log(d)
    return pw / exp(log_piid)
def extrap(ns, vals):
    k = len(ns)
    M = np.array([[1.0 / m ** j for j in range(k)] for m in ns])
    return np.linalg.solve(M, np.array(vals))[0]
C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
R5 = [[1 if (j - i) % 5 in (1, 2) else 0 for j in range(5)] for i in range(5)]
R7 = [[1 if (j - i) % 7 in (1, 2, 3) else 0 for j in range(7)] for i in range(7)]
P7 = [[1 if (j - i) % 7 in (1, 2, 4) else 0 for j in range(7)] for i in range(7)]
for name, A, ns, target in [("C3", C3, list(range(20, 41, 4)), Fraction(2)), ("R5", R5, list(range(6, 15, 2)), Fraction(44, 19)),
                            ("R7", R7, [3, 4, 5, 6, 7], Fraction(3032, 1231)), ("P7", P7, [3, 4, 5, 6, 7], Fraction(64, 27))]:
    vals = [ratio(A, n) for n in ns]
    print(f"{name}: n={ns} ratio={['%.6f' % v for v in vals]}")
    for k in (2, 3, 4):
        print(f"    extrapolation from last {k} points: {extrap(ns[-k:], vals[-k:]):.6f}   (predicted rho = {float(target):.6f})")
