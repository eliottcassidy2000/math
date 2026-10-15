#!/usr/bin/env python3
"""Audit B (2026-10-14): independent rebuild of the Voronov-Neopryatnaya-Dergachev spindle-free graph (series 2,
case 10) and of the session's sub-pieces, with an exact representation and a neighbour search that share no code
with the session's hn_udg.py / spindle_free_build.py.

Representation.  Every point is a complex number in K = Q(zeta)(sqrt5), zeta = exp(i pi/12) a primitive 24th root of
unity, stored as 16 integer numerators over one global denominator D:
    z = (sum_j c_j zeta^j + sqrt5 * sum_j c'_j zeta^j) / D,   j = 0..7,
reduced modulo Phi_24(x) = x^8 - x^4 + 1 (so the representation is canonical: {zeta^j} is a Q-basis of Q(zeta),
and sqrt5 is not in Q(zeta)).  The session uses real coordinates x, y in Q(sqrt2, sqrt3, sqrt5) instead.

Unit distance.  |z - w|^2 = (z - w) * conj(z - w) computed exactly in K (conj: zeta -> zeta^-1, sqrt5 -> sqrt5),
compared with 1.

Neighbour search.  Brute force over ALL pairs, in blocks, in float64: a pair is a candidate if |d^2 - 1| < 1e-6
(float error is ~1e-15, so no exact unit pair can be missed); every candidate is then decided exactly.  The session
uses a KD-tree ball query with a 1e-9 window instead.

Usage: python3 -u hn_spindle_free_rebuild_zeta24.py OUTDIR
Writes OUTDIR/<name>.npz with float coordinates, exact numerators and the edge list for each graph.
"""
import sys
import time
from fractions import Fraction as Fr
from math import lcm
import cmath

import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
T0 = time.time()


def log(*a):
    print(*a, f"[{time.time() - T0:.1f}s]", flush=True)


# ------------------------------------------------------------------ Q(zeta_24) arithmetic on Fractions
ZETA = cmath.exp(1j * cmath.pi / 12)


def zred(c):
    """reduce a coefficient list (any length) modulo x^8 - x^4 + 1 -> length 8"""
    c = list(c) + [Fr(0)] * max(0, 8 - len(c))
    for d in range(len(c) - 1, 7, -1):
        a = c[d]
        if a:
            c[d] = Fr(0)
            c[d - 4] += a      # x^d = x^(d-4) * x^4 ... using x^8 = x^4 - 1
            c[d - 8] -= a
    return c[:8]


def zmul(a, b):
    out = [Fr(0)] * 15
    for i in range(8):
        if a[i]:
            for j in range(8):
                if b[j]:
                    out[i + j] += a[i] * b[j]
    return zred(out)


def zpow_elem(j):
    """zeta^j for any integer j as a reduced vector"""
    j %= 24
    c = [Fr(0)] * (j + 1)
    c[j] = Fr(1)
    return zred(c)


def zconj(a):
    out = [Fr(0)] * 8
    for j in range(8):
        if a[j]:
            e = zpow_elem(-j)
            for k in range(8):
                out[k] += a[j] * e[k]
    return out


def zval(a):
    return sum(complex(float(a[j])) * ZETA ** j for j in range(8))


# K elements: (a, b) meaning a + b sqrt5, a, b in Q(zeta)
def K(a, b=None):
    return (list(a), list(b) if b is not None else [Fr(0)] * 8)


def kadd(x, y):
    return ([u + v for u, v in zip(x[0], y[0])], [u + v for u, v in zip(x[1], y[1])])


def kmul(x, y):
    a0 = zmul(x[0], y[0])
    a1 = zmul(x[1], y[1])
    b0 = zmul(x[0], y[1])
    b1 = zmul(x[1], y[0])
    return ([u + 5 * v for u, v in zip(a0, a1)], [u + v for u, v in zip(b0, b1)])


def kconj(x):
    return (zconj(x[0]), zconj(x[1]))


def kval(x):
    return zval(x[0]) + 5 ** 0.5 * zval(x[1])


def kscal(q, x):
    return ([q * u for u in x[0]], [q * u for u in x[1]])


ONE = K(zpow_elem(0))
ZERO = K([Fr(0)] * 8)
zeta = K(zpow_elem(1))
i_ = K(zpow_elem(6))
sqrt2 = K([Fr(0), Fr(1), Fr(0), Fr(1), Fr(0), Fr(-1), Fr(0), Fr(0)])        # zeta + zeta^3 - zeta^5
sqrt3 = K([Fr(0), Fr(0), Fr(2), Fr(0), Fr(0), Fr(0), Fr(-1), Fr(0)])        # 2 zeta^2 - zeta^6
sqrt5 = K([Fr(0)] * 8, zpow_elem(0))
for name, e, val in (('sqrt2', sqrt2, 2 ** 0.5), ('sqrt3', sqrt3, 3 ** 0.5), ('sqrt5', sqrt5, 5 ** 0.5),
                     ('i', i_, 1j)):
    assert abs(kval(e) - val) < 1e-12, name
assert kmul(sqrt2, sqrt2) == kscal(2, ONE) and kmul(sqrt3, sqrt3) == kscal(3, ONE)
assert kmul(sqrt5, sqrt5) == kscal(5, ONE) and kmul(i_, i_) == kscal(-1, ONE)

phi0 = zeta
phi1 = kscal(Fr(1, 3), kadd(kmul(sqrt2, sqrt3), kmul(i_, sqrt3)))               # (sqrt6 + i sqrt3)/3
rho = kadd(kscal(Fr(7, 8), ONE), kscal(Fr(1, 8), kmul(i_, kmul(sqrt3, sqrt5))))  # (7 + i sqrt15)/8
psi = kmul(kconj(phi0), rho)
for name, u in (('phi0', phi0), ('phi1', phi1), ('rho', rho), ('psi', psi)):
    assert kmul(u, kconj(u)) == ONE, name
assert abs(kval(phi1) - (6 ** 0.5 + 1j * 3 ** 0.5) / 3) < 1e-12
assert abs(kval(rho) - (7 + 1j * 15 ** 0.5) / 8) < 1e-12
assert abs(kval(psi) - cmath.exp(-1j * cmath.pi / 12) * (7 + 1j * 15 ** 0.5) / 8) < 1e-12


def kpow(x, n):
    if n < 0:
        x, n = kconj(x), -n
    out = ONE
    for _ in range(n):
        out = kmul(out, x)
    return out


assert kpow(phi0, 24) == ONE and kpow(phi0, 12) == kscal(-1, ONE)

# ------------------------------------------------------------------ point sets as integer arrays over denominator D
gens_clock = [kpow(phi0, a) for a in range(24)]
gens_M1 = [kmul(c, kpow(phi1, b)) for c in gens_clock for b in (-1, 0, 1)]
psi_K = psi


def flat(x):
    return x[0] + x[1]


allfr = [c for g in gens_M1 + [psi_K] for c in flat(g)]
D0 = 1
for c in allfr:
    D0 = lcm(D0, c.denominator)
log(f"denominator of generators and psi: {D0}")

# integer multiplication tensor for Q(zeta) and conjugation matrix
TENS = np.zeros((8, 8, 8), dtype=np.int64)
for a_ in range(8):
    for b_ in range(8):
        e = zred([Fr(0)] * (a_ + b_) + [Fr(1)])
        TENS[a_, b_] = [int(v) for v in e]
CONJ = np.zeros((8, 8), dtype=np.int64)       # conj(e_j) = sum_k CONJ[j, k] e_k
for j in range(8):
    CONJ[j] = [int(v) for v in zpow_elem(-j)]


def to_int(points, D):
    arr = np.array([[int(c * D) for c in flat(p)] for p in points], dtype=np.int64)
    assert all((c * D).denominator == 1 for p in points for c in flat(p))
    return arr


def zmul_arr(A, B):
    """rowwise product in Q(zeta) of integer arrays (n,8)"""
    out = np.zeros_like(A)
    for a_ in range(8):
        for b_ in range(8):
            t = TENS[a_, b_]
            nz = np.nonzero(t)[0]
            if len(nz) == 0:
                continue
            prod = A[:, a_] * B[:, b_]
            for k in nz:
                out[:, k] += t[k] * prod
    return out


def kmul_arr(X, Y):
    a0, a1 = X[:, :8], X[:, 8:]
    b0, b1 = Y[:, :8], Y[:, 8:]
    return np.hstack([zmul_arr(a0, b0) + 5 * zmul_arr(a1, b1), zmul_arr(a0, b1) + zmul_arr(a1, b0)])


def kconj_arr(X):
    return np.hstack([X[:, :8] @ CONJ, X[:, 8:] @ CONJ])


ZPOW = np.array([ZETA ** j for j in range(8)])


def floats(Xint, D):
    z = (Xint[:, :8] @ ZPOW + 5 ** 0.5 * (Xint[:, 8:] @ ZPOW)) / D
    return np.stack([z.real, z.imag], axis=1)


def norm_is(Xint, D, num, den):
    """exact test |z|^2 == num/den for rows of Xint (over denominator D)"""
    N = kmul_arr(Xint, kconj_arr(Xint))          # over D^2
    target = np.zeros(16, dtype=np.int64)
    assert (num * D * D) % den == 0
    target[0] = num * D * D // den
    return np.all(N == target, axis=1)


def uniq(A):
    return np.unique(A, axis=0)


def minkowski(A, B):
    S = (A[:, None, :] + B[None, :, :]).reshape(-1, A.shape[1])
    return uniq(S)


def rot_set(A, u, D):
    """multiply each row (over D) by the K element u (Fractions); result over D * Du"""
    Du = 1
    for c in flat(u):
        Du = lcm(Du, c.denominator)
    U = np.array([[int(c * Du) for c in flat(u)]], dtype=np.int64)
    R = kmul_arr(A, np.repeat(U, len(A), axis=0))
    return R, D * Du


D = D0
M1 = uniq(np.vstack([np.zeros((1, 16), dtype=np.int64), to_int(gens_M1, D)]))
M1c = uniq(np.vstack([np.zeros((1, 16), dtype=np.int64), to_int(gens_clock, D)]))
log(f"|M1| = {len(M1)}, |M1clock| = {len(M1c)}, D = {D}")
assert np.all(norm_is(M1[np.any(M1 != 0, axis=1)], D, 1, 1))


def clip1(A, D):
    xy = floats(A, D)
    r2 = (xy ** 2).sum(1)
    near = np.abs(r2 - 1) < 1e-9
    ex = norm_is(A[near], D, 1, 1)
    assert ex.all(), "near-unit radius that is not exactly 1"
    return A[(r2 < 1) | near]


M11 = minkowski(M1, M1)
M2 = clip1(M11, D)
M3 = minkowski(M2, M1)
P1 = minkowski(M2, M1c)
M2c = clip1(minkowski(M1c, M1c), D)
log(f"|M1+M1| = {len(M11)}, |M2| = {len(M2)}, |M3| = {len(M3)}, |P1 = M2+M1clock| = {len(P1)}, |M2c| = {len(M2c)}")

# phi0-invariance of the pieces (so psi B = rho B up to the rotation phi0^-1 of B onto itself)
for name, A in (('M3', M3), ('P1', P1), ('M1+M1', M11)):
    R, DR = rot_set(A, phi0, D)
    assert DR == D
    assert set(map(tuple, uniq(R).tolist())) == set(map(tuple, A.tolist())), name
log("M1+M1, M3 and P1 are invariant under multiplication by phi0 = zeta (so psi B = rho B)")


def all_pairs_edges(Xint, D, block=96):
    xy = floats(Xint, D)
    n = len(xy)
    x, y = xy[:, 0], xy[:, 1]
    I, J = [], []
    for s in range(0, n, block):
        e = min(n, s + block)
        dx = x[s:e, None] - x[None, :]
        dy = y[s:e, None] - y[None, :]
        d2 = dx * dx + dy * dy
        m = np.abs(d2 - 1.0) < 1e-6
        ii, jj = np.nonzero(m)
        ii = ii + s
        keep = jj > ii
        I.append(ii[keep])
        J.append(jj[keep])
    I = np.concatenate(I)
    J = np.concatenate(J)
    ncand = len(I)
    ok = np.zeros(ncand, dtype=bool)
    for s in range(0, ncand, 200000):
        dd = Xint[I[s:s + 200000]] - Xint[J[s:s + 200000]]
        ok[s:s + 200000] = norm_is(dd, D, 1, 1)
    return I[ok], J[ok], ncand, xy


def build_graph(name, A, B=None, rotB=None):
    """A u rotB*B (rotB a K element); returns integer array over common denominator, edges"""
    if B is None:
        X, DX = A, D
    else:
        R, DR = rot_set(B, rotB, D)
        f = DR // D
        X = uniq(np.vstack([A * f, R]))
        DX = DR
    I, J, ncand, xy = all_pairs_edges(X, DX)
    log(f"{name}: |V| = {len(X)}, candidates {ncand}, exact unit edges {len(I)}, near-miss non-edges {ncand - len(I)}")
    np.savez_compressed(f"{OUT}/{name}.npz", X=X, D=np.array([DX]), I=I, J=J, xy=xy)
    return X, DX, I, J, xy


if __name__ == '__main__':
    which = sys.argv[2].split(',') if len(sys.argv) > 2 else ['M11', 'P1', 'M3', 'G48', 'W']
    if 'M11' in which:
        build_graph('M1+M1', M11)
    if 'P1' in which:
        build_graph('P1', P1)
    if 'M3' in which:
        build_graph('M3', M3)
    if 'G48' in which:
        build_graph('M3_u_psiP1', M3, P1, psi)
    if 'W' in which:
        build_graph('W_M3_u_psiM3', M3, M3, psi)
    if 'M3M11' in which:
        build_graph('M3_u_psiM11', M3, M11, psi)
    if 'P1P1' in which:
        build_graph('P1_u_psiP1', P1, P1, psi)
    if 'P1M3' in which:
        build_graph('P1_u_psiM3', P1, M3, psi)
    if 'orbits' in which:
        for nm, A in (('M1+M1', M11), ('P1', P1), ('M3', M3)):
            log(f"{nm}: radius-2 points {int(norm_is(A, D, 4, 1).sum())}, |z|^2 = 2/3 points {int(norm_is(A, D, 2, 3).sum())}, "
                f"|z|^2 = 8/3 points {int(norm_is(A, D, 8, 3).sum())}")
    if 'small' in which:
        build_graph('M11_u_psiM11', M11, M11, psi)
        build_graph('M2cM1_u_psi', minkowski(M2c, M1), minkowski(M2c, M1), psi)
