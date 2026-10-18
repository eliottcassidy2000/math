#!/usr/bin/env python3
"""Independent q-series checks for HYP-9296 / note section 5 (audit2).
1. Schur index of 4d N=2 SU(2) with N_f = 4 (8 half-hypers in 2 x 8_v), unflavoured and so(8)-zero-weight part,
   from the contour integral  (1/2) CT_z (1-z^2)(1-z^-2) PE[-2q/(1-q) chi_adj(z)] PE[q^{1/2}/(1-q) chi_2(z) chi_8v(a)].
2. Universal V^k(so8) zero-weight dims (PE[q/(1-q) chi_28]) at weights 5, 6.
3. sl3: CT_z of PE[sum_{u not | n} q^n chi_8(z)]  (= Dhat(q^u,z)/Dhat(q,z)) for u = 2 and u = 4.
4. sl2: Dhat(q^3,z)/Dhat(q,z) vs PE[q/(1-q^2) chi_3] at low order (sanity)."""
from collections import defaultdict
import itertools

N = 6          # q-order
S = 2 * N      # order in s = q^{1/2}


def mul(A, B, maxs=S):
    out = defaultdict(int)
    for ka, va in A.items():
        for kb, vb in B.items():
            s = ka[0] + kb[0]
            if s > maxs:
                continue
            k = (s,) + tuple(x + y for x, y in zip(ka[1:], kb[1:]))
            out[k] += va * vb
    return {k: v for k, v in out.items() if v}


def geom(mono, nvars, maxs=S):
    """1/(1 - mono), mono = (s_exp, exps...) with s_exp > 0"""
    out = {(0,) + (0,) * nvars: 1}
    k = 1
    while mono[0] * k <= maxs:
        out[(mono[0] * k,) + tuple(e * k for e in mono[1:])] = 1
        k += 1
    return out


def one_minus(mono, nvars):
    return {(0,) + (0,) * nvars: 1, mono: -1}


# ---------------- 1. SU(2) N_f = 4 ----------------
def su2_nf4(zero_weight):
    # variables: (s, z, a)   [a only used for zero-weight via factorisation]
    # G(z) = CT_a prod_n 1/((1-s^{2n-1} z a)(1-s^{2n-1} z/a)(1-s^{2n-1} a/z)(1-s^{2n-1}/(z a)))
    if zero_weight:
        F = {(0, 0, 0): 1}
        for n in range(1, N + 1):
            e = 2 * n - 1
            if e > S:
                break
            for zz, aa in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                F = mul(F, geom((e, zz, aa), 2))
        G = {(k[0], k[1]): v for k, v in F.items() if k[2] == 0}
        H = {(0, 0): 1}
        for _ in range(4):
            H = mul(H, G)
    else:
        H = {(0, 0): 1}
        for n in range(1, N + 1):
            e = 2 * n - 1
            if e > S:
                break
            for zz in (1, -1):
                g = geom((e, zz), 1)
                for _ in range(8):
                    H = mul(H, g)
    # vector multiplet: prod_n (1-q^n)^2 (1-z^2 q^n)^2 (1-z^-2 q^n)^2 ; Haar (1/2)(1-z^2)(1-z^-2)
    V = {(0, 0): 1}
    for n in range(1, N + 1):
        for zz in (0, 2, -2):
            f = one_minus((2 * n, zz), 1)
            V = mul(V, f)
            V = mul(V, f)
    haar = mul(one_minus((0, 2), 1), one_minus((0, -2), 1))
    tot = mul(mul(H, V), haar)
    res = [0] * (N + 1)
    for (s, z), v in tot.items():
        if z == 0:
            assert s % 2 == 0 or v == 0
            if s % 2 == 0:
                res[s // 2] += v
    assert all(r % 2 == 0 for r in res)
    return [r // 2 for r in res]


print('SU(2) N_f=4 Schur index (unflavoured):', su2_nf4(False))
print('SU(2) N_f=4 Schur index, so(8) zero-weight part:', su2_nf4(True))

# ---------------- 2. universal V^k(so8) zero-weight dims ----------------
roots_so8 = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            r = [0] * 4
            r[i] = si; r[j] = sj
            roots_so8.append(tuple(r))
assert len(roots_so8) == 24


def universal_zero_weight(roots, rank, Nq):
    # PE[q/(1-q) (rank + sum_roots e^alpha)] constant term, as series in q (s variable used as q here)
    F = {(0,) + (0,) * len(roots[0]): 1}
    for n in range(1, Nq + 1):
        for r in roots:
            F = mul(F, geom((n,) + r, len(r), Nq), Nq)
        g = geom((n,) + (0,) * len(roots[0]), len(roots[0]), Nq)
        for _ in range(rank):
            F = mul(F, g, Nq)
    out = [0] * (Nq + 1)
    for k, v in F.items():
        if all(x == 0 for x in k[1:]):
            out[k[0]] += v
    return out


print('universal V^k(so8) zero-weight dims, weights 0..6:', universal_zero_weight(roots_so8, 4, 6))

# ---------------- 3. sl3 boundary products ----------------
# roots of sl3 in the basis (e1-e2, e2-e3) coordinates
roots_sl3 = [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]


def sl3_product_zero_weight(u, Nq):
    F = {(0, 0, 0): 1}
    for n in range(1, Nq + 1):
        if n % u == 0:
            continue
        for r in roots_sl3:
            F = mul(F, geom((n,) + r, 2, Nq), Nq)
        g = geom((n, 0, 0), 2, Nq)
        F = mul(F, g, Nq); F = mul(F, g, Nq)
    out = [0] * (Nq + 1)
    for k, v in F.items():
        if k[1] == 0 and k[2] == 0:
            out[k[0]] += v
    return out


print('sl3 u=2: CT_z Dhat(q^2,z)/Dhat(q,z), weights 0..8:', sl3_product_zero_weight(2, 8))
print('sl3 u=4: CT_z Dhat(q^4,z)/Dhat(q,z), weights 0..8:', sl3_product_zero_weight(4, 8))
print('universal V^k(sl3) zero-weight dims 0..8:', universal_zero_weight([(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)], 2, 8))
