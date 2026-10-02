#!/usr/bin/env python3
"""procgen_tbij_20261001_twist.py -- Theorem G (twisted Mallows-Sloane for every S_n-submodule U of F_2^E)
checked by two independent computations:
  (a) Burnside on the torsor 'tournaments mod U' by F_2 linear algebra per permutation (one permutation per
      cycle type, weighted by the class size);
  (b) direct orbit count of graphs X in U^perp having no odd automorphism (nauty geng + dreadnaut).
Also: the per-permutation identity  #Fix_{T/U}(g) = sum_{X in (U^perp)^g} sgn_X(g)  for small n."""
import itertools
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_tbij_20261001_lib as L  # noqa: E402
import procgen_tbij_20261001_nat as N  # noqa: E402


def perm_edge_map(g, n):
    idx = L.edge_index(n)
    return [idx[tuple(sorted((g[i], g[j])))] for (i, j) in L.edge_list(n)]


def apply_edges(emap, x):
    y = 0
    k = 0
    while x:
        if x & 1:
            y |= 1 << emap[k]
        x >>= 1
        k += 1
    return y


def cocycle(g, n):
    """pairs where g(T0) differs from T0 (T0: i -> j for i < j)"""
    T0 = L.from_flip(0, n)
    return L.flip_vector(L.apply_perm_t(T0, g))


def span_basis(vecs):
    piv = {}
    for v in vecs:
        w = v
        while w:
            hb = w.bit_length() - 1
            if hb in piv:
                w ^= piv[hb]
            else:
                piv[hb] = w
                break
    return list(piv.values())


def in_span(basis_piv, v):
    w = v
    while w:
        hb = w.bit_length() - 1
        if hb in basis_piv:
            w ^= basis_piv[hb]
        else:
            return False
    return True


def piv_dict(basis):
    piv = {}
    for v in basis:
        w = v
        while w:
            hb = w.bit_length() - 1
            if hb in piv:
                w ^= piv[hb]
            else:
                piv[hb] = w
                break
    return piv


def submodules(n):
    """named S_n-submodules U of V = F_2^E as spanning lists"""
    E = L.edge_list(n)
    m = len(E)
    J = (1 << m) - 1
    cuts = [L.cut_mask(1 << v, n) for v in range(n)]
    C = span_basis(cuts)
    # cycle space Z = C^perp: triangles {0,i,j} span it
    tri = []
    for i in range(1, n):
        for j in range(i + 1, n):
            tri.append((1 << L.eid(n, 0, i)) | (1 << L.eid(n, 0, j)) | (1 << L.eid(n, i, j)))
    Z = span_basis(tri)
    Cpiv = piv_dict(C)
    Zpiv = piv_dict(Z)
    CcapZ = [c for c in (L.cut_mask(W, n) for W in range(1 << n)) if in_span(Zpiv, c)]
    mods = {
        '0': [],
        'J': [J],
        'C': C,
        'C+J': span_basis(C + [J]),
        'Z': Z,
        'Z+J': span_basis(Z + [J]),
        'CcapZ': span_basis(CcapZ),
        'C+Z': span_basis(C + Z),
    }
    return mods


def perp_basis(Ub, m):
    """basis of U^perp in F_2^m"""
    # solve <u, x> = 0 for all u in Ub
    piv = {}
    rows = []
    for u in Ub:
        w = u
        for hb in sorted(piv, reverse=True):
            if (w >> hb) & 1:
                w ^= piv[hb]
        if w:
            hb = w.bit_length() - 1
            # reduce existing
            for k in list(piv):
                if (piv[k] >> hb) & 1:
                    piv[k] ^= w
            piv[hb] = w
    pivots = set(piv)
    free = [b for b in range(m) if b not in pivots]
    basis = []
    for f in free:
        x = 1 << f
        for hb, row in piv.items():
            # row has leading bit hb; condition <row, x> = 0 determines x_hb
            if bin(row & x).count('1') % 2:
                x |= 1 << hb
        basis.append(x)
    for u in Ub:
        for b in basis:
            assert bin(u & b).count('1') % 2 == 0
    return basis


def fix_count_torsor(g, n, Ub, Uperp):
    """#fixed points of g on (tournaments mod U):  #{x : (P_g + 1) x + c_g in U} / |U|"""
    m = len(L.edge_list(n))
    emap = perm_edge_map(g, n)
    c = cocycle(g, n)
    # condition: for each q in Uperp basis: <q, P_g x + x + c> = 0  <=>  <P_g^{-1} q + q, x> = <q, c>
    ginv = L.inverse(g)
    emap_inv = perm_edge_map(ginv, n)
    rows = []
    for q in Uperp:
        r = apply_edges(emap_inv, q) ^ q
        rhs = bin(q & c).count('1') % 2
        rows.append((r, rhs))
    # Gaussian elimination with rhs
    piv = {}
    for r, b in rows:
        w, bb = r, b
        while w:
            hb = w.bit_length() - 1
            if hb in piv:
                w ^= piv[hb][0]
                bb ^= piv[hb][1]
            else:
                piv[hb] = (w, bb)
                break
        if w == 0 and bb == 1:
            return 0
    rank = len(piv)
    sols = 2 ** (m - rank)
    dimU = len(Ub)
    assert sols % (2 ** dimU) == 0
    return sols // (2 ** dimU)


def torsor_orbits(n, Ub, Uperp):
    tot = Fraction(0)
    for mu in L.partitions(n):
        g = N.standard_perm(mu, n)
        tot += Fraction(fix_count_torsor(g, n, Ub, Uperp), L.zee(mu))
    assert tot.denominator == 1
    return int(tot)


def even_graph_orbits_in(n, Uperp, graphs=None):
    """count iso classes of even graphs lying in U^perp (U^perp is S_n-stable so membership is iso-invariant);
    also return counts split by edge parity"""
    piv = piv_dict(Uperp)
    if graphs is None:
        graphs = N.graph_orbits(n)
    cnt = 0
    for o in graphs:
        if not o.info['even']:
            continue
        x = L.mask_from_adj(o.rep)
        if in_span(piv, x):
            cnt += 1
    return cnt


def self_converse_tournaments(n, Ts=None):
    if Ts is None:
        Ts = L.gentourng(n)
    full = (1 << n) - 1
    conv = [tuple(full & ~t[i] & ~(1 << i) for i in range(n)) for t in Ts]
    a = L.labelg_batch([L.d6_encode(t) for t in Ts])
    b = L.labelg_batch([L.d6_encode(t) for t in conv])
    return sum(1 for x, y in zip(a, b) if x == y)
