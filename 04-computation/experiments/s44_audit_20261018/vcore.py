#!/usr/bin/env python3
"""Independent Virasoro vacuum-module engine for the S44 audit (written from scratch; does NOT import vir_ss.py).

Conventions (deliberately different from vir_ss.py, which uses DESCENDING PBW order):
  * A PBW monomial is an ASCENDING tuple (k1 <= k2 <= ... <= kr), k_i >= 2, meaning L_{-k1} L_{-k2} ... L_{-kr}|0>.
  * The commutative symbol of a monomial is the same tuple read as x_{k1-2} ... x_{kr-2} in C[x_0, x_1, ...]
    (length r = polynomial degree; weight = sum k_i; Li degree = weight - 2r).
  * d on C[x]: x_j -> (j+1) x_{j+1}, i.e. part k -> (k-1) * part (k+1).
Exact arithmetic: integer structure constants; central charge enters only through [L_n, L_{-n}].
"""
from functools import lru_cache
from fractions import Fraction


def asc(t):
    return tuple(sorted(t))


@lru_cache(maxsize=None)
def neg_mono(m, mono):
    """L_{-m} (m >= 1) applied on the left of the ascending PBW monomial `mono`.  Returns tuple of (mono, int)."""
    if m == 1:
        if not mono:
            return ()
        k1, rest = mono[0], mono[1:]
        out = {}
        # L_{-1} L_{-k1} rest = L_{-k1} (L_{-1} rest) + (k1 - 1) L_{-(k1+1)} rest
        for mo, co in neg_mono(1, rest):
            # every part of mo is >= k1, so prepending keeps ascending order
            assert not mo or mo[0] >= k1
            key = (k1,) + mo
            out[key] = out.get(key, 0) + co
        for mo, co in neg_mono(k1 + 1, rest):
            out[mo] = out.get(mo, 0) + (k1 - 1) * co
        return tuple((k, v) for k, v in out.items() if v)
    if not mono or m <= mono[0]:
        return (((m,) + mono, 1),)
    k1, rest = mono[0], mono[1:]
    out = {}
    # m > k1:  L_{-m} L_{-k1} rest = L_{-k1} (L_{-m} rest) + (k1 - m) L_{-(m+k1)} rest
    for mo, co in neg_mono(m, rest):
        assert not mo or mo[0] >= k1
        key = (k1,) + mo
        out[key] = out.get(key, 0) + co
    for mo, co in neg_mono(m + k1, rest):
        out[mo] = out.get(mo, 0) + (k1 - m) * co
    return tuple((k, v) for k, v in out.items() if v)


@lru_cache(maxsize=None)
def pos_mono_cached(n, mono, c):
    return tuple(pos_mono(n, mono, c).items())


def pos_mono(n, mono, c):
    """L_n (n >= 0) applied to ascending monomial; c = central charge (Fraction or int mod-P residue handled by caller).
    Returns dict mono -> coefficient (Fraction/int)."""
    if n == 0:
        return {mono: sum(mono)} if True else {}
    if not mono:
        return {}
    k1, rest = mono[0], mono[1:]
    out = {}
    # L_n L_{-k1} rest = L_{-k1} (L_n rest) + (n + k1) L_{n-k1} rest + delta_{n,k1} c/12 (n^3 - n) rest
    for mo, co in pos_mono_cached(n, rest, c):
        for mo2, co2 in neg_mono(k1, mo):
            out[mo2] = out.get(mo2, 0) + co * co2
    j = n - k1
    if j < 0:
        for mo, co in neg_mono(-j, rest):
            out[mo] = out.get(mo, 0) + (n + k1) * co
    elif j == 0:
        out[rest] = out.get(rest, 0) + (n + k1) * sum(rest)
        out[rest] = out.get(rest, 0) + c * Fraction(n ** 3 - n, 12)
    else:
        for mo, co in pos_mono_cached(j, rest, c):
            out[mo] = out.get(mo, 0) + (n + k1) * co
    return {k: v for k, v in out.items() if v}


@lru_cache(maxsize=None)
def parts(D, minpart=2):
    """All ascending tuples of parts >= minpart summing to D."""
    if D == 0:
        return ((),)
    res = []
    for k in range(minpart, D + 1):
        for rest in parts(D - k, k):
            res.append((k,) + rest)
    return tuple(res)


@lru_cache(maxsize=None)
def basis_by_length(D):
    """Basis of V_D ordered by length DESCENDING (longest first), ties lexicographic."""
    return tuple(sorted(parts(D), key=lambda t: (-len(t), t)))


@lru_cache(maxsize=None)
def index_of(D):
    return {t: i for i, t in enumerate(basis_by_length(D))}


@lru_cache(maxsize=None)
def block(D, r):
    return tuple(sorted(t for t in parts(D) if len(t) == r))


@lru_cache(maxsize=None)
def block_index(D, r):
    return {t: i for i, t in enumerate(block(D, r))}


def central_charge(p, q):
    return 1 - Fraction(6 * (p - q) ** 2, p * q)


def apply_neg_vec(m, vec, P=None):
    out = {}
    for mono, co in vec.items():
        for mo, c2 in neg_mono(m, mono):
            out[mo] = out.get(mo, 0) + co * c2
    if P is not None:
        return {k: v % P for k, v in out.items() if v % P}
    return {k: v for k, v in out.items() if v}


def apply_pos_vec(n, vec, c):
    out = {}
    for mono, co in vec.items():
        for mo, c2 in pos_mono_cached(n, mono, c):
            out[mo] = out.get(mo, 0) + co * c2
    return {k: v for k, v in out.items() if v}


def dpoly(poly, P=None):
    """d on C[x] (commutative), monomials = ascending tuples of parts."""
    out = {}
    for mono, co in poly.items():
        for i, k in enumerate(mono):
            new = asc(mono[:i] + (k + 1,) + mono[i + 1:])
            out[new] = out.get(new, 0) + co * (k - 1)
    if P is not None:
        return {k: v % P for k, v in out.items() if v % P}
    return {k: v for k, v in out.items() if v}


def polymul_mono(poly, m):
    return {asc(k + m): v for k, v in poly.items()}


def by_len(vec):
    out = {}
    for mono, co in vec.items():
        out.setdefault(len(mono), {})[mono] = co
    return out
