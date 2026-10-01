#!/usr/bin/env python3
"""procgen_localglobal_20260930_duality.py -- the local/global duality between LRC and Collatz gates.

  bad_classes(k, ell)        all residue speed vectors u in (F_ell^*)^k (up to sign and order) with NO
                             lonely time m/ell at level 1/(k+1); exact density beta_ell(k)
  min_relation_l1(u, ell, B) least ||m||_1 over m != 0 with m.u = 0 (mod ell), searched to norm B
  lonely_count_fourier(V,q)  the Poisson/relation-lattice identity
                             N_lone(V,q) = q * sum_{xi in Lambda_q(V)} prod_i ghat(xi_i)   (exact check)
  single_relation_blocks     a single relation xi of Lambda_q(V) forces 'no lonely m' by itself iff
                             ||xi||_1 = 1, i.e. q | v_i (coverage) -- checked against brute force
stdout only.
"""
import cmath
import itertools
import math

import numpy as np


def is_bad(u, ell, k=None):
    if k is None:
        k = len(u)
    ua = np.array(u, dtype=np.int64) % ell
    if np.any(ua == 0):
        return True
    m = np.arange(1, ell, dtype=np.int64)
    r = (m[:, None] * ua[None, :]) % ell
    d = np.minimum(r, ell - r)
    return not bool(((k + 1) * d >= ell).all(axis=1).any())


def _box(k, B):
    """all integer vectors of length k with l1 norm <= B (as an array)."""
    pts = [()]
    for _ in range(k):
        new = []
        for p in pts:
            s = sum(abs(x) for x in p)
            for x in range(-(B - s), B - s + 1):
                new.append(p + (x,))
        pts = new
    return np.array(pts, dtype=np.int64)


_BOX = {}


def min_relation_l1(u, ell, B):
    """least l1 norm of a nonzero integer m with m.u = 0 mod ell (None if > B)."""
    k = len(u)
    key = (k - 1, B)
    if key not in _BOX:
        _BOX[key] = _box(k - 1, B)
    R = _BOX[key]
    ua = np.array(u, dtype=np.int64) % ell
    inv = pow(int(ua[0]), -1, ell)
    rest = (R @ ua[1:]) % ell
    m1 = (-rest * inv) % ell
    m1c = np.where(m1 > ell // 2, m1 - ell, m1)
    norms = np.abs(m1c) + np.abs(R).sum(axis=1)
    zero = np.abs(R).sum(axis=1) == 0
    norms = np.where(zero & (m1c == 0), 10 ** 9, norms)  # exclude m = 0
    norms = np.where(zero & (m1c != 0), 10 ** 9, norms)
    best = int(norms.min())
    best = min(best, ell)  # m = (ell, 0, ..., 0)
    return best if best <= B else None


def bad_classes(k, ell):
    """canonical classes (sorted |u_i| in 1..(ell-1)/2) of bad vectors, with the exact density beta of
    bad vectors among all of (F_ell^*)^k."""
    half = (ell - 1) // 2
    bad = []
    nbad = 0
    for c in itertools.combinations_with_replacement(range(1, half + 1), k):
        if is_bad(c, ell, k):
            bad.append(c)
            # number of ordered signed vectors with this multiset of absolute values
            cnt = math.factorial(k)
            for v in set(c):
                cnt //= math.factorial(c.count(v))
            nbad += cnt * 2 ** k
    return bad, nbad / (ell - 1) ** k


def lonely_count_brute(V, q):
    k = len(V)
    n = 0
    for m in range(q):
        if all((k + 1) * min((m * v) % q, q - (m * v) % q) >= q for v in V):
            n += 1
    return n


def lonely_count_fourier(V, q):
    """q * sum over the relation group Lambda_q(V) = {xi in (Z/q)^k : xi.V = 0} of prod ghat(xi_i)."""
    k = len(V)
    g = [1 if (k + 1) * min(x, q - x) >= q else 0 for x in range(q)]
    ghat = [sum(g[x] * cmath.exp(-2j * math.pi * xi * x / q) for x in range(q)) / q for xi in range(q)]
    tot = 0
    for xi in itertools.product(range(q), repeat=k):
        if sum(a * b for a, b in zip(xi, V)) % q == 0:
            pr = 1
            for t in xi:
                pr *= ghat[t]
            tot += pr
    return q * tot


def single_relation_blocks(xi, k):
    """does the single relation xi (integers) force <xi, y> in Z to be impossible on the box
    [delta, 1-delta]^k, delta = 1/(k+1)?  (interval of <xi,y> avoids Z)"""
    from fractions import Fraction
    d = Fraction(1, k + 1)
    lo = sum(Fraction(x) * (d if x > 0 else 1 - d) for x in xi)
    hi = sum(Fraction(x) * (1 - d if x > 0 else d) for x in xi)
    return math.floor(hi) < math.ceil(lo)  # no integer in [lo, hi]


if __name__ == "__main__":
    # Poisson identity on small cases
    for V, q in [((1, 2, 3), 7), ((1, 3, 4), 11), ((2, 5, 7), 13), ((1, 2, 3, 4), 7)]:
        a = lonely_count_brute(V, q)
        b = lonely_count_fourier(V, q)
        assert abs(a - b) < 1e-8, (V, q, a, b)
    print("duality helpers OK")


def short_relation_rate(k, ell, B=4):
    """fraction (weighted by the number of ordered signed vectors) of ALL u in (F_ell^*)^k that carry a
    relation m.u = 0 mod ell with 0 < ||m||_1 <= B -- the base rate against which 'every bad vector has
    one' must be compared."""
    half = (ell - 1) // 2
    tot = hit = 0
    for c in itertools.combinations_with_replacement(range(1, half + 1), k):
        cnt = math.factorial(k)
        for v in set(c):
            cnt //= math.factorial(c.count(v))
        tot += cnt
        if min_relation_l1(c, ell, B) is not None:
            hit += cnt
    return hit / tot


TIGHT = {3: [(1, 2, 3)], 4: [(1, 2, 3, 4), (1, 3, 4, 7)], 5: [(1, 2, 3, 4, 5), (1, 3, 4, 5, 9)],
         6: [(1, 2, 3, 4, 5, 6)],
         7: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]}


def canon(u, ell):
    return tuple(sorted(min(x % ell, ell - x % ell) for x in u))


def on_tight_line(u, ell, k):
    """is u (mod ell, up to signs and order) a scalar multiple of the reduction of a known tight k-set?"""
    targets = {canon(T, ell) for T in TIGHT.get(k, [])}
    for x in u:
        lam = pow(int(x) % ell, -1, ell)
        if canon([lam * y for y in u], ell) in targets:
            return True
    return False


def bad_classes_fast(k, ell):
    """same as bad_classes (canonical classes + exact density), vectorised over the last coordinate."""
    half = (ell - 1) // 2
    m = np.arange(1, ell, dtype=np.int64)
    out, nb = [], 0
    for pre in itertools.combinations_with_replacement(range(1, half + 1), k - 1):
        last = np.arange(pre[-1], half + 1, dtype=np.int64)
        P = np.array(pre, dtype=np.int64)
        rp = (m[:, None] * P[None, :]) % ell
        okp = ((k + 1) * np.minimum(rp, ell - rp) >= ell).all(axis=1)
        rl = (m[:, None] * last[None, :]) % ell
        okl = (k + 1) * np.minimum(rl, ell - rl) >= ell
        good = (okp[:, None] & okl).any(axis=0)
        for x in last[~good]:
            c = pre + (int(x),)
            out.append(c)
            cnt = math.factorial(k)
            for v in set(c):
                cnt //= math.factorial(c.count(v))
            nb += cnt * 2 ** k
    return out, nb / (ell - 1) ** k


def tight_line_prediction(k, ell):
    """density of the union of the tight lines (each: k! 2^(k-1) (ell-1) ordered signed vectors)."""
    return len(TIGHT[k]) * math.factorial(k) * 2 ** (k - 1) * (ell - 1) / (ell - 1) ** k
