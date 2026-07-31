#!/usr/bin/env python3
"""Build the ILP-found fast block rules explicitly and verify them from
scratch: exact composition balance on the whole shell, exact fairness as a
polynomial identity, causality, and the deadline profile.

If these pass, stitching the blocks [m,2m) over m = 2^r gives an exactly
fair extractor with T(n) <= rho * n, rho = max_r rho(2^r).
"""
from fractions import Fraction
from itertools import product
from math import comb
import sys

from amm12592_block_profile_ilp_thm3006 import feasible_ilp


def build_rule(m, a, wit):
    """S_k subsets realizing the witness counts p_{k,i}."""
    S = []
    for k in range(m):
        chosen = set()
        for i in range(a[k] + 1):
            pool = sorted(y for y in product((0, 1), repeat=a[k]) if sum(y) == i)
            cnt = wit[(k, i)]
            assert 0 <= cnt <= len(pool) == comb(a[k], i)
            chosen.update(pool[:cnt])
        S.append(chosen)
    return S


def heads_rel(m, a, S, z):
    """Verdict on the relative tail z (branch b=0).  z != 0."""
    k = 0
    while z[k] == 0:
        k += 1                       # k = index of leading 1 (0-based)
    w = z[k + 1:]
    return tuple(w[:a[k]]) in S[k]


def verify(m, a, wit, do_poly=True):
    S = build_rule(m, a, wit)

    # 1. relative layers 1..m-1 bisected
    layer = [0] * (m + 1)
    for z in product((0, 1), repeat=m):
        if not any(z):
            continue
        if heads_rel(m, a, S, z):
            layer[sum(z)] += 1
    for j in range(1, m):
        assert 2 * layer[j] == comb(m, j), ("layer", m, j, layer[j], comb(m, j))

    # 2. full shell: composition classes + fairness + causality + deadline
    def head(word):
        b = word[0]
        z = tuple(x ^ b for x in word[m:])
        hd = heads_rel(m, a, S, z)
        return hd if b == 0 else not hd

    def crit(w):
        n = 1
        while n < len(w) and w[n] == w[0]:
            n += 1
        return n

    cls, by_stop = {}, {}
    for b in (0, 1):
        for tail in product((0, 1), repeat=m):
            if all(x == b for x in tail):
                continue
            word = (b,) * m + tail
            n = crit(word)
            k = n - m
            tau = m + k + 1 + a[k]
            assert tau == n + 1 + a[k]
            assert tau <= 2 * m, (m, n, tau)
            cls.setdefault(sum(word), []).append(word)
            by_stop.setdefault(word[:tau], set()).add(head(word))

    for kk, ws in cls.items():
        hh = sum(1 for w in ws if head(w))
        assert 2 * hh == len(ws), ("class", m, kk, hh, len(ws))
    for pre, vals in by_stop.items():
        assert len(vals) == 1, ("not causal", m, pre)

    if do_poly:
        for p in (Fraction(1, 3), Fraction(5, 9)):
            q = 1 - p
            tot = hd = Fraction(0)
            for ws in cls.values():
                for w in ws:
                    mass = p ** (2 * m - sum(w)) * q ** sum(w)
                    tot += mass
                    if head(w):
                        hd += mass
            assert 2 * hd == tot, ("fairness", m, p)

    T = [m + k + 1 + a[k] for k in range(m)]
    rho = max(Fraction(T[k], m + k) for k in range(m))
    return T, rho


if __name__ == "__main__":
    PROFILES = {
        4: [1, 1, 1, 0],
        8: [3, 4, 4, 4, 3, 2, 1, 0],
        16: [8, 8, 9, 9, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
    }
    for m, a in PROFILES.items():
        ok, wit = feasible_ilp(m, a)
        assert ok, ("ILP lost feasibility", m)
        T, rho = verify(m, a, wit, do_poly=(m <= 8))
        print(f"m={m:3d}  VERIFIED   rho = {rho} = {float(rho):.5f}")
        print(f"       a = {a}")
        print(f"       T = {T}")
        print(f"       ratios = {[f'{float(Fraction(T[k], m+k)):.3f}' for k in range(m)]}")
        print()
    print("all explicit rules verified: layers bisected, classes bisected,")
    print("causal, deadlines as claimed (and exact fairness for m <= 8)")
