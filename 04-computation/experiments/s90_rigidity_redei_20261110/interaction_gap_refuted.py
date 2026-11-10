#!/usr/bin/env python3
"""A tested and REFUTED strengthening of THM-5124 (recorded in HYP-9657 (C)).

For two vertices a, b that are both active over R = V - {a, b} put c = H(R), d_a = H(R+a) - c, d_b = H(R+b) - c and
mu = H(R+a+b) - c - d_a - d_b (the interaction, >= 0).  THM-5124 (1) classifies mu = 0.

Candidate "second-order gap":   mu = 0   or   2 mu >= max(d_a, d_b).
  * TRUE on every labelled tournament of orders 4-7 (pair {0, 1} over the full rest; exhaustive at orders 4-6 in
    this script; the order-7 line printed below is a constant recorded from a separate numba run and was confirmed
    by the independent audit on every pair), with equality on 14 profiles, among them (c, d_a, d_b, mu) =
    (1, 4, 2, 2), (1, 4, 4, 2), (5, 4, 4, 2), (11, 8, 8, 4);
  * FALSE at order 8: this script finds an explicit witness with (c, d_a, d_b, mu) = (15, 14, 14, 6) by a seeded
    random search and re-counts its four path numbers by brute force over permutations.
Run:  python3 interaction_gap_refuted.py > interaction_gap_refuted.out
"""
import itertools
import random
import sys


def hs_table(n, beats):
    size = 1 << n
    cnt = [[0] * n for _ in range(size)]
    for v in range(n):
        cnt[1 << v][v] = 1
    Hs = [0] * size
    Hs[0] = 1
    for mask in range(1, size):
        row = cnt[mask]
        tot = 0
        for v in range(n):
            c = row[v]
            if c:
                tot += c
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
        Hs[mask] = tot
    return Hs


def brute(vertices, beats):
    total = 0
    for perm in itertools.permutations(vertices):
        ok = True
        for i in range(len(perm) - 1):
            if not (beats[perm[i]] >> perm[i + 1]) & 1:
                ok = False
                break
        if ok:
            total += 1
    return total


def exhaustive_small(n):
    """min of mu / max(d) over all labelled tournaments of order n, pair {0, 1}, base = the rest."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    full = (1 << n) - 1
    best = None
    for bits in range(1 << len(pairs)):
        beats = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                beats[i] |= 1 << j
            else:
                beats[j] |= 1 << i
        Hs = hs_table(n, beats)
        R = full ^ 3
        c = Hs[R]
        da = Hs[R | 1] - c
        db = Hs[R | 2] - c
        mu = Hs[full] - c - da - db
        if da > 0 and db > 0 and mu > 0:
            key = (mu, max(da, db))
            if best is None or key[0] * best[0][1] < best[0][0] * key[1]:
                best = (key, (c, da, db, mu))
    return best


def main():
    ok = True
    for n in (4, 5, 6):
        best = exhaustive_small(n)
        if best is None:
            print(f"order {n}: no active interacting pair")
            continue
        (mu, hi), rec = best
        good = (2 * mu >= hi)
        ok &= good
        print(f"order {n} (exhaustive): min mu/max(d) = {mu}/{hi} at (c, d_a, d_b, mu) = {rec}; 2 mu >= max(d) everywhere:", good)
    print("order 7 (recorded constant from a separate exhaustive numba run; confirmed by the audit): min mu/max(d) = 2/4, e.g. at (1, 4, 4, 2)")
    random.seed(20261110)
    n = 8
    full = (1 << n) - 1
    witness = None
    trials = 0
    while witness is None and trials < 200000:
        trials += 1
        beats = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.85:
                    beats[i] |= 1 << j
                else:
                    beats[j] |= 1 << i
        Hs = hs_table(n, beats)
        for a in range(n):
            for b in range(a + 1, n):
                R = full ^ (1 << a) ^ (1 << b)
                c = Hs[R]
                da = Hs[R | (1 << a)] - c
                db = Hs[R | (1 << b)] - c
                mu = Hs[full] - c - da - db
                if da > 0 and db > 0 and mu > 0 and 2 * mu < max(da, db):
                    witness = (list(beats), a, b, c, da, db, mu)
                    break
            if witness:
                break
    if witness is None:
        print("no witness found")
        sys.exit(1)
    beats, a, b, c, da, db, mu = witness
    verts = list(range(n))
    R = [v for v in verts if v not in (a, b)]
    c2 = brute(R, beats)
    ha = brute(R + [a], beats)
    hb = brute(R + [b], beats)
    hab = brute(verts, beats)
    good = (c2 == c and ha - c2 == da and hb - c2 == db and hab - ha - hb + c2 == mu and 2 * mu < max(da, db))
    ok &= good
    print(f"order 8 witness after {trials} seeded trials: out-neighbour masks {beats}, a = {a}, b = {b}")
    print(f"   brute-force counts: H(R) = {c2}, H(R+a) = {ha}, H(R+b) = {hb}, H(R+a+b) = {hab}")
    print(f"   (c, d_a, d_b, mu) = ({c2}, {ha - c2}, {hb - c2}, {hab - ha - hb + c2}); 2 mu < max(d): candidate gap REFUTED:", good)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
