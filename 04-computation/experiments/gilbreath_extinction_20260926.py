#!/usr/bin/env python3
"""gilbreath_extinction_20260926.py -- the extinction function of Gilbreath's automaton (session
gilbreath-fermat-platonic-20260926, opus, 2026-09-26).

p_d(F) = probability that a single defect of size d, at distance F from the leading 1 (columns 1..F-1 and
F+1..F+R uniformly random in {0,2}), ever destroys the leading 1 (some entry >= 4 reaches column 1). Exact
enumeration over all 2^(F-1+R) configurations for small F (the truncation after F+R rows misses only losses
caused by fronts emitted while the stationary copy is alive beyond row R, probability O(2^-R) times a survival
probability), Monte Carlo beyond. Compared with the front-only tail sum_(t<=j-2) C(F-1,t) 2^-(F-1) (d = 2j),
which counts the first front alone (proved exact for the first front by the unit-triangular path map).
In a zero sea the defect spreads as Pascal mod 2 in units of d (|d-0| = d, |d-d| = 0); twos absorb it.
Then the risk sum over the actual prime frontier (primes below P): sum over fresh fronts of p_d(F).
Usage: python3 gilbreath_extinction_20260926.py [P=200000]
"""
import math, sys
import numpy as np


def evolve_loss(row, rows):
    lost = np.zeros(row.shape[0], dtype=bool)
    r = row
    for s in range(rows):
        if r.shape[1] < 2:
            break
        r = np.abs(r[:, :-1].astype(np.int16) - r[:, 1:].astype(np.int16)).astype(np.uint8)
        lost |= (r[:, 0] != 1)
    return lost


def build(bits, d, F, R):
    N = bits.shape[0]
    row = np.empty((N, F + R + 1), dtype=np.uint8)
    row[:, 0] = 1; row[:, 1:F] = bits[:, :F - 1]; row[:, F] = d; row[:, F + 1:] = bits[:, F - 1:]
    return row


def exact(d, F, R):
    nb = F - 1 + R
    N = 1 << nb
    idx = np.arange(N, dtype=np.uint32)
    bits = (((idx[:, None] >> np.arange(nb, dtype=np.uint32)[None, :]) & 1) * 2).astype(np.uint8)
    return evolve_loss(build(bits, d, F, R), F + R).mean()


def mc(d, F, R, samples, rng):
    bits = (rng.integers(0, 2, size=(samples, F - 1 + R)) * 2).astype(np.uint8)
    return evolve_loss(build(bits, d, F, R), F + R).mean()


def front_only(d, F):
    j = d // 2
    return sum(math.comb(F - 1, t) for t in range(0, j - 1)) / 2 ** (F - 1)


def primes_below(P):
    s = np.ones(P, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def main():
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    rng = np.random.default_rng(20260926)
    print("== exact extinction probabilities p_d(F) (R right bits; truncation error O(2^-R)) vs front-only tail ==")
    print(" d   F   R   exact p_d(F)     front-only      ratio")
    table = {}
    for d in (4, 6, 8):
        for F in range(3, 13):
            R = 10 if F <= 11 else 9
            p = exact(d, F, R)
            fo = front_only(d, F)
            table[(d, F)] = p
            print(" %d  %2d  %2d   %.9f    %.9f    %.4f" % (d, F, R, p, fo, p / fo if fo else float('nan')))
    print("== wall test for size 4: minimum column ever holding a 4 versus F - z, z = trailing zeros of the sea left of the defect ==")
    for F in range(3, 13):
        R = 10 if F <= 11 else 9
        nb = F - 1 + R
        N = 1 << nb
        idx = np.arange(N, dtype=np.uint32)
        bits = (((idx[:, None] >> np.arange(nb, dtype=np.uint32)[None, :]) & 1) * 2).astype(np.uint8)
        left = bits[:, :F - 1]
        z = np.zeros(N, dtype=np.int64)
        alive = np.ones(N, dtype=bool)
        for c in range(F - 1, 0, -1):  # columns F-1 down to 1
            alive &= (left[:, c - 1] == 0)
            z += alive
        r = build(bits, 4, F, R)
        mincol = np.full(N, F, dtype=np.int64)
        for s in range(F + R):
            if r.shape[1] < 2:
                break
            r = np.abs(r[:, :-1].astype(np.int16) - r[:, 1:].astype(np.int16)).astype(np.uint8)
            has4 = (r[:, :F + 1] == 4)
            cols = np.where(has4, np.arange(r.shape[1])[None, :F + 1][:, :has4.shape[1]], F + 1)
            mincol = np.minimum(mincol, cols.min(1))
        ok = np.array_equal(mincol, F - z)
        print(" F=%2d R=%2d: min column of any 4 == F - z for all %d configurations: %s" % (F, R, N, ok))
        assert ok
    print(" WALL: a size-4 defect never passes the first sea 2 to its left, whatever lies to its right (all contexts, F <= 12)")
    print(" larger Monte Carlo for size 4 (4e6 samples, R = 16):", ["F=%d: %.3e vs 2^-(F-1) = %.3e" % (F, mc(4, F, 16, 4000000, rng), 2.0 ** (1 - F)) for F in (14, 16)])
    print(" stabilization in R (d=4, F=8):", ["R=%d: %.9f" % (R, exact(4, 8, R)) for R in (6, 8, 10, 12)])
    print(" stabilization in R (d=8, F=8):", ["R=%d: %.9f" % (R, exact(8, 8, R)) for R in (6, 8, 10, 12)])
    print("== Monte Carlo for larger F (400000 samples, R = 16) ==")
    print(" d   F   MC p_d(F)      front-only     ratio")
    for d in (4, 6, 8):
        for F in (12, 14, 16, 18, 20, 24, 28):
            p = mc(d, F, 16, 400000, rng)
            fo = front_only(d, F)
            print(" %d  %2d   %.3e      %.3e     %.3f" % (d, F, p, fo, p / fo if fo else float('nan')))
    print("== prime frontier, primes below %d ==" % P)
    row = primes_below(P)
    fronts = []  # (row, F, d) fresh fronts
    prevF = None
    risk_rows = 0.0
    nrows = 0
    for r in range(1, len(row) - 1):
        row = np.abs(np.diff(row))
        assert row[0] == 1
        big = np.nonzero(row[1:] >= 4)[0]
        if len(big) == 0:
            break
        nrows += 1
        F = int(big[0]) + 1; d = int(row[F])
        risk_rows += front_only(d, F)
        if prevF is None or F != prevF - 1:
            fronts.append((r, F, d))
        prevF = F
    print(" rows with a defect: %d (1..%d); fresh fronts (F(r) != F(r-1) - 1): %d" % (nrows, nrows, len(fronts)))
    print(" fresh fronts (row, distance F, size d):", fronts[:40])
    risk_fresh = sum(front_only(d, F) for _, F, d in fronts)
    risk_fresh_exact = 0.0
    for _, F, d in fronts:
        if (min(d, 8), F) in table and d <= 8:
            risk_fresh_exact += table[(d, F)]
        else:
            risk_fresh_exact += front_only(d, F) * 1.35
    print(" risk sum over fresh fronts, front-only tail: %.6f ; with exact/estimated re-emission: %.6f ; sum over all rows (each front counted once per row): %.6f" % (risk_fresh, risk_fresh_exact, risk_rows))
    print(" largest single terms:", sorted(((front_only(d, F), r, F, d) for r, F, d in fronts), reverse=True)[:6])


if __name__ == '__main__':
    main()
