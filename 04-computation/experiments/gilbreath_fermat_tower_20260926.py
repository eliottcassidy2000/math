#!/usr/bin/env python3
"""gilbreath_fermat_tower_20260926.py -- the 0/2 sea of Gilbreath's triangle as the Frobenius tower of the
Fermat numbers (session gilbreath-fermat-platonic-20260926, opus, 2026-09-26).

Facts checked exactly:
 (1) Lucas: C(n, j) is odd iff the bits of j are a subset of the bits of n; hence row n of Pascal's triangle
     mod 2, read as a binary number, equals prod_(i in bits(n)) F_i with F_i = 2^(2^i) + 1 (Fermat numbers).
 (2) The 32 rows n < 32 are exactly the 32 products of distinct KNOWN Fermat primes (the odd orders of
     constructible regular polygons, Gauss-Wantzel); row 32 is F_5 = 641 * 6700417 (composite).
 (3) In the single-seed diagram (Gilbreath orientation new[i] = a[i] XOR a[i+1], seed with room on both sides)
     every inverted zero triangle has side 2^m - 1 (Mersenne), topped by a run of 2^m ones; the count of side
     2^m - 1 among top rows n < 2^K equals the sum over n < 2^K with exactly m trailing ones of
     2^(popcount(n) - m). The all-ones rows 2^m - 1 have value prod_(i<m) F_i = 2^(2^m) - 1.
 (4) Sea kernel theorem on the actual prime triangle (primes below P): wherever the cone of a cell lies in the
     0/2 region of row r, entry(r + t, i) = 2 * XOR_(j subset of t) b_r(i + j), b = a/2; checked from the first
     all-0/2 row (65) for t up to 2047, and inside the sea left of the frontier for r = 10.
 (5) The front-only survival count is a unit-triangular linear map: the diagonal cells met by a speed-1 front
     are L_s(b) = XOR_(j subset s) b(F-1-s+j), s = 0..F-2, an invertible F_2-linear image of the F-1 sea bits;
     hence exactly C(F-1, w) sea words give a front path of weight w (checked by enumeration for F <= 12).
 (6) Side probes: drift log_2(F_k) - 2 of the Conway maps (F_k n + 1)/2; leading column of the absolute
     difference triangle of the constructible sequence s_n = prod_(i in bits n) F_i.
Usage: python3 gilbreath_fermat_tower_20260926.py [P=200000]
"""
import math, sys, itertools
import numpy as np

FERMAT = [2 ** (2 ** i) + 1 for i in range(8)]
KNOWN = FERMAT[:5]


def fermat_product(n):
    v = 1
    for i in range(n.bit_length()):
        if n >> i & 1:
            v *= FERMAT[i]
    return v


def row_value(n):
    return sum(1 << j for j in range(n + 1) if (j & ~n) == 0)


def part1_2():
    print("== (1) Lucas + Fermat products ==")
    for n in range(0, 128):
        for j in range(n + 1):
            assert (math.comb(n, j) % 2 == 1) == ((j & ~n) == 0)
        assert row_value(n) == fermat_product(n)
    print("rows 0..127: C(n,j) odd iff j subset n; row value = prod F_i over bits of n: OK")
    print("rows 2^k (value F_k):", [(2 ** k, row_value(2 ** k)) for k in range(6)])
    print("rows 2^k-1 (value 2^(2^k)-1 = prod_(i<k) F_i):", [(2 ** k - 1, row_value(2 ** k - 1)) for k in range(1, 6)])
    print("== (2) constructible odd polygon orders ==")
    subset_products = sorted(set(math.prod(c) for r in range(6) for c in itertools.combinations(KNOWN, r)))
    rows = sorted(row_value(n) for n in range(32))
    assert subset_products == rows and len(rows) == 32
    print("rows 0..31 as binary numbers = the 32 products of distinct known Fermat primes:", rows[:12], "...", rows[-2:])
    assert row_value(32) == FERMAT[5] == 641 * 6700417
    print("row 32 = F_5 =", FERMAT[5], "= 641 * 6700417 (composite): the first non-constructible row")
    print("every row n with 32 <= n < 2^33 contains a Fermat factor F_i with 5 <= i <= 32, all known composite")


def single_seed_triangles(K):
    p = 2 ** K  # seed with a zero at column 0 to its left and zeros to its right
    a = np.zeros(2 ** (K + 1) + 2, dtype=np.uint8); a[p] = 1
    rows = [a.copy()]
    for r in range(2 ** (K + 1)):
        a = a[:-1] ^ a[1:]
        rows.append(a.copy())
    hist = {}
    for r in range(2 ** K):
        row = rows[r]
        n = len(row)
        i = 0
        while i < n:
            if row[i] == 1:
                j = i
                while j < n and row[j] == 1:
                    j += 1
                L = j - i
                assert i > 0 and j < n and row[i - 1] == 0 and row[j] == 0, (r, i, j)
                t = L - 1
                if t >= 1:
                    for s in range(1, t + 1):
                        seg = rows[r + s][i:i + t - s + 1]
                        assert len(seg) == t - s + 1 and not seg.any(), (r, i, t, s)
                        assert rows[r + s][i - 1] == 1 and rows[r + s][i + t - s + 1] == 1, (r, i, t, s)
                    hist[t] = hist.get(t, 0) + 1
                i = j
            else:
                i += 1
    return hist


def part3():
    print("== (3) single-seed diagram: zero triangles ==")
    for K in (5, 7, 9):
        hist = single_seed_triangles(K)
        sides = sorted(hist)
        assert all((t + 1) & t == 0 for t in sides), sides
        pred = {}
        for n in range(2 ** K):
            m = 0
            while n >> m & 1:
                m += 1
            if m >= 1:
                t = 2 ** m - 1
                pred[t] = pred.get(t, 0) + 2 ** (bin(n).count('1') - m)
        print(" top rows n < %d: sides %s counts %s predicted %s" % (2 ** K, sides, [hist[t] for t in sides], [pred.get(t, 0) for t in sides]))
        assert hist == pred, (hist, pred)
    print(" every zero triangle is edged by ones (vertical left edge, diagonal right edge), every side is a Mersenne number 2^m - 1, and the count formula sum_(n: m trailing ones) 2^(popcount n - m) is exact")


def primes_below(P):
    s = np.ones(P, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def kernel(w, t):
    """XOR_(j subset t) w[i + j] for all i with i + t < len(w)"""
    n = len(w) - t
    k = np.zeros(n, dtype=np.uint8)
    sub = t
    while True:
        k ^= w[sub:sub + n]
        if sub == 0:
            break
        sub = (sub - 1) & t
    return k


def part4(P):
    print("== (4) sea kernel theorem on the prime triangle, primes below %d ==" % P)
    row = primes_below(P)
    want = set(range(1, 66)) | {65 + t for t in list(range(1, 65)) + [127, 255, 511, 1023, 2047]} | set(range(10, 70))
    keep = {}
    frontier = {}
    first_all = None
    for r in range(1, 65 + 2047 + 1):
        row = np.abs(np.diff(row))
        assert row[0] == 1
        big = np.nonzero(row[1:] >= 4)[0]
        frontier[r] = (int(big[0]) + 1) if len(big) else None
        if first_all is None and len(big) == 0:
            first_all = r
        if r in want:
            keep[r] = row
    print(" first row whose entries after the leading 1 are all in {0,2}: %d" % first_all)
    assert first_all == 65 and all(frontier[r] is not None for r in range(1, 65))
    print(" frontier F(r) for r = 1, 2, 5, 10, 20, 30, 40, 50, 60, 64:", [frontier[r] for r in (1, 2, 5, 10, 20, 30, 40, 50, 60, 64)])
    r0 = 65
    w = (keep[r0][1:] // 2).astype(np.uint8)
    assert set(np.unique(keep[r0][1:]).tolist()) <= {0, 2}
    checked = 0
    for t in list(range(1, 65)) + [127, 255, 511, 1023, 2047]:
        act = (keep[r0 + t][1:] // 2).astype(np.uint8)
        k = kernel(w, t)
        assert len(act) == len(k), (len(act), len(k))
        assert np.array_equal(act, k), t
        checked += 1
    print(" rows 65+t equal the Lucas-kernel image of row 65 for t in 1..64 and t = 127, 255, 511, 1023, 2047 (%d checks); leading 1 throughout" % checked)
    r = 10; F = frontier[r]
    w = (keep[r][1:F] // 2).astype(np.uint8)
    assert set(np.unique(keep[r][1:F]).tolist()) <= {0, 2} and keep[r][F] >= 4
    tot_in = 0
    for t in range(1, F - 1):
        k = kernel(w, t)
        act = keep[r + t][1:1 + len(k)]
        assert np.array_equal(act, 2 * k), t
        tot_in += len(k)
    print(" row 10 (frontier F=%d, defect %d): all %d in-cone cells over t = 1..%d equal the kernel image of the sea word" % (F, int(keep[r][F]), tot_in, F - 2))
    return frontier


def part5():
    print("== (5) front path = unit-triangular image of the sea bits ==")
    for F in range(3, 13):
        counts = {}
        for bits in range(1 << (F - 1)):
            b = [(bits >> i) & 1 for i in range(F - 1)]
            wgt = 0
            for s in range(F - 1):
                col = F - 1 - s
                v = 0
                sub = s
                while True:
                    v ^= b[col - 1 + sub]
                    if sub == 0:
                        break
                    sub = (sub - 1) & s
                wgt += v
            counts[wgt] = counts.get(wgt, 0) + 1
        assert all(counts.get(w_, 0) == math.comb(F - 1, w_) for w_ in range(F)), (F, counts)
    print(" for F = 3..12 the number of sea words whose front path carries exactly w twos is C(F-1, w): the path map is a bijection of F_2^(F-1)")
    print(" hence, front alone, a defect of size 2j at distance F survives to the edge with probability exactly sum_(t<=j-2) C(F-1,t) 2^-(F-1)")


def part6():
    print("== (6) side probes ==")
    print(" Conway maps (F_k n + 1)/2, drift per odd step log_2 F_k - 2:", ["%d: %+.4f" % (p, math.log2(p) - 2) for p in KNOWN])
    s = [row_value(n) for n in range(64)]
    lead = []
    cur = s
    for r in range(40):
        cur = [abs(cur[i + 1] - cur[i]) for i in range(len(cur) - 1)]
        lead.append(cur[0])
    print(" constructible sequence s_n (n<64): leading column of its absolute difference triangle, rows 1..40:", lead)


def main():
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    part1_2(); part3(); part4(P); part5(); part6()


if __name__ == '__main__':
    main()
