#!/usr/bin/env python3
"""Artin-type correction factors for Collatz: the orbit visit measure against
counting measure (the audit's unexplained 0.497), effective sample sizes, the
arrival bias above the start range, the joint 2-adic/3-adic residue law along
orbits; and the seeds of the session decoded exactly (consecutive-sum identities
with pivots 2T_n and 4T_n, square triangular numbers as the Pell unit 3 + 2 sqrt 2,
Schur S(2) = 4, the Heule-Kullmann-Marek threshold 7825, Artin's 20/19).

Session: opus, collatz-poset-dag-20260927 (S18), 2026-09-27.
Run: python 04-computation/experiments/collatz_artin_corrections_20260927.py
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from fractions import Fraction

from sympy import primerange, factorint, n_order


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1


def U(m: int) -> tuple[int, int]:
    x = 3 * m + 1
    v = v2(x)
    return x >> v, v


# ----------------------------------------------------------------------------
# P1: orbit visit measure versus counting measure
# ----------------------------------------------------------------------------

def part1(N: int = 10 ** 6):
    print(f"== P1: Syracuse orbits of all odd n <= {N}: visit weights w(m) = number of starts visiting m ==")
    w = Counter()
    nxt = {}
    for n in range(1, N + 1, 2):
        m = n
        while True:
            w[m] += 1
            if m == 1:
                break
            if m in nxt:
                m2 = nxt[m]
            else:
                m2, _ = U(m)
                nxt[m] = m2
            m = m2
    print(f"   distinct odd values visited: {len(w)}; total visits: {sum(w.values())}; max weight {max(w.values())} at m = {max(w, key=w.get)}")

    def val(m):
        return v2(3 * m + 1)

    # per dyadic band: orbit-weighted and distinct-value persistence P(v(U m) = 1 | v(m) = 1), effective sample size
    print("   band [2^i, 2^(i+1)): R_w = weighted P(next v=1 | v=1), R_1 = distinct-value version, n_eff = (sum w)^2/sum w^2 over v=1 values, #distinct v=1 values")
    bands = defaultdict(lambda: [0, 0, 0, 0, 0.0])  # sum w (v=1), sum w (v=1 and next v=1), count (v=1), count (both), sum w^2
    for m, wt in w.items():
        if m == 1:
            continue
        if m % 4 != 3:
            continue
        i = m.bit_length() - 1
        b = bands[i]
        nx = nxt[m]
        both = (nx % 4 == 3)
        b[0] += wt
        b[1] += wt * both
        b[2] += 1
        b[3] += both
        b[4] += wt * wt
    for i in sorted(bands):
        b = bands[i]
        if b[2] < 200:
            continue
        neff = b[0] ** 2 / b[4]
        tag = "below N" if 2 ** (i + 1) <= N else ("straddles N" if 2 ** i <= N else "above N")
        print(f"   i={i:>2} ({tag:>11}): R_w = {b[1]/b[0]:.4f}  R_1 = {b[3]/b[2]:.4f}  n_eff = {neff:9.0f}  distinct = {b[2]:8d}  (2 s.e. of R_w ~ {1/math.sqrt(neff):.4f})")
    # aggregate above and below N
    for tag, cond in (("m <= N", lambda m: m <= N), ("m > N", lambda m: m > N)):
        sw = sb = c = cb = 0
        for m, wt in w.items():
            if m == 1 or m % 4 != 3 or not cond(m):
                continue
            both = (nxt[m] % 4 == 3)
            sw += wt
            sb += wt * both
            c += 1
            cb += both
        print(f"   {tag}: R_w = {sb/sw:.5f}, R_1 = {cb/c:.5f} over {c} distinct v=1 values")
    # orbit-weighted and distinct valuation law above N versus 2^-k
    lawW = Counter()
    law1 = Counter()
    for m, wt in w.items():
        if m <= N or m == 1:
            continue
        k = val(m)
        lawW[min(k, 8)] += wt
        law1[min(k, 8)] += 1
    tw = sum(lawW.values())
    t1 = sum(law1.values())
    print("   valuation law of visited m > N:  k: weighted / distinct / 2^-k")
    for k in range(1, 9):
        print(f"      {k}: {lawW[k]/tw:.4f} / {law1[k]/t1:.4f} / {2.0**-k:.4f}")
    driftW = sum(wt * (math.log2(3) - val(m)) for m, wt in w.items() if m > N) / tw
    drift1 = sum((math.log2(3) - val(m)) for m in w if m > N) / t1
    print(f"   mean step drift log2(3) - v on visited m > N: weighted {driftW:+.4f}, distinct {drift1:+.4f} (Haar {math.log2(3)-2:+.4f})")
    # joint residue law (m mod 8, m mod 9) of distinct visited values above N: mutual information in bits
    joint = Counter()
    for m in w:
        if m > N:
            joint[(m % 8, m % 9)] += 1
    tot = sum(joint.values())
    p8 = Counter()
    p9 = Counter()
    for (a, b), c in joint.items():
        p8[a] += c
        p9[b] += c
    mi = sum(c / tot * math.log2((c / tot) / ((p8[a] / tot) * (p9[b] / tot))) for (a, b), c in joint.items())
    print(f"   distinct visited m > N: marginal mod 8 (odd classes 1,3,5,7): {[round(p8[a]/tot, 4) for a in (1, 3, 5, 7)]}; mod 3 classes 0,1,2: "
          f"{[round(sum(p9[b] for b in range(9) if b % 3 == r)/tot, 4) for r in range(3)]}; I(m mod 8; m mod 9) = {mi:.5f} bits")
    # the arrival bias: residue mod 8 of m > N given the valuation of the step INTO m along the visiting orbits
    into = Counter()
    for m, nx in nxt.items():
        if nx > N and nx != 1:
            into[(val(m), nx % 8)] += 1
    for k in (1, 2, 3):
        row = [into[(k, r)] for r in (1, 3, 5, 7)]
        s = sum(row)
        if s:
            print(f"   values m > N entered by a step of valuation {k}: residues mod 8 (1,3,5,7) = {[round(x/s, 3) for x in row]}")
    return w, nxt


# ----------------------------------------------------------------------------
# P2: the seeds decoded
# ----------------------------------------------------------------------------

def T(n: int) -> int:
    return n * (n + 1) // 2


def part2():
    print("\n== P2: consecutive-sum identities and their pivots ==")
    for n in range(1, 41):
        c = n * n
        left = sum(range(c, c + n + 1))
        right = sum(range(c + n + 1, c + 2 * n + 1))
        assert left == right and c + n == 2 * T(n)
        c2 = n * (2 * n + 1)
        left2 = sum(k * k for k in range(c2, c2 + n + 1))
        right2 = sum(k * k for k in range(c2 + n + 1, c2 + 2 * n + 1))
        assert left2 == right2 and c2 + n == 4 * T(n) and c2 == T(2 * n)
    print("   n^2 + ... + (n^2+n) = (n^2+n+1) + ... + (n^2+2n): pivot n^2 + n = 2 T_n, checked n <= 40")
    print("   (T_(2n))^2 + ... + (4T_n)^2 = (4T_n+1)^2 + ... + (4T_n+n)^2: start T_(2n) = n(2n+1), pivot 4 T_n, checked n <= 40; n = 1 is 3^2 + 4^2 = 5^2")
    # cubes: no such family for small n
    nocube = []
    for n in range(1, 7):
        found = [c for c in range(1, 20000) if sum((c + k) ** 3 for k in range(n + 1)) == sum((c + n + k) ** 3 for k in range(1, n + 1))]
        nocube.append((n, found))
    print(f"   cubes: c with (c)^3 + ... + (c+n)^3 = (c+n+1)^3 + ... + (c+2n)^3, c < 20000: {nocube}")
    # square triangular numbers and the Pell unit
    sq = []
    x, y = 3, 2  # x^2 - 2 y^2 = 1
    while len(sq) < 8:
        n = (x - 1) // 2
        m = y // 2
        assert T(n) == m * m
        sq.append((n, m, T(n)))
        x, y = 3 * x + 4 * y, 2 * x + 3 * y
    print(f"   square triangular numbers T_n = m^2 from the Pell unit 3 + 2 sqrt 2 (THM-4505's law 3^2 - 2*2^2 = 1): {sq}")
    # Schur S(2) = 4: sum-free 2-colourings of [1, n]
    def schur_ok(n):
        from itertools import product
        good = []
        for col in product((0, 1), repeat=n):
            ok = True
            for a in range(1, n + 1):
                for b in range(a, n + 1):
                    c = a + b
                    if c <= n and col[a - 1] == col[b - 1] == col[c - 1]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                good.append(col)
        return good
    g4 = schur_ok(4)
    g5 = schur_ok(5)
    print(f"   Schur: sum-free 2-colourings of [1,4]: {len(g4)} ({[tuple(i+1 for i in range(4) if c[i] == c[0]) for c in g4]}), of [1,5]: {len(g5)} -> S(2) = 4; the colouring is squares {{1,4}} against {{2,3}}")
    # HKM: Pythagorean triples within [1, 7825]
    M = 7825
    in_triple = set()
    triples_with = 0
    sq = {k * k: k for k in range(1, M + 1)}
    for a in range(1, M + 1):
        for b in range(a, M + 1):
            c2 = a * a + b * b
            if c2 > M * M:
                break
            if c2 in sq:
                c = sq[c2]
                in_triple.update((a, b, c))
                if M in (a, b, c):
                    triples_with += 1
    print(f"   Pythagorean triples within [1, {M}]: numbers in no triple: {M - len(in_triple)}; triples containing {M}: {triples_with}; 7825 = {factorint(7825)}")
    # Artin: the density of primes with 2 as a primitive root, and 5 with the 20/19 correction
    A = 1.0
    for p in primerange(2, 10 ** 6):
        A *= 1 - 1 / (p * (p - 1))
    cnt2 = cnt5 = tot = 0
    for p in primerange(3, 3 * 10 ** 5):
        if p == 5:
            continue
        tot += 1
        if n_order(2, p) == p - 1:
            cnt2 += 1
        if n_order(5, p) == p - 1:
            cnt5 += 1
    print(f"   Artin's constant A = {A:.6f}; primes < 3*10^5 with primitive root 2: {cnt2/tot:.4f}; with primitive root 5: {cnt5/tot:.4f}; A*20/19 = {A*20/19:.4f} (the entanglement of p = 1 mod 5 with (5/p) = 1)")


if __name__ == "__main__":
    part1()
    part2()
    print("\nDONE")
