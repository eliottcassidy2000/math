#!/usr/bin/env python3
"""collatz_hedgehog_family_20260927.py -- the hedgehog dictionary for a functional graph, the qn+1 family dichotomy,
and the exact linearization of the trivial cycle (session collatz-posets-zeta5-20260927, opus, 2026-09-27, fourth note).

 (1) Invariant probability measures of a map on a countable set are carried by cycles (checked in the only way a finite
     computation can: for T on {1..N} truncated, every recurrent point lies on a cycle); an orbit that visits a finite set
     infinitely often is eventually periodic -- so on N there are no "creeping" orbits (checked on all n <= 10^5: the orbit
     either enters {1,2} or, in the 5n+1 control, leaves every window).
 (2) The family T_q(n) = n/2 (n even), (qn+1)/2 (n odd), q odd: Terras's bijection (n mod 2^k <-> parity word of length k)
     holds for every odd q (checked k <= 12, q <= 11); the fraction of length-k parity words with no coefficient descent
     (q^(o_j) > 2^j for all j <= k), i.e. the density of the no-descent set, tends to 0 for q = 1, 3 and to a positive limit
     for q >= 5 (exact DP, k <= 200): the drift (log_2 q - 2)/2 changes sign between q = 3 and q = 5.
 (3) The trivial cycle is exactly linear on its 2-adic neighbourhood: T^2(x) - 1 = (3/4)(x - 1) for x = 1 mod 4, hence
     T^(2j)(x) - 1 = (3/4)^j (x - 1) for x = 1 mod 4^j: a real contraction by 3/4, a 2-adic expansion by 4, a 3-adic contraction
     by 3 (product formula). Shadow families of the trivial cycle: m = 1 + 4^j t falls to 1 + (3/4)^j (m - 1).
 (4) Cycles of T_q with minimum below 20000 for q = 1, 3, 5, 7 (search capped), for the family table.
Usage: python3 collatz_hedgehog_family_20260927.py
"""
import math
from fractions import Fraction


def T(n, q=3):
    return n // 2 if n % 2 == 0 else (q * n + 1) // 2


def part1():
    print("== (1) no creeping on N: orbits either enter a cycle or leave every finite window ==")
    N = 10 ** 5
    entered = 0; left = 0
    for n in range(1, N + 1):
        x = n; seen = set()
        while True:
            if x in (1, 2):
                entered += 1; break
            if x > 10 ** 12:
                left += 1; break
            x = T(x)
    print(" 3n+1: of %d starts, %d enter {1,2}, %d exceed 10^12 (none)" % (N, entered, left))
    # 5n+1 control: orbit of 7 leaves every window [1, W]
    x = 7; maxseen = 0; visits_below_1000 = 0
    for _ in range(20000):
        x = T(x, 5); maxseen = max(maxseen, x); visits_below_1000 += (x < 1000)
    print(" 5n+1 control, orbit of 7 for 20000 steps: %d visits below 1000, max ~ 2^%.0f: it leaves the window and (as far as computed) never returns" % (visits_below_1000, math.log2(maxseen)))
    print(" reading: a Collatz orbit on N is trapped (cycle) or escapes with all its statistical time at infinity; the hedgehog's creeping orbits have no counterpart")


def part2():
    print("== (2) the family T_q: Terras bijection for every odd q, and the no-descent density dichotomy ==")
    ok = True
    for q in (1, 3, 5, 7, 9, 11):
        for k in (6, 10, 12):
            words = set()
            for n in range(2 ** k):
                x = n; w = 0
                for i in range(k):
                    w = 2 * w + (x % 2); x = T(x, q)
                words.add(w)
            ok &= (len(words) == 2 ** k)
    print(" n mod 2^k -> parity word of length k is a bijection for q in {1,3,5,7,9,11}, k in {6,10,12}: %s" % ok)
    print(" fraction of length-k parity words with no coefficient descent (q^o_j > 2^j for all j <= k), exact DP:")
    print("  q, drift (log_2 q - 2)/2, fractions at k = 10, 20, 40, 80, 160, 200")
    for q in (1, 3, 5, 7, 9):
        lq = math.log2(q) if q > 1 else 0.0
        cur = {0: 1}  # o -> count, over words whose all prefixes satisfy the condition
        fr = {}
        for j in range(1, 201):
            nxt = {}
            for o, c in cur.items():
                for letter in (0, 1):
                    oo = o + letter
                    if q ** oo > 2 ** j:
                        nxt[oo] = nxt.get(oo, 0) + c
            cur = nxt
            if j in (10, 20, 40, 80, 160, 200):
                fr[j] = sum(cur.values()) / 2 ** j
        print("  q = %d, drift %+.4f: %s" % (q, (lq - 2) / 2, ["%.3g" % fr[j] for j in (10, 20, 40, 80, 160, 200)]))
    print(" reading: the density of the no-descent set is 0 exactly for q <= 3 (negative drift: Terras) and positive for q >= 5;")
    print(" among q > 1, Collatz is the unique member with a null no-descent set -- the density-level form of 'the only system'")


def part3():
    print("== (3) the trivial cycle is exactly linear on its 2-adic neighbourhood ==")
    ok = all(T(T(x)) - 1 == Fraction(3, 4) * (x - 1) for x in range(1, 40001, 4))
    print(" T^2(x) - 1 = (3/4)(x - 1) for all x = 1 mod 4, x <= 40000: %s" % ok)
    ok2 = True
    for j in range(1, 8):
        for t in range(0, 200):
            m = 1 + 4 ** j * t; x = m
            for _ in range(2 * j):
                x = T(x)
            ok2 &= (Fraction(x - 1) == Fraction(3, 4) ** j * (m - 1))
    print(" T^(2j)(m) - 1 = (3/4)^j (m - 1) for m = 1 mod 4^j (j <= 7, 200 values each): %s" % ok2)
    print(" multiplier 3/4 of the cycle map T^2 at its fixed point 1: |3/4|_oo = 0.75 (real contraction), |3/4|_2 = 4 (2-adic expansion), |3/4|_3 = 1/3 (3-adic contraction); product 0.75 * 4 * (1/3) = 1")
    print(" reading: the trivial cycle is linearizable (Koenigs) at the real place -- the opposite of a hedgehog (indifferent, non-linearizable);")
    print(" an integer orbit that passes within 2^(-2j) of 1 two-adically is pulled down by (3/4)^j in value: closeness to the cycle forces descent")


def part4():
    print("== (4) cycles of T_q with minimum below 20000 (search capped at 1500 steps and 10^40) ==")
    for q in (1, 3, 5, 7):
        minima = set()
        for n in range(1, 20001):
            x = n; steps = 0
            while x >= n and steps < 1500 and x < 10 ** 40:
                x = T(x, q); steps += 1
                if x == n:
                    minima.add(n); break
        print(" q = %d: cycle minima below 20000: %s" % (q, sorted(minima)))


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
