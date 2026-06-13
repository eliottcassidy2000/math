#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S702
THE LOWER-HALF DELETION LAW for AP_n.

CONJECTURED THEOREM (S702):  for 1 <= a < n/2,
    M(AP_n \\ {a})  =  2 / D*(n,a),
    where  D*(n,a) = min { D >= n+a : gcd(D, a) = 1 }.

Special case gcd(a,n)=1:  D* = n+a, so M = 2/(n+a).
This completes S701's profile: upper half a>=n/2 gives M=1/a (proved, S701).

Fast exact M: every optimum of  max_t min_v ||v t||  sits at t=p/D with
D = b1+b2 <= 2*max(V) (two opposite-slope binders), gap = c/D.  So we search
all D in [2, 2*maxV], all residues p, in pure integer arithmetic, then reduce.
"""
from fractions import Fraction
from math import gcd
from itertools import combinations


def M_exact_int(V):
    """max_t min_{v in V} ||v t||, exact, via integer denominator search."""
    mx = max(V)
    best_num, best_den = 0, 1          # best gap = best_num/best_den
    best_t = None
    for D in range(2, 2 * mx + 1):
        for p in range(1, D):
            c = D  # min residue distance over V
            for v in V:
                r = (v * p) % D
                d = r if r <= D - r else D - r
                if d < c:
                    c = d
                    if c == 0:
                        break
            if c == 0:
                continue
            # gap c/D vs best_num/best_den
            if c * best_den > best_num * D:
                best_num, best_den = c, D
                best_t = Fraction(p, D)
    g = Fraction(best_num, best_den)
    return g, best_t


def Dstar(n, a):
    D = n + a
    while gcd(D, a) != 1:
        D += 1
    return D


if __name__ == "__main__":
    print("=" * 78)
    print("S702: LOWER-HALF DELETION LAW  M(AP_n\\{a}) = 2/D*,  "
          "D*=min{D>=n+a: gcd(D,a)=1}")
    print("=" * 78)
    fails = 0
    wit_fails = 0
    total = 0
    for n in range(4, 49):
        for a in range(1, n):
            if 2 * a >= n:                 # only strict lower half a<n/2
                continue
            R = [v for v in range(1, n) if v != a]
            M, t = M_exact_int(R)
            D = Dstar(n, a)
            pred = Fraction(2, D)
            total += 1
            if M != pred:
                fails += 1
                print(f"  FAIL n={n} a={a}: M={M} pred=2/{D}={pred}")
            # lower-bound witness t=p/D*, p ≡ a^{-1} mod D*
            p = pow(a, -1, D)
            tw = Fraction(p, D)
            gw = min(min((v * p) % D, D - (v * p) % D) for v in R)
            if Fraction(gw, D) != Fraction(2, D):
                wit_fails += 1
                print(f"  WITNESS FAIL n={n} a={a}: realized {gw}/{D} != 2/{D}")
    print(f"\n  law:     {total-fails}/{total} exact matches "
          f"({'ALL PASS' if fails==0 else f'{fails} FAILS'})  [n=4..48]")
    print(f"  witness: {total-wit_fails}/{total} lower-bound witnesses realize 2/D* "
          f"({'ALL PASS' if wit_fails==0 else f'{wit_fails} FAILS'})")

    print("\n" + "=" * 78)
    print("Minimal UPPER-BOUND witness subset  W <= AP_n\\{a}  with M(W)=2/D*")
    print("=" * 78)
    for n in range(6, 19):
        for a in range(1, n):
            if 2 * a >= n:
                continue
            D = Dstar(n, a)
            target = Fraction(2, D)
            R = [v for v in range(1, n) if v != a]
            found = None
            for k in range(2, 5):
                for W in combinations(R, k):
                    mw, _ = M_exact_int(list(W))
                    if mw == target:
                        found = W
                        break
                if found:
                    break
            g = gcd(a, n)
            tag = f"gcd={g}" if g > 1 else "coprime"
            print(f"  n={n:>2} a={a:>2} ({tag:>7}): D*={D:>2}  min witness W={found}")
