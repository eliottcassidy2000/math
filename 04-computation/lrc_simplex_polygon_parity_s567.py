#!/usr/bin/env python3
"""Tournaments as simplex (mesh) and polygon (dihedral); the regular/dihedral face
is odd-m only; the LRC even-n tightness is the observer-marked regular n-gon.
opus-2026-06-02-S567."""
from fractions import Fraction
from sympy import totient, divisors
import collections

def A000016(m):  # round/polygon (the 'outside')
    return sum(int(totient(d)) * 2**(m // d) for d in divisors(m) if d % 2 == 1) // (2 * m)

def regular_exists(m):                    # out-degree (m-1)/2 integer
    return m % 2 == 1

def ap_block_outdegrees(n):               # AP {1..n-1} at the tight time 1/n
    V = list(range(1, n)); t = Fraction(1, n); od = []
    for i in V:
        od.append(sum(1 for j in V if i != j and 0 < ((i - j) * t) % 1 < Fraction(1, 2)))
    return collections.Counter(od)

if __name__ == "__main__":
    print("m   regular/dihedral D_m exists (out-deg (m-1)/2 integer)?   round=A000016")
    for m in range(2, 14):
        print(f"{m:2d}  {str(regular_exists(m)):>5s}  {'D_'+str(m) if regular_exists(m) else '--':>6}"
              f"   {A000016(m) if regular_exists(m) else '--'}")
    print("\n=> regular/dihedral face exists EXACTLY at ODD m ('every other n').\n")
    print("LRC AP runner-block out-degrees at t=1/n (n-gon MINUS observer => near-regular):")
    for n in [7, 9, 11, 13, 14, 15]:
        m = n - 1
        print(f"  n={n:2d} (m={m:2d}, {'odd m / EVEN n -> polygon face' if m % 2 else 'even m / ODD n -> simplex'}):"
              f" out-degrees {dict(sorted(ap_block_outdegrees(n).items()))}")
