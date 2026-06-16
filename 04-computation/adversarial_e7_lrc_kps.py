#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adversarial_e7_lrc_kps.py

Sanity-check the number facts in the SHARED CONTEXT that the worker's task touches
(the verifier is asked to judge 7/21 facts harshly for genuine-vs-coincidence):

 E7 claims: rank 7, dim 133 = 7*19, 126 roots, 63 positive roots = 9*7,
            Coxeter h=18, exponents {1,5,7,9,11,13,17}, |W(E7)| = 2^10*3^4*5*7.
 Forbidden: 7 = Phi_3(2) = 2^2+2+1 = 2^3-1 = M_3 ; 21 = 3*7.
 LRC dictionary: s(t)=sin(pi t/7)/(pi t); 7-vanishing iff 7|t; sign pattern.

These are arithmetic facts I can check directly.
"""
import math
from fractions import Fraction


def main():
    print("=== E7 numerology checks ===")
    # dim 133 = 7*19
    print(f"  133 == 7*19 ? {133 == 7*19}")
    # 126 roots, 63 positive = 9*7
    print(f"  63 == 9*7 ? {63 == 9*7}; 126 == 2*63 ? {126 == 2*63}")
    # exponents of E7: 1,5,7,9,11,13,17 ; sum should = #positive roots = 63
    exps = [1,5,7,9,11,13,17]
    print(f"  E7 exponents {exps}: count={len(exps)} (=rank 7? {len(exps)==7}), "
          f"sum={sum(exps)} (=#pos roots 63? {sum(exps)==63})")
    print(f"  7 is an exponent of E7? {7 in exps}")
    # Coxeter number h = (#roots)/rank = 126/7 = 18
    print(f"  Coxeter h = 126/7 = {126//7} (claim 18? {126//7==18})")
    # max exponent = h-1 = 17
    print(f"  max exponent {max(exps)} == h-1 = {18-1}? {max(exps)==17}")
    # |W(E7)| = 2^10 * 3^4 * 5 * 7
    WE7 = 2**10 * 3**4 * 5 * 7
    print(f"  |W(E7)| = 2^10*3^4*5*7 = {WE7} = {2903040}? {WE7==2903040}")
    # known value of |W(E7)| is 2903040 (standard). 7 divides it.
    print(f"  7 | |W(E7)| ? {WE7 % 7 == 0}")

    print("\n=== Forbidden 7/21 arithmetic ===")
    Phi3_2 = 2**2 + 2 + 1
    print(f"  Phi_3(2) = 2^2+2+1 = {Phi3_2} (==7? {Phi3_2==7})")
    print(f"  2^3 - 1 = {2**3-1} (Mersenne M_3, ==7? {2**3-1==7})")
    print(f"  21 = 3*7 = {3*7}")
    print(f"  Fano plane PG(2,2): points = (2^3-1)/(2-1) = {(2**3-1)//(2-1)} (==7? {(2**3-1)//1==7})")

    print("\n=== LRC s(t)=sin(pi t/7)/(pi t) sign/vanishing pattern ===")
    print("  t :   s(t)         sign   7|t?")
    for t in range(1, 15):
        s = math.sin(math.pi * t / 7) / (math.pi * t)
        sgn = "0" if abs(s) < 1e-12 else ("+" if s > 0 else "-")
        print(f"  {t:2d}: {s: .6f}    {sgn:>3s}   {t%7==0}")
    print("  Claimed: + on {1..6}, 0 at 7, - on {8..13}, 0 at 14. 7-vanish iff 7|t.")

    print("\n=== JUDGEMENT NOTES ===")
    print("  All E7 numerology facts are STANDARD Lie-theory values (verifiable in any")
    print("  reference: Bourbaki / Carter). 7 appears genuinely (rank, exponent, |W| factor).")
    print("  These are TRUE facts but their RELEVANCE to the OCF forbidden-7 is an")
    print("  ANALOGY, not a derived bridge: nothing connects W(E7)'s factor-7 to H=7's")
    print("  non-realizability. 7=Phi_3(2) is the load-bearing fact for the gas; E7 is")
    print("  decorative coincidence unless a representation-theoretic map is exhibited.")


if __name__ == "__main__":
    main()
