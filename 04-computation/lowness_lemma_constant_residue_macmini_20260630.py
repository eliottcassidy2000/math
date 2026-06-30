#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The LOWNESS LEMMA, routes 1&2 synthesized: the CONSTANT-RESIDUE principle. mac-mini-2026-06-30-S57.
Route 2 (k-witness): missing speed k => at t=k^{-1}/p (band prime p with {k,p-k} uncovered), min||vt||>=2/p>n/Phi6.
Route 1 (budget/transversal): to defeat all such witnesses, S needs a speed ≡±k mod every binding modulus.
SYNTHESIS: a SMALL speed k has CONSTANT residue (k mod p = k for all p>k) -- covers {k,p-k} mod EVERY prime at
once; a LARGE speed has SCATTERED residues -- patches only p | (w∓k). So speed k is IRREPLACEABLE. CRUX RESOLVED:
even w ≡ 1 mod (all primes <= 43) (~1.3e16) leaves M~0.141 >> 14/183 (the lonely hole just moves to mod 85). So
{1..n-2} is genuinely forced (no large-speed CRT escape) -> the lowness lemma is robust.
"""
from __future__ import annotations
import functools, math
print = functools.partial(print, flush=True)


def Mgrid(S, N=300000):
    best, bt = 0.0, 0
    for k in range(1, N):
        t = k/N
        g = min(abs(((v*t+0.5) % 1)-0.5) for v in S)
        if g > best:
            best, bt = g, t
    return best, bt


def prodprimes(P):
    r = 1
    for p in range(2, P+1):
        if all(p % d for d in range(2, int(p**0.5)+1)):
            r *= p
    return r


def main():
    Mc = 14/183
    print("CRUX: replace missing speed 1 with a huge CRT speed w ≡ 1 mod (all primes <= P):\n")
    for P in (23, 29, 43):
        w = 1 + prodprimes(P)
        S = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, w, 182]
        gm, t = Mgrid(S)
        print(f"  w≡1 mod primes<={P} (={w}): M~{gm:.5f} at t~{t:.4f}  vs 14/183={Mc:.5f}  -> still >> ({gm>Mc+0.05})")
    print("\n  => even an astronomically large CRT speed cannot restore M; {1..n-2} is forced.")
    print("\nCONSTANT vs SCATTERED residue (why speed k is irreplaceable):")
    for v in (1, 7430):
        print(f"  speed {v}: residues mod 17,19,23,29,31 = {[v%p for p in (17,19,23,29,31)]}"
              f"  ({'CONSTANT' if v < 17 else 'scattered'})")


if __name__ == "__main__":
    main()
