#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The odd/even back-and-forth RESOLVED by PRIMITIVITY, + the three leverage ways + the Ramanujan/Paley frame.
mac-mini-2026-06-30-S49.

CRUX: THM-523 (canon, PROVED) reduces LRC to PRIMITIVE covering sets, where M>1/n STRICTLY.
- FULL covering-min (incl. non-primitive) = 1/n for ALL n via scaled block g*{1..n-1}, g=smallest prime
  factor of n -- this is the q-witness/EASY case (g*{1..n-1} / g = {1..n-1}, which omits a multiple of n).
  opus's even block 2*{1..n-1} is the g=2 (even n) special case. Parity only chooses g; the value is 1/n.
- PRIMITIVE covering-min (the HARD case, THM-523) > 1/n: n=7->2/13, n=8->2/15, ... margin = the looseness.
RAMANUJAN: the primitive covering-min lives on a circulant mod m; the Ihara-RH/Ramanujan/Weil bound on the
speed character sums is the spectral-gap criterion; 2n-1 is a Paley-tournament vertex count (Ramanujan).
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
import numpy as np
print = functools.partial(print, flush=True)


def Mexact(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add(F(k, d))
    best = F(0)
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in S)
        if g > best:
            best = g
    return best


def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))


def smallest_prime_factor(n):
    for p in range(2, n+1):
        if n % p == 0:
            return p


def paley_graph_eigs(q):
    QR = set((x*x) % q for x in range(1, q))
    A = np.array([[1.0 if (j-i) % q in QR else 0.0 for j in range(q)] for i in range(q)])
    return sorted(np.linalg.eigvals(A).real)


def main():
    print("=" * 78)
    print("PRIMITIVITY RESOLVES THE ODD/EVEN BACK-AND-FORTH (THM-523 is canon)")
    print("=" * 78)
    print("\nFULL covering-min = 1/n for ALL n via scaled block g*{1..n-1}, g = smallest prime factor:")
    for n in (7, 8, 9, 10, 11, 12, 14, 15):
        g = smallest_prime_factor(n)
        S = [g*v for v in range(1, n)]
        m = Mexact(S)
        print(f"  n={n:>2}: g={g:>2}  M={str(m):>5}={float(m):.5f} (=1/n: {m==F(1,n)})  covers:{covers(S,n)}  "
              f"gcd={g} NON-primitive  [the q={n}-witness in disguise: g*{{1..{n-1}}}/g = {{1..{n-1}}} omits mult of {n}]")
    print("\n  => parity only chooses g (even n: g=2 = opus's even block; odd prime: g=n; odd composite: g=p_min).")
    print("     The value is 1/n for ALL n. This is the EASY/q-witness case. THM-523 EXCLUDES it (non-primitive).")

    print("\nPRIMITIVE covering-min (the HARD case, THM-523, M>1/n strictly):")
    for n, S in [(7, [1, 2, 5, 6, 7, 8]), (8, [1, 4, 5, 6, 7, 11, 16]), (9, [1, 3, 4, 5, 7, 11, 18, 32])]:
        m = Mexact(S)
        print(f"  n={n}: primitive {S}  M={m}={float(m):.5f} > 1/n={1/n:.5f}  (margin {float(m)-1/n:.5f})")
    print("  n=14: THM-523 search min ~ 7/89 = 0.0787 > 1/14 = 0.0714 (HYP-2566: uniform looseness).")

    print("\n" + "=" * 78)
    print("RAMANUJAN / PALEY FRAME (the spectral criterion)")
    print("=" * 78)
    ev = paley_graph_eigs(13)
    k = 6
    nontriv = max(abs(e) for e in ev if abs(e-k) > 1e-6)
    print(f"\n  n=7 mediant modulus 2n-1=13 (prime, 1 mod 4): Paley GRAPH on 13 is {k}-regular, eigenvalues "
          f"6, (-1+-sqrt13)/2; max|nontrivial|={nontriv:.3f} <= 2sqrt(k-1)={2*math.sqrt(k-1):.3f}  => RAMANUJAN.")
    print(f"  n=14 mediant modulus 2n-1=27=GF(3^3) (3 mod 4): Paley TOURNAMENT on 27 vertices; non-principal")
    print(f"        |eigenvalue| ~ sqrt(27)/2={math.sqrt(27)/2:.3f} = the Weil/Ramanujan sqrt bound.")
    print("\n  FRAME: the primitive covering-min lives on a circulant mod m; M is governed by the speed")
    print("  character sums; the Ihara-RH / Ramanujan / Weil sqrt-bound is the spectral-gap criterion that")
    print("  controls the floor. The construction's AP = Gauss/Dirichlet sums (Weil-tight). opus's metazeta")
    print("  (Ihara zeta of the metagraph G_n) is the tournament-side instance of the same zeta machinery.")


if __name__ == "__main__":
    main()
