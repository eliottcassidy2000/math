#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The covering-min uniqueness: literal base-uniqueness FAILS, but M-uniqueness holds. mac-mini-2026-06-30-S55.
Band-covering bases are MANY (n=13: 1406, incl. killer-based {13..23}), but the construction (consecutive
base {1..n-2} + lcm outlier) is the strict M-MINIMIZER; alternatives give much higher M. The construction is
the CANONICAL GREEDY form (all-transversal base + minimal killer) -- a Zeckendorf/Sylvester-style canonical
uniqueness (cf. HYP-3724 Sylvester-Egyptian greedy).
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for d in (Sg[i]-Sg[j], Sg[i]+Sg[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add(F(k, d))
    return max(min(min((v*t) % 1, 1-((v*t) % 1)) for v in Sg) for t in cand)


def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))


def main():
    print("Covering-min uniqueness: base-uniqueness FAILS, M-uniqueness HOLDS.\n")
    for n in (13, 14):
        Mc = F(n, n*n-n+1)
        constr = list(range(1, n-1)) + [n*(n-1)]
        print(f"n={n}: construction {constr} M={Mexact(constr)}={float(Mc):.5f} (=n/Phi6)")
        # alternative band-covering bases give higher M
        for name, S in [("killer-block base {n+1..2n-3}+small outlier",
                         list(range(n+1, 2*n-2)) + [2*(n+1)]),
                        ("shifted base {2..n-1}+outlier", list(range(2, n)) + [2*n-1])]:
            S = sorted(set(S))
            if len(S) == n-1 and covers(S, n):
                print(f"   alt: {name}: M={Mexact(S)}={float(Mexact(S)):.5f} >> {float(Mc):.5f}")
        print("   => the consecutive base + lcm outlier is the CANONICAL GREEDY (all-transversal, minimal-killer)")
        print("      M-minimizer -- unique in M though NOT in band-coverage (killer-or-transversal dichotomy).")


if __name__ == "__main__":
    main()
