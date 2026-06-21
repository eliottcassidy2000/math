#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: the GENUINE-WIDE finite-window check (III) -- generalized doublet
over (base, gap, M in [15, M*]), is_gw-filtered, confirm p0 < cap.

Completes the leg-C closure outline:
  genuine-wide E = generalized doublet B u {M,M+g} (HYP-2807), p0 = Phi_frozen(B,g)+g(M)/M.
  (I) frozen room Phi_frozen(B,g) < cap [tail M>=M* auto].
  (II) Tornheim-uniform R-tail => G=period-max+sup|R| ~ 3.7-7, M*=ceil(G/(cap-Phi)) ~ tens [HYP-2808].
  (III) finite window 15<=M<M*: p0(B u {M,M+g}) < cap EXACTLY.   <-- THIS SCRIPT.

CRITICAL (this-session hygiene): apply the is_gw filter. Dilated configs (e.g. even-AP base +
even doublet => all-even, reduces to single-far) are BINDING (THM-563's job), NOT genuine-wide;
they reach higher p0 but must be excluded here.

Checks all bounded bases B (size k-2, primitive), gaps g=1..4, M in [15, MSTAR] with
is_gw(B u {M,M+g}); reports the WORST genuine-wide p0 and cap - worst. MSTAR=60 (well past
M*~44; the maximizer sits at small M). k=10,11,12 (binding).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP
from lrc14_wide_branch_ridge_codex_s47 import primitive

MSTAR = 60


def reprim(E):
    E = tuple(sorted(set(int(x) for x in E)))
    g = reduce(gcd, E)
    return tuple(x // g for x in E) if g > 1 else E


def is_gw(E):
    E = tuple(sorted(set(E)))
    if 0 not in E or max(E) - min(E) <= 14 or not primitive(E):
        return False
    for e in E:
        sub = reprim(tuple(x for x in E if x != e))
        if len(sub) < 2 or max(sub) - min(sub) <= 14:
            return False
    return True


def bounded_bases(size):
    for S in combinations(range(1, 15), size - 1):
        yield (0,) + S


def main():
    print("=" * 78)
    print("GENUINE-WIDE finite-window check (III): generalized doublet, is_gw-filtered  0622")
    print(f"   all bounded bases (size k-2) x gap g=1..4 x M in [15,{MSTAR}]")
    print("=" * 78)
    for k in (10, 11, 12):
        cap = CAP[k]
        size = k - 2
        worst = F(0)
        worstE = None
        nfail = 0
        nchk = 0
        nbases = 0
        for B in bounded_bases(size):
            nbases += 1
            for g in (1, 2, 3, 4):
                for M in range(15, MSTAR + 1):
                    E = tuple(sorted(B + (M, M + g)))
                    if len(set(E)) != k:
                        continue
                    if not is_gw(E):
                        continue
                    nchk += 1
                    v = p0_fast(E)
                    if v >= cap:
                        nfail += 1
                    if v > worst:
                        worst, worstE = v, E
        print(f"\nk={k}  cap={cap}={float(cap):.6f}  bases={nbases}  gw-configs checked={nchk}")
        print(f"  WORST genuine-wide p0 = {float(worst):.6f}  at {worstE}")
        print(f"  cap - worst = {float(cap - worst):+.6f}   FAILS(>=cap): {nfail}")
        print(f"  => genuine-wide finite window {'PASSES (all < cap)' if nfail == 0 else 'FAILS'}")
    print("\n" + "=" * 78)
    print("If all PASS: with the Tornheim R-tail (II) + frozen room (I), the genuine-wide leg")
    print("CLOSES = p0 < cap for all genuine-wide E. (Binding/dilated = THM-563, separate.)")


if __name__ == "__main__":
    main()
