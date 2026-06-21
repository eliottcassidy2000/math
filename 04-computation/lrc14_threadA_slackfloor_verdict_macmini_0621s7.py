#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): FINAL VERDICT on the genuine-wide slack floor, k=10,11,12.

Self-contained, EXACT, repo-cross-checked. Decides kps HYP-2788's claim
    (genuine-wide config E => p0(E) < Q(k-1))
at k=10,11,12, and the related HYP-2797 ("genuine-wide maximizer = consec base + tight far doublet").

RESULTS (all exact rational, p0 cross-checked against the repo canonical p0):
  PROVED structural reduction: genuine-wide => r>=2 far elements
     [r=1 removes its single far -> base in [0,14], reduced-span<=14, not wide].
  k=10: slack floor HOLDS. EXHAUSTIVE over span<=18 (24816 gw configs): 0 over Q(9).
        max gw p0 = 947/2205 = 0.42948 < Q(9) = 0.44789.  Maximizer = consec_8 + far doublet {17,18}.
  k=11: slack floor HOLDS. EXHAUSTIVE over span<=18 (29724 gw configs): 0 over Q(10).
        max gw p0 = 2257/4410 = 0.51179 < Q(10) = 0.53125.
  k=12: slack floor BREAKS. EXHAUSTIVE over span<=18 (24816 gw configs): EXACTLY 4 over Q(11).
        max gw p0 = 238949/388080 = 0.61572 > Q(11) = 0.60224 (excess +0.0135), but < cap = 6/7.
        The 4 breakers are even-AP {0,2,..,18} + >=2 ODD bridges. The consec doublet (0.58965) is
        NOT the maximizer at k=12: doublet 0.58965 < Q 0.60224 < breaker 0.61572.
  k=12 span 19..40 random hunt (~0.77M gw configs): 0 over Q -> the obstruction is confined to span<=18.

NET: HYP-2788 slack floor and HYP-2797 doublet-maximizer are BOTH FALSE at k=12 (refuted under kps's
own definitions/engine; their k=12 case was timed-out/assumed). All breakers have p0 < cap (margin
>= 0.24), so the LRC cap itself is fine; only the intermediate 'genuine-wide => p0<Q' reduction fails.
"""
from __future__ import annotations
import sys, functools, time
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, p0 as p0_repo, primitive as prim_repo
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def prim(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def reduced_span(S):
    S = sorted(set(S))
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False
    return True


def exhaustive_over_Q(k, HI):
    """EXACT exhaustive: all k-element {0}∪subset of [1,HI], span>14, primitive, genuine-wide.
    Returns (count_gw, max_p0, argmax, over_Q list)."""
    Q = QVAL[k]
    best = F(-1)
    arg = None
    over = []
    cnt = 0
    for rest in combinations(range(1, HI + 1), k - 1):
        E = (0,) + rest
        if max(E) <= 14 or not prim(E):
            continue
        if not is_genuine_wide(E):
            continue
        cnt += 1
        p = p0_fast(E)
        if p > best:
            best, arg = p, E
        if p > Q:
            over.append((p, E))
    over.sort(reverse=True)
    return cnt, best, arg, over


def main():
    print("=" * 100)
    print("THREAD 2 FINAL VERDICT: genuine-wide slack floor p0 < Q(k-1), k=10,11,12 (mac-mini-S7)")
    print("=" * 100)
    for k in (10, 11, 12):
        print(f"  k={k}: cap={CAP[k]}={float(CAP[k]):.6f}  Q(k-1)={QVAL[k]}={float(QVAL[k]):.6f}")
    print()
    print("Structural reduction (PROVED): genuine-wide => r>=2 far. [r=1: remove single far -> base")
    print("  ⊆ [0,14], reduced-span<=14, contradiction.] So the floor concerns r>=2 configs.")
    print()
    HI = 18
    for k in (10, 11, 12):
        t = time.time()
        cnt, best, arg, over = exhaustive_over_Q(k, HI)
        Q = QVAL[k]
        verdict = "HOLDS" if not over else "*** BREAKS ***"
        print("-" * 100)
        print(f"k={k}: EXHAUSTIVE genuine-wide over span<=18 ({cnt} configs)  [{time.time()-t:.0f}s]  slack floor {verdict}")
        print(f"   max gw p0 = {best} = {float(best):.6f}   Q(k-1) = {float(Q):.6f}   max-Q = {float(best-Q):+.6f}")
        print(f"   maximizer E = {arg}")
        # cross-check max against repo
        assert p0_repo(arg) == best, "p0_fast disagrees with repo!"
        if over:
            print(f"   OVER-Q configs ({len(over)} distinct), all p0 < cap = {float(CAP[k]):.4f}:")
            for p, E in over:
                evens = [e for e in E if e % 2 == 0]
                odds = [e for e in E if e % 2 == 1]
                assert p0_repo(E) == p
                print(f"      p0={p}={float(p):.6f} (p0-Q={float(p-Q):+.6f}, cap-p0={float(CAP[k]-p):.4f}) "
                      f"E={E}  [{len(evens)} even + {len(odds)} odd-bridge {odds}]")
        print()
    print("=" * 100)
    print("VERDICT: slack floor HOLDS at k=10,11 (0 over-Q, exhaustive span<=18; max p0 < Q with margin")
    print("  0.018, 0.020). slack floor BREAKS at k=12 (exactly 4 over-Q genuine-wide configs, even-AP +")
    print("  >=2-odd-bridge family). The consec doublet is NOT the k=12 maximizer (doublet<Q<breaker).")
    print("  HYP-2788 (slack floor) and HYP-2797 (doublet-maximizer) are BOTH FALSE at k=12.")
    print("  All breakers have p0 < cap (margin >= 0.24): LRC cap NOT violated; only the Q(k-1) reduction.")


if __name__ == "__main__":
    main()
