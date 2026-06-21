#!/usr/bin/env python3
"""
ROUTE 5, part 7 -- the FUGACITY contrast: the deep structural reason the
two extremalities are NOT one theorem.

Both H and measS7 are signed/weighted level-sums over a subset complex:
  H       = sum_j  (+2)^j  alpha_j      (odd-cycle independence poly, fugacity +2)
  measS7  = sum_k  (-1)^k  MISS_k       (sector-miss inclusion-exclusion, fugacity -1)

CLAIM (this script tests):
  - The tournament AP (interval) advantage lives in HIGH levels; positive
    fugacity +2 AMPLIFIES high levels => AP wins at large p.
  - The LRC consec advantage lives in the LOW level (MISS_1); the alternating
    fugacity -1 DAMPS high levels => low-order dominance => consec wins.
  - So the SAME "AP extremal" surface, but driven by OPPOSITE ends of the
    level spectrum under OPPOSITE-sign fugacities. NOT one c_k mechanism.

TESTS:
  (T1) tournament: define H_x = sum_j x^j alpha_j (general fugacity). At p=7,11
       does Paley win for SMALL x and Interval for LARGE x? Find the crossover
       fugacity x*. (At x=2 it's H.)  Shows AP-win is a HIGH-fugacity effect.
  (T2) LRC: define G_x = sum_k x^k MISS_k (general fugacity; x=-1 is measS7).
       Does consec maximize/extremize for x near -1 but lose for other x?
       I.e. consec-optimality is a LOW-order (x->0) effect.

Author: mac-mini-2026-06-21-S20 (ROUTE 5)
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.path.insert(0, '04-computation')
from lrc14_route1_conductance_minimax_opus_0621 import (
    sec, breakpoints, is_full_residue, primitive, consec, measS7,
)
sys.stdout.reconfigure(line_buffering=True)
HALF = F(1, 14)


# ---------- tournament side ----------
def is_qr(a, p):
    return a % p != 0 and pow(a, (p - 1) // 2, p) == 1


def paley_set(p):
    return frozenset(j for j in range(1, p) if is_qr(j, p))


def interval_set(p):
    return frozenset(range(1, (p - 1) // 2 + 1))


def adj(p, S):
    Sset = set(S)
    return [[1 if (j != i and (j - i) % p in Sset) else 0 for j in range(p)]
            for i in range(p)]


def all_odd_cycles(A, p):
    nbrs = [[j for j in range(p) if A[i][j]] for i in range(p)]
    cyc = []
    for ell in range(3, p + 1, 2):
        for start in range(p):
            def dfs(cur, mask, path):
                if len(path) == ell:
                    if start in nbrs[cur]:
                        cyc.append(frozenset(path))
                    return
                for nx in nbrs[cur]:
                    if nx <= start or (mask & (1 << nx)):
                        continue
                    dfs(nx, mask | (1 << nx), path + [nx])
            dfs(start, 1 << start, [start])
    return cyc


def alpha_levels_from_cycles(cyc, p):
    masks = [sum(1 << v for v in s) for s in cyc]
    N = len(masks)
    res = defaultdict(int)
    def rec(i, used, depth):
        res[depth] += 1
        for j in range(i, N):
            if masks[j] & used:
                continue
            rec(j + 1, used | masks[j], depth + 1)
    rec(0, 0, 0)
    return dict(res)


def Hx(alpha, x):
    return sum(F(x) ** j * cnt for j, cnt in alpha.items())


# ---------- LRC side ----------
def occ_set(E, a, y):
    return set(sec(e, a, y) for e in E)


def miss_levels(E):
    exact = defaultdict(F)
    for a in range(1, 7):
        bps = breakpoints(E, a, HALF)
        for lo, hi in zip(bps, bps[1:]):
            if hi <= lo:
                continue
            mid = (lo + hi) / 2
            occ = occ_set(E, a, mid)
            empty = frozenset(s for s in range(1, 7) if s not in occ)
            exact[empty] += hi - lo
    miss = defaultdict(F)
    allsec = list(range(1, 7))
    for k in range(0, 7):
        tot = F(0)
        for S in itertools.combinations(allsec, k):
            Sset = frozenset(S)
            for E_set, v in exact.items():
                if Sset <= E_set:
                    tot += v
        miss[k] = tot
    return miss


def Gx(miss, x):
    return sum(F(x) ** k * miss[k] for k in range(0, 7))


def main():
    print("ROUTE 5 part 7: FUGACITY CONTRAST")
    print("=" * 70)

    # ---- T1 tournament fugacity sweep ----
    print("\n(T1) TOURNAMENT H_x = sum_j x^j alpha_j  (x=2 is H)")
    for p in [7, 11]:
        aP = alpha_levels_from_cycles(all_odd_cycles(adj(p, paley_set(p)), p), p)
        aI = alpha_levels_from_cycles(all_odd_cycles(adj(p, interval_set(p)), p), p)
        print(f"\n  p={p}: alpha(Paley)={dict(sorted(aP.items()))}")
        print(f"        alpha(Interval)={dict(sorted(aI.items()))}")
        print(f"  {'x':>5} {'H_x(P)':>18} {'H_x(I)':>18} {'winner':>8}")
        prev = None
        for xnum in range(1, 9):
            x = F(xnum)
            hp, hi = Hx(aP, x), Hx(aI, x)
            w = "PALEY" if hp > hi else ("INT" if hi > hp else "tie")
            print(f"  {xnum:>5} {str(hp):>18} {str(hi):>18} {w:>8}")
        # find crossover fugacity between Paley and Interval (real root of H_x(P)=H_x(I))
        # diff(x) = sum_j (aP_j - aI_j) x^j ; level1>0 (Paley), level>=2 <0 (Int)
        diffs = {j: aP.get(j, 0) - aI.get(j, 0) for j in set(aP) | set(aI)}
        print(f"  level diffs (aP-aI): {dict(sorted(diffs.items()))}")
        # numerically find x where diff=0 (besides x=0)
        def diff(xf):
            return sum(d * xf ** j for j, d in diffs.items())
        lo, hiX = 0.001, 100.0
        if diff(lo) > 0 and diff(hiX) < 0:
            for _ in range(80):
                mid = (lo + hiX) / 2
                if diff(mid) > 0:
                    lo = mid
                else:
                    hiX = mid
            print(f"  >>> crossover fugacity x* = {lo:.5f}  "
                  f"(Paley wins x<x*, Interval wins x>x*)")
        else:
            print(f"  (no sign change in (0,100): diff(0+)={diff(lo):.2f}, diff(100)={diff(hiX):.2e})")

    # ---- T2 LRC fugacity sweep ----
    print("\n(T2) LRC G_x = sum_k x^k MISS_k  (x=-1 is measS7)")
    shapes = {
        "consec": consec(8),
        "0,2..8": [0,2,3,4,5,6,7,8],
        "0..6,9": [0,1,2,3,4,5,6,9],
        "0..6,14": [0,1,2,3,4,5,6,14],
    }
    ml = {name: miss_levels(E) for name, E in shapes.items()}
    print(f"  {'x':>6} " + " ".join(f"{n:>12}" for n in shapes) + "   argmax")
    for xt in [F(-1), F(-1,2), F(0), F(1,2), F(1), F(2)]:
        vals = {n: Gx(ml[n], xt) for n in shapes}
        am = max(vals, key=lambda n: vals[n])
        print(f"  {str(xt):>6} " + " ".join(f"{float(vals[n]):>12.6f}" for n in shapes)
              + f"   {am}")
    print("\n  At x=-1 (measS7) consec is argmax. Check whether consec stays")
    print("  argmax for x>0 (positive fugacity) -- if NOT, consec-optimality is")
    print("  specific to the alternating (-1) inclusion-exclusion.")

    # explicit: at x=-1 is consec max? compare full sweep small set
    print("\n  measS7 (=G_{-1}) values:")
    for n in shapes:
        print(f"    {n:>10}: {float(Gx(ml[n], F(-1))):.6f}")


if __name__ == "__main__":
    main()
