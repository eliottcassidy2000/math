#!/usr/bin/env python3
"""
lrc_n14_seven_impossibility_s552.py    oracle-2026-06-01-S552o

"THE 7 IMPOSSIBILITY" as the key to LRC@14.  n=14 = 2*7, n*=7 prime (S546). S524:
   LRC@14  <=>  the 7 mod-7 CRT classes are simultaneously safe
   (6 pairs {i,i+7} for residues 1..6, + the singleton {7} = the multiples of 7).
It is a coupon-collector over 7 classes; the open part is the 7-way correlation.
Per S551, the AP core is measure-zero and needs CONSTRUCTION, not measure.

THE CONSTRUCTIVE ATTACK -- the 7-gon-vertex windows. At t = j/7 (j=1..6):
  - every PAIR-class runner v = 7u+c (c != 0 mod 7) has ||v*(j/7)|| = ||cj/7|| >= 1/7
    = 2/14 (SAFE with margin); perturbing t = j/7 + s keeps it safe for |s| <= 1/(14 V).
  - every SINGLETON runner 7w has 7w*(j/7) = wj in Z (BLOCKED at the vertex), but
    ||7w*(j/7 + s)|| = ||7 w s|| -- so near the vertex the whole problem reduces to
    clearing ONLY the multiples-of-7 sub-system {w_k} in a small s-window.
So: a lonely time exists in the window near t=j/7  <=>  the singleton sub-system {w_k}
clears in that window. With r = #multiples-of-7, if r is small the windows clear easily.

THE 7 IMPOSSIBILITY (precise): for a counterexample, ALL 7 classes must block at every
t -- in particular ALL SIX windows near t=1/7..6/7 must fail to clear the singleton
sub-system, AND the non-window region must be covered too. We TEST whether the 6
windows + the sieve already force loneliness for every n=14 set (leaving only the AP
wall), i.e. whether the 7-gon-window construction is the key.
"""
from itertools import combinations
from functools import reduce
from math import gcd
import random

def frac(x): return x - int(x // 1)
def dist0(p):
    p = p % 1.0
    return min(p, 1 - p)

N = 14

def lonely_at(speeds, t):
    return all(dist0(s*t) >= 1.0/N - 1e-12 for s in speeds)

def lonely_anywhere(speeds, G=400000):
    for i in range(G):
        t = (i+0.5)/G
        if all(1.0/N < (s*t) % 1.0 < 1.0 - 1.0/N for s in speeds): return True, t
    return False, None

def lonely_in_7gon_windows(speeds, sub=4000):
    """search ONLY the windows near t=j/7 (j=1..6), each of half-width 1/(14*Vmax),
    for a lonely time -- the constructive 7-impossibility witness."""
    V = max(speeds); delta = 1.0/(14*V)
    for j in range(1, 7):
        c = j/7.0
        for k in range(sub):
            t = c - delta + 2*delta*(k+0.5)/sub
            if all(1.0/N < (s*t) % 1.0 < 1.0 - 1.0/N for s in speeds):
                return True, t, j
    return False, None, None

def crt_class_safe_measures(speeds, G=200000):
    """measure each mod-7 class is safe (all its runners >= 1/14 from 0)."""
    cls = {c: [s for s in speeds if s % 7 == c] for c in range(7)}
    safe = {c: 0 for c in range(7)}
    for i in range(G):
        t = (i+0.5)/G
        for c in range(7):
            if all(dist0(s*t) >= 1.0/N for s in cls[c]): safe[c] += 1
    return {c: safe[c]/G for c in cls if cls[c]}

def main():
    print("="*74)
    print("LRC@14 and THE 7 IMPOSSIBILITY: the 7 mod-7 CRT classes + the 7-gon windows")
    print("="*74)
    print("  (S524: LRC@14 <=> all 7 mod-7 classes simultaneously safe.)")
    print()
    print("(1) class-safe measures for the AP (initial segment) -- ~ (6/7)^2 = 0.735")
    ap = tuple(range(1, 14))
    m = crt_class_safe_measures(ap)
    for c in sorted(m):
        print(f"     class {c} (speeds {sorted([s for s in ap if s%7==c])}): safe {m[c]:.3f}")
    print()
    print("(2) THE 7-GON-WINDOW CONSTRUCTION: lonely time in some window near t=j/7?")
    print("    (vs lonely anywhere). r = #multiples of 7.")
    rnd = random.Random(14)
    tot = 0; win_ok = 0; any_ok = 0; fails = []
    sets = [ap]  # include AP (the wall/core)
    while len(sets) < 25:
        v = tuple(sorted(rnd.sample(range(2, 60), 13)))
        if reduce(gcd, v) == 1: sets.append(v)
    for v in sets:
        r = sum(1 for s in v if s % 7 == 0)
        aok, at = lonely_anywhere(v)
        wok, wt, wj = lonely_in_7gon_windows(v)
        tot += 1
        if aok: any_ok += 1
        if wok: win_ok += 1
        elif aok: fails.append((v, r))
        tag = "AP/core" if v == ap else ""
        print(f"     v={v if v!=ap else 'AP 1..13'}: r(mult7)={r:2d}  lonely_anywhere={aok}  "
              f"7gon-window={wok}{(' (j='+str(wj)+')') if wok else ''} {tag}")
    print()
    print(f"  => 7-gon-window construction found a lonely time for {win_ok}/{tot} sets; "
          f"lonely-anywhere {any_ok}/{tot}.")
    print(f"     sets lonely-but-window-missed: {fails[:6]}")
    print()
    print("="*74)
    print("ASSESSMENT of 'the 7 impossibility'")
    print("="*74)
    print("  LRC@14 = the 7 mod-7 classes never all-block. Near t=j/7 the 6 pair-classes are")
    print("  SAFE automatically (||cj/7||>=1/7); the whole problem collapses to clearing the")
    print("  SINGLETON {multiples of 7} in a small window. So the '7 impossibility' = the")
    print("  multiples-of-7 sub-system cannot block ALL SIX windows AND the bulk simultaneously.")
    print("  This is the constructive (S551) form of LRC@14: not a measure bound but explicit")
    print("  7-gon-vertex witnesses, with the residual = the small multiples-of-7 sub-LRC.")

if __name__ == "__main__":
    main()
