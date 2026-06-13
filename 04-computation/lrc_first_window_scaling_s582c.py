#!/usr/bin/env python3
r"""
lrc_first_window_scaling_s582c.py    oracle-2026-06-03-S582o (overnight cycle 4)

Attacks the remaining gap via the FIRST-WINDOW conjecture (S556o): every off-wall config is
lonely in (0, c/n) for an absolute c. If max over configs of (first-lonely-time * n) is
bounded, LRC@n localizes to a bounded structured search in (0, c/n).

For many primitive (n-1)-sets, find the SMALLEST t>0 with all ||v_i t|| >= 1/n (lonely), via
a fine grid, and record t0*n. Report the distribution and the max (the empirical c). Flag
configs whose first lonely time is large (near-wall) or absent in (0, C/n).
"""
from functools import reduce
from math import gcd
import random

def d0(p):
    p %= 1.0
    return min(p, 1 - p)

def first_lonely_scaled(S, n, C=6.0, G=400000):
    """smallest t>0 (grid over (0, C/n]) with min_i ||v_i t|| >= 1/n; return t*n or None."""
    th = 1.0 / n
    hi = C / n
    steps = int(G * C / n) if False else 60000
    for i in range(1, steps + 1):
        t = hi * i / steps
        if all(d0(v * t) >= th - 1e-12 for v in S):
            return t * n
    return None

def main():
    print("=" * 72)
    print("FIRST-WINDOW scaling: smallest lonely t0, as t0*n, over configs (S556o)")
    print("=" * 72)
    rnd = random.Random(5823)
    for n in (8, 10, 12, 14):
        k = n - 1
        AP = tuple(range(1, n))
        configs = [AP]
        while len(configs) < 150:
            S = tuple(sorted(rnd.sample(range(1, 7 * n), k)))
            if reduce(gcd, S) == 1:
                configs.append(S)
        vals = []
        none_in_window = 0
        ap_val = None
        for S in configs:
            v = first_lonely_scaled(S, n)
            if v is None:
                none_in_window += 1
            else:
                vals.append(v)
                if S == AP:
                    ap_val = v
        vals.sort()
        if vals:
            mx = max(vals); md = vals[len(vals)//2]
            frac_lt1 = sum(1 for x in vals if x < 1.0 - 1e-9) / len(vals)
            print(f"  n={n:2d}: configs={len(configs)}  first-lonely t0*n:  median={md:.3f}  "
                  f"max={mx:.3f}  frac(<1)={frac_lt1:.2f}  AP={ap_val}  none-in-(0,6/n)={none_in_window}")
    print("\n" + "=" * 72)
    print("READING")
    print("=" * 72)
    print("""  The AP (tight wall) is lonely first at t0*n = 1 (=t=1/n). If every off-wall config has
  first-lonely t0*n bounded by an absolute c (and none are 'none-in-window'), the LRC lonely
  time always lives in (0, c/n) -- a bounded window. max(t0*n) is the empirical c; a small,
  n-stable c supports localizing the proof to a finite structured search in the first window.
  'none-in-window' configs (if any) would be the hard near-wall cases needing a wider C.""")

if __name__ == "__main__":
    main()
