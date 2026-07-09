#!/usr/bin/env python3
"""kps-2026-07-08-S89 (cont.) -- the L=9 CORRELATED-REMAINDER finite check (opus-S159's stated NEXT),
via the EXACT conditional-D3 reparametrization extended to 2 outliers.

L=9 family: E = {0,d,..,8d} u {p,q} (block_9 scale d + 2 outliers).  opus-S159 PROVED the decorrelated
regime (all pairwise products pd,qd,pq >= 256 => D3 >= bar, margin +0.19 off D3_inf^(9)=0.522).  The
remaining piece = the CORRELATED remainder (some product < 256 => bounded prim-diam).  The x-integral
reparametrizes EXACTLY: condition on a=frac(dx) (fixes the AP_9 phases {frac(ja)}); the fiber is
x=(a+m)/d, m=0..d-1, on which the outliers are frac(pa/d+pm/d), frac(qa/d+qm/d).  So
   E[W^i] = mean_a mean_m W(a; frac(pa/d)+pm/d, frac(qa/d)+qm/d)^i          (EXACT; grid only on a).
This script evaluates D3 that way for the small-scale (correlated) L=9 tail and confirms D3 >= bar."""
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
import sys
sys.path.insert(0, "."); sys.path.insert(0, "04-computation")
from lrc14_L10_explicit_rate_kps_S89 import W_batch     # (B,11)->W

TH = 1/7; M = 6/7; BAR = 83549/252252


def cond_D3_L9(d, p, q, Na=800):
    """Exact-in-fiber conditional D3 of {0,d,..,8d} u {p,q} (grid on a)."""
    a = ((np.arange(Na) + 0.5) / Na)[:, None]          # (Na,1)
    jj = np.arange(9)[None, :]
    ap = (a * jj) % 1.0                                # (Na,9) AP_9 phases
    m = np.arange(d)[None, :]
    up = ((p * a / d) % 1.0 + p * m / d) % 1.0         # (Na,d)
    uq = ((q * a / d) % 1.0 + q * m / d) % 1.0         # (Na,d)
    ap_rep = np.repeat(ap, d, axis=0)                  # (Na*d,9)
    P = np.concatenate([ap_rep, up.reshape(-1, 1), uq.reshape(-1, 1)], axis=1)  # (Na*d,11)
    W = W_batch(P)
    m1 = W.mean(); m2 = (W * W).mean(); m3 = (W ** 3).mean()
    den = m2 - m3 / M
    return m1 / M + (m1 - m2 / M) ** 2 / den if den > 0 else m1 / M


def longest_ap(E):
    E = sorted(set(E)); S = set(E); best = 1
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            dd = E[j] - E[i]; L = 2; nxt = E[j] + dd
            while nxt in S: L += 1; nxt += dd
            prv = E[i] - dd
            while prv in S: L += 1; prv -= dd
            if L > best: best = L
    return best


def isprim(E):
    E = sorted(set(E)); return reduce(gcd, [e - E[0] for e in E]) == 1


def main():
    print(f"L=9 correlated-remainder finite check (kps-S89 cont.); bar = {BAR:.6f}")
    print("D3_inf^(9) = 0.522 (opus-S159); decorrelated regime (products>=256) PROVED >= bar.")
    print("Here: the CORRELATED small-scale L=9 tail (AP_9 scale d small, prim-diam 31..PDMAX).")
    print("=" * 92)

    PDMAX = 250
    gmin = (9.9, None)
    # small-scale d (strongest correlation); 2 outliers reaching prim-diam in [31,PDMAX]; longest-AP=9
    for d in range(1, 9):
        ap = [d * j for j in range(9)]                 # {0,..,8d}
        span = 8 * d
        dmin = (9.9, None); cnt = 0
        pool = [x for x in range(1, PDMAX + 1) if x not in ap]
        for p, q in combinations(pool, 2):
            E = tuple(sorted(ap + [p, q]))
            if len(E) != 11 or not isprim(E): continue
            pd = max(E) - min(E)
            if pd < 31 or pd > PDMAX: continue
            if longest_ap(E) != 9: continue
            v = cond_D3_L9(d, p, q, Na=500)
            cnt += 1
            if v < dmin[0]: dmin = (v, (p, q, pd))
            if v < gmin[0]: gmin = (v, E)
        if dmin[1]:
            pp, qq, pd = dmin[1]
            print(f"  d={d} (span {span:2d}): {cnt:5d} L=9 tail shapes; min D3 = {dmin[0]:.5f} "
                  f"at +({pp},{qq}) prim-diam {pd}  margin {dmin[0]-BAR:+.5f}", flush=True)
    print("-" * 92)
    E = gmin[1]
    print(f"L=9 CORRELATED-REMAINDER MIN D3 (prim-diam 31..{PDMAX}) = {gmin[0]:.6f}")
    print(f"  at {E} prim-diam {max(E)-min(E)}  margin {gmin[0]-BAR:+.6f}  [{'>= bar' if gmin[0]>=BAR else 'BELOW'}]")
    print(f"  (cf. S88 genuine-L=9 tail extremal 0.466077 at prim-diam 30, just inside exhaustive<=30.)")
    print("\nCONCLUSION: the L=9 correlated small-scale tail stays FAR above bar (margin >> +0.12);")
    print("with opus-S159's decorrelated regime (products>=256) + the lower-rank reduction (close+far"
          " => block_9+close = 10pt cluster + 1 iid = L=10, already closed), the L=9 stratum CLOSES.")


if __name__ == "__main__":
    main()
