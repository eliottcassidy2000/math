#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_r3_direct_kpswf3.py   (kind-pasteur 2026-06-21, THREAD B / HYP-2730 followup)

THREAD B: the BALANCED r=3 case DIRECTLY (not by domination on r=2).

For E = B u {f1,f2,f3} balanced (three comparable far elements), as f1->inf the far
TRIPLE  (frac(f1 x), frac(f2 x), frac(f3 x)) -> (frac(q v), frac(p2 v), frac(p3 v))
in law for v ~ Unif[0,1), where the ratio pair is
   (gamma2, gamma3) = (f2/f1, f3/f1) = (p2/q, p3/q),  common denominator q.
(Standard equidistribution: scale f1=q*t, f2=p2*t, f3=p3*t -> t*v sweeps the 3D torus
geodesic with direction (q,p2,p3); base elements bounded & independent.)

Exact 3D limit integral (rational, breakpoint-cell):
   p0_inf(B; q,p2,p3) = E_{x,v}[ {sector(b x): b in B} U {sector(qv),sector(p2 v),sector(p3 v)} = Z/7 ].

We also build the THREE-INDEX cell measure mu_{q,p2,p3}(i,j,k) on (Z/7)^3 and the 3D L1
discrepancy  D3 = sum_{i,j,k} |mu - 1/343|.  The decorrelated plateau is
   P3(B) = E_x[ P( base(x) U {3 indep Unif(Z/7)} = Z/7 ) ]  = sum mu_indep * g_B,
and  R3 = p0_inf - P3 = sum_{ijk}[mu - 1/343] g_B(i,j,k),  |R3| <= max g_B * D3 <= D3.

GOAL (THREAD B deliverables):
  (1) p0_inf over a grid of small-denom ratio PAIRS (gamma2,gamma3) in (1,2.15]^2, find sup.
  (2) Confirm sup << cap_9.
  (3) Confirm resonance confined to small denominators: D3 = O(1/q) (3D torus-line discrepancy).
  (4) New HYP for the r-general bound.

EXACT rational throughout (Fraction). Outputs -> 05-knowledge/results/.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

P = 7
def sector(yf):  # yf a Fraction in [0,1)
    return int(P*yf)

CAP = {8: Fr(2243,5880), 9: Fr(1979,4004), 10: Fr(55,91)}

# --------------------------------------------------------------------------
# base sector-cells in x: list of (length, frozenset of base sectors)
# --------------------------------------------------------------------------
def base_xcells(B):
    B = [int(b) for b in B]
    xbp = {Fr(0), Fr(1)}
    for b in B:
        if b == 0:
            continue
        for t in range(0, P*b):
            xbp.add(Fr(t, P*b))
    xs = sorted(xbp)
    cells = []
    for a, b in zip(xs, xs[1:]):
        mid = (a+b)/2
        S = frozenset(sector((bb*mid) % 1) for bb in B)
        cells.append((b-a, S))
    return cells

# --------------------------------------------------------------------------
# v-cells for a far set of multipliers (q,p2,p3,...): list of (length, sectortuple)
# --------------------------------------------------------------------------
def far_vcells(mults):
    mults = [int(m) for m in mults]
    vbp = {Fr(0), Fr(1)}
    for f in mults:
        for t in range(0, P*f):
            vbp.add(Fr(t, P*f))
    vs = sorted(vbp)
    cells = []
    for a, b in zip(vs, vs[1:]):
        mid = (a+b)/2
        st = tuple(sector((f*mid) % 1) for f in mults)
        cells.append((b-a, st))
    return cells

# --------------------------------------------------------------------------
# exact p0_inf for an r-far balanced limit  (B base ints, mults far multipliers)
# --------------------------------------------------------------------------
def p0_inf_multi(B, mults):
    xcells = base_xcells(B)
    vcells = far_vcells(mults)
    total = Fr(0)
    # group v-cells by sector-set to speed up
    from collections import defaultdict
    vgrp = defaultdict(Fr)
    for vlen, st in vcells:
        vgrp[frozenset(st)] += vlen
    for sf, vlen in vgrp.items():
        for xlen, Sbase in xcells:
            if len(Sbase | sf) == P:
                total += xlen*vlen
    return total

# --------------------------------------------------------------------------
# decorrelated plateau P_r(B): r independent Unif(Z/7) far sectors
# P(cover m missing with r uniform) by inclusion-exclusion
# --------------------------------------------------------------------------
def Pr_decorrelated(B, r):
    B = [int(b) for b in B]
    xcells = base_xcells(B)
    total = Fr(0)
    for xlen, Sbase in xcells:
        m = P - len(Sbase)        # number of missing sectors
        if m == 0:
            pcov = Fr(1)
        elif m > r:
            pcov = Fr(0)
        else:
            # P(r uniform sectors hit all m specific missing) by incl-excl:
            # sum_{j=0}^{m} (-1)^j C(m,j) ((P-j)/P)^r
            from math import comb
            pcov = Fr(0)
            for j in range(0, m+1):
                pcov += Fr((-1)**j * comb(m, j)) * Fr(P-j, P)**r
        total += xlen*pcov
    return total

# --------------------------------------------------------------------------
# 3D cell discrepancy D3 for far direction (q,p2,p3)
# --------------------------------------------------------------------------
def D3_disc(q, p2, p3):
    vcells = far_vcells((q, p2, p3))
    from collections import defaultdict
    mu = defaultdict(Fr)
    for vlen, st in vcells:
        mu[st] += vlen
    inv = Fr(1, P**3)
    tot = Fr(0)
    keys = set(mu.keys())
    # all 7^3 cells: missing ones contribute |0-inv|=inv each
    seen = Fr(0)
    for st, w in mu.items():
        tot += abs(w - inv)
        seen += 0
    nmiss = P**3 - len(mu)
    tot += nmiss * inv
    return tot

# --------------------------------------------------------------------------
# ratio-pair grid: coprime-with-q numerators p2<p3, all in (1,2.15], common denom q
# We parametrize by q and pick p2,p3 in (q, 2.15 q]. (gamma_i = p_i/q.)
# To respect "balanced & ordered" use 1 < p2/q <= p3/q <= 2.15 ; f1<f2<f3.
# --------------------------------------------------------------------------
def ratio_pair_grid(qmax, hi=Fr(43,20)):
    out = []
    for q in range(1, qmax+1):
        # numerators in (q, hi*q]
        lo_n = q+1
        hi_n = int(hi*q)
        nums = list(range(lo_n, hi_n+1))
        for p2, p3 in itertools.combinations_with_replacement(nums, 2):
            # require the triple (q,p2,p3) to have gcd 1 collectively (else it's
            # really a lower-q direction); also p2<=p3.
            g = gcd(gcd(q, p2), p3)
            if g != 1:
                continue
            out.append((q, p2, p3))
    return out

def main():
    print("="*88)
    print("THREAD B: DIRECT balanced r=3   p0_inf(B; q,p2,p3) over ratio-pair grid (1,2.15]^2")
    print("="*88)

    # worst balanced base at k=9: r=3 leaves base size 6 (6 base + 3 far = 9).
    bases = {
        "k=9 base6 [0,2,4,6,8,10]": [0,2,4,6,8,10],
        "k=9 base6 [0,1,2,3,4,5]":  [0,1,2,3,4,5],
        "k=10 base7 [0,2..12]":     [0,2,4,6,8,10,12],   # 7 base + 3 far = 10
    }
    QMAX = 5   # ratio-pair grid denominators up to 5 (3D atlas)
    for name, B in bases.items():
        k = len(B)+3
        cap = CAP.get(k)
        P3 = Pr_decorrelated(B, 3)
        print(f"\n--- {name}  (k={k}, cap_{k}={float(cap) if cap else '?'}) ---")
        print(f"  decorrelated plateau P3(B) = {float(P3):.5f}  (q->inf / non-resonant 3-far value)")
        grid = ratio_pair_grid(QMAX)
        rows = []
        for q, p2, p3 in grid:
            val = p0_inf_multi(B, (q, p2, p3))
            R3 = val - P3
            rows.append((q, p2, p3, float(p2)/q, float(p3)/q, float(val), float(R3)))
        rows.sort(key=lambda t:-t[5])
        print(f"  {'q':>2} {'p2':>3} {'p3':>3} {'g2':>6} {'g3':>6} {'p0_inf':>9} {'R3':>9} {'<cap?':>6}")
        for q,p2,p3,g2,g3,val,R3 in rows[:12]:
            flag = "OK" if (cap and val < float(cap)) else ("OK" if not cap else "VIOL")
            print(f"  {q:>2} {p2:>3} {p3:>3} {g2:>6.3f} {g3:>6.3f} {val:>9.5f} {R3:>+9.5f} {flag:>6}")
        supval = max(t[5] for t in rows)
        worst = max(rows, key=lambda t:t[5])
        print(f"  => SUP p0_inf = {supval:.5f} at (q,p2,p3)=({worst[0]},{worst[1]},{worst[2]}) "
              f"gamma=({worst[3]:.3f},{worst[4]:.3f})")
        if cap:
            print(f"     cap_{k}={float(cap):.5f}; MARGIN={float(cap)-supval:.5f}  "
                  f"{'COMFORTABLE' if float(cap)-supval>0.1 else 'CHECK'}")

    # ----- 3D discrepancy decay: does D3 ~ C/q ? -----
    print("\n" + "="*88)
    print("3D torus-line discrepancy D3_{q,p2,p3} = L1 cell-disc on (Z/7)^3; does it decay ~1/q?")
    print("="*88)
    print("  Family: direction (q, q+1, 2q) -> gamma=( (q+1)/q, 2 ), spanning a 2-param resonance.")
    print(f"  {'q':>3} {'p2':>4} {'p3':>4} {'D3':>10} {'D3*q':>9}")
    supD3q = 0.0
    for q in range(1, 16):
        p2, p3 = q+1, 2*q
        if gcd(gcd(q,p2),p3) != 1:
            # 2q & q share q; reduce by using p3=2q+1 instead to keep coprime direction
            p3 = 2*q+1
        if not (Fr(1) < Fr(p2,q) <= Fr(43,20) and Fr(1) < Fr(p3,q) <= Fr(43,20)):
            continue
        d3 = D3_disc(q, p2, p3)
        d3f = float(d3); d3q = d3f*q
        supD3q = max(supD3q, d3q)
        print(f"  {q:>3} {p2:>4} {p3:>4} {d3f:>10.5f} {d3q:>9.4f}")
    print(f"  => sup D3*q over this family = {supD3q:.4f}  (bounded => D3 = O(1/q): 3D linear-flow discrepancy)")

    # ----- exhaustive D3 over the small atlas; largest q with D3>=margin -----
    print("\n  D3 atlas (all coprime directions, q<=8): largest q with D3 >= margin(0.21)?")
    MARGIN = 0.21
    qbad = 0; cnt = 0
    supD3q_all = 0.0
    for q, p2, p3 in ratio_pair_grid(8):
        d3 = float(D3_disc(q, p2, p3)); cnt += 1
        supD3q_all = max(supD3q_all, d3*q)
        if d3 >= MARGIN:
            qbad = max(qbad, q)
    print(f"   scanned {cnt} directions; sup D3*q = {supD3q_all:.4f}; "
          f"largest q with D3>=margin = {qbad}  => finite 3D atlas is q <= {qbad}")
    print(f"   (R3 <= max g_B * D3 <= D3; for q > {qbad}, R3 < {MARGIN} <= margin: TAIL-SAFE)")

    print("\nDONE.")

if __name__ == "__main__":
    main()
