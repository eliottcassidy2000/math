#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_resonance_atlas_kps.py   (kind-pasteur 2026-06-21, HYP-2729)

DIRECT ATTACK on L7 = the SOLE remaining LRC(14) lemma (the balanced multi-cluster
correction; OPEN-Q-108). Per the kps-S23 ledger, everything else is PROVED/feasible.

STRUCTURE (HYP-2729): for E = B u {f1,f2}, balanced (gamma=f2/f1 in (1,2.15)), the
f1->inf limit p0_inf(B,gamma):
  * gamma=p/q coprime  -> far pair (frac(f1 x),frac(f2 x)) -> (frac(q v),frac(p v)),
                          the (q,p) torus geodesic (CORRELATED) -> RESONANT.
  * gamma irrational / large q -> (q,p) curve FILLS torus -> far pair INDEPENDENT
                          -> p0_inf = the decorrelated plateau P2(B) (R=0).
So R(p/q)=p0_inf(B,p/q)-P2(B) is nonzero only at small-denom p/q, decaying ~1/q.

This script computes EXACTLY (rational):
  (1) P2(B) = decorrelated 2-far plateau (far pair independent uniform on the torus).
  (2) p0_inf(B,p,q) for the resonance atlas { p/q in (1,2.15), q <= QMAX }.
  (3) R(p/q) and its decay; sup over the atlas vs cap_k.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

P = 7
def sector(yfrac):  # yfrac in [0,1) as Fraction
    return int(P*yfrac)

CAP = {8: Fr(2243,5880), 9: Fr(1979,4004), 10: Fr(55,91)}

# ---- exact 2D limit integral p0_inf(B, p, q) over (x, v) in [0,1)^2 -------------
def p0_inf(B, p, q):
    """E_{x,v}[ {sector(b x): b in B} U {sector(q v), sector(p v)} covers all 7 ].
       Exact rational via breakpoint cells in x (base) and v (far (q,p)-curve)."""
    B = [int(b) for b in B]
    xbp = {Fr(0), Fr(1)}
    for b in B:
        if b == 0: continue
        for t in range(0, P*b): xbp.add(Fr(t, P*b))
    vbp = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(0, P*f): vbp.add(Fr(t, P*f))
    xs = sorted(xbp); vs = sorted(vbp)
    # precompute base sector-set per x-cell
    xcells = []
    for a, b in zip(xs, xs[1:]):
        mid = (a+b)/2
        S = set()
        for bb in B:
            S.add(sector((bb*mid) % 1))
        xcells.append((b-a, frozenset(S)))
    total = Fr(0)
    for va, vb in zip(vs, vs[1:]):
        midv = (va+vb)/2
        sf = {sector((q*midv) % 1), sector((p*midv) % 1)}
        vlen = vb-va
        for xlen, Sbase in xcells:
            if len(Sbase | sf) == P:
                total += xlen*vlen
    return total

# ---- decorrelated plateau P2(B): far pair independent uniform on the 7x7 torus --
def P2_decorrelated(B):
    """E_x[ P_{w1,w2 ~ Unif(Z/7) indep}( base(x) U {w1,w2} = all 7 ) ].
       For each x: base covers Sbase; need 2 indep uniform sectors to cover the
       m=7-|Sbase| missing. P(cover m missing with 2 uniform) = 0 if m>2, else
       inclusion-exclusion over the missing set."""
    B = [int(b) for b in B]
    xbp = {Fr(0), Fr(1)}
    for b in B:
        if b == 0: continue
        for t in range(0, P*b): xbp.add(Fr(t, P*b))
    xs = sorted(xbp)
    total = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a+b)/2
        S = set(sector((bb*mid) % 1) for bb in B)
        m = P - len(S)            # missing sectors
        if m == 0:
            pcov = Fr(1)
        elif m == 1:
            # need at least one of 2 uniform hits the 1 missing: 1-(6/7)^2
            pcov = Fr(1) - Fr(6,7)**2
        elif m == 2:
            # 2 uniform sectors must hit both missing: only if {w1,w2}={the 2}: 2/49
            pcov = Fr(2, 49)
        else:
            pcov = Fr(0)
        total += (b-a)*pcov
    return total

def farey_window(lo, hi, qmax):
    """coprime p/q in (lo,hi], q<=qmax, p/q>1."""
    out = []
    for q in range(1, qmax+1):
        for p in range(q+1, int(hi*q)+1):
            if gcd(p, q) != 1: continue
            r = Fr(p, q)
            if lo < r <= hi:
                out.append((p, q, r))
    out.sort(key=lambda t: t[2])
    return out

def main():
    # worst balanced base per ledger: k=9 base [0,2,4,6,8,10,12] (7 base + 2 far = 9)
    bases = {
        "k=9 base[0,2..12]": [0,2,4,6,8,10,12],
        "k=8 base[0,2..10]": [0,2,4,6,8,10],
        "k=9 base[0..6]":    [0,1,2,3,4,5,6],
    }
    print("="*84)
    print("L7 RESONANCE ATLAS  (HYP-2729): balanced two-far p0_inf(B, p/q), gamma in (1,2.15)")
    print("="*84)
    for name, B in bases.items():
        k = len(B)+2
        cap = CAP.get(k, None)
        P2 = P2_decorrelated(B)
        print(f"\n--- {name}  (k={k}, cap_{k}={float(cap) if cap else '?'}) ---")
        print(f"  decorrelated plateau P2(B) = {float(P2):.5f}   (the q->inf / non-resonant value)")
        atlas = farey_window(Fr(1), Fr(43,20), 12)   # (1, 2.15], q<=12
        rows = []
        for p, q, r in atlas:
            val = p0_inf(B, p, q)
            R = val - P2
            rows.append((float(r), p, q, float(val), float(R)))
        rows.sort(key=lambda t:-t[3])
        print(f"  {'gamma':>7} {'p/q':>7} {'p0_inf':>9} {'R=corr':>9}  {'<cap?':>6}")
        for r,p,q,val,R in rows[:14]:
            flag = "OK" if (cap and val < float(cap)) else ("OK" if not cap else "VIOLATION")
            print(f"  {r:>7.4f} {f'{p}/{q}':>7} {val:>9.5f} {R:>+9.5f}  {flag:>6}")
        supval = max(v for _,_,_,v,_ in rows)
        print(f"  => sup p0_inf over atlas = {supval:.5f}; cap_{k}={float(cap):.5f}; "
              f"margin = {float(cap)-supval:.5f}  {'COMFORTABLE' if float(cap)-supval>0.1 else 'CHECK'}")

    # decay of R(p/q) in q at fixed gamma~2 (resonant family p/q near 2)
    print("\n" + "="*84)
    print("DECAY OF R(p/q) IN q  (does the resonance correction -> 0 ? the tail bound)")
    print("="*84)
    B = [0,2,4,6,8,10,12]; P2 = P2_decorrelated(B)
    print(f"  base k=9 [0,2..12], P2={float(P2):.5f}; ratios near 2 with growing q:")
    for q in range(1, 16):
        p = 2*q+1   # ratio (2q+1)/q -> 2+, coprime since gcd(2q+1,q)=1
        if not (Fr(1) < Fr(p,q) <= Fr(43,20)): continue
        val = p0_inf(B, p, q); R = val-P2
        print(f"   {p}/{q}={float(Fr(p,q)):.4f}:  p0_inf={float(val):.5f}  R={float(R):+.5f}")
    print("\nDONE.")

if __name__ == "__main__":
    main()
