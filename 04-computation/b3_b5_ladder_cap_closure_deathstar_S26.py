#!/usr/bin/env python3
"""death-star-2026-07-16-S26 (HYP-7065): THE B3/B5 BERNOULLI RUNGS + THE THM-905 CAP TABLE.

RUNG DERIVATIONS (this session, referee below):
  B5 (five runners, full-support k, singleton boxes [c_i/7,(c_i+1)/7)):
    c_B(k) = (1/(120 prod k_i)) * sum_{eps in {0,1}^5} (-1)^{sum eps} B5({sum k_i (c_i+eps_i)/7})
    [i-powers: (2 pi i j)^5 = i (2pi)^5 j^5;  sum_{j!=0} e(-j x)/j^5 = i (2pi)^5 B5({x})/120]
  B3 (support-3 term inside a 4/5-tuple lattice sum; omitted runners contribute |B|/7):
    term(k) = (1/(6 prod_3 k_i)) * sum_{eps in {0,1}^3} (-1)^{sum eps} B3({sum k_i(c_i+eps_i)/7}) * prod_omitted (|B|/7)
    [(2 pi i j)^3 = -i (2pi)^3 j^3;  sum_{j!=0} e(-jx)/j^3 = -i (2pi)^3 B3({x})/6]
  (Signs fixed by the numeric referee; B3(x) = x^3 - 3x^2/2 + x/2, B5(x) = x^5 - 5x^4/2 + 5x^3/3 - x/6.)

CAP TABLE (THM-905 (3)): limit values of P_V(S_a)+P_V(S_b), P_V(S_c), P_E({0} u S_ab)
over relation-lattice strata = 24/2401-type independent + B2 pair-limits + B3/B4(/B5)
lattice mass; verify < 1/12, 5/42, 40/441 with margins; dilation-invariance makes exact
scans cover the fully-pinned strata.
"""
from fractions import Fraction as Fr
from itertools import permutations, product as iproduct, combinations
from math import gcd
import sys, time

def B2(x): return x*x - x + Fr(1,6)
def B3(x): return x**3 - Fr(3,2)*x**2 + Fr(1,2)*x
def B4(x): return x**4 - 2*x**3 + x**2 - Fr(1,30)
def B5(x): return x**5 - Fr(5,2)*x**4 + Fr(5,3)*x**3 - Fr(1,6)*x
def fr(x): return x - (x.numerator // x.denominator)

def corner_sum(k, secs, Bfun, denom_c):
    """(1/(denom_c * prod k)) * sum_eps (-1)^{#} Bfun({sum k_i(c_i+eps)/7}) — sign referee'd."""
    n = len(k)
    prodk = 1
    for ki in k: prodk *= ki
    tot = Fr(0)
    for eps in iproduct([0,1], repeat=n):
        arg = Fr(0)
        for ki, ci, ei in zip(k, secs, eps):
            arg += Fr(ki * (ci + ei), 7)
        tot += (-1)**sum(eps) * Bfun(fr(arg))
    return tot / (denom_c * prodk)

def hit_measure(vs, boxes):
    bps = sorted(set(Fr(kk, 7*v) for v in vs for kk in range(7*v+1)))
    tot = Fr(0)
    for i in range(len(bps)-1):
        mid = (bps[i]+bps[i+1])/2
        if all(int((v*mid % 1)*7) in B for v, B in zip(vs, boxes)):
            tot += bps[i+1]-bps[i]
    return tot

def referee_rungs():
    print("RUNG REFEREE")
    # B5: planted 5-runner relation k=(1,1,1,-1,-2): v s.t. v1+v2+v3-v4-2v5=0
    k5 = (1,1,1,-1,-2)
    secs = (0,2,4,6,1)
    pred = corner_sum(k5, secs, B5, 120)
    print(f"  B5 k={k5} secs={secs}: closed = {pred} ≈ {float(pred):.6e}")
    for (a,b,c,e) in [(11,14,17,9),(23,30,35,19),(47,62,71,39)]:
        d5 = a+b+c-e
        if d5 % 2: continue
        v5 = d5//2
        v = [a,b,c,e,v5]
        if len(set(v))<5 or min(v)<1: continue
        # remainder after independent+pairs (pairs ~ 1/49-limits, no pair relations if generic)
        full = hit_measure(v, [[s] for s in secs])
        indep = Fr(1,7**5)
        p2 = Fr(0)
        for i in range(5):
            for j in range(i+1,5):
                pij = hit_measure([v[i],v[j]], [[secs[i]],[secs[j]]])
                p2 += (pij - Fr(1,49)) * Fr(1,7**3)
        R = full - indep - p2
        print(f"    v={v}: R = {float(R):.6e}  R/closed = {float(R/pred) if pred else float('nan'):.3f}")
        sys.stdout.flush()
    # B3 support-3 inside a 4-tuple: k=(1,1,-2,0): v1+v2=2v3, v4 free
    k3 = (1,1,-2)
    secs4 = (0,2,4,6)
    pred3 = corner_sum(k3, secs4[:3], B3, 6) * Fr(1,7)
    print(f"  B3 k={k3}+free secs={secs4}: closed = {pred3} ≈ {float(pred3):.6e}")
    for (a,b,w) in [(9,13,15),(19,27,23),(37,55,41),(75,109,83)]:
        c = (a+b)//2
        if (a+b)%2: continue
        v = [a,b,c,w]
        if len(set(v))<4: continue
        full = hit_measure(v, [[s] for s in secs4])
        indep = Fr(1,7**4)
        p2 = Fr(0)
        for i in range(4):
            for j in range(i+1,4):
                pij = hit_measure([v[i],v[j]], [[secs4[i]],[secs4[j]]])
                p2 += (pij - Fr(1,49)) * Fr(1,49)
        R = full - indep - p2
        print(f"    v={v}: R = {float(R):.6e}  R/closed = {float(R/pred3) if pred3 else float('nan'):.3f}")
        sys.stdout.flush()

SA = [1,3,5,6]; SB = [2,3,4,6]; SC = [1,2,4,5]

def PV(v, S):
    tot = Fr(0)
    for pi in permutations(S):
        tot += hit_measure(v, [[s] for s in pi])
    return tot

def cap_table():
    print("\nCAP TABLE: exact P-values on dangerous families vs candidates")
    capA, capB, capC = Fr(1,12), Fr(5,42), Fr(40,441)
    print(f"  candidates: A = {capA} ≈ {float(capA):.5f}, B = {capB} ≈ {float(capB):.5f}, C = {capC} ≈ {float(capC):.5f}")
    print("  -- cap A/B on dilate ladders (dilation-invariance check + base values)")
    for base in [(2,3,4,6),(1,2,3,4)]:
        for c in [1, 5, 11]:
            v = [c*x for x in base]
            a = PV(v, SA) + PV(v, SB); b = PV(v, SC)
            print(f"    {tuple(v)}: A-val = {a} ≈ {float(a):.5f} ({'OK' if a<=capA else 'VIOLATION?'})  "
                  f"B-val = {b} ≈ {float(b):.5f} ({'OK' if b<=capB else 'VIOLATION?'})")
        sys.stdout.flush()
    print("  -- large one-relation family (limit ≈ independent + tiny mass)")
    for (a,b,c) in [(61,108,82)]:
        v = [a,b,c,a+b-c]
        A = PV(v, SA)+PV(v, SB); Bv = PV(v, SC)
        print(f"    {tuple(v)}: A-val ≈ {float(A):.5f}  B-val ≈ {float(Bv):.5f}  "
              f"(indep 48/2401 = {float(Fr(48,2401)):.5f}, 24/2401 = {float(Fr(24,2401)):.5f})")
    print("  -- commensurate-mixed (small ratio pair + large): v = (n,2n,a,b)")
    for (n,a,b) in [(3,50,77),(5,81,123)]:
        v = [n,2*n,a,b]
        A = PV(v, SA)+PV(v, SB); Bv = PV(v, SC)
        print(f"    {tuple(v)}: A-val ≈ {float(A):.5f}  B-val ≈ {float(Bv):.5f}")
    print("  -- cap C on the quintuple maximizer ladder")
    for c in [1, 3]:
        v = [c*x for x in (3,5,7,9,11)]
        Cv = PV(v, [0]+SA) + PV(v, [0]+SB)
        print(f"    {tuple(v)}: C-val = {Cv} ≈ {float(Cv):.5f} ({'OK' if Cv<=capC else 'VIOLATION?'})")
    sys.stdout.flush()

if __name__ == "__main__":
    t0 = time.time()
    referee_rungs()
    cap_table()
    print(f"[total {time.time()-t0:.1f}s]")
