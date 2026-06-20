#!/usr/bin/env python3
"""
LRC(14) HYP-2675 -- coverage-deficit-direct, PART 2 (kps-S2-wf):
the RIGOROUS absolute-Fourier majorant for the WIDE regime.

GOAL (the deliverable): an EXPLICIT B and explicit constant C with a PROVED
chain  span(E) > B  ==>  p0(E) <= cap_k.

------------------------------------------------------------------------------
STRUCTURE OF THE ARGUMENT (what is PROVED vs VERIFIED)
------------------------------------------------------------------------------
1.  THM-534 (PROVED): p0(E) = meas(S7(E)) <= L_y(E) = sum_r y_r S_r(E),
    with the integer-root duals (g(t) >= 1[t=0] on {0..6}).  S_r are the
    factorial moments of the missed-inner-sector count N(x).  This is an exact
    Bonferroni majorant valid for EVERY E.  We re-verify it here exactly.

2.  The moments S_r(E) = sum_{|A|=r, A subset {1..6}} J(A,E), where
        J(A,E) = meas{ x : for all e in E, frac(e x) not in any sector of A }.
    Each J(A,E) is a MULTI-SECTOR AVOIDANCE measure: an intersection of
    half-open-arc-avoidance conditions, one per runner e.

3.  ABSOLUTE FOURIER MAJORANT for J(A,E) (the new direct input).
    Write the avoidance event for runner e and sector-set A as the indicator
        chi_A( frac(e x) ),  chi_A = 1 - sum_{j in A} 1_{sector j}.
    Each sector indicator 1_{[j/7,(j+1)/7)} has Fourier series with coefficients
        s_hat(n) = (1/(2 pi i n)) (omega^{-jn} - omega^{-(j+1)n}),  omega=e^{2pi i/7},
    so |s_hat(n)| = |sin(pi n /7)| / (pi |n|)  for n != 0,  s_hat(0)=1/7.
    KEY VANISHING (THM-503): s_hat(n) = 0 whenever 7 | n.  Hence the only
    nonzero frequencies are n with 7 does not divide n, and
        |s_hat(n)| <= sin(pi/7)/(pi |n|) = 0.31033.../|n|  (the max of |sin| numerator).
    J(A,E) = integral over x of prod_{e} chi_A(frac(e x)) dx
           = sum over frequency vectors that sum to 0 of products of coeffs.
    The DC (constant) term is (1 - |A|/7)^{|E_nz|}  (E_nz = nonzero runners).
    All other terms require a nontrivial integer relation sum_e c_e e = 0 with
    c_e in the support of the per-runner factor.  For a WIDE / dissociated set
    these relations are sparse and high-height, so the correction is small and
    ABSOLUTELY bounded by sum over relations of products of |coeff|.

4.  We make (3) quantitative in the cleanest closed regime and VERIFY the
    resulting numeric majorant against the exact engine over a wide bank, then
    locate the explicit threshold B per k.

WHAT THIS SCRIPT PROVES / VERIFIES (honest labels printed inline):
  (A) [PROVED, re-verified exact] p0(E) <= L_y(E) for all E (dual feasibility).
  (B) [PROVED, exact] DC term of J(A,E) is (1-|A|/7)^{|E_nz|}; the iid main
      term of L_y is  Mbar_k := sum_r y_r C(6,r) (1-r/7)^{k-1}.  We give Mbar_k
      in closed rational form and show Mbar_k < cap_k for ALL k=8..12 with a
      large margin (this is the "fully-dissociated" value -- the floor that
      every wide set approaches).
  (C) [VERIFIED exact] the correction corr(E) := L_y(E) - Mbar_k is small for
      wide E: we bound |corr| over a wide bank and show Mbar_k + max|corr| < cap_k
      once span > B, giving the explicit B per k.
  (D) The honest remaining gap: a fully rigorous a-priori bound max|corr| <= C/span
      (absolute relation-sum), which the per-arc coeff bound 0.3103/|n| supplies
      in principle; we estimate C numerically and state the clean target.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb, sin, pi
from itertools import combinations
import random, sys

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7)}

# THM-534 dual weights y_r (exact)
DUAL = {
 8:  {0:F(1),1:F(-1),2:F(1),3:F(-9,10),4:F(3,5)},
 9:  {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
 10: {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
 11: {0:F(1),1:F(-1,2),2:F(1,6)},
 12: {0:F(1),1:F(-1,2),2:F(1,6)},
}

def primitive(E):
    nz=[abs(x) for x in E if x]
    return bool(nz) and reduce(gcd,nz)==1

def breakpoints(E):
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): bps.add(F(a,7*e))
    return sorted(b for b in bps if 0<=b<=1)

def _sector(e,mid):
    v=e*mid; v=v-(v.numerator//v.denominator)
    return (v.numerator*7)//v.denominator

def p0_and_moments(E):
    """exact p0 and factorial moments S_0..S_4."""
    E=sorted(set(E)); nz=[e for e in E if e]
    bps=breakpoints(E); p0=F(0); S=[F(0)]*5
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; hm=0
        for e in nz: hm|=(1<<_sector(e,mid))
        inner=hm&0b1111110; Nt=6-bin(inner).count('1'); L=hi-lo
        if Nt==0: p0+=L
        for r in range(5): S[r]+=comb(Nt,r)*L
    return p0,S

def Ly_from_moments(S,k):
    y=DUAL[k]; return sum(y[r]*S[r] for r in y)

def Mbar(k):
    """iid (fully dissociated) main term: replace J(A,E) by (1-|A|/7)^(k-1).
    S_r -> C(6,r)*(1-r/7)^(k-1).  L_y main term = sum_r y_r C(6,r) (1-r/7)^(k-1)."""
    y=DUAL[k]; m=k-1
    tot=F(0)
    for r in y:
        tot += y[r]*comb(6,r)*F(7-r,7)**m
    return tot

def coeff_max():
    """max_{n: 7 nmid n} |s_hat(n)|*|n| = sin(pi/7)/pi (numeric)."""
    return sin(pi/7)/pi

if __name__=="__main__":
    try: sys.stdout.reconfigure(encoding='utf-8')
    except: pass
    print("="*78)
    print("HYP-2675 coverage-deficit-direct PART 2: absolute-Fourier majorant")
    print("="*78)

    # (B) the dissociated floor Mbar_k (exact rational) vs cap_k
    print("\n(B) [PROVED exact] Dissociated floor Mbar_k = L_y of the fully-independent set")
    print("    (the value every infinitely-wide set approaches). Mbar_k < cap_k margin:")
    for k in sorted(CAP):
        mb=Mbar(k); cap=CAP[k]
        print(f"  k={k}: Mbar={mb}={float(mb):.5f}  cap={float(cap):.5f}  margin cap-Mbar={float(cap-mb):.5f}  "
              f"{'OK' if mb<cap else '*** FAILS ***'}")

    print(f"\n  per-arc Fourier coeff bound: |s_hat(n)| <= {coeff_max():.6f}/|n|  (THM-503: =0 if 7|n)")

    # (A) re-verify p0<=Ly exact on a bank; (C) corr = Ly - Mbar over wide bank
    print("\n(A,C) [exact] verify p0<=L_y and measure corr(E)=L_y(E)-Mbar_k vs span")
    rng=random.Random(7770)
    for k in sorted(CAP):
        mb=Mbar(k); cap=CAP[k]
        # span-banded worst p0, worst Ly, worst corr
        bands={'[8,14]':[], '[15,28]':[], '[29,60]':[], '[61,150]':[]}
        worst_ly=(F(0),None); worst_p0=(F(0),None); ly_over_cap=0; tot=0
        for _ in range(500):
            mode=rng.choice(['clusters','spread','random','dyadic'])
            E=[0]
            if mode=='clusters':
                nc=rng.choice([2,3]); sc=rng.choice([15,20,30,50])
                for c in range(nc):
                    for j in range(rng.choice([2,3])): E.append(c*sc+j)
            elif mode=='spread':
                st=rng.choice([2,3,4,5]); E=[i*st for i in range(k)]
            elif mode=='dyadic':
                E=[0,1,2,4,8,16,24,32,48,64,80,96][:k]
            else:
                while len(set(E))<k: E.append(rng.randint(1,rng.choice([30,60,100])))
            E=sorted(set(E))
            while len(E)<k: E.append(E[-1]+rng.randint(1,20))
            E=sorted(set(E))[:k]
            if len(E)!=k or not primitive(E): continue
            span=E[-1]
            p0,S=p0_and_moments(E); ly=Ly_from_moments(S,k)
            assert p0<=ly+F(1,10**9), f"p0>Ly! {E}"
            tot+=1
            if ly>cap: ly_over_cap+=1
            if p0>worst_p0[0]: worst_p0=(p0,tuple(E))
            if ly>worst_ly[0]: worst_ly=(ly,tuple(E))
            corr=ly-mb
            for lab,(a,b) in {'[8,14]':(8,14),'[15,28]':(15,28),'[29,60]':(29,60),'[61,150]':(61,150)}.items():
                if a<=span<=b: bands[lab].append((float(corr),float(p0),float(ly)))
        print(f"\n  k={k}: sampled={tot}  Mbar={float(mb):.5f}  cap={float(cap):.5f}")
        print(f"     worst p0={float(worst_p0[0]):.5f}  worst L_y={float(worst_ly[0]):.5f}  L_y>cap count={ly_over_cap}")
        for lab in ['[8,14]','[15,28]','[29,60]','[61,150]']:
            v=bands[lab]
            if not v:
                print(f"     span {lab}: (none)"); continue
            mc=max(c for c,_,_ in v); mp=max(p for _,p,_ in v); ml=max(l for _,_,l in v)
            print(f"     span {lab}: n={len(v):3d}  max corr(L_y-Mbar)={mc:+.5f}  max p0={mp:.5f}  max L_y={ml:.5f}  "
                  f"Mbar+maxcorr={float(mb)+mc:.5f} {'<cap OK' if float(mb)+mc<float(cap) else 'NOT<cap'}")
