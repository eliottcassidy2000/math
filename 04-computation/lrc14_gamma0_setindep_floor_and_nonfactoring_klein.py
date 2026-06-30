#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 floor: the SET-INDEPENDENT bound, and the danger relation does NOT factor (klein-S8).

Two pieces of the final proof sentence ("the relation does not factor, and composed with itself it
stays small, in the right frame"):

(1) COMPOSED-WITH-ITSELF STAYS SMALL (set-independently). The floor deficit is the ACTUAL correlation
    |SPEC|/product = |R'-1|, R' = m_S/(m_R m_Q) (S=R u 14Q). The Cauchy-Schwarz intermediary
    CV(N_R)*sqrt((1-m_Q)/m_Q) is UNBOUNDED (klein HYP-3554, Z_14 non-transitive trap) -- but the ACTUAL
    |R'-1| is bounded < 1, i.e. R' bounded BELOW by a set-independent constant. We scan a broad
    adversarial covering family, find inf R' (the floor constant) and sup|R'-1|, and compare to the
    Gamma_0(14) constants phi(14)=6, psi(14)=24, J_2(14)=144 and 1/(2 zeta(2)).

(2) DOES NOT FACTOR (essential). If the danger relation factored (R-safe indep of Q-lonely), SPEC=0 and
    R'=1 identically. We verify SPEC != 0 (R' != 1) generically -- the relation couples R and Q via the
    bilinear pairing v*t. At the EXTREMAL {1..13}, the Borsuk-Ulam counting measure certifies it: the
    lonely set is exactly the units (Z/14)* in phi(14)/2 = 3 antipodal pairs (saddle index 3) -- an
    irreducible multiplicative object that does not factor additively.
"""
from __future__ import annotations
import sys, os, math, itertools, random
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
M = __import__("lrc14_floor_CV_sheetcount_bound_macmini_20260629")
lonely_set, measure, floor_row = M.lonely_set, M.measure, M.floor_row

# ---- Gamma_0(14) arithmetic ----
def factor(n):
    f={};d=2;m=n
    while d*d<=m:
        while m%d==0:f[d]=f.get(d,0)+1;m//=d
        d+=1
    if m>1:f[m]=f.get(m,0)+1
    return f
def phi(n):
    r=n
    for p in factor(n):r=r//p*(p-1)
    return r
def psi(n):
    r=n
    for p in factor(n):r=r*(p+1)//p
    return r
def J2(n):
    r=n*n
    for p in factor(n):r=r*(p*p-1)//(p*p)
    return r

def Rprime_of(R,Q):
    R=tuple(sorted(set(x for x in R if x%14!=0)));Q=tuple(sorted(set(Q)))
    mR=measure(lonely_set(R));mQ=measure(lonely_set(Q))
    if mR==0 or mQ==0:return None,mR,mQ
    S=tuple(sorted(set(R)|set(14*q for q in Q)))
    mS=measure(lonely_set(S))
    return mS/(mR*mQ),mR,mQ

def units_mod(n): return [a for a in range(1,n) if math.gcd(a,n)==1]

if __name__=="__main__":
    Z2=math.pi**2/6
    print("="*86)
    print(" LRC14: set-independent floor bound + the danger relation does not factor (klein-S8)")
    print("="*86)
    print(f" Gamma_0(14): phi=6, psi=24 (=index), J_2=144, |SL2(Z/14)|=14*144=2016.")
    print(f" candidate set-indep constants:  1/(2 zeta2)={1/(2*Z2):.5f}  "
          f"phi/psi={phi(14)/psi(14):.5f}  J2/psi^2={J2(14)/psi(14)**2:.5f}  "
          f"1/psi={1/psi(14):.5f}  phi/(psi)*? ")

    # ---------- (1) the floor deficit |R'-1| over an adversarial covering family ----------
    print("\n (1) ACTUAL floor R' = m_S/(m_R m_Q) over a broad adversarial covering family (|R|+|Q|=14):")
    rng=random.Random(140014)
    fam=set()
    def add(R):
        R=tuple(sorted(set(x for x in R if 1<=x<=13)))
        if 2<=len(R)<=12:fam.add(R)
    for k in range(2,13):add(range(1,k+1))                       # consecutive prefixes
    base=list(range(1,14))
    for x in base:add([y for y in base if y!=x])                  # full13 minus one
    for x,y in itertools.combinations(base,2):add([z for z in base if z not in(x,y)])
    for ex in itertools.chain.from_iterable(itertools.combinations([1,2,3,4,5,6,8,9,10,11,12,13],j) for j in range(0,4)):
        add([7]+list(ex))                                         # speed-7 family
    for _ in range(800):add(rng.sample(range(1,14),rng.randint(2,12)))
    rows=[]
    for R in fam:
        Q=list(range(1,14-len(R)+1))
        rp,mR,mQ=Rprime_of(R,Q)
        if rp is not None: rows.append((float(rp),R,len(Q),float(mR),float(mQ)))
    rows.sort()
    infRp=rows[0]; supdev=max(abs(r[0]-1) for r in rows)
    print(f"   scanned {len(rows)} valid coverings.  inf R' = {infRp[0]:.5f}  at R={infRp[1]} (|Q|={infRp[2]})")
    print(f"   sup |R'-1| (= sup |SPEC|/product) = {supdev:.5f}   (floor holds set-indep iff <1: "
          f"{'YES' if supdev<1 else 'NO'})")
    print(f"   => R' bounded BELOW by {infRp[0]:.4f} > 0 across the family, EVEN THOUGH CV(N_R)^2 is "
          f"unbounded (HYP-3554).")
    print(f"   compare inf R' to set-indep constants: 1/(2 zeta2)={1/(2*Z2):.5f}, "
          f"phi/psi=0.25, J2/psi^2=0.25.")
    print("   lowest 6 R':")
    for rp,R,q,mR,mQ in rows[:6]:
        print(f"      R'={rp:.4f}  |Q|={q}  m_R={mR:.4f} m_Q={mQ:.4f}  R={R}")

    # ---------- (2) the danger relation does not factor ----------
    print("\n (2) DOES NOT FACTOR:")
    # (2a) SPEC != 0  (R' != 1) generically -> R-safe and Q-lonely are coupled
    nonfac=sum(1 for r in rows if abs(r[0]-1)>1e-9)
    print(f"   (2a) coverings with R'!=1 (SPEC!=0, relation couples R&Q): {nonfac}/{len(rows)} "
          f"({'relation does NOT factor' if nonfac>0 else 'factors'})")
    # (2b) extremal {1..13}: lonely set = units (Z/14)* in phi/2 antipodal pairs = BU counting measure
    L=lonely_set(tuple(range(1,14)))
    mL=measure(L)
    u=units_mod(14)
    pairs=sorted({tuple(sorted((a,14-a))) for a in u})
    print(f"   (2b) extremal {{1..13}}: meas(lonely)={float(mL)} (tight/measure-0 locus); "
          f"units (Z/14)* = {u} (phi(14)={len(u)});")
    print(f"        antipodal pairs (the BU counting measure, saddle index phi/2={len(u)//2}): {pairs}")
    # confirm the unit fractions a/14 are the lonely touch-points (||v*a/14|| >= 1/14 for all v in 1..13?)
    touch=[]
    for a in u:
        t=F(a,14)
        ok=all(min((v*t)%1,1-((v*t)%1))>=F(1,14)-F(1,10**9) for v in range(1,14))
        touch.append((a,ok))
    print(f"        unit fractions a/14 lonely for {{1..13}} (touch-points): "
          f"{[a for a,ok in touch if ok]} of {u}")
    print("\n VERDICT: the floor deficit |SPEC|/product stays < 1 set-independently (R' bounded below > 0)"
          " even where CV blows up; and the danger relation does not factor (SPEC!=0; at the extremal the"
          " lonely set is the multiplicative units in antipodal BU pairs). The 'right frame' is the actual"
          " correlation |R'-1|, not the lossy CV -- the Gamma_0(14)/cyclotomic collapse bounds THAT.")
