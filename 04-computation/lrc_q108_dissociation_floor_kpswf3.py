#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_dissociation_floor_kpswf3.py  (kind-pasteur 2026-06-21, THREAD D lead)

The one genuinely under-probed observation from the diverse sweep:
   MULTIPLICATIVE / additive DISSOCIATION dramatically LOWERS measS7 toward the
   decorrelated coupon floor C(k).  consec (max additive energy / Schur triples)
   is the MAX; Sidon/dissociated sets (min additive energy) approach the floor.

This is the COMPLEMENT of the consec-maximizer (HYP-2635(2)/HYP-2709): the
"anti-maximizer" side.  Relevant because the LRC bound only needs the MAX, but a
clean *monotone* relationship  measS7 ~ increasing function of additive energy /
Schur-triple count  would give a PROOF route for consec-max: maximize additive
energy <=> consec (Freiman), and measS7 monotone in it.  HYP-2719 already located the
binding carrier error at SUPPORT-3 = Schur triples (a+b=c).  So test directly:

  Is measS7(E) monotone in the number of Schur triples (support-3 additive relations)?

If YES (monotone), then consec (max Schur triples among k-sets) => max measS7, giving
the consec-max crux a structural proof handle via additive energy.  If NO (non-monotone),
record the counterexample (this matches HYP-2708's 'measS7 non-concave / mixed-sign'
warning and would be an honest DEAD rating for this specific monotonicity).
"""
import itertools, random
from fractions import Fraction as Fr
from math import gcd, comb

P=7
def sector(yf): return int(P*yf)
def measS7(E):
    E=[int(e) for e in E if int(e)!=0]
    bp={Fr(0),Fr(1)}
    for e in E:
        for t in range(0,P*e): bp.add(Fr(t,P*e))
    xs=sorted(bp); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        if len(set(sector((e*mid)%1) for e in E))==P: tot+=(b-a)
    return tot

def schur_triples(E):
    """count ordered (a,b,c) in E^3 with a+b=c (a<=b), the support-3 additive relations."""
    S=set(E); cnt=0
    for a in E:
        for b in E:
            if a<=b and (a+b) in S: cnt+=1
    return cnt

def additive_energy(E):
    """|{(a,b,c,d): a+b=c+d}| = sum_s r(s)^2, r(s)=#pairs summing to s."""
    from collections import Counter
    r=Counter()
    for a in E:
        for b in E: r[a+b]+=1
    return sum(v*v for v in r.values())

def coupon(k):
    return sum((-1)**j*comb(P,j)*Fr((P-j)**k,P**k) for j in range(P+1))

def main():
    random.seed(20260621)
    print("#"*78)
    print("# DISSOCIATION FLOOR: is measS7 monotone in additive energy / Schur triples?")
    print("#"*78)

    for k in [8,9]:
        floor=coupon(k)
        consec=list(range(1,k+1))
        print(f"\n=== k={k}: coupon floor C(k)={float(floor):.5f}, consec measS7={float(measS7(consec)):.5f} ===")
        # sample many primitive k-sets, record (additive_energy, schur, measS7)
        data=[]
        seen=set()
        # include consec, AP's, Sidon-ish, primes, geometric
        specials={
            "consec[1..k]":consec,
            "AP d=2":[1+2*i for i in range(k)],
            "primes":[2,3,5,7,11,13,17,19,23][:k],
            "powers2":[2**i for i in range(k)],
            "Sidon-ish":[1,2,5,11,22,33,44,55,66][:k] if k<=9 else None,
        }
        for name,E in specials.items():
            if E is None: continue
            ae=additive_energy(E); st=schur_triples(E); m=measS7(E)
            print(f"   {name:14s} AE={ae:5d} Schur={st:3d} measS7={float(m):.5f}")
        # random sample to test monotonicity in Schur count
        for _ in range(250):
            E=sorted(random.sample(range(1,26),k))
            if gcd(*E) if k==2 else 0: pass
            key=tuple(E)
            if key in seen: continue
            seen.add(key)
            data.append((schur_triples(E),additive_energy(E),float(measS7(E))))
        # correlation of measS7 with schur count and AE
        import statistics
        sc=[d[0] for d in data]; ae=[d[1] for d in data]; ms=[d[2] for d in data]
        def corr(xs,ys):
            n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
            num=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
            dx=sum((x-mx)**2 for x in xs)**.5; dy=sum((y-my)**2 for y in ys)**.5
            return num/(dx*dy) if dx*dy else 0
        print(f"   [{len(data)} random k-sets] corr(measS7, Schur#)={corr(sc,ms):+.3f}"
              f"   corr(measS7, AddEnergy)={corr(ae,ms):+.3f}")
        # strict monotonicity test: any pair with more schur but LOWER measS7?
        viol=0; checks=0
        for i in range(len(data)):
            for j in range(len(data)):
                if data[i][0]>data[j][0]:
                    checks+=1
                    if data[i][2]<data[j][2]-1e-12: viol+=1
        print(f"   strict monotone (Schur# up => measS7 up)? violations {viol}/{checks}"
              f"  => {'MONOTONE' if viol==0 else 'NON-MONOTONE (Schur# alone does not order measS7)'}")

    print("\nDONE.")

if __name__=="__main__":
    main()
