#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the HYP-2675 "coverage-deficit-direct" claim
(kps-Sx-wf). Independent re-derivation + counterexample hunt.

The claimed proof skeleton (S1..S4):
  LINK 1  p0(E) <= L_y(E)        [PROVED, THM-534 integer-root Bonferroni dual]
  LINK 2  Mbar_k < cap_k         [PROVED, dissociated floor]
  LINK 4a span<=14 finite check  [PROVED upstream]
  LINK 4b 15<=span<=B exhaustive L_y<=cap  [VERIFIED-EXACT]
  LINK 3  span>B => corr<=cap-Mbar [CONJECTURE: |corr|<=C*W2, C not pinned]

This script INDEPENDENTLY:
  (1) re-derives LINK 1's pointwise inequality g(Nt) >= [Nt==0] exactly;
  (2) HUNTS for a wide primitive k-set (span>14, k=8,9,10) with either
        p0(E) > cap_k         (a TRUE LRC counterexample), or
        L_y(E) > cap_k        (a CERTIFIER failure -- the tail mechanism breaks);
      using adversarial families: near-APs (primitive-perturbed dilated APs),
      balanced multi-cluster sets, all-residues-mod-7 sets, resonant sets with
      small step keeping W2 large;
  (3) tests whether W2(E) actually decays with span on PRIMITIVE sets (the
      load-bearing assumption of LINK 3), by exhibiting wide primitive sets
      with LARGE W2;
  (4) checks the glue: is there a gap between the exhaustive window (S3 caps
      span at 18/19) and the asymptotic regime where corr is provably small?
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb
from itertools import combinations
import random, sys

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7)}
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

def p0_Ly(E,k):
    E=sorted(set(E)); nz=[e for e in E if e]
    bps=breakpoints(E); p0=F(0); S=[F(0)]*5
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; hm=0
        for e in nz: hm|=(1<<_sector(e,mid))
        inner=hm&0b1111110; Nt=6-bin(inner).count('1'); L=hi-lo
        if Nt==0: p0+=L
        for r in range(5): S[r]+=comb(Nt,r)*L
    y=DUAL[k]; ly=sum(y[r]*S[r] for r in y)
    return p0,ly

def W2(E):
    nz=[e for e in E if e]; tot=F(0)
    for i in range(len(nz)):
        for j in range(i+1,len(nz)):
            a,b=nz[i],nz[j]; tot+=F(gcd(a,b),a*b)
    return tot

def part1_dual():
    print("="*70)
    print("(1) LINK 1 re-derivation: g(Nt)=sum_r y_r C(Nt,r) >= [Nt==0]?")
    print("="*70)
    allok=True
    for k,y in DUAL.items():
        bad=[]
        for Nt in range(7):
            g=sum(y[r]*comb(Nt,r) for r in y)
            tgt=1 if Nt==0 else 0
            if g<tgt: bad.append((Nt,g))
            if g<0: bad.append(("NEG",Nt,g))
        status="OK (valid majorant, all g>=0)" if not bad else f"FAIL {bad}"
        if bad: allok=False
        print(f"  k={k}: {status}")
    print(f"  => LINK 1 majorant valid for all k: {allok}")
    return allok

def hunt(ks=(8,9,10), trials=40000, seed=2026):
    print("="*70)
    print("(2)+(3) HUNT wide primitive sets: p0>cap, L_y>cap, max W2")
    print("="*70)
    rng=random.Random(seed)
    bestp0={k:(F(0),None) for k in ks}
    bestly={k:(F(0),None) for k in ks}
    bestw2={k:(F(0),None) for k in ks}  # max W2 among WIDE primitive sets
    viol_p0=[]; viol_ly=[]

    def consider(k,E):
        E=sorted(set(E))
        if len(E)!=k or 0 not in E or not primitive(E): return
        if max(E)-min(E)<=14: return
        p0,ly=p0_Ly(E,k)
        if p0>bestp0[k][0]: bestp0[k]=(p0,tuple(E))
        if ly>bestly[k][0]: bestly[k]=(ly,tuple(E))
        w=W2(E)
        if w>bestw2[k][0]: bestw2[k]=(w,tuple(E))
        if p0>CAP[k]: viol_p0.append((k,tuple(E),float(p0)))
        if ly>CAP[k]: viol_ly.append((k,tuple(E),float(ly)))

    for k in ks:
        # near-AP: dilated AP step d, perturb one element to make primitive (keeps resonance)
        for d in range(2,13):
            for pos in range(1,k):
                for delta in (1,-1,2):
                    E=[d*i for i in range(k)]
                    E[pos]+=delta
                    E=sorted(set(E))
                    if len(E)!=k: continue
                    m=min(E); E=[x-m for x in E]
                    consider(k,E)
        # 2- and 3-cluster balanced
        for _ in range(trials):
            nclust=rng.choice([2,2,3])
            if nclust==2:
                s1=rng.randint(2,k-2); sizes=[s1,k-s1]
            else:
                s1=rng.randint(2,k-4); s2=rng.randint(2,k-2-s1); s3=k-s1-s2
                if s3<2: continue
                sizes=[s1,s2,s3]
            scales=sorted(rng.sample(range(0,45),nclust))
            E=set()
            ok=True
            for sz,b in zip(sizes,scales):
                el=set(); t=0
                while len(el)<sz and t<80: el.add(b+rng.randint(0,3)); t+=1
                if len(el)<sz: ok=False; break
                E|=el
            if not ok: continue
            E=sorted(E)
            if len(E)!=k: continue
            m=min(E); E=[x-m for x in E]
            consider(k,E)
        # all residues mod 7 present + wide filler
        for _ in range(trials//2):
            E={0}
            for r in range(7): E.add(r+7*rng.randint(0,3))
            while len(E)<k: E.add(rng.randint(1,45))
            E=sorted(E)
            if len(E)!=k: continue
            m=min(E); E=[x-m for x in E]
            consider(k,E)
        # small-step resonant wide: {0, d, 2d,..., (k-2)d, +1 primitivizer at far end}
        for d in range(2,8):
            for far in range(d*(k-1), d*(k-1)+8):
                E=[d*i for i in range(k-1)]+[far+1]
                E=sorted(set(E))
                if len(E)!=k: continue
                consider(k,E)

    print(f"  p0>cap violations: {len(viol_p0)}")
    for v in viol_p0[:10]: print("    P0VIOL", v)
    print(f"  L_y>cap violations: {len(viol_ly)}")
    for v in viol_ly[:10]: print("    LYVIOL", v)
    for k in ks:
        c=float(CAP[k])
        print(f"  k={k} cap={c:.4f}")
        print(f"     max p0  = {float(bestp0[k][0]):.4f}  E={bestp0[k][1]}  span={max(bestp0[k][1])-min(bestp0[k][1]) if bestp0[k][1] else '-'}")
        print(f"     max L_y = {float(bestly[k][0]):.4f}  E={bestly[k][1]}  (cap margin {c-float(bestly[k][0]):.4f})")
        print(f"     max W2  = {float(bestw2[k][0]):.4f}  E={bestw2[k][1]}  span={max(bestw2[k][1])-min(bestw2[k][1]) if bestw2[k][1] else '-'}")
    return viol_p0, viol_ly, bestw2

def part4_w2_vs_span(k=8):
    print("="*70)
    print(f"(4) Does W2 decay with span on PRIMITIVE sets? (k={k})")
    print("    Show wide primitive sets with LARGE W2 (LINK 3 assumes W2->0).")
    print("="*70)
    # near-AP step 2, span grows; perturb last to primitivize
    for d in (2,3,4,5):
        E=[d*i for i in range(k-1)]+[d*(k-1)+1]
        if not primitive(E): continue
        span=max(E)-min(E)
        p0,ly=p0_Ly(E,k)
        print(f"  step-{d} near-AP {E}: span={span} W2={float(W2(E)):.4f} p0={float(p0):.4f} L_y={float(ly):.4f} cap={float(CAP[k]):.4f}")

if __name__=="__main__":
    try: sys.stdout.reconfigure(encoding='utf-8')
    except: pass
    ok=part1_dual()
    part4_w2_vs_span(8)
    part4_w2_vs_span(10)
    vp,vl,bw=hunt(ks=(8,9,10), trials=15000)
    print("="*70)
    print("VERDICT INPUTS:")
    print(f"  LINK1 majorant valid: {ok}")
    print(f"  p0>cap counterexamples found: {len(vp)}")
    print(f"  L_y>cap (certifier failures on wide sets) found: {len(vl)}")
    print("="*70)
