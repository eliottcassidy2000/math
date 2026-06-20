import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

def meas_S7(E):
    E=sorted(set(e for e in E if e!=0))
    if not E: return Fraction(0)
    bps={Fraction(0),Fraction(1)}
    for e in E:
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(bps); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        x=(lo+hi)/2; sectors={0}
        for e in E:
            v=e*x; v=v-(v.numerator//v.denominator); sectors.add((v.numerator*7)//v.denominator)
            if len(sectors)==7: break
        if len(sectors)==7: tot+=hi-lo
    return tot

def dyadic_decomp(E):
    d={}
    for e in E:
        b,a=e,0
        while b%2==0: b//=2; a+=1
        d.setdefault(b,[]).append(a)
    for b in d: d[b].sort()
    return d
def c2(E): return sum(len(v)*(len(v)-1)//2 for v in dyadic_decomp(E).values())
def max_consec_run(E):
    E=sorted(set(E)); run=1; mx=1
    for i in range(1,len(E)):
        if E[i]==E[i-1]+1: run+=1; mx=max(mx,run)
        else: run=1
    return mx

print("="*78)
print("k=9 doubling-monotonicity: replace one elt of consec{1..9} by larger value")
print("="*78)
base=list(range(1,10)); base_p0=meas_S7(base)
print(f"base {base}: p0={float(base_p0):.5f}")
down=tot=0; worst=None
for i in range(9):
    for new in range(10,14):
        if new in base: continue
        E=base[:i]+[new]+base[i+1:]; p=meas_S7(E); tot+=1
        if p<base_p0: down+=1
        else: worst=(E,p)
print(f"   {down}/{tot} strictly decreased. counterexample(>=): {worst}")

print()
print("="*78)
print("ISOLATE: fix max-consec-run, vary dyadic richness. Does richness still help")
print("WITHIN a fixed contiguity stratum? (k=8, group by max_consec_run)")
print("="*78)
from collections import defaultdict
strata=defaultdict(list)
for E in itertools.combinations(range(1,14),8):
    strata[max_consec_run(E)].append((float(meas_S7(E)), c2(E), E))
def pearson(xs,ys):
    n=len(xs);
    if n<2: return float('nan')
    mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy) if sx*sy else float('nan')
for run in sorted(strata):
    rows=strata[run]; p0s=[r[0] for r in rows]; r2s=[r[1] for r in rows]
    print(f"   run={run}: n={len(rows):4d}  mean(p0)={sum(p0s)/len(p0s):.4f}  max(p0)={max(p0s):.4f}  corr(p0,richness)={pearson(p0s,r2s):+.4f}")
print("=> mean p0 should INCREASE with contiguity-run; richness-corr within strata ~0 if contiguity is the true driver.")

print()
print("="*78)
print("ERGODIC / EQUIDISTRIBUTION ANGLE (quantitative):")
print("The dyadic doubling-orbit of b under x->2x equidistributes. For consec, the")
print("orbit {frac(e x)} is the densest 'multi-scale' net. Measure the WORST-CASE")
print("uncovered: max over x of (#sectors MISSED). Lower worst-miss => higher p0.")
print("="*78)
def worst_miss(E):
    """max over breakpoint-cells of (7 - #sectors hit). 0 means always full cover."""
    E=sorted(set(e for e in E if e!=0))
    bps={Fraction(0),Fraction(1)}
    for e in E:
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(bps); wm=0; argx=None
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        x=(lo+hi)/2; sectors={0}
        for e in E:
            v=e*x; v=v-(v.numerator//v.denominator); sectors.add((v.numerator*7)//v.denominator)
        miss=7-len(sectors)
        if miss>wm: wm=miss; argx=x
    return wm,argx
for lab,E in [("consec{1..8}",list(range(1,9))),
              ("consec{1..9}",list(range(1,10))),
              ("dyadic-rich {1,2,3,4,6,8,12,5}",[1,2,3,4,6,8,12,5]),
              ("apex {1..7,13}",list(range(1,8))+[13])]:
    wm,ax=worst_miss(E)
    print(f"   {lab:34s} p0={float(meas_S7(E)):.5f}  worst-miss={wm} sectors @x={ax}")

print()
print("="*78)
print("VERDICT TABLE: top-3 controlling variables for p0 (k=8, full enumeration)")
print("="*78)
rows=[(float(meas_S7(E)),max_consec_run(E),c2(E),E) for E in itertools.combinations(range(1,14),8)]
p0s=[r[0] for r in rows]
print(f"   corr(p0, max-consec-run) = {pearson(p0s,[r[1] for r in rows]):+.4f}  <- contiguity")
print(f"   corr(p0, dyadic-richness)= {pearson(p0s,[r[2] for r in rows]):+.4f}  <- dyadic (weaker, partly confounded)")
# partial: among run=1 sets only, is richness predictive?
run1=[(r[0],r[2]) for r in rows if r[1]==1]
print(f"   Within non-contiguous (run=1) sets: corr(p0,richness)={pearson([a for a,_ in run1],[b for _,b in run1]):+.4f}")
