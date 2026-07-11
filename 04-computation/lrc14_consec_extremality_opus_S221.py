"""
opus-2026-07-11-S221: attacking the JOINT Phi-consec-extremality (the single remaining LRC(14)-S3 lemma).

Targets: (A) consec minimizes 4m1 - m2 = 5 E[N] - E[N^2] over bounded k-cores (k=9,10) -- the degree-2 (k=9,10)
piece; (B) the JOINT Phi = p0 + (1/3)p1 extremality (kps: p0 maxed at consec but p1 NOT; joint via lambda=1/3
< lambda*). N = # missed inner sectors {1..6}; sector 0 always covered by the stationary offset 0.

Approach: EXHAUSTIVE over primitive 9-/10-subsets of {0..D} (0 in E) -- confirm consec is the argmin of 4m1-m2
AND the argmax of Phi and of (p0 + lambda p1) for lambda up to lambda*; measure the STABILITY gap (2nd best)
and WHICH perturbations move away from consec (the structure a proof must exploit).
"""
from fractions import Fraction as F
from itertools import combinations

def sector(e,x): return int(((e*x)%1)*7)
def brk(E):
    Es=[abs(e) for e in E if e!=0]; bps={F(0),F(1)}
    for e in Es:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(b for b in bps if 0<=b<1)
def missdist(E):
    pts=brk(E)+[F(1)]; p=[F(0)]*7
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        occ=set(sector(e,(a+b)/2) for e in E)
        p[len(set(range(1,7))-occ)]+=(b-a)
    return p
def stats(E):
    p=missdist(E)
    m1=sum(F(n)*p[n] for n in range(7)); m2=sum(F(n*(n-1))*p[n] for n in range(7))
    return p[0], p[1], m1, m2   # p0, p1, m1, m2

def consec(k): return list(range(k))

print("=== consec exact m1,m2,4m1-m2, p0,p1 for k=9,10 (closed-form target) ===")
for k in [9,10]:
    p0,p1,m1,m2=stats(consec(k))
    print(f"  consec_{k}: m1={m1}={float(m1):.5f}  m2={m2}={float(m2):.5f}  4m1-m2={4*m1-m2}={float(4*m1-m2):.5f}  p0={p0}  p1={p1}")

for k, D in [(9,12),(10,12)]:
    print(f"\n=== EXHAUSTIVE k={k}, diameter<=D={D}: argmin(4m1-m2) and argmax(Phi) ===")
    base=consec(k); bp0,bp1,bm1,bm2=stats(base); bval=4*bm1-bm2; bphi=bp0+F(1,3)*bp1
    results=[]  # (4m1-m2, Phi, E)
    n=0
    for rest in combinations(range(1,D+1), k-1):
        E=[0]+list(rest); n+=1
        p0,p1,m1,m2=stats(E)
        results.append((4*m1-m2, p0+F(1,3)*p1, p0, p1, E))
    results_val=sorted(results, key=lambda r: r[0])            # by 4m1-m2 ascending (consec should be min)
    results_phi=sorted(results, key=lambda r: -r[1])           # by Phi descending (consec should be max)
    print(f"  scanned {n} cores. consec_{k}: 4m1-m2={float(bval):.5f}, Phi={float(bphi):.5f}")
    print(f"  --- smallest 4m1-m2 (consec should be #1): ---")
    for r in results_val[:4]:
        tag=" <-- consec" if r[4]==base else ""
        print(f"     4m1-m2={float(r[0]):.5f}  Phi={float(r[1]):.5f}  {r[4]}{tag}")
    print(f"  --- largest Phi (consec should be #1): ---")
    for r in results_phi[:4]:
        tag=" <-- consec" if r[4]==base else ""
        print(f"     Phi={float(r[1]):.5f}  p0={float(r[2]):.5f}  p1={float(r[3]):.5f}  {r[4]}{tag}")
    is_min = results_val[0][4]==base
    is_max = results_phi[0][4]==base
    print(f"  consec is argmin(4m1-m2)? {is_min}   argmax(Phi)? {is_max}")
    if is_min:
        gap = results_val[1][0]-bval
        print(f"  stability gap (2nd - consec) in 4m1-m2 = {float(gap):.6f}; 2nd = {results_val[1][4]}")
