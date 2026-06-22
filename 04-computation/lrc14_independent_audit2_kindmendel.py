#!/usr/bin/env python3
"""CHECK 3 (positivity linchpin) + CHECK 4 (cap-bound search). Fast integer p0.
kind-mendel-2026-06-22-S1."""
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
from itertools import combinations
def lcm(a,b): return a*b//gcd(a,b)
def gl(xs): return reduce(lcm,xs,1)

def p0_fast(E):
    "exact meas(S7(E)) as Fraction, integer breakpoints."
    pos=[e for e in E if e!=0]
    if not pos: return F(0)
    D=7*gl(pos)
    bset=set([0,D])
    for e in pos:
        step=D//(7*e)        # breakpoints at multiples of 1/(7e) = step/D
        x=0
        while x<=D:
            bset.add(x); x+=step
    bps=sorted(bset)
    num=0
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b              # 2*midpoint*D
        # frac(e*mid) sector = floor(7*frac); mid = mid2/(2D)
        sec=set()
        for e in E:
            r=(e*mid2)%(2*D)          # 2D * frac(e*mid)
            sec.add((7*r)//(2*D))
            if len(sec)==7: break
        if len(sec)==7:
            num+=b-a
    return F(num,D)

def meas_GP(P):
    if not P: return F(1)
    D=14*gl(P)
    bset=set([0,D])
    for p in P:
        for n in range(p):
            c=n*D//p
            bset.add((c-D//(14*p))%D); bset.add((c+D//(14*p))%D)
    bps=sorted(bset)
    num=0
    for a,b in zip(bps,bps[1:]):
        mid2=a+b
        ok=True
        for p in P:
            r=(p*mid2)%(2*D)         # 2D*frac(p*mid)
            # ||p*mid|| = min(r,2D-r)/(2D) >= 1/14  <=> min(r,2D-r) >= 2D/14 = D/7
            if min(r,2*D-r) < 2*D//14: ok=False; break
        if ok: num+=b-a
    return F(num,D)

assert p0_fast(list(range(8)))==F(481,1470), "p0_fast self-check failed"
m_P=F(14249,252252)
# cap_k = min meas(G_P) over |P|=13-k (published minima)
capvals={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
print("=== CHECK 3: positivity linchpin  cap_k vs p0(consec_k) ===")
print(f"m_P = {m_P} = {float(m_P):.6f}")
for k in range(8,14):
    pc=p0_fast(list(range(k)))
    cap=capvals[k]
    marg=cap-pc
    print(f"k={k:2d}: cap={float(cap):.5f}  p0(consec)={float(pc):.5f}  cap-p0={float(marg):+.5f}  "
          f">0:{marg>0}  >=m_P:{marg>=m_P}")

print("\n=== CHECK 4: cap-bound search  max_E p0(E) vs cap_k (E=0+7-subset, spread<=S) ===")
for k,S in [(8,16),(9,16),(10,16)]:
    cap=capvals[k]
    best=(F(-1),None); over=0; cnt=0
    for sub in combinations(range(1,S+1),k-1):
        E=[0]+list(sub)
        if gl(E)!=1: continue          # primitive
        cnt+=1
        v=p0_fast(E)
        if v>best[0]: best=(v,E)
        if v>=cap: over+=1
    consec=p0_fast(list(range(k)))
    print(f"k={k} spread<={S}: primitive sets={cnt}  max p0={float(best[0]):.5f} at {best[1]}  "
          f"cap_{k}={float(cap):.5f}  #(p0>=cap)={over}  consec_max:{best[1]==list(range(k))} (consec p0={float(consec):.5f})")
