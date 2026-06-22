#!/usr/bin/env python3
"""CHECK 5: true witnessG2 = meas(GOOD ∩ G_P) for consec_8 (P=cap-achiever {1,5,7,8,9})
vs the lossy Bonferroni bound (cap-p0) vs m_P. Exact + pure-python grid. kind-mendel-S1."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def lcm(a,b): return a*b//gcd(a,b)
def gl(xs): return reduce(lcm,[x for x in xs if x],1)
E=list(range(8)); P=[1,5,7,8,9]
m_P=F(14249,252252); cap8=F(2243,5880); p0c=F(481,1470)
D=14*gl(E+P); bset=set([0,D])
for e in E:
    if e==0: continue
    step=D//(7*e); x=0
    while x<=D: bset.add(x); x+=step
for p in P:
    for n in range(p):
        c=n*D//p; bset.add((c-D//(14*p))%D); bset.add((c+D//(14*p))%D)
bps=sorted(bset); ng=cg=0
for a,b in zip(bps,bps[1:]):
    if b<=a: continue
    mid2=a+b
    if not all(min((p*mid2)%(2*D),2*D-(p*mid2)%(2*D))>=2*D//14 for p in P): continue
    ng+=b-a
    sec={ (7*((e*mid2)%(2*D)))//(2*D) for e in E }
    if len(sec)==7: cg+=b-a
measGP=F(ng,D); covInGP=F(cg,D); complInGP=measGP-covInGP
print("=== consec_8 cluster, P={1,5,7,8,9} (achieves cap_8) ===")
print(f"meas(G_P)              = {complInGP+covInGP}  = {float(measGP):.6f}  (=cap_8 {float(cap8):.6f})")
print(f"meas(coverSet ∩ G_P)   = {covInGP} = {float(covInGP):.6f}   (Bonferroni uses p0={float(p0c):.4f} here -> lossy by {float(p0c-covInGP):.4f})")
print(f"meas(coverSet^c ∩ G_P) = {complInGP} = {float(complInGP):.6f}  [RIGOROUS lower bound on true G2, since coverSet^c ⊆ GOOD a.e.]")
print(f"Bonferroni floor cap-p0= {cap8-p0c} = {float(cap8-p0c):.6f}  (< m_P!)")
print(f"m_P                    = {float(m_P):.6f}")
print(f"--> true G2 >= {float(complInGP):.4f} = {float(complInGP/m_P):.2f}x m_P. Bonferroni bound is ~{float(complInGP/(cap8-p0c)):.1f}x too small.")
# pure-python grid confirmation of true maxgap-based G2
N=300000; good_gp=0
for i in range(N):
    x=(i+0.5)/N
    if not all(min((p*x)%1,1-(p*x)%1)>=1/14 for p in P): continue
    ph=sorted((e*x)%1 for e in E)
    mg=max([ph[0]+1-ph[-1]]+[ph[j+1]-ph[j] for j in range(len(ph)-1)])
    if mg>1/7: good_gp+=1
print(f"true G2 (grid {N}, maxgap>1/7) ~ {good_gp/N:.6f}  (in [{float(complInGP):.4f}, {float(measGP):.4f}], confirms >> m_P)")
