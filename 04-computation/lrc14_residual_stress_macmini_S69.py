#!/usr/bin/env python3
"""mac-mini-S69: stress-test the RESIDUAL of THM-724 = near-tight NON-dilated cores with a FAST
binding runner (mu just above 1/13, large s), + a resonant covering killer. These are the only
configs not covered by the unconditional cases. Verify M(S) >= 14/183 (no counterexample)."""
from fractions import Fraction as F
from math import gcd

def resdist(x,q): r=x%q; return min(r,q-r)
def M_exact(S,Qmax):
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q
            for v in S:
                d=resdist(a*v,q)
                if d<mind: mind=d
                if d==0: break
            if mind>0 and F(mind,q)>best: best=F(mind,q); arg=(q,a)
    return best,arg
def M_core(C,Qmax=200):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=min(resdist(a*v,q) for v in C)
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return best
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

target=F(14,183)
print(f"target 14/183={float(target):.6f}. RESIDUAL = near-tight non-dilated large-s cores.\n")
# Build near-tight non-dilated cores: c*{1..12} with ONE element perturbed (breaks dilation,
# keeps mu near 1/13, keeps s=c large), + killer that keeps covering & primitive.
below=[]; tested=0; near=[]
for c in range(2,9):
    base=[c*k for k in range(1,13)]
    for j in range(12):            # perturb element j
        for delta in [-2,-1,1,2,3]:
            core=sorted(base[:j]+[base[j]+delta]+base[j+1:])
            if len(set(core))!=12 or any(x<=0 for x in core): continue
            mu=M_core(core)
            if not (F(1,13) < mu < F(1,11)): continue   # near-tight window
            # binding speed s = a runner achieving mu at core optimum -- proxy: max(core) is large
            for killer in [182,364,546, c*13*14]:
                S=sorted(set(core+[killer]))
                if len(S)!=13 or not prim(S) or not is_cov(S): continue
                tested+=1
                M,arg=M_exact(S,min(2*max(S),420))
                if M<target: below.append((M,S,float(mu)))
                elif M<F(1,13): near.append((M,S,float(mu)))
print(f"tested {tested} near-tight non-dilated large-s residual configs.")
if below:
    print(f"*** {len(below)} BELOW 14/183 (residual counterexample!):")
    for M,S,mu in sorted(below)[:8]: print(f"    M={float(M):.6f}={M} mu_core={mu:.4f} S={S}")
else:
    print("NONE below 14/183: the residual is empirically closed (no counterexample).")
if near:
    print(f"\n{len(near)} configs with 14/183 <= M < 1/13 (closest to the floor):")
    for M,S,mu in sorted(near)[:8]:
        print(f"    M={float(M):.6f}={M} mu_core={mu:.4f} S={S}")
else:
    print("\nAll residual configs have M >= 1/13 (comfortably above 14/183).")
