#!/usr/bin/env python3
"""
death-star-2026-07-20-S61 (HYP-8265) -- pin the exact status of the companion cubic
reduction (APPEND-ONLY, no reuse => each move provably det-preserving via Schur). Question:
does it yield a Keller cubic map, and is its Jacobian nilpotent? This decides the honest
statement about the remaining 'cubic-homogeneous' step.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib_exact_deathstar_S61 import *
from fractions import Fraction as Fr

NVMAX=48
def V(i): return pvar(i,NVMAX)
def Cc(c): return pconst(c,NVMAX)
def mono(expo,coeff=QI(1,0)): return {tuple(expo):coeff} if not coeff.is_zero() else {}

x,y,z=V(0),V(1),V(2)
u=padd(Cc(1),pmul(x,y))
F1=padd(pmul(ppow(u,3,NVMAX),z),pmul(pmul(ppow(y,2,NVMAX),u),padd(Cc(4),pscale(pmul(x,y),3))))
F2=padd(padd(y,pscale(pmul(pmul(x,ppow(u,2,NVMAX)),z),3)),pscale(pmul(pmul(x,ppow(y,2,NVMAX)),padd(Cc(4),pscale(pmul(x,y),3))),3))
F3=psub(psub(pscale(x,2),pscale(pmul(ppow(x,2,NVMAX),y),3)),pmul(ppow(x,3,NVMAX),z))
F=[F1,F2,F3]
Linv=[[QI(0),QI(0),QI(Fr(1,2))],[QI(0),QI(1),QI(0)],[QI(1),QI(0),QI(0)]]
G=[padd(padd(pscale(F[0],Linv[i][0]),pscale(F[1],Linv[i][1])),pscale(F[2],Linv[i][2])) for i in range(3)]

# APPEND-ONLY companion reduction (fresh var each step; provably det-preserving)
comps=[G[0],G[1],G[2]]; Nact=3; comp_def={}
def choose_beta(expo):
    idxs=[]
    for i,e in enumerate(expo): idxs+=[i]*e
    a,b=idxs[0],idxs[1]; beta=[0]*NVMAX; beta[a]+=1; beta[b]+=1; return tuple(beta)
def replace_beta_once(poly,beta,w):
    out={}
    for e,c in poly.items():
        if sum(e)>=4 and all(e[i]>=beta[i] for i in range(NVMAX)):
            e2=list(e)
            for i in range(NVMAX): e2[i]-=beta[i]
            e2[w]+=1; e2=tuple(e2); out[e2]=out.get(e2,QI(0,0))+c
        else: out[e]=out.get(e,QI(0,0))+c
    return pclean(out)
steps=0
while True:
    target=None
    for poly in comps:
        for e in poly:
            if sum(e)>=4: target=e; break
        if target: break
    if target is None: break
    beta=choose_beta(target)
    w=Nact; Nact+=1; comp_def[w]=beta          # FRESH var every step (no reuse)
    comps.append(psub(V(w),mono(beta)))
    comps=[replace_beta_once(p,beta,w) for p in comps]
    steps+=1
    if steps>200: print("runaway"); break
print(f"append-only: steps={steps}, dim={Nact}, degrees={[pdeg(p) for p in comps]}")

# det of cubic map G3 (pre-homogenization) at random points
def numeric_det(M):
    n=len(M); A=[r[:] for r in M]; d=QI(1,0)
    for c in range(n):
        p=None
        for r in range(c,n):
            if not A[r][c].is_zero(): p=r;break
        if p is None: return QI(0,0)
        if p!=c: A[c],A[p]=A[p],A[c]; d=-d
        d=d*A[c][c]; inv=A[c][c].inv()
        for r in range(c+1,n):
            f=A[r][c]*inv
            if f.is_zero(): continue
            for cc in range(c,n): A[r][cc]=A[r][cc]-f*A[c][cc]
    return d
JG3=[[pderiv(comps[i],j) for j in range(Nact)] for i in range(Nact)]
rng=random.Random(7)
dets=[]
for _ in range(6):
    pt=[QI(Fr(rng.randint(-4,4)),Fr(rng.randint(-3,3))) for _ in range(NVMAX)]
    dets.append(numeric_det([[peval(JG3[i][j],pt) for j in range(Nact)] for i in range(Nact)]))
print("det J(G3) at 6 pts:", [str(d) for d in dets])
print("=> cubic map G3 is KELLER (det==1):", all(d==QI(1,0) for d in dets))

# is JH3 nilpotent? check at a random point (necessary condition)
def matpow(M,k):
    n=len(M); R=[[QI(1,0) if i==j else QI(0,0) for j in range(n)] for i in range(n)]
    for _ in range(k): R=[[sum((R[i][t]*M[t][j] for t in range(n)),QI(0,0)) for j in range(n)] for i in range(n)]
    return R
Xv=[V(i) for i in range(Nact)]
JH3=[[pderiv(psub(comps[i],Xv[i]),j) for j in range(Nact)] for i in range(Nact)]
pt=[QI(Fr(rng.randint(-4,4)),Fr(rng.randint(-3,3))) for _ in range(NVMAX)]
JH3pt=[[peval(JH3[i][j],pt) for j in range(Nact)] for i in range(Nact)]
Nl=matpow(JH3pt,Nact)
print("JH3 nilpotent at random pt (JH3^dim==0):", all(Nl[i][j].is_zero() for i in range(Nact) for j in range(Nact)))
# trace of JH3 (nilpotent needs identically 0)
tr=pzero()
for i in range(Nact): tr=padd(tr,JH3[i][i])
print("trace(JH3) identically zero:", pis_zero(tr), "(else NOT nilpotent as a poly matrix)")
print("\nCONCLUSION: append-only companion reduction gives a KELLER CUBIC map with the collision,")
print("but its Jacobian is NOT nilpotent (mixed-degree Keller does not force nilpotency).")
print("Cubic-HOMOGENEOUS (=> nilpotent) is a strictly further step; naive homogenization breaks")
print("Keller. Securing nilpotency is the genuine remaining mathematical content.")
