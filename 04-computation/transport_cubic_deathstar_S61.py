#!/usr/bin/env python3
"""
death-star-2026-07-20-S61 (HYP-8265) -- STAGE 2-3 of the VC-witness transport:
reduce the Yagzhev form G = X+H (deg 7) to a CUBIC-HOMOGENEOUS Keller map G_c via
companion variables (each move Phi_hat=(Phi with X^b->W, W-X^b) preserves det EXACTLY --
Schur complement = graph-restricted Jacobian), then homogenize with x0. Verify:
(a) G_c cubic homogeneous, (b) det J(G_c) == 1 (Keller => JH_c nilpotent), (c) the triple
collision transports to explicit points A,B,C with G_c(A)=G_c(B)=G_c(C). Credits mac-mini
S323 cotangent-lift idea (Poisson version); this is the Hessian-nilpotent-quartic branch.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib_exact_deathstar_S61 import *
from fractions import Fraction as Fr

NVMAX = 48
def V(i): return pvar(i, NVMAX)
def C(c): return pconst(c, NVMAX)
def monomial(expo, coeff=QI(1,0)):
    return {tuple(expo): coeff} if not coeff.is_zero() else {}

# ---------- build F, Yagzhev G, H  (in NVMAX vars; use 0,1,2 = x,y,z) ----------
x,y,z = V(0),V(1),V(2)
u = padd(C(1), pmul(x,y))
F1 = padd(pmul(ppow(u,3,NVMAX), z), pmul(pmul(ppow(y,2,NVMAX),u), padd(C(4), pscale(pmul(x,y),3))))
F2 = padd(padd(y, pscale(pmul(pmul(x,ppow(u,2,NVMAX)),z),3)),
          pscale(pmul(pmul(x,ppow(y,2,NVMAX)), padd(C(4), pscale(pmul(x,y),3))),3))
F3 = psub(psub(pscale(x,2), pscale(pmul(ppow(x,2,NVMAX),y),3)), pmul(ppow(x,3,NVMAX),z))
F = [F1,F2,F3]
Linv = [[QI(0),QI(0),QI(Fr(1,2))],[QI(0),QI(1),QI(0)],[QI(1),QI(0),QI(0)]]
G = []
for i in range(3):
    gi = pzero()
    for k in range(3): gi = padd(gi, pscale(F[k], Linv[i][k]))
    G.append(gi)
# collision points (x,y,z) for F (== for G, since G=Linv*F)
coll_orig = [[QI(0),QI(0),QI(Fr(-1,4))],
             [QI(1),QI(Fr(-3,2)),QI(Fr(13,2))],
             [QI(-1),QI(Fr(3,2)),QI(Fr(13,2))]]
print("G degrees:", [pdeg(g) for g in G])

# ---------- reduce to cubic via companion variables (with reuse) ----------
comps = [G[0],G[1],G[2]]         # map components (outputs); comps[i] = X_i + H_i
Nact = 3                          # active variable count; new companions appended at Nact
comp_def = {}                     # var index -> beta exponent tuple (its defining monomial X^beta)
beta_to_var = {}                  # beta tuple -> var index (reuse)

def choose_beta(expo):
    """pick a degree-2 sub-monomial of expo (>=4): two lowest-index units."""
    idxs = []
    for i,e in enumerate(expo):
        idxs += [i]*e
    # take two lowest-index factors
    a,b = idxs[0], idxs[1]
    beta = [0]*NVMAX; beta[a]+=1; beta[b]+=1
    return tuple(beta)

def replace_beta(poly, beta, w):
    """replace one factor X^beta by var w in every monomial of degree>=4 divisible by beta."""
    out = {}
    for e,c in poly.items():
        if sum(e) >= 4 and all(e[i] >= beta[i] for i in range(NVMAX)):
            e2 = list(e)
            for i in range(NVMAX): e2[i] -= beta[i]
            e2[w] += 1
            e2 = tuple(e2)
            out[e2] = out.get(e2, QI(0,0)) + c
        else:
            out[e] = out.get(e, QI(0,0)) + c
    return pclean(out)

steps = 0
while True:
    # find a monomial of degree >=4 in any comp
    target = None
    for ci,poly in enumerate(comps):
        for e in poly:
            if sum(e) >= 4:
                target = e; break
        if target: break
    if target is None: break
    beta = choose_beta(target)
    if beta in beta_to_var:
        w = beta_to_var[beta]
    else:
        w = Nact; Nact += 1
        beta_to_var[beta] = w; comp_def[w] = beta
        # append component for W: output = W - X^beta   (i.e. w-th coord maps to W - monomial(beta))
        comps.append(psub(V(w), monomial(beta)))
    # replace beta->w in ALL comps (degree>=4 monomials)
    comps = [replace_beta(p, beta, w) for p in comps]
    steps += 1
    if steps > 400: print("TOO MANY STEPS"); break

print(f"reduction steps={steps}, companion vars added={Nact-3}, cubic-map dim Nact={Nact}")
print("cubic map degrees:", [pdeg(p) for p in comps], "(expect all <=3)")

# ---------- section: evaluate all companion vars at an original point ----------
def section(orig_xyz):
    val = [QI(0)]*NVMAX
    val[0],val[1],val[2] = orig_xyz
    for w in sorted(comp_def):
        beta = comp_def[w]
        t = QI(1,0)
        for i in range(NVMAX):
            for _ in range(beta[i]): t = t*val[i]
        val[w] = t
    return val

# verify collision on the cubic (pre-homogenization) map
print("\ncubic map at transported collision pts (expect equal):")
cub_imgs=[]
for a in coll_orig:
    pa = section(a)
    img = [peval(p, pa) for p in comps]
    cub_imgs.append(img)
    print("  ", [str(t) for t in img[:6]], "...")
print("  collision holds:", all(all(cub_imgs[0][j]==cub_imgs[k][j] for j in range(len(comps))) for k in [1,2]))

# ---------- homogenize to cubic-homogeneous with x0 at index Nact ----------
h0 = Nact                       # x0 index
x0 = V(h0)
def homogenize(poly):
    """nonlinear part deg2 -> *x0 (deg3), deg3 stays, linear/const stay."""
    out = {}
    for e,c in poly.items():
        d = sum(e)
        if d == 2:
            e2 = list(e); e2[h0]+=1; out[tuple(e2)] = out.get(tuple(e2),QI(0,0))+c
        else:
            out[e] = out.get(e,QI(0,0))+c
    return pclean(out)
Gc = [homogenize(p) for p in comps]
Gc.append(V(h0))               # x0 -> x0 (identity)
Ndim = Nact+1
print(f"\ncubic-HOMOGENEOUS map G_c: dim N={Ndim}")
# check each nonlinear part homogeneous degree 3
def nonlin_degrees(poly, ident_var):
    ds=set()
    for e,c in poly.items():
        if sum(e)==1 and e[ident_var]==1:   # the identity linear term
            continue
        ds.add(sum(e))
    return sorted(ds)
allcubic = True
for i,p in enumerate(Gc):
    nd = nonlin_degrees(p, i)
    if nd and nd != [3]: allcubic=False
print("  every nonlinear part homogeneous deg 3:", allcubic)

# ---------- verify Keller: det J(G_c) == 1 at random points ----------
def numeric_det(M):
    n=len(M); A=[row[:] for row in M]; det=QI(1,0)
    for col in range(n):
        piv=None
        for r in range(col,n):
            if not A[r][col].is_zero(): piv=r; break
        if piv is None: return QI(0,0)
        if piv!=col: A[col],A[piv]=A[piv],A[col]; det=-det
        det = det*A[col][col]
        inv = A[col][col].inv()
        for r in range(col+1,n):
            f = A[r][col]*inv
            if f.is_zero(): continue
            for cc in range(col,n): A[r][cc]=A[r][cc]-f*A[col][cc]
    return det
JGc = jacobian(Gc, Ndim if False else NVMAX)   # NVMAX cols but only first Ndim vars used
# restrict to first Ndim vars for square Jacobian
JGc_sq = [[pderiv(Gc[i], j) for j in range(Ndim)] for i in range(Ndim)]
import random
rng = random.Random(12345)
dets=[]
for _ in range(8):
    pt=[QI(Fr(rng.randint(-5,5)), Fr(rng.randint(-5,5))) for _ in range(NVMAX)]
    M=[[peval(JGc_sq[i][j], pt) for j in range(Ndim)] for i in range(Ndim)]
    dets.append(numeric_det(M))
print("  det J(G_c) at 8 random pts:", [str(d) for d in dets], "=> Keller/det==1:", all(d==QI(1,0) for d in dets))

# spot-check JH_c nilpotent at a random point (necessary check)
def mat_pow(M,k):
    n=len(M); R=[[QI(1,0) if i==j else QI(0,0) for j in range(n)] for i in range(n)]
    for _ in range(k):
        R=[[sum((R[i][t]*M[t][j] for t in range(n)),QI(0,0)) for j in range(n)] for i in range(n)]
    return R
pt=[QI(Fr(rng.randint(-4,4)),Fr(rng.randint(-4,4))) for _ in range(NVMAX)]
JHc_pt=[[peval(pderiv(psub(Gc[i], V(i) if i<Ndim else pzero()), j), pt) for j in range(Ndim)] for i in range(Ndim)]
Nilp = mat_pow(JHc_pt, Ndim)
print("  JH_c nilpotent at random pt (JH_c^N==0):", all(Nilp[i][j].is_zero() for i in range(Ndim) for j in range(Ndim)))

# ---------- collision on G_c (with x0=1) ----------
print("\nG_c at transported collision pts A,B,C (x0=1):")
def sectionH(a):
    v=section(a); v[h0]=QI(1,0); return v
imgs=[]
for a in coll_orig:
    pa=sectionH(a); imgs.append([peval(p,pa) for p in Gc])
coll_ok = all(all(imgs[0][j]==imgs[k][j] for j in range(Ndim)) for k in [1,2])
print("  G_c(A)=G_c(B)=G_c(C):", coll_ok)
for k,a in enumerate(coll_orig):
    pa=sectionH(a)
    print(f"  A{k} first coords:", [str(t) for t in pa[:6]], "...")

# save the cubic-homog map + collision points for stage 4 (lift+rotation)
import pickle
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"transport_Gc_S61.pkl"),"wb") as f:
    pickle.dump({"Gc":Gc,"Ndim":Ndim,"coll_pts":[sectionH(a) for a in coll_orig],
                 "NVMAX":NVMAX}, f)
print(f"\nSTAGE 2-3 done. N={Ndim}, lift target dim 2N={2*Ndim}. Saved Gc for stage 4.")
