#!/usr/bin/env python3
"""
death-star-2026-07-20-S61b (HYP-8300) -- FIX the S61 global-det bug for the BCW reduction.
S61's companion move Phi=(Phi with X^b->W, W-X^b) was Keller only ON the section {W=X^b}.
The CORRECT BCW move is a SHEAR COMPOSITION E1.(F(+)id).E2 with E1,E2 elementary shears
(det J = 1 identically), so det is preserved GLOBALLY. Concrete gadget on (vars,u):
   (vars,u) -> (F(vars) + h(u + g(vars)),  u + g(vars))
= [shear vars by h(u)] . (F(+)id) . [shear u by g(vars)]. Demonstrate: (A) det is globally
constant = det JF for ANY g,h (unlike S61); (B) with a perfect-power term, the shear
REDUCES its degree while keeping det constant -- the mechanism BCW iterates to cubic.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib_exact_deathstar_S61 import *
from fractions import Fraction as Fr

NV=8
def V(i): return pvar(i,NV)
def Cc(c): return pconst(c,NV)

# F (Alpoge), det JF = -2
x,y,z=V(0),V(1),V(2)
u_=padd(Cc(1),pmul(x,y))
F1=padd(pmul(ppow(u_,3,NV),z),pmul(pmul(ppow(y,2,NV),u_),padd(Cc(4),pscale(pmul(x,y),3))))
F2=padd(padd(y,pscale(pmul(pmul(x,ppow(u_,2,NV)),z),3)),pscale(pmul(pmul(x,ppow(y,2,NV)),padd(Cc(4),pscale(pmul(x,y),3))),3))
F3=psub(psub(pscale(x,2),pscale(pmul(ppow(x,2,NV),y),3)),pmul(ppow(x,3,NV),z))
F=[F1,F2,F3]

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
def psubst(poly, subs):
    res=pzero()
    for e,c in poly.items():
        term={tuple([0]*NV):c}
        for i,k in enumerate(e):
            if k: term=pmul(term, ppow(subs.get(i,pvar(i,NV)),k,NV))
        res=padd(res,term)
    return res

def shear_gadget(comps, g, h_of_u, uidx):
    """(vars,u) -> (comps + h(u+g), u+g). h_of_u: poly in var uidx. g: poly in old vars."""
    ug = padd(V(uidx), g)                      # u + g(vars)
    h_at = psubst(h_of_u, {uidx: ug})          # h(u+g)
    new = [padd(comps[i], h_at) for i in range(len(comps))]  # NOTE: h maps into ALL comps? no-- into chosen targets
    return new, ug

print("=== (A) shear gadget preserves det GLOBALLY (contrast S61's on-section-only) ===")
# gadget: shear only the X-components; here put h into component 0 only (shear x by h(u))
# build map on (x,y,z,u): comp = [F1 + h(u+g), F2, F3, u+g], with g=xy, h(u)=u^2
g = pmul(x,y)
uidx=3
hh = ppow(V(uidx),2,NV)                        # h(u)=u^2
ug = padd(V(uidx), g)
h_at = psubst(hh, {uidx: ug})                  # (u+xy)^2
mp = [padd(F[0], h_at), F[1], F[2], ug]        # 4-var map
rng=random.Random(3)
dets=[]
for _ in range(6):
    pt=[QI(Fr(rng.randint(-4,4)),Fr(rng.randint(-3,3))) for _ in range(NV)]
    J=[[peval(pderiv(mp[i],j),pt) for j in range(4)] for i in range(4)]
    dets.append(numeric_det(J))
print("  det at 6 random pts:", [str(d) for d in dets])
print("  GLOBALLY CONSTANT = det JF = -2:", all(d==QI(-2,0) for d in dets), "  <-- S61 bug FIXED")

print("\n=== (B) a shear REDUCES a perfect-power term's degree, det still constant ===")
# craft F' with an added degree-6 term (xy)^3 in component 0 (still det -2? no -- adding changes det).
# Instead: demonstrate on a Yagzhev-like map M0 = X + (xy)^3 e_2 (JH nilpotent, det 1), reduce (xy)^3.
M0 = [x, y, padd(z, ppow(pmul(x,y),3,NV))]     # z -> z + (xy)^3 ; JH strictly upper => nilpotent, det 1
print("  M0 degree:", [pdeg(p) for p in M0], " det(M0):", end=" ")
d0=[]
for _ in range(4):
    pt=[QI(Fr(rng.randint(-3,3)),Fr(rng.randint(-2,2))) for _ in range(NV)]
    d0.append(numeric_det([[peval(pderiv(M0[i],j),pt) for j in range(3)] for i in range(3)]))
print([str(v) for v in d0])
# shear: g=xy (deg2) in u; h(u) = -u^3 ; put h into component 2 (the z-component)
uidx=3
g=pmul(x,y); hh=pscale(ppow(V(uidx),3,NV),-1)  # h(u) = -u^3
ug=padd(V(uidx), g)
h_at=psubst(hh,{uidx:ug})                       # -(u+xy)^3
M1=[M0[0], M0[1], padd(M0[2], h_at), ug]        # (x,y,z+(xy)^3-(u+xy)^3, u+xy)
M1=[pclean(p) for p in M1]
print("  after shear: degrees", [pdeg(p) for p in M1], " (deg-6 (xy)^3 term cancels -> top deg drops)")
d1=[]
for _ in range(5):
    pt=[QI(Fr(rng.randint(-3,3)),Fr(rng.randint(-2,2))) for _ in range(NV)]
    d1.append(numeric_det([[peval(pderiv(M1[i],j),pt) for j in range(4)] for i in range(4)]))
print("  det after shear at 5 pts:", [str(v) for v in d1], " GLOBALLY CONSTANT = det M0 = 1:", all(v==QI(1,0) for v in d1))
print("\n  => shear composition = the correct BCW move: det preserved GLOBALLY (not on-section),")
print("     and it reduces degree. Iterating to cubic-homogeneous is the exact construction the")
print("     research agent is retrieving; the framework fix (S61 bug) is DEMONSTRATED here.")
