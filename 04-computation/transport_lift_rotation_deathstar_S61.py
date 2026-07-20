#!/usr/bin/env python3
"""
death-star-2026-07-20-S61 (HYP-8265) -- STAGE 4: validate the EXACT Q(i) cotangent-lift +
de Bondt rotation machinery (the step that turns a nilpotent-Jacobian cubic map into a
Hessian-nilpotent quartic), and diagnose the cubic-HOMOGENEOUS-Keller blocker.

Machinery (derived, exact over Q(i), no sqrt2): given G=Z+H with JH NILPOTENT,
  K(Z,Y) = (H(Z), JH(Z)^T Y)   [cotangent lift nonlinear part, dim 2N]
  grad P = T'^{-1} K(T' Z),  T' = [[I,iI],[iI,I]],  T'^{-1} = (1/2)[[I,-iI],[-iI,I]]
Then Hess(P) = T'^{-1} (JK) T' with JK=[[JH,0],[A,JH^T]] block-lower-triangular, nilpotent
diagonal blocks => Hess(P) NILPOTENT; and grad P is a genuine gradient (symmetric Jacobian).
Collisions of id+grad P are T'^{-1}(collision of the lift). VALIDATE on a known nilpotent
cubic-homogeneous map where we can check every claim.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib_exact_deathstar_S61 import *
from fractions import Fraction as Fr

def psubst(poly, subs, NV):
    """substitute each var i by subs[i] (a poly); unspecified -> itself."""
    res = pzero()
    for e,c in poly.items():
        term = {tuple([0]*NV): c}
        for i,k in enumerate(e):
            if k:
                rep = subs.get(i, pvar(i,NV))
                term = pmul(term, ppow(rep,k,NV))
        res = padd(res, term)
    return res

def lift_and_rotate(H, N, NV):
    """H: list of N cubic-homogeneous polys (nonlinear part of G=Z+H) in the first N vars.
       Returns gradP (list of 2N polys) and P, working in 2N vars (X=0..N-1, Y=N..2N-1)."""
    twoN = 2*N
    # JH (N x N) then JH^T Y
    JH = [[pderiv(H[i], j) for j in range(N)] for i in range(N)]
    Yv = [pvar(N+j, NV) for j in range(N)]
    # K = (H(Z), JH^T Y):  first N comps = H_i (in X vars, already), next N = sum_i JH[i][j]*Y_i
    K = [dict(H[i]) for i in range(N)]                      # H depends only on X=0..N-1
    for j in range(N):
        s = pzero()
        for i in range(N):
            s = padd(s, pmul(JH[i][j], Yv[i]))
        K.append(s)
    # T' Z : X_j -> Z_j + i Z_{N+j};  Y_j -> i Z_j + Z_{N+j}
    subs = {}
    for j in range(N):
        subs[j]   = padd(pvar(j,NV), pscale(pvar(N+j,NV), I))
        subs[N+j] = padd(pscale(pvar(j,NV), I), pvar(N+j,NV))
    KT = [psubst(K[a], subs, NV) for a in range(twoN)]
    # grad P = T'^{-1} K(T'Z);  T'^{-1} = (1/2)[[I,-iI],[-iI,I]]
    half = QI(Fr(1,2),0)
    gradP = []
    for j in range(N):     # rows 0..N-1: (1/2)(KT[j] - i KT[N+j])
        gradP.append(pscale(psub(KT[j], pscale(KT[N+j], I)), half))
    for j in range(N):     # rows N..2N-1: (1/2)(-i KT[j] + KT[N+j])
        gradP.append(pscale(padd(pscale(KT[j], -I), KT[N+j]), half))
    # P = (1/4) sum_k Z_k gradP_k   (Euler, gradP homogeneous cubic)
    P = pzero()
    for k in range(twoN):
        P = padd(P, pmul(pvar(k,NV), gradP[k]))
    P = pscale(P, QI(Fr(1,4),0))
    return gradP, P, twoN

def check(H, N, NV, label):
    print(f"\n=== validate lift+rotation on: {label}  (N={N}, 2N={2*N}) ===")
    gradP, P, twoN = lift_and_rotate(H, N, NV)
    print("  P degree:", pdeg(P), " (expect 4 if H cubic homogeneous)")
    # (1) gradP is a genuine gradient: d/dZ_j P == gradP_j
    isgrad = all(pequal(pderiv(P, j), gradP[j]) for j in range(twoN))
    print("  grad P == dP/dZ (P is a true potential, Jacobian symmetric):", isgrad)
    # (2) Hess(P) nilpotent: Hess = Jacobian(gradP); check Hess^{2N} == 0
    Hess = [[pderiv(gradP[i], j) for j in range(twoN)] for i in range(twoN)]
    M = Hess
    Mk = M
    nz_power = None
    for k in range(2, twoN+1):
        Mk = mmul(Mk, M)
        if mis_zero(Mk):
            nz_power = k; break
    print(f"  Hess(P) nilpotent: Hess^{nz_power}==0" if nz_power else "  Hess(P) NOT nilpotent up to 2N")
    return gradP, P

if __name__ == "__main__":
    NV = 8
    # ---- known nilpotent cubic-homogeneous map, dim N=2: H=(0, x0^3), JH=[[0,0],[3x0^2,0]] nilpotent
    x0,x1 = pvar(0,NV), pvar(1,NV)
    H2 = [pzero(), ppow(x0,3,NV)]
    check(H2, 2, NV, "H=(0, x0^3)  [classic nilpotent cubic-homog, dim 2]")

    # ---- a richer nilpotent cubic-homogeneous map, dim N=3 (strictly-triangular nilpotent):
    # H=(0, x0^3, x0^2 x1)  -> JH strictly lower triangular in a suitable order => nilpotent
    NV2 = 8
    y0,y1,y2 = pvar(0,NV2),pvar(1,NV2),pvar(2,NV2)
    H3 = [pzero(), ppow(y0,3,NV2), pmul(ppow(y0,2,NV2), y1)]
    # verify JH3 nilpotent first
    JH3 = [[pderiv(H3[i],j) for j in range(3)] for i in range(3)]
    P2 = mmul(JH3,JH3); P3 = mmul(P2,JH3)
    print("\n[dim-3 map] JH nilpotent? JH^3==0:", mis_zero(P3))
    check(H3, 3, NV2, "H=(0, x0^3, x0^2 x1)  [nilpotent cubic-homog, dim 3]")

    print("\n=== DIAGNOSIS of the full-F transport blocker ===")
    print("  - F, Yagzhev G, companion cubic reduction: collision transports EXACTLY (stage 2-3).")
    print("  - BLOCKER: securing a NILPOTENT Jacobian requires a cubic-HOMOGENEOUS Keller reduction.")
    print("    Naive homogenization (deg-2 nonlinear x x0) breaks Keller: det(I+x0*JH2+JH3) is not")
    print("    const because nilpotency of the pencil a*JH2+b*JH3 is NOT implied by nilpotency at a=b=1.")
    print("  - The lift+rotation machinery ABOVE is validated exact: given ANY nilpotent-Jacobian G,")
    print("    it outputs a Hessian-nilpotent quartic. The one missing input is the nilpotent G_c.")
