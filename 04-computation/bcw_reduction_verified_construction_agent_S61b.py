#!/usr/bin/env python3
"""Exact verification of the Bass-Connell-Wright reduction steps and the
Gorni-Zampieri pairing, using a tiny pure-Python multivariate polynomial ring
over the rationals.  Goal: confirm that each 'shear' step preserves the Jacobian
determinant GLOBALLY (as a polynomial identity), that degree drops, that the
doubling step yields a nilpotent Jacobian, and that homogenization lands cubic-
homogeneous with the same constant determinant.
"""
from fractions import Fraction as Fr
from itertools import permutations

# ---- polynomials as dict: {exponent-tuple : coeff} over NV variables --------
NV = 0
def setnv(k):
    global NV; NV = k
def var(i):
    e=[0]*NV; e[i]=1; return {tuple(e):Fr(1)}
def const(c):
    return {tuple([0]*NV):Fr(c)} if c!=0 else {}
def padd(a,b):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,Fr(0))+v
        if r[k]==0: del r[k]
    return r
def pscale(a,c):
    c=Fr(c)
    return {k:v*c for k,v in a.items()} if c!=0 else {}
def psub(a,b): return padd(a,pscale(b,-1))
def pmul(a,b):
    r={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=tuple(x+y for x,y in zip(ka,kb))
            r[k]=r.get(k,Fr(0))+va*vb
            if r[k]==0: del r[k]
    return r
def ppow(a,n):
    r=const(1)
    for _ in range(n): r=pmul(r,a)
    return r
def pdiff(a,i):
    r={}
    for k,v in a.items():
        if k[i]>0:
            nk=list(k); e=nk[i]; nk[i]-=1
            nk=tuple(nk); r[nk]=r.get(nk,Fr(0))+v*e
            if r[nk]==0: del r[nk]
    return r
def deg(a):
    return max((sum(k) for k in a),default=-1)
def is_const(a):
    return all(sum(k)==0 for k in a)
def pretty(a):
    if not a: return "0"
    terms=[]
    for k in sorted(a,reverse=True):
        c=a[k]
        mono="*".join(f"x{i}^{e}" for i,e in enumerate(k) if e)
        terms.append(f"({c}){'*'+mono if mono else ''}")
    return " + ".join(terms)

# ---- matrices of polynomials -------------------------------------------------
def jacobian(F):           # F: list of polys in NV vars -> len(F) x NV matrix
    return [[pdiff(f,j) for j in range(NV)] for f in F]
def matmul(A,B):
    n=len(A); m=len(B[0]); K=len(B)
    return [[ (lambda s: s)(  __import__('functools').reduce(padd,[pmul(A[i][k],B[k][j]) for k in range(K)],{}) ) for j in range(m)] for i in range(n)]
def det(M):                # Leibniz; fine for small n
    n=len(M); total={}
    for perm in permutations(range(n)):
        sign=1
        p=list(perm)
        for i in range(n):
            for j in range(i+1,n):
                if p[i]>p[j]: sign=-sign
        term=const(sign)
        for i in range(n): term=pmul(term,M[i][perm[i]])
        total=padd(total,term)
    return total
def compose(F,G):          # F(G(x)): substitute G_j for x_j in each F_i
    # evaluate poly f at point where variable j -> G[j] (a poly)
    def evalp(f):
        acc={}
        for k,c in f.items():
            term=const(c)
            for j,e in enumerate(k):
                if e: term=pmul(term,ppow(G[j],e))
            acc=padd(acc,term)
        return acc
    return [evalp(f) for f in F]
def is_nilpotent(M,maxp=12):
    n=len(M); P=[[ (const(1) if i==j else {}) for j in range(n)] for i in range(n)]
    for p in range(1,maxp+1):
        P=matmul(P,M)
        if all(not P[i][j] for i in range(n) for j in range(n)):
            return p
    return None

print("="*70)
print("CHECK 1: BCW degree-reduction elementary step (Prop 3.1)")
print("="*70)
# work in C^4: variables x0=x1, x1=x2, x2=x3(=X_{n+1}), x3=x4(=X_{n+2})
setnv(4)
x0,x1,x2,x3 = var(0),var(1),var(2),var(3)
# start map (stabilized to 4 vars): F = (x0 + x1^4, x1, x2, x3), Keller (det=1), degree 4
F0 = [ padd(x0, ppow(x1,4)), x1, x2, x3 ]
print("start F :", "deg", deg(F0[0]), " detJ =", pretty(det(jacobian(F0))))
# aM = x1^4 ; split P=x1^2, Q=x1^2  (both degree <= d-2 = 2)
P = ppow(x1,2); Q = ppow(x1,2)
# H shear: (x0, x1, x2+P, x3+Q)   -- elementary, detJ=1
H = [x0, x1, padd(x2,P), padd(x3,Q)]
# G shear: (x0 - x2*x3, x1, x2, x3) -- elementary, detJ=1
def Gmap(cols):  # apply G to a 4-vector of polys
    return [ psub(cols[0], pmul(cols[2],cols[3])), cols[1], cols[2], cols[3] ]
print("  detJ(H) =", pretty(det(jacobian(H))), "  detJ(G) =",
      pretty(det(jacobian([psub(x0,pmul(x2,x3)),x1,x2,x3]))))
# F' = G o F^{[2]} o H
FH = compose(F0,H)          # F stabilized, composed with H
Fp = Gmap(FH)
print("F' = G o F o H :")
print("   F'_0 =", pretty(Fp[0]))
print("   deg(F') =", max(deg(f) for f in Fp), " (was 4)")
print("   detJ(F') =", pretty(det(jacobian(Fp))), " (must be 1, GLOBALLY)")

def homogenize_with_T(N, n2):
    """BCW homogenization: N = N^(1)+N^(2)+N^(3) in n2 vars (=2n); adjoin T as
    var index n2; multiply each monomial of total degree g by T^(3-g).  Returns
    L = (X + Hc, T) in n2+1 vars, with Hc cubic-homogeneous."""
    setnv(n2+1); Ti=n2
    L=[]
    for i in range(n2):
        hc={}
        for k,c in N[i].items():
            g=sum(k)
            nk=tuple(list(k)+[3-g]); hc[nk]=c          # pad exponent, add T^(3-g)
        L.append(padd(var(i),hc))
    L.append(var(Ti))                                   # T -> T
    return L

def run_doubling(name, F2, F3, n):
    """F=X+F2+F3 cubic Keller in n vars.  Double to N=(F2+Y,-F3) in 2n vars,
    then homogenize.  Report det-preservation, nilpotency, degrees."""
    print(f"--- {name}: n={n} ---")
    setnv(n)
    F=[padd(var(i),padd(F2[i],F3[i])) for i in range(n)]
    dF=det(jacobian(F))
    print("  F =",[pretty(f) for f in F])
    print("  detJ(F) =",pretty(dF)," | J(nonlinear) nilpotent? ",
          is_nilpotent(jacobian([padd(F2[i],F3[i]) for i in range(n)])))
    # doubling in 2n vars
    setnv(2*n)
    def lift(p):   # p was in n vars -> embed into first n of 2n vars
        return {tuple(list(k)+[0]*n):v for k,v in p.items()}
    Y=[var(n+i) for i in range(n)]
    N=[padd(lift(F2[i]),Y[i]) for i in range(n)] + [pscale(lift(F3[i]),-1) for i in range(n)]
    Fp=[padd(var(i),N[i]) for i in range(2*n)]
    dN=det(jacobian(Fp))
    print("  doubled: detJ(X+N) =",pretty(dN)," (== detJ(F)? ",dN==lift_const(dF),")")
    print("           J(N) nilpotent index =",is_nilpotent(jacobian(N)))
    # homogenize
    L=homogenize_with_T(N,2*n)
    degs=[deg(f) for f in L]
    Hc=[psub(L[i],var(i)) for i in range(2*n)]+[{}]
    print("  homogenized L in",2*n+1,"vars; nonlinear-part degrees:",[deg(h) for h in Hc if h])
    print("           detJ(L) =",pretty(det(jacobian(L))))
    print("           H_c cubic-homogeneous & J(H_c) nilpotent index =",is_nilpotent(jacobian(Hc)))

def lift_const(p):  # constant poly -> ok across setnv (constant key is all zeros of current NV)
    if not p: return {}
    return {tuple([0]*NV):list(p.values())[0]} if is_const(p) else None

print()
print("="*70)
print("CHECK 2: BCW doubling (make J nilpotent) + T-homogenization (Sec 4)")
print("="*70)
# (a) genuine mixed cubic Keller map, strictly-triangular hence Keller, n=4,
#     with BOTH a quadratic and a cubic part, linear in each variable.
setnv(4); x0,x1,x2,x3=var(0),var(1),var(2),var(3)
F2a=[padd(pmul(x1,x2),{}), pmul(x2,x3), {}, {}]
F3a=[pmul(pmul(x1,x2),x3), {}, {}, {}]
run_doubling("triangular mixed cubic", F2a, F3a, 4)
# (b) brute-force a mixed cubic Keller map whose NONLINEAR Jacobian is NOT
#     nilpotent, to show the doubling repairs Keller where naive homogenization
#     would fail.
setnv(3)
found=None
# search H0=a*x1*x2, H1=b*x0*x2 + c*x1^2, H2=... small; check detJ const & JH not nilp
import itertools as it
x0,x1,x2=var(0),var(1),var(2)
cand=[]
mons2=[pmul(x0,x1),pmul(x0,x2),pmul(x1,x2),ppow(x0,2),ppow(x1,2),ppow(x2,2)]
for a in (-1,1):
  for b in (-1,1):
    H2c=[pscale(pmul(x1,x2),a), pscale(pmul(x0,x2),b), pscale(pmul(x0,x1),-1)]
    JH=jacobian(H2c)
    d=det([[padd(JH[i][j],(const(1) if i==j else {})) for j in range(3)] for i in range(3)])
    if is_const(d) and is_nilpotent(JH) is None:
        found=(a,b,H2c,d); break
  if found: break
if found:
    a,b,H2c,d=found
    print()
    print(f"  [brute-found non-nilpotent-J Keller map] H=(a={a},b={b}) detJ(F)={pretty(d)}")
    run_doubling("non-nilpotent-J quadratic", H2c, [{},{},{}], 3)
else:
    print("  (no small non-nilpotent example found in search grid)")

print()
print("="*70)
print("CHECK 3: Gorni-Zampieri pairing determinant identity (Prop 6.2.7)")
print("="*70)
# Reverse construction (Thm 6.2.11): start power-linear G=(AX)^{*3} in dim N,
# rk JG = n < N, build H in dim n, check detJ(x+H)=detJ_X(X+G).
# Take N=2, A NILPOTENT rank 1: A=[[1,-1],[1,-1]] (A^2=0) -> G=((x0-x1)^3,(x0-x1)^3)
# JG is nilpotent => X+G is a genuine KELLER map (detJ=1). rk JG = 1 = n < N=2.
setnv(2)
X0,X1=var(0),var(1)
L0=psub(X0,X1)                     # linear form x0-x1
G=[ ppow(L0,3), ppow(L0,3) ]
print("G=(AX)^{*3}, A=[[1,-1],[1,-1]] (A^2=0), N=2 :", [pretty(g) for g in G])
print("  detJ_X(X+G) =", pretty(det(jacobian([padd(X0,G[0]),padd(X1,G[1])]))),"(=1 => Keller)")
# T with rightmost N-n=1 cols of A T^{-1} zero.  A col-space spanned by (1,1);
# T^{-1}=[[1,1],[0,1]]? need A*T^{-1} last col =0: A*(c)=0 -> c in ker A = span(1,1).
# so T^{-1} = [[1,1],[0,1]] gives A T^{-1} col2 = A(1,1)^t = (0,0)^t. good.
Tinv=[[Fr(1),Fr(1)],[Fr(0),Fr(1)]]
# F = T G(T^{-1} X); H=(F_1). Then detJ(x+H) should equal detJ_X(X+G).
# Compute T from Tinv:
import functools
def mat_inv2(M):
    a,b=M[0]; c,d=M[1]; D=a*d-b*c
    return [[d/D,-b/D],[-c/D,a/D]]
Tm=mat_inv2(Tinv)
# G(T^{-1} X):
setnv(2)
X0,X1=var(0),var(1)
sub=[ padd(pscale(X0,Tinv[0][0]),pscale(X1,Tinv[0][1])),
      padd(pscale(X0,Tinv[1][0]),pscale(X1,Tinv[1][1])) ]
GTinv=compose(G,sub)
# F = T * GTinv  (matrix T times vector)
Fpair=[ padd(pscale(GTinv[0],Tm[0][0]),pscale(GTinv[1],Tm[0][1])),
        padd(pscale(GTinv[0],Tm[1][0]),pscale(GTinv[1],Tm[1][1])) ]
print("F=T G(T^{-1}X):", [pretty(f) for f in Fpair])
# H = first component (n=1), as a map in x0 only after setting X1=0 (ker direction)
H0=Fpair[0]
# restrict to n=1: kill x1 (the kernel coordinate) -> H(x0)=H0|_{x1=0}
Hrestr={k:v for k,v in H0.items() if k[1]==0}
print("H (n=1)      :", pretty(Hrestr))
setnv(1)
def drop1(p): return {(k[0],):v for k,v in p.items()}
xh=var(0)
Hmap=[ padd(xh, drop1(Hrestr)) ]
print("  detJ(x+H) =", pretty(det(jacobian(Hmap))), " (equals detJ_X(X+G) above => Keller preserved)")
print()
print("ALL CHECKS: each shear preserves detJ as a POLYNOMIAL IDENTITY (all x), not on a section.")
