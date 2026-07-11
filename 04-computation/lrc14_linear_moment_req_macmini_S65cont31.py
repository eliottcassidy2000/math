#!/usr/bin/env python3
"""THM-705: the linear moment requirement for the deg-2 rows (k=9,10).
(a) OPTIMAL quadratic majorant by exact vertex enumeration: q = a + b N + c N(N-1) >= g
on {0..7}; at moment point (m1,m2) the optimal value = min over feasible (a,b,c) of
a + b m1 + c m2 -- LP; vertices = triples of active constraints. We enumerate all C(8,3)
triples, solve 3x3 exactly, keep feasible ones; the bound is the min over vertices
(valid for ANY (m1,m2) with a genuine N-distribution).
(b) THE LINEAR REQUIREMENT per row: with the optimal vertex (a*,b*,c*) at consec,
row k closes iff  a* + b* m1(F) + c* m2(F) <= cap_{k+1} for all k-cores.
(c) uniform margin sweep: exact (m1,m2) over structured + random cores."""
from fractions import Fraction as F
from itertools import combinations
import random
CAP={9:F(1979,4004),10:F(55,91),11:F(66,91)}
g=[F(1),F(1,3)]+[F(0)]*6
def pvec(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(p for p in pts if 0<=p<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
def moments(E):
    p=pvec(E); m1=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    Phi=p[0]+p[1]/3; return m1,m2,Phi
def vertices():
    V=[]
    for trip in combinations(range(8),3):
        # solve a + b N + c N(N-1) = g(N) at the triple
        (x,y,z)=trip
        import fractions
        # 3x3 solve by hand (Cramer) over Fractions
        M=[[F(1),F(x),F(x*(x-1))],[F(1),F(y),F(y*(y-1))],[F(1),F(z),F(z*(z-1))]]
        R=[g[x],g[y],g[z]]
        det=lambda A: A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
        D=det(M)
        if D==0: continue
        import copy
        sol=[]
        for c in range(3):
            Mc=[row[:] for row in M]
            for r in range(3): Mc[r][c]=R[r]
            sol.append(det(Mc)/D)
        a,b,c=sol
        if all(a+b*N+c*N*(N-1)>=g[N]-F(0) for N in range(8)):
            V.append((a,b,c,trip))
    return V
VS=vertices()
print(f"feasible majorant vertices: {len(VS)}")
def best_bound(m1,m2):
    return min(a+b*m1+c*m2 for a,b,c,_ in VS)
print()
print("row requirements (optimal deg-2 bound at exact moments) vs caps:")
for k,Es in [(9,list(range(9))),(10,list(range(10)))]:
    m1,m2,Phi=moments(Es)
    bb=best_bound(m1,m2)
    # which vertex is active
    a,b,c,trip=min(VS,key=lambda v: v[0]+v[1]*m1+v[2]*m2)
    print(f"  k={k} consec: Phi={float(Phi):.4f} optimal-deg2={float(bb):.4f} cap={float(CAP[k+1]):.4f} margin={float(CAP[k+1]-bb):+.4f} active-vertex g-touch at N={trip} coeffs=({a},{b},{c})")
print()
print("uniform margin sweep (worst optimal-deg2 bound over adversarial cores):")
random.seed(5)
for k in [9,10]:
    worst=None
    fams=[list(range(k)),list(range(1,k+1)),[0,1,2,3,4,5,6,7,k+5][:k]]
    for _ in range(30):
        fams.append(sorted(random.sample(range(0,26),k)))
    for Es in fams:
        m1,m2,_=moments(Es)
        bb=best_bound(m1,m2)
        if worst is None or bb>worst[0]: worst=(bb,Es)
    print(f"  k={k}: WORST deg-2 bound = {float(worst[0]):.4f} at {worst[1]} vs cap {float(CAP[k+1]):.4f} -> uniform margin {float(CAP[k+1]-worst[0]):+.4f}")
