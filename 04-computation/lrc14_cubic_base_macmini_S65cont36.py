#!/usr/bin/env python3
"""THM-712: the k=8 cubic base form. Optimal degree-3 majorant of g on {0..7} by exact
vertex enumeration (4 active constraints, C(8,4)=70 candidates, 4x4 Cramer in Fractions);
the k=8 requirement in product form; adversarial hunt for unconditional validity."""
from fractions import Fraction as F
from itertools import combinations
import random
cap9=F(1979,4004)
g=[F(1),F(1,3)]+[F(0)]*6
def det4(M):
    from itertools import permutations
    s=F(0)
    for perm in permutations(range(4)):
        sign=1; seen=list(perm)
        # parity
        inv=sum(1 for i in range(4) for j in range(i+1,4) if perm[i]>perm[j])
        sign=(-1)**inv
        prod=F(1)
        for r in range(4): prod*=M[r][perm[r]]
        s+=sign*prod
    return s
def basis(N): return [F(1),F(N),F(N*(N-1)),F(N*(N-1)*(N-2))]
VS=[]
for quad in combinations(range(8),4):
    M=[basis(N) for N in quad]; R=[g[N] for N in quad]
    D=det4(M)
    if D==0: continue
    sol=[]
    for c in range(4):
        Mc=[row[:] for row in M]
        for r in range(4): Mc[r][c]=R[r]
        sol.append(det4(Mc)/D)
    a,b,c3,d=sol
    if all(a+b*N+c3*N*(N-1)+d*N*(N-1)*(N-2) >= g[N] for N in range(8)):
        VS.append((a,b,c3,d,quad))
print(f"feasible deg-3 vertices: {len(VS)}")
def pvec(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
def falling(j,r):
    v=1
    for i in range(r): v*=(j-i)
    return v
def mom(E):
    p=pvec(E); return [sum(falling(j,r)*p[j] for j in range(8)) for r in range(4)]
def bound(E):
    m=mom(E)
    return min(a+b*m[1]+c3*m[2]+d*m[3] for a,b,c3,d,_ in VS)
E8=list(range(8))
m=mom(E8)
vals=[(a+b*m[1]+c3*m[2]+d*m[3],(a,b,c3,d,quad)) for a,b,c3,d,quad in VS]
best=min(vals)
a,b,c3,d,quad=best[1]
print(f"consec_8: optimal deg-3 bound = {best[0]} = {float(best[0]):.4f} vs cap_9 = {float(cap9):.4f} margin {float(cap9-best[0]):+.4f}")
print(f"active vertex touches N={quad}, coeffs a={a} b={b} c={c3} d={d}")
# requirement form: a + b m1 + c m2 + d m3 <= cap_9 <=> -(b m1 + c m2 + d m3) >= a - cap_9
print()
print("adversarial hunt at k=8 (worst deg-3 bound):")
random.seed(4)
worst=(None,None)
fams=[list(range(8)),list(range(1,9)),[1,2,3,4,5,6,7,18],[1,2,3,4,5,6,7,26],[1,8,15,22,29,36,43,50]]
for _ in range(35): fams.append(sorted(random.sample(range(0,40),8)))
for _ in range(10):
    base=sorted(random.sample(range(0,11),7)); fams.append(base+[random.randint(60,300)])
for E in fams:
    bb=bound(E)
    if worst[0] is None or bb>worst[0]: worst=(bb,E)
print(f"  {len(fams)} families: WORST = {float(worst[0]):.4f} at {worst[1]}")
print(f"  UNCONDITIONAL? worst <= cap_9: {'YES margin '+format(float(cap9-worst[0]),'+.4f') if worst[0]<=cap9 else '*** NO ***'}")
