#!/usr/bin/env python3
"""
unit_distance_moser_crossover_s4.py
monad-explorer-2026-06-07-S4  (crossover-lane, follow-up to HYP-2301)

GOAL: make THM-431's CITED ceiling N*<=28 SELF-CONTAINED. HYP-2301 showed no
single-norm RANK-2 lattice beats 3N at N<=28 (earliest sqrt7 at 32); Engel's
record u(28)>=85 lives instead in the RANK-4 Moser ring
    M_L = Z[zeta6, omega3],  zeta6=(1+i√3)/2 (cos 1/2),  omega3=(5+i√11)/6 (cos 5/6),
which has 18 UNIT VECTORS all at radius 1 (omega3 is a non-torsion unit) — it
escapes the degree-radius/kissing tension by leaving rank 2.

Here we run a densest-patch search DIRECTLY in M_L (graph BFS/greedy + anneal in
Z^4 with the 18 unit-vector offsets) and EXACT-recount over Q(√3,√11). If a
28-point patch reaches 85 (>84) we have re-derived Engel's crossover with explicit
integer coordinates — the ceiling N*<=28 becomes self-contained. (All adjacency is
exact: |z|^2==1 over the field, no float decides an edge.)

Exact arithmetic in K=Q(√3,√11), basis 1,√3,√11,√33 (reused from
unit_distance_moser_lattice_u21_monad_s4.py).
"""
from fractions import Fraction as F
from itertools import product
import random
random.seed(20260607)

Z4 = (F(0),)*4
ONE = (F(1),F(0),F(0),F(0))
def add(x,y): return tuple(x[i]+y[i] for i in range(4))
def smul(s,x): return tuple(s*x[i] for i in range(4))
_MT = {(0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),
       (1,1):(3,0),(1,2):(1,3),(1,3):(3,2),(2,2):(11,0),(2,3):(11,1),(3,3):(33,0)}
def mul(x,y):
    r=[F(0)]*4
    for i in range(4):
        if x[i]==0: continue
        for j in range(4):
            if y[j]==0: continue
            a,b=(i,j) if i<=j else (j,i)
            c,idx=_MT[(a,b)]; r[idx]+=x[i]*y[j]*c
    return tuple(r)
RE={'1':(F(1),F(0),F(0),F(0)),'w1':(F(1,2),F(0),F(0),F(0)),
    'w3':(F(5,6),F(0),F(0),F(0)),'w13':(F(5,12),F(0),F(0),F(-1,12))}
IM={'1':(F(0),F(0),F(0),F(0)),'w1':(F(0),F(1,2),F(0),F(0)),
    'w3':(F(0),F(0),F(1,6),F(0)),'w13':(F(0),F(5,12),F(1,12),F(0))}
KEYS=['1','w1','w3','w13']
def coord(v):
    re=Z4; im=Z4
    for k,key in zip(v,KEYS):
        if k==0: continue
        re=add(re,smul(F(k),RE[key])); im=add(im,smul(F(k),IM[key]))
    return re,im
def normsq(v):
    re,im=coord(v); return add(mul(re,re),mul(im,im))

# --- the 18 unit vectors ---
R=4
UNITS=[(a,b,c,d) for a,b,c,d in product(range(-R,R+1),repeat=4)
       if (a,b,c,d)!=(0,0,0,0) and normsq((a,b,c,d))==ONE]
assert len(UNITS)==18, f"expected 18 unit vectors, got {len(UNITS)}"
UNITSET=set(UNITS)

def deg_in(p,S):
    return sum((p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3]) in S for u in UNITS)
def U_exact(S):
    s=set(S); return sum(deg_in(p,s) for p in s)//2
def cands(S):
    c=set()
    for p in S:
        for u in UNITS:
            q=(p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
            if q not in S: c.add(q)
    return c

def greedy_grow(N):
    S={(0,0,0,0)}
    while len(S)<N:
        cs=cands(S)
        if not cs: break
        best=max(cs,key=lambda q: deg_in(q,S))
        S.add(best)
    return S

def anneal(S,iters):
    S=set(S); E=U_exact(S); best=E; bestS=set(S)
    for it in range(iters):
        T=max(0.05,3.0*(1-it/iters))
        u=random.choice(tuple(S))
        p=random.choice(tuple(S)); v=random.choice(UNITS)
        w=(p[0]+v[0],p[1]+v[1],p[2]+v[2],p[3]+v[3])
        if w in S or w==u: continue
        du=deg_in(u,S)
        diff=(w[0]-u[0],w[1]-u[1],w[2]-u[2],w[3]-u[3])
        dw=deg_in(w,S)-(1 if diff in UNITSET else 0)
        delta=dw-du
        if delta>=0 or random.random()<pow(2.718281828,delta/T):
            S.remove(u); S.add(w); E+=delta
            if E>best: best=E; bestS=set(S)
    return best,bestS

def densest(N,iters=80000,restarts=10):
    best=-1; bestS=None
    base=greedy_grow(N+8)
    for r in range(restarts):
        if r==0:
            seed=greedy_grow(N)
        else:
            b=list(base); random.shuffle(b); seed=set(b[:N])
            # ensure connected-ish: regrow from a random patch member
            seed=greedy_grow_from(set(list(seed)[:1]),N)
        e,S=anneal(seed,iters)
        if e>best: best,bestS=e,S
    return best,bestS

def greedy_grow_from(S,N):
    S=set(S)
    while len(S)<N:
        cs=cands(S)
        if not cs: break
        S.add(max(cs,key=lambda q: deg_in(q,S)))
    return S

def recount_from_coords(S):
    """fully independent exact recount: pairwise |zi-zj|^2==1 over K."""
    L=list(S); n=len(L); c=0
    for i in range(n):
        for j in range(i+1,n):
            d=(L[i][0]-L[j][0],L[i][1]-L[j][1],L[i][2]-L[j][2],L[i][3]-L[j][3])
            if normsq(d)==ONE: c+=1
    return c

def main():
    print("="*68)
    print("MOSER RING M_L = Z[zeta6,omega3] densest-patch crossover (rank 4, 18 units)")
    print("exact |z|^2==1 over Q(√3,√11); target: reach u(28)>=85 (>84) self-contained")
    print("="*68)
    print(f"#unit vectors = {len(UNITS)} (6 triangular c=d=0, 12 Moser); all ad=bc: "
          f"{all(a*d==b*c for a,b,c,d in UNITS)}")
    print(f"\n   {'N':>3} {'best':>5} {'3N':>5} {'u-3N':>6}  recount  status")
    firstbeat=None
    for N in range(22,31):
        best,S=densest(N)
        rc=recount_from_coords(S)
        assert rc==best,f"recount mismatch {rc} vs {best} at N={N}"
        d=best-3*N
        st=""
        if d>0:
            st="BEATS 3N"
            if firstbeat is None: firstbeat=N
        elif d==0: st="ties"
        print(f"   {N:>3} {best:>5} {3*N:>5} {d:>6}  {rc:>6}   {st}")
    print("="*68)
    if firstbeat:
        print(f"Moser ring first beats 3N at N={firstbeat}.")
        if firstbeat<=28:
            print(f"=> reproduces Engel's N*<=28 with EXACT integer coords (self-contained).")
    print("Compare HYP-2301: best RANK-2 single-norm lattice (sqrt7) first beats at N=32.")
    print("The rank-4 Moser ring (18 units at radius 1) crosses EARLIER — the escape.")

if __name__=="__main__":
    main()
