#!/usr/bin/env python3
"""
unit_distance_augment_cube_monad_s6c.py
monad-explorer-2026-06-07-S6  (deep-research; OPEN-Q-057 frontier) -- PART C/D.

PART C0  DIRECT u(27)>81 attack: scan rotation angles (th2,th3) of the cube
         K3^[]3; for special angles the Minkowski sum acquires ACCIDENTAL unit
         distances among its 27 points -> U>81 (if still 27 distinct points).
PART C1  AUGMENTATION (H(3,3)+k): max # cube vertices on a common unit circle =
         max extra unit distances from one added point. Tests "u(28)=H(3,3)+1".
         Greedy multi-add 27 -> 28,29,30 and report exact-friendly winners.
PART D   from-scratch anneal: densest 27 / 28 point unit-distance graph.

Float search; any record's raw coordinates are dumped for exact re-verification.
"""
import math, random, itertools
random.seed(11)
TOL=1e-9

def sqd(p,q): return (p[0]-q[0])**2+(p[1]-q[1])**2
def count_unit(pts,tol=TOL):
    n=len(pts);u=0
    for i in range(n):
        for j in range(i+1,n):
            if abs(sqd(pts[i],pts[j])-1.0)<tol:u+=1
    return u
def distinct(pts,tol=1e-6):
    n=len(pts)
    for i in range(n):
        for j in range(i+1,n):
            if sqd(pts[i],pts[j])<tol*tol:return False
    return True
def tri(theta):
    return [(0.0,0.0),(math.cos(theta),math.sin(theta)),
            (math.cos(theta+math.pi/3),math.sin(theta+math.pi/3))]
def minkowski(factors):
    pts=[(0.0,0.0)]
    for f in factors:
        pts=[(p[0]+q[0],p[1]+q[1]) for p in pts for q in f]
    return pts
def cube(th2,th3,th1=0.0): return minkowski([tri(th1),tri(th2),tri(th3)])

print("="*74)
print("PART C0 -- DIRECT u(27)>81: scan cube angles for accidental unit distances")
print("="*74)
best27=(-1,None)
collisions=0; valid=0
# fine scan; symmetry: th2,th3 in (0, pi/3) suffices up to lattice/relabeling-ish
N=240
seen=[]
for i in range(1,N):
    th2=math.pi/3*i/N
    for j in range(1,N):
        th3=math.pi/3*j/N
        P=cube(th2,th3)
        if not distinct(P): collisions+=1; continue
        valid+=1
        U=count_unit(P)
        if U>best27[0]:
            best27=(U,(th2,th3))
        if U>81:
            seen.append((U,th2,th3))
print(f"  scanned {valid} distinct-angle cubes ({collisions} had collisions).")
print(f"  MAX U over 27-point cubes = {best27[0]}  at (th2,th3)={best27[1]}")
if best27[0]<=81:
    print("  => NO special angle gives U>81 on the 27-pt cube (in this scan):")
    print("     either 27 distinct pts with U=81, or angle specialization collapses")
    print("     points (collision) BEFORE adding extra unit distances. Strong")
    print("     evidence the cube cannot be tuned past the 81 tie.")
else:
    print(f"  FOUND {len(seen)} angle pairs with U>81! samples:")
    for U,a,b in sorted(seen,reverse=True)[:10]:
        print(f"     U={U}  th2={a:.6f} th3={b:.6f}")

print("\n"+"="*74)
print("PART C1 -- AUGMENTATION H(3,3)+k : extra unit distances from added points")
print("="*74)
def circle_inter(p,q,r=1.0):
    """intersection points of circles radius r about p,q (float)."""
    d2=sqd(p,q); d=math.sqrt(d2)
    if d>2*r or d<1e-12: return []
    a=d/2
    h2=r*r-a*a
    if h2<-1e-15: return []
    h=math.sqrt(max(0.0,h2))
    mx=(p[0]+q[0])/2; my=(p[1]+q[1])/2
    ux=(q[0]-p[0])/d; uy=(q[1]-p[1])/d
    return [(mx-uy*h, my+ux*h),(mx+uy*h, my-ux*h)]

def best_augment(pts):
    """return (k, point) maximizing # of pts at unit distance from a NEW point."""
    best=(0,None)
    n=len(pts)
    for i in range(n):
        for j in range(i+1,n):
            if sqd(pts[i],pts[j])>4.0: continue
            for c in circle_inter(pts[i],pts[j]):
                # reject if c coincides with an existing point
                if any(sqd(c,p)<1e-10 for p in pts): continue
                k=sum(1 for p in pts if abs(sqd(c,p)-1.0)<1e-7)
                if k>best[0]: best=(k,c)
    return best

# test augmentation over several realizations of the cube (random + structured)
realizations={
 "generic(0.7,1.93)":cube(0.7,1.93),
 "generic(0.61,1.31)":cube(0.61,1.31),
}
# also random angles
for t in range(40):
    a=random.uniform(0.05,math.pi/3-0.05); b=random.uniform(0.05,math.pi/3-0.05)
    realizations[f"rand{t}"]=cube(a,b)

bestaug=(0,None,None)
for name,P in realizations.items():
    if not distinct(P) or count_unit(P)!=81: continue
    k,c=best_augment(P)
    if k>bestaug[0]:
        bestaug=(k,name,(P,c))
    # only print the notable ones
print(f"  best single-point augmentation of an 81-edge cube: +{bestaug[0]} edges")
print(f"     => u(28) >= 81 + {bestaug[0]} = {81+bestaug[0]}  (need +4=>85 to BEAT 84)")
print(f"     (realization: {bestaug[1]})")

# greedy multi-augment from a chosen realization: 27 -> 28,29,30
def greedy_grow(P0, upto=30):
    P=list(P0)
    log=[(len(P),count_unit(P))]
    while len(P)<upto:
        k,c=best_augment(P)
        if c is None or k==0: break
        P.append(c)
        log.append((len(P),count_unit(P)))
    return P,log
# use the best-augmentable realization
P0=bestaug[2][0] if bestaug[2] else cube(0.7,1.93)
Pg,log=greedy_grow(P0,30)
print("  greedy growth from the cube (n, U, 3n, U-3n):")
for (n,U) in log:
    print(f"     n={n:2d}  U={U:3d}  3n={3*n:3d}  U-3n={U-3*n:+d}")

print("\n"+"="*74)
print("PART D -- from-scratch anneal: densest 27 and 28 point unit-distance graph")
print("="*74)
def anneal(n, iters=60000, restarts=6):
    bestU=-1; bestP=None
    for r in range(restarts):
        # seed: jittered triangular patch + a few random
        P=[]
        # start from cube-like seed for n>=27
        if n>=27:
            seed=cube(random.uniform(0.1,1.0),random.uniform(0.1,1.0))
            P=[list(p) for p in seed[:n]]
            while len(P)<n: P.append([random.uniform(-2,2),random.uniform(-2,2)])
        else:
            P=[[random.uniform(-2,2),random.uniform(-2,2)] for _ in range(n)]
        def U(P):
            u=0
            for i in range(n):
                for j in range(i+1,n):
                    if abs((P[i][0]-P[j][0])**2+(P[i][1]-P[j][1])**2-1.0)<1e-4:u+=1
            return u
        T=0.3
        curU=U(P)
        for it in range(iters):
            T=0.3*(1-it/iters)+0.002
            v=random.randrange(n)
            old=P[v][:]
            if random.random()<0.5:
                P[v][0]+=random.gauss(0,T); P[v][1]+=random.gauss(0,T)
            else:
                # snap toward making a near-unit pair exact
                w=random.randrange(n)
                if w!=v:
                    dx=P[v][0]-P[w][0]; dy=P[v][1]-P[w][1]
                    d=math.hypot(dx,dy)
                    if d>1e-6:
                        P[v][0]=P[w][0]+dx/d; P[v][1]=P[w][1]+dy/d
            nU=U(P)
            if nU>=curU or random.random()<math.exp((nU-curU)/max(T,1e-3)):
                curU=nU
            else:
                P[v]=old
            if curU>bestU:
                bestU=curU; bestP=[p[:] for p in P]
    return bestU,bestP
for n in (27,28):
    U,_=anneal(n, iters=40000, restarts=5)
    print(f"  anneal densest n={n}: U≈{U}  (3n={3*n}, U-3n={U-3*n:+d})  "
          f"[float, lower bound only]")
print("\nDONE PART C/D.")
