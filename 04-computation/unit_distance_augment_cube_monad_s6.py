#!/usr/bin/env python3
"""
unit_distance_augment_cube_monad_s6.py
monad-explorer-2026-06-07-S6  (deep-research; OPEN-Q-057 frontier)

NUMERIC LANDSCAPE MAP for the 3N unit-distance crossover at N in {27,28}.

Background (OPEN-Q-057, THM-431/432/433): u(N) = max # unit distances among N
planar points; N* = smallest N with u(N) > 3N. Proven N* in [25,28].
  - u(27) >= 81 = 3*27 (TIE), realized EXACTLY by the Hamming cube H(3,3)=K3^[]3
    (6-regular, generic-angle Minkowski sum of 3 unit triangles).  -- THM-432
  - u(28) >= 85 (BEAT), AMP arXiv:2412.11914 lower bound, only CITED in repo,
    NOT exactly realized.  The S711 handoff asks: is the u(28)=85 crosser
    "H(3,3)+1" -- a 28th point unit-distant from 4 of the cube's 27 vertices?

THIS SCRIPT (numeric, float; an EXACT verifier follows for any winner):
  PART A  validate: H(3,3) -> 81/27/6-reg; K3[]W7 -> 57 (proven u(21)).
  PART B  Harborth check: best triangular(Eisenstein)-lattice patch n=22..28,
          to resolve the S4 reflection's "triangular blob ~78 @ 27" vs Harborth=63.
  PART C  AUGMENTATION: for many rotation-angle realizations of H(3,3), find a
          28th point lying on the max number of the 27 unit circles (= max extra
          unit distances). Tests H(3,3)+1 -> can a single point add >=4 (=>85)?
          Also multi-point greedy augmentation 27 -> 28,29,30.
  PART D  Direct anneal: densest 27- and 28-point unit-distance graph from scratch
          (independent check of u(27)=81, u(28)=85 lower bounds).

Float only here; coordinates of any record config are dumped for exact re-verify.
"""
import math, random, itertools

random.seed(20260607)
TOL = 1e-9

# ---------------- basic float geometry ----------------
def sqd(p, q):
    return (p[0]-q[0])**2 + (p[1]-q[1])**2

def count_unit(pts, tol=TOL):
    n = len(pts); u = 0
    for i in range(n):
        for j in range(i+1, n):
            if abs(sqd(pts[i], pts[j]) - 1.0) < tol:
                u += 1
    return u

def degseq(pts, tol=TOL):
    n = len(pts); deg = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if abs(sqd(pts[i], pts[j]) - 1.0) < tol:
                deg[i]+=1; deg[j]+=1
    return deg

def distinct(pts, tol=1e-7):
    n=len(pts)
    for i in range(n):
        for j in range(i+1,n):
            if sqd(pts[i],pts[j]) < tol*tol:
                return False
    return True

def tri(theta):
    """unit equilateral triangle rotated by theta: vertices 0, e^{i th}, e^{i(th+60)}."""
    return [(0.0,0.0),
            (math.cos(theta), math.sin(theta)),
            (math.cos(theta+math.pi/3), math.sin(theta+math.pi/3))]

def minkowski(factors):
    pts=[(0.0,0.0)]
    for f in factors:
        pts=[(p[0]+q[0], p[1]+q[1]) for p in pts for q in f]
    return pts

def cube(th2, th3, th1=0.0):
    return minkowski([tri(th1), tri(th2), tri(th3)])

# ---------------- PART A: validation ----------------
print("="*74)
print("PART A -- validate machinery against PROVEN values")
print("="*74)
# generic angles -> H(3,3)
H = cube(0.7, 1.93)
print(f"H(3,3) generic: pts={len(H)} distinct={distinct(H)} U={count_unit(H)} "
      f"(expect 27, True, 81); degset={sorted(set(degseq(H)))}")

# K3 [] W7 : W7 = hub + unit hexagon ; should give u=57 (proven u(21))
def wheel7(theta=0.0):
    pts=[(0.0,0.0)]
    for k in range(6):
        a=theta+k*math.pi/3
        pts.append((math.cos(a), math.sin(a)))
    return pts
K3W7 = minkowski([tri(0.0), [ (math.cos(0.61)*x-math.sin(0.61)*y, math.sin(0.61)*x+math.cos(0.61)*y) for (x,y) in wheel7()] ])
print(f"K3[]W7 generic: pts={len(K3W7)} distinct={distinct(K3W7)} U={count_unit(K3W7)} "
      f"(expect 21, True, 57)")

# ---------------- PART B: Harborth triangular patch ----------------
print("\n"+"="*74)
print("PART B -- best triangular(Eisenstein) lattice patch, n=22..28")
print("         resolve S4 'triangular blob ~78 @ 27' vs Harborth floor(3n-sqrt(12n-3))")
print("="*74)
# Eisenstein integer coords: point (a,b) -> a*(1,0)+b*(1/2, sqrt3/2)
def eis(a,b): return (a+b*0.5, b*math.sqrt(3)/2)
R=8
lat=[(a,b) for a in range(-R,R+1) for b in range(-R,R+1)]
latpts={(a,b):eis(a,b) for (a,b) in lat}
# spiral/greedy: pick n points closest to centroid that maximize compactness.
# Use distance-from-origin ordering (centered hex grows optimally for hex numbers).
def harborth(n): return math.floor(3*n - math.sqrt(12*n-3))
def best_tri_patch(n, trials=4000):
    best=-1; bestset=None
    # candidate centers near origin to break ties of the hex shells
    centers=[(0,0),(0.5,0.0),(0.33,0.33),(0.0,0.5),(1/3,2/3)]
    for cx,cy in centers:
        order=sorted(lat, key=lambda ab:(eis(*ab)[0]-cx)**2+(eis(*ab)[1]-cy)**2)
        sel=order[:n]
        S=set(sel)
        # count edges exactly via integer norm: (a,b)-(c,d) is unit iff
        # (da)^2+(da)(db)+(db)^2==1  (Eisenstein norm with this basis)
        e=0
        sl=list(S)
        for i in range(len(sl)):
            ai,bi=sl[i]
            for j in range(i+1,len(sl)):
                aj,bj=sl[j]
                da=ai-aj; db=bi-bj
                if da*da+da*db+db*db==1: e+=1
        if e>best: best=e; bestset=sl
    # local improvement: try swapping a boundary point for an outside neighbor
    improved=True
    S=set(bestset)
    def edges(S):
        e=0; sl=list(S)
        for i in range(len(sl)):
            ai,bi=sl[i]
            for j in range(i+1,len(sl)):
                aj,bj=sl[j]
                da=ai-aj; db=bi-bj
                if da*da+da*db+db*db==1: e+=1
        return e
    cur=edges(S)
    nbr=[(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
    for _ in range(6):
        improved=False
        frontier=set()
        for (a,b) in S:
            for da,db in nbr:
                if (a+da,b+db) not in S and -R<=a+da<=R and -R<=b+db<=R:
                    frontier.add((a+da,b+db))
        for inn in list(S):
            for out in list(frontier):
                S2=set(S); S2.discard(inn); S2.add(out)
                e2=edges(S2)
                if e2>cur:
                    cur=e2; S=S2; improved=True
        if not improved: break
    return cur
print(f"  n |  greedy tri-patch U | Harborth floor(3n-sqrt(12n-3)) | 3n | AMP u>=")
AMP={22:60,23:64,24:68,25:72,26:76,27:81,28:85}
for n in range(22,29):
    g=best_tri_patch(n)
    print(f" {n:2d} |        {g:4d}         |          {harborth(n):4d}                | {3*n:2d} | {AMP[n]}")

print("\nDONE PART A/B. (PART C/D in next block.)")
