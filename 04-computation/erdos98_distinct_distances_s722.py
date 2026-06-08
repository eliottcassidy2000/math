#!/usr/bin/env python3
"""
S722 — Erdos problem 98 (OPEN): h(n) = min #distinct distances over n points with NO 3 collinear,
NO 4 concyclic. Does h(n)/n -> infinity? Erdos couldn't prove h(n)>=n. Upper bounds (constructions):
Pach h(n)<n^(log2 3); Erdos-Furedi-Pach h(n)<n*exp(c sqrt(log n)).

RELATION TO THE REPO (the bridge):
  - Problem 98 is the DUAL of our unit-distance work (S711-S719): unit distances MAXIMIZE one distance's
    multiplicity (the autocorrelation PEAK, u(n)); problem 98 MINIMIZES the #distinct distances (the
    autocorrelation SUPPORT) under a general-position constraint.
  - THE SHARED "3": no 4 concyclic => no config point has 4 others at one distance (4 equidistant points
    are concyclic on a circle centered at that point) => every distance-graph has MAX DEGREE <= 3 =>
    each distance multiplicity m(r) <= 3n/2 => D = #distinct distances >= (n-1)/3. This "3" is exactly
    THM-440/S719's "max unit-cocyclic extension degree = 3" that caps u(22) at 60: problem 98 is the
    GLOBAL (all-radii) form of the local obstruction that caps u(22).
  - The CM/twistor TENSION (S719): few distinct distances wants a CM lattice (few norms = peaked support),
    but a CM lattice has r_Q(D) >= 4 points per circle (12 for the sqrt(13) layer) = MASSIVELY 4-concyclic
    = FORBIDDEN. So problem 98 forbids precisely the CM structure our unit-distance optima exploit; its
    extremal configs must be "CM-broken" (the Erdos-Furedi-Pach perturbation).
  - The COLORING reframe (repo's coloring-unification): minimizing D = decomposing K_n into the fewest
    <=3-regular distance-classes that are simultaneously realizable in the plane with no 3 collinear.

This session: (1) verify no-4-concyclic => per-point distance-multiplicity <= 3 (the bound mechanism),
(2) the lower bound D >= (n-1)/3, (3) construction search: general-position integer configs minimizing
distinct distances for small n (upper bounds on h(n)), (4) quantify the CM tension.

Exact integer arithmetic (integer coords).
"""
from itertools import combinations
import random

def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def collinear3(a,b,c): return cross(a,b,c)==0
def concyclic4(p):
    # 4 points concyclic or collinear iff det[[x^2+y^2,x,y,1]]=0 (4x4). With no-3-collinear, =0 => concyclic
    M=[[x*x+y*y,x,y,1] for (x,y) in p]
    return det4(M)==0
def det4(M):
    # exact 4x4 determinant by cofactor (integers)
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
               -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
               +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    d=0
    for j in range(4):
        minor=[[M[i][k] for k in range(4) if k!=j] for i in range(1,4)]
        d+=((-1)**j)*M[0][j]*det3(minor)
    return d

def has_3_collinear(S):
    for a,b,c in combinations(S,3):
        if collinear3(a,b,c): return True
    return False
def has_4_concyclic(S):
    for q in combinations(S,4):
        if concyclic4(q): return True
    return False
def general_position(S):
    return (not has_3_collinear(S)) and (not has_4_concyclic(S))
def distinct_sq_dists(S):
    return len({(a[0]-b[0])**2+(a[1]-b[1])**2 for a,b in combinations(S,2)})
def max_distance_graph_degree(S):
    """max over points p, over distances r, of #other points at distance r from p."""
    best=0
    for p in S:
        cnt={}
        for q in S:
            if q==p: continue
            d=(p[0]-q[0])**2+(p[1]-q[1])**2
            cnt[d]=cnt.get(d,0)+1
        best=max(best, max(cnt.values()) if cnt else 0)
    return best

# ---- Eisenstein CM (for the tension) ----
def r_Q(D):
    return sum(1 for a in range(-D,D+1) for b in range(-D,D+1) if a*a+a*b+b*b==D)

if __name__=="__main__":
    rng=random.Random(0)
    print("="*88)
    print("S722 — Erdos problem 98: distinct distances under no-3-collinear, no-4-concyclic")
    print("="*88)

    # (1) no 4 concyclic => per-point distance-multiplicity <= 3 (=> distance-graph max degree <= 3)
    print("\n(1) MECHANISM: no 4 concyclic => every distance-graph has max degree <= 3")
    ok=True
    for _ in range(400):
        S=set()
        while len(S)<7: S.add((rng.randint(-6,6),rng.randint(-6,6)))
        S=list(S)
        if has_4_concyclic(S): continue
        if max_distance_graph_degree(S)>3: ok=False; break
    print(f"  verified on random 7-pt no-4-concyclic configs: max distance-graph degree <= 3 : {ok}")
    print("  (a point with 4 others at one distance => those 4 lie on a circle centered at it => 4 concyclic)")

    # (2) lower bound D >= (n-1)/3
    print("\n(2) LOWER BOUND: each point sees each distance <=3 times => >= ceil((n-1)/3) distinct distances")
    for n in (5,8,11,15,22,50):
        print(f"   n={n:3d}: h(n) >= ceil((n-1)/3) = {-(-(n-1)//3)}   (Erdos wants h(n)>=n, conj h(n)/n->inf; OPEN)")

    # (3) construction search: general-position integer configs minimizing distinct distances
    print("\n(3) CONSTRUCTIONS (upper bounds on h(n)): general-position integer configs, fewest distinct distances")
    for k in (5,6,7,8):
        best=None; bestS=None
        # random + greedy over a grid
        for _ in range(60000):
            R=3 if k<=6 else 4
            S=set()
            while len(S)<k: S.add((rng.randint(-R,R),rng.randint(-R,R)))
            S=list(S)
            if has_3_collinear(S) or has_4_concyclic(S): continue
            d=distinct_sq_dists(S)
            if best is None or d<best: best=d; bestS=S
        lb=-(-(k-1)//3)
        print(f"   n={k}: fewest distinct distances found = {best}  (lower bound ceil((n-1)/3)={lb}); config={sorted(bestS) if bestS else None}")

    # (4) the CM tension: CM lattices violate no-4-concyclic massively (the structure our optima use)
    print("\n(4) CM TENSION: a CM lattice point has r_Q(D) others per circle => >=4 concyclic = FORBIDDEN by 98")
    for D in (1,3,7,13):
        print(f"   Eisenstein norm-{D}: r_Q({D})={r_Q(D)} points on each such circle "
              + ("(<=3 OK)" if r_Q(D)<=3 else "(>=4 => 4-CONCYCLIC, forbidden by problem 98)"))
    print("   => few-distinct-distance CM/lattice configs (S719) are EXACTLY what problem 98 forbids:")
    print("      the unit-distance OPTIMUM (peaked, concyclic) and the problem-98 OPTIMUM (spread, general")
    print("      position) are the two opposite extremizations of the SAME radial autocorrelation.")
    print("="*88)
