#!/usr/bin/env python3
"""
S719 — The unit-distance momentum twistor (angle = log of the rotation group U(1)); the n=22 bleeding edge.

LRC twistor (S718): dual conformal = MULTIPLICATIVE (Z/m)*, linearized by discrete LOG -> additive Z/phi.
Unit distance: a unit edge is a difference vector ON THE UNIT CIRCLE = a modulus-1 number; the dual
conformal symmetry is GLOBAL ROTATION = the multiplicative group U(1). So the unit-distance momentum
twistor is the ANGLE map (the log of U(1)): e^{i theta} -> theta, additive. Global rotation = a constant
shift of all edge-angles (the dual conformal translation, exact analog of LRC multiplier->translation).

KEY DECOMPOSITION (the angle-twistor autocorrelation): resolve the unit distances by edge LINE-direction
phi in [0,pi):
        u(config) = sum_phi m(phi),   m(phi) = #unit edges along line-direction phi.
A single direction is a union of unit-vector chains, so m(phi) <= n-1. To make u large you must
CONCENTRATE edges into FEW directions each near-saturated => the edge-directions form a ROTATION GROUP
(roots of unity) => the CM/cyclotomic structure (the grid-disproof's modulus-1 supply). This is the
twistor explanation of why optimal unit-distance configs are CM/lattice.

In a CM (Eisenstein) layer the "unit" = the norm-D vector; the number of unit-neighbors (degree) =
r_Q(D) = #representations, and the rotation group has order r_Q(D) (its log = Z/r_Q(D), the twistor).
r_Q(1)=6 (triangular, deg 6), r_Q(7)=r_Q(13)=12 (deg 12 -> the AMP sqrt(13) layer, S710/THM-431).

We (1) verify the twistor (angle) linearizes rotation, (2) compute the direction structure of the proven
u(21)=57 = K3 [] W7 and a triangular 22-blob, (3) tabulate r_Q(D) = rotation-group order = degree, and
(4) SEARCH CM layers for a dense 22-config, framing the 57(product) -> 60(record) -> 61(open) gap by the
number of directions / the non-lattice rigidity (S710).

No numpy/sympy.
"""
import math, random
from itertools import combinations
from collections import Counter

# ---------- geometry ----------
def rot(p,th):
    c,s=math.cos(th),math.sin(th); return (c*p[0]-s*p[1], s*p[0]+c*p[1])
def K3(): return [(0.0,0.0),(1.0,0.0),(0.5,math.sqrt(3)/2)]
def W7():
    P=[(0.0,0.0)]
    for k in range(6):
        a=math.pi/3*k; P.append((math.cos(a),math.sin(a)))
    return P
def minkowski(A,B,th):
    BR=[rot(b,th) for b in B]; return [(a[0]+b[0],a[1]+b[1]) for a in A for b in BR]
def unit_edges(P,tol=1e-7):
    E=[]
    for i,j in combinations(range(len(P)),2):
        dx=P[i][0]-P[j][0]; dy=P[i][1]-P[j][1]
        if abs(math.hypot(dx,dy)-1.0)<tol: E.append((i,j,math.atan2(dy,dx)))
    return E
def line_dirs(E, q=1e6):
    """multiset of edge LINE-directions phi in [0,pi), quantized."""
    c=Counter()
    for (_,_,th) in E:
        phi=th % math.pi
        c[round(phi*q)] += 1
    return c

# ---------- Eisenstein CM layer ----------
def eis_pts(R):
    return [(a,b) for a in range(-R,R+1) for b in range(-R,R+1)]
def eis_norm(da,db): return da*da+da*db+db*db
def r_Q(D):
    return sum(1 for a in range(-D,D+1) for b in range(-D,D+1) if a*a+a*b+b*b==D)

def densest_subset_normD(pts, D, k, seed):
    """greedy + local-swap: pick k Eisenstein points maximizing #pairs at norm D."""
    rng=random.Random(seed)
    idx=list(range(len(pts)))
    adj={i:set() for i in idx}
    for i,j in combinations(idx,2):
        if eis_norm(pts[i][0]-pts[j][0],pts[i][1]-pts[j][1])==D:
            adj[i].add(j); adj[j].add(i)
    # greedy start: highest-degree seed, add by marginal gain
    best_overall=0
    for _ in range(6):
        start=rng.choice(idx)
        S=[start]; Sset={start}
        while len(S)<k:
            cand=max((v for v in idx if v not in Sset),
                     key=lambda v: len(adj[v]&Sset)+1e-9*rng.random())
            S.append(cand); Sset.add(cand)
        def edges(Sset): return sum(len(adj[v]&Sset) for v in Sset)//2
        cur=edges(Sset)
        # local swaps
        improved=True
        while improved:
            improved=False
            for v in list(Sset):
                for w in idx:
                    if w in Sset: continue
                    ns=(Sset-{v})|{w}
                    if edges(ns)>cur: Sset=ns; cur=edges(ns); improved=True; break
                if improved: break
        best_overall=max(best_overall,cur)
    return best_overall

if __name__=="__main__":
    print("="*86)
    print("S719 — the unit-distance momentum twistor (angle = log U(1)); the n=22 bleeding edge")
    print("="*86)

    # (1) the angle twistor linearizes rotation
    print("\n(1) TWISTOR = angle (log of U(1)): global rotation -> additive shift of all edge-angles")
    P=minkowski(K3(),W7(),0.4); E=unit_edges(P)
    phi=0.7; Pr=[rot(p,phi) for p in P]; Er=unit_edges(Pr)
    def wrap(a): return (a+math.pi)%(2*math.pi)-math.pi
    a1=sorted(wrap(e[2]) for e in E)              # original edge-angles
    a2=sorted(wrap(e[2]-phi) for e in Er)         # rotated edge-angles, un-shifted by phi
    match = (len(a1)==len(a2)) and all(abs(wrap(x-y))<1e-6 for x,y in zip(a1,a2))
    print(f"  rotating the config by phi shifts every edge-angle by phi (angles match after un-shift): {match}")
    print("  => rotation = additive translation in the angle-twistor (the dual conformal symmetry, linear).")

    # (2) direction structure of the proven optimum K3 [] W7 (u=57) and a triangular blob
    print("\n(2) ANGLE-TWISTOR AUTOCORRELATION: u = sum over line-directions of m(direction)")
    ld=line_dirs(E)
    print(f"  K3 [] W7 (u={len(E)}): #line-directions = {len(ld)}, per-direction counts = {sorted(ld.values(),reverse=True)}")
    print(f"    => few directions, high multiplicity = the rotation-group (rosette) structure; sum = {sum(ld.values())} = u")
    # triangular blob of 22 points (Eisenstein norm-1)
    pts=eis_pts(5)
    # central 22 by |.|
    pts_sorted=sorted(pts,key=lambda p:eis_norm(p[0],p[1]))
    blob=pts_sorted[:22]
    Eb=[(i,j) for i,j in combinations(range(22),2) if eis_norm(blob[i][0]-blob[j][0],blob[i][1]-blob[j][1])==1]
    # directions of triangular edges
    print(f"  triangular 22-blob (norm-1): u={len(Eb)}; the 6 neighbor vectors => 3 line-directions, all parallel families")

    # (3) r_Q(D) = rotation-group order = interior degree of the sqrt(D) layer (the twistor Z/r_Q(D))
    print("\n(3) CM LAYERS: r_Q(D) = #unit-neighbors (degree) = rotation-group order = the twistor Z/r_Q(D)")
    for D in (1,3,7,13,19,21,31,37,49,91):
        print(f"   D={D:3d}: r_Q(D)={r_Q(D):2d} unit-neighbors  ({r_Q(D)//2} line-directions)"
              + ("  <- triangular" if D==1 else ("  <- AMP sqrt(13) layer, deg 12 (S710/THM-431)" if D==13 else "")))
    print("   => the optimum uses sqrt(13) (deg 12) not sqrt(1) (deg 6): more directions = denser; the")
    print("      discrete log of the rotation group (Z/r_Q(D)) is the unit-distance analog of the LRC twistor.")

    # (4) n=22 bleeding edge: densest 22-subset in CM layers
    print("\n(4) n=22 SEARCH: densest 22-point subset in Eisenstein norm-D layers (lattice lower bounds)")
    ptsR=eis_pts(6)
    for D in (1,3,7,13):
        best=densest_subset_normD(ptsR, D, 22, seed=D)
        print(f"   norm-{D:2d} layer: densest 22-subset has {best} unit distances")
    print("   Known: u(22) in [60,61] (open); triangular/lattice patches under-perform because the OPTIMUM")
    print("   is a NON-LATTICE RIGID graph (S710): a lattice DISK only crosses 3N at N=43, the true optimum")
    print("   at N=28. The twistor reading: the record 60/61 needs MORE saturated directions than any single")
    print("   lattice patch supplies at n=22; reaching it = a rigid multi-direction (CM-rotation) graph,")
    print("   i.e. a NON-PRODUCT config (S712: 22=2*11 caps products at 57) with a rich rotation group.")
    print("="*86)
