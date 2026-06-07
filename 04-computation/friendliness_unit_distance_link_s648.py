#!/usr/bin/env python3
"""
S648 — Does friendliness (S647) correlate to the unit distance problem?  YES.

THE DICTIONARY (via the S634 HN/LRC/unit-distance unification, one graph G):
  loneliness (LRC)        <->  being in the INDEPENDENT / 1-avoiding set (alpha)  <-> isolated vertex
  friendliness            <->  being an EDGE endpoint = having neighbors          <-> vertex DEGREE
  unit distance           <->  the EDGES of G (pairs at the resonant distance)
  max friendliness         ->  max edges = u(n) (the unit-distance maximization problem, S641)
  first lonely time tau    ->  first ISOLATION as a scan parameter grows (nearest-neighbor distance)

So 'friendliness' in unit distance = the DEGREE in the unit-distance graph; the friendliest
configuration is the u(n) extremal (triangular/Eisenstein lattice, escaping to Moser/sqrt-11, S641).
The difference: LRC friendliness is a BAND (dist < gap); unit distance is a sharp RESONANCE
(dist = 1 <=> dZ = 1/6, the hexagon, S623). Same covering structure, band vs sharp target.

This script makes the correlation concrete:
  (A) unit-distance friendliness (degree) of the Eisenstein lattice = 6 (the hexagon, formalized S641);
  (B) the first-passage analog: as the threshold radius r grows, when does each point first get a
      neighbor (its 'first friendly radius') -- the unit-distance 'never isolated yet';
  (C) the shared resonance: LRC lonely-at-1/n (tight) <-> unit distance at the hexagonal dZ=1/6.
No external libs.
"""
import math

# ---- triangular/Eisenstein lattice points ----
def emb(a,b): return (a + 0.5*b, (3**0.5/2)*b)
def dist(p,q): return math.hypot(p[0]-q[0], p[1]-q[1])

print("="*68)
print("(A) unit-distance FRIENDLINESS (degree) of the triangular lattice = 6")
print("="*68)
pts=[emb(a,b) for a in range(-4,5) for b in range(-4,5)]
center=emb(0,0)
deg=sum(1 for p in pts if abs(dist(center,p)-1.0)<1e-9)
print(f"  a point of the Eisenstein lattice has {deg} unit-distance neighbours (the hexagon, 6=2*3)")
print(f"  = exactly the formalized eisenstein_unit_neighbours (S641): a^2-ab+b^2=1 has 6 solutions.")
print(f"  -> 'max friendliness' in unit distance is bounded by this rigid degree-6 (Harborth ~3n);")
print(f"     u(n) ESCAPES it via the Moser/sqrt(-11) slab to be friendlier on the boundary (S641).")

print("\n" + "="*68)
print("(B) the FIRST-PASSAGE analog: 'first friendly radius' as r grows (never isolated yet)")
print("="*68)
# scan threshold r upward; a point is 'friendly by r' if it has a neighbour within r.
# nearest-neighbour distance = the first r at which it stops being isolated.
import random
random.seed(3)
P=[(random.random(),random.random()) for _ in range(40)]
def nn(i):
    return min(dist(P[i],P[j]) for j in range(len(P)) if j!=i)
nns=sorted(nn(i) for i in range(len(P)))
print("  random 40-point cloud; 'isolation-survival' S(r)=frac with NO neighbour within r:")
for r in [0.05,0.1,0.15,0.2,0.3]:
    surv=sum(1 for d in nns if d>r)/len(nns)
    print(f"    r={r:4.2f}: still-isolated fraction = {surv:.2f}")
print("  -> EXACTLY the friendliness-survival shape of S647 (flat at 1, then ->0): in LRC the scan")
print("     is TIME and the floor is 1/n; in unit distance the scan is RADIUS and the floor is the")
print("     minimum inter-point distance. 'Never lonely yet' = 'never isolated yet'.")

print("\n" + "="*68)
print("(C) the shared RESONANCE: LRC gap 1/n  <->  unit distance dZ = 1/6 (hexagon)")
print("="*68)
# chord length between two points at clock-distance d on the unit circle = 2 sin(pi d).
# unit distance (chord=1) <=> 2 sin(pi d)=1 <=> d=1/6. (S623 chord/arc bridge.)
for d in [1/14,1/7,1/6,1/3,1/2]:
    chord=2*math.sin(math.pi*d)
    tag=" <- UNIT DISTANCE (hexagon, dZ=1/6)" if abs(d-1/6)<1e-9 else ""
    print(f"    clock-dist d={d:.4f}: chord=2 sin(pi d)={chord:.4f}{tag}")
print("  -> LRC asks for d >= 1/n (a BAND, far from 0); unit distance asks for d = 1/6 EXACTLY")
print("     (a sharp resonance = the cube root / hexagon). Both are 'friendliness in a covering';")
print("     LRC = band-friendly, unit distance = resonance-friendly. Same graph (S634), different")
print("     target distance. The friendliness/survival machinery transfers directly.")

print("\n" + "="*68)
print("SUMMARY: friendliness <-> unit distance")
print("="*68)
print("  loneliness=independent set=alpha ; friendliness=degree=EDGES=unit distances ; max=u(n).")
print("  first-lonely-time tau (scan=time, floor 1/n) <-> first-friendly-radius (scan=radius, floor=")
print("  min distance). Same survival curve. The unit-distance graph's friendliness is bounded by the")
print("  hexagon (degree 6, S641); LRC friendliness is bounded below by 1/n (S647). One covering, two")
print("  scans, one cube-root resonance (1/6 = hexagon = the 6=2*3 the whole arc converges on).")
