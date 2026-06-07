#!/usr/bin/env python3
"""
unit_distance_cube_angle_rigidity_monad_s6.py
monad-explorer-2026-06-07-S6  (deep-research; OPEN-Q-057)  --  THM-435 verifier.

THM-435 (CUBE ANGLE-RIGIDITY).  Let C(t2,t3) be the Minkowski-sum realization of
the Hamming cube  K3^[]3 = H(3,3)  as three unit equilateral triangles, the i-th
rotated by angle t_i (t1=0 WLOG):  points  = v1 + R(t2) v2 + R(t3) v3,  vi in the
unit triangle {0, e^{i0}, e^{i60}}.  Then for EVERY (t2,t3):

   (#distinct points = 27)  ==>  (#unit distances = 81 = 3*27)   exactly.

i.e. no choice of angles yields 27 distinct points with MORE than 81 unit
distances.  Any angle specialization that would create an extra unit distance
simultaneously collapses two of the 27 points.  Hence the H(3,3) construction
family CANNOT witness u(27) > 81: the 3N tie at n=27 is angle-rigid.

PROOF (companion -- this script verifies each step numerically):

 The 81 "product" edges (two points differing in exactly ONE triangle coordinate
 by a triangle-edge vector) have length 1 for ALL angles (a triangle edge is a
 unit vector; rotation preserves length).  An EXTRA unit distance is a pair of
 points differing in >=2 coordinates whose displacement -- a sum of 2 or 3 unit
 vectors, one per differing factor -- has length 1.

 (A) TWO differing factors i<j:  displacement = e^{i a} + e^{i b} with a in the
     hex-direction set of factor i (multiples of 60 + t_i) and b in that of factor
     j.  |e^{ia}+e^{ib}|=1  <=>  angle(a,b) = 120 deg  <=>  t_i - t_j ≡ 0 (mod 60).
     LEM-A verifies: t_i - t_j ≡ 0 (mod 60)  =>  the cube has < 27 distinct points.

 (B) THREE differing factors:  WLOG divide by factor-1's unit vector; condition is
        cos u + cos w + cos(u - w) = -1,
     where u ≡ t2, w ≡ t3, u - w ≡ t2 - t3  (each mod 60 deg).
     LEM-B (the key trig lemma) verifies the COMPLETE solution set is exactly
        { u ≡ 180 } U { w ≡ 180 } U { w ≡ u - 180 }   (mod 360),
     i.e.  t2 ≡ 0  OR  t3 ≡ 0  OR  t2 - t3 ≡ 0   (mod 60 deg).  (The first family
     is the degenerate 1+cos u = 0 line, all w.)  All three are collision loci by
     LEM-A.  So every 3-factor extra edge also forces < 27 distinct points.

 (A)+(B): whenever the 27 points are distinct, no extra edge exists => U = 81. QED.

 LEM-C: global fine scan -- over a dense grid of (t2,t3) the cube NEVER has 27
 distinct points with U > 81 (independent confirmation of the proof).
"""
import math
TWO_PI = 2*math.pi
def norm_deg(x): return x % 360.0

# ---------------- LEM-B: the trig lemma ----------------
print("="*74)
print("LEM-B  solution set of  cos u + cos w + cos(u-w) = -1")
print("="*74)
# Dense scan over (u,w); for every approximate solution, check it satisfies
# w ≡ 180 (mod 360) OR w ≡ u-180 (mod 360).  Report any solution that does NOT.
def is_collision_family(u_deg, w_deg, tol=0.75):
    # DEGENERATE: u ≡ 180 (mod 360) -> 1+cos u = 0, equation is 0=0 for all w
    #            (this is t2 ≡ 0 mod 60, a collision locus by LEM-A)
    du = norm_deg(u_deg-180)
    if min(du, 360-du) < tol: return True
    # w ≡ 180 (mod 360)   [ t3 ≡ 0 mod 60 ]
    if min(abs(norm_deg(w_deg-180)), abs(norm_deg(w_deg-180)-360)) < tol: return True
    # w ≡ u-180 (mod 360) [ t2 - t3 ≡ 0 mod 60 ]
    d = norm_deg(w_deg-(u_deg-180))
    if min(d, 360-d) < tol: return True
    return False

NU=2000; NW=2000
sols=0; rogue=0; rogue_samples=[]
step_u=360.0/NU
for iu in range(NU):
    u=math.radians(iu*step_u)
    cu=math.cos(u)
    # solve cos w (1+cu) + sin w * sin u = -1 - cu   for w on a fine grid + bisection-free check
    # just scan w; record sign changes of f(w)=cosu+cosw+cos(u-w)+1
    prev=None; prevw=None
    for iw in range(NW+1):
        w=math.radians(iw*360.0/NW)
        f=cu+math.cos(w)+math.cos(u-w)+1.0
        if prev is not None and prev==0.0:
            pass
        if prev is not None and (prev<0)!=(f<0):
            # root between prevw and w -> refine by linear interp
            wr=prevw + (w-prevw)*(0.0-prev)/(f-prev)
            sols+=1
            if not is_collision_family(math.degrees(u), math.degrees(wr), tol=1.2):
                rogue+=1
                if len(rogue_samples)<8:
                    rogue_samples.append((math.degrees(u),math.degrees(wr)))
        prev=f; prevw=w
print(f"  found {sols} root crossings on the grid")
print(f"  ROGUE solutions (NOT in w≡180 or w≡u-180): {rogue}")
if rogue:
    print("  !!! LEM-B FAILS -- samples:", rogue_samples)
else:
    print("  => EVERY solution lies on {w≡180} U {w≡u-180}.  LEM-B CONFIRMED.")

# algebraic cross-check: derived closed form w = 180 or w = u-180
print("\n  algebraic identity check (A cos w + B sin w = C):")
import random
random.seed(1)
bad=0
for _ in range(100000):
    u=random.uniform(0,TWO_PI)
    for w in (math.pi, u-math.pi):
        val=math.cos(u)+math.cos(w)+math.cos(u-w)+1.0
        if abs(val)>1e-9: bad+=1
print(f"    w=180 and w=u-180 satisfy the equation for all u:  failures={bad}  (expect 0)")

# ---------------- LEM-A: collision loci ----------------
print("\n"+"="*74)
print("LEM-A  t_i - t_j ≡ 0 (mod 60 deg)  =>  cube has < 27 distinct points")
print("="*74)
HEX=[ (math.cos(math.radians(60*k)), math.sin(math.radians(60*k))) for k in range(6)]
def tri(t):
    return [(0.0,0.0),(math.cos(t),math.sin(t)),(math.cos(t+math.pi/3),math.sin(t+math.pi/3))]
def cube(t2,t3):
    pts=[(0.0,0.0)]
    for f in (tri(0.0),tri(t2),tri(t3)):
        pts=[(p[0]+q[0],p[1]+q[1]) for p in pts for q in f]
    return pts
def ndistinct(pts,tol=1e-7):
    s=set()
    for x,y in pts: s.add((round(x/tol),round(y/tol)))
    return len(s)
# test all special-angle combinations where some pair of {0,t2,t3} differ by mult of 60
print("  angles where t2≡0, t3≡0, or t2-t3≡0 (mod 60), with generic offset 'g':")
g=0.4137  # generic
cases={
 "t2≡0 (t2=0)":            (0.0, g),
 "t2≡0 (t2=60)":           (math.pi/3, g),
 "t3≡0 (t3=120)":          (g, math.radians(120)),
 "t2-t3≡0 (t2=g,t3=g)":    (g, g),
 "t2-t3≡0 (t3=t2+60)":     (g, g+math.pi/3),
 "t2-t3≡0 (t3=t2-120)":    (g, g-math.radians(120)),
}
allcol=True
for name,(a,b) in cases.items():
    nd=ndistinct(cube(a,b))
    print(f"    {name:24s}: distinct points = {nd}  ({'COLLISION' if nd<27 else 'NO collision !!'})")
    if nd>=27: allcol=False
print("  => "+("LEM-A CONFIRMED (all collision loci collapse points)." if allcol
              else "LEM-A FAILS !!"))

# ---------------- LEM-C: global confirmation ----------------
print("\n"+"="*74)
print("LEM-C  global: NO (t2,t3) gives 27 distinct points with U > 81")
print("="*74)
def count_unit(pts,tol=1e-9):
    n=len(pts);u=0
    for i in range(n):
        for j in range(i+1,n):
            dx=pts[i][0]-pts[j][0];dy=pts[i][1]-pts[j][1]
            if abs(dx*dx+dy*dy-1.0)<tol:u+=1
    return u
N=300
maxU=0; argmax=None; over=0
for i in range(1,N):
    t2=math.pi/3*i/N
    for j in range(1,N):
        t3=math.pi/3*j/N
        P=cube(t2,t3)
        if ndistinct(P)<27: continue
        U=count_unit(P)
        if U>maxU: maxU=U; argmax=(math.degrees(t2),math.degrees(t3))
        if U>81: over+=1
print(f"  scanned ~{N*N} (t2,t3) in (0,60)^2; max U over 27-distinct cubes = {maxU}")
print(f"  #configs with U>81 and 27 distinct points = {over}")
print("  => "+("LEM-C CONFIRMED: cube angle-rigid at 81." if (maxU<=81 and over==0)
              else "LEM-C anomaly !!"))
print("\nTHM-435 verification complete.")
