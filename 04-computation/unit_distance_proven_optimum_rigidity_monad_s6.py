#!/usr/bin/env python3
"""
unit_distance_proven_optimum_rigidity_monad_s6.py
monad-explorer-2026-06-07-S6  (OPEN-Q-057) -- does THM-437's angle-rigidity
mechanism also govern the PROVEN optimum u(21)=57?

THM-437 proved the n=27 cube K3^[]3 is angle-rigid at 81 (accidental edge <=>
collision, via 60-deg quantization of the equilateral triangle's directions).
Question: the PROVEN-optimal u(21)=57 extremal is K3 [] W7 (W7 = hub + unit
hexagon; AMP arXiv:2412.11914 + THM-431/432). Is it also angle-rigid -- i.e. does
tuning the relative angle of the K3 vs the W7 ever give 21 distinct points with
U > 57?  W7 has NON-unit internal distances (sqrt3, 2), so unlike the pure cube an
accidental edge could be (triangle unit) + (rotated W7 sqrt3-or-2 vector) of total
length 1 at a NON-aligned angle -- a richer mechanism. We scan exhaustively.

Also: K3 [] (hub+hexagon) vs alternative wheels. Float scan; report any U>57 with
21 distinct points (would be a NEW realization detail, though u(21)=57 is PROVEN so
no config can EXCEED 57 -- a U>57 with 21 distinct would signal a counting bug).
"""
import math
TOL=1e-9
def sqd(p,q): return (p[0]-q[0])**2+(p[1]-q[1])**2
def count_unit(pts,tol=TOL):
    n=len(pts);u=0
    for i in range(n):
        for j in range(i+1,n):
            if abs(sqd(pts[i],pts[j])-1.0)<tol:u+=1
    return u
def ndistinct(pts,tol=1e-7):
    s=set()
    for x,y in pts: s.add((round(x/tol),round(y/tol)))
    return len(s)
def rot(p,th):
    c,s=math.cos(th),math.sin(th)
    return (c*p[0]-s*p[1], s*p[0]+c*p[1])
def tri(): return [(0.0,0.0),(1.0,0.0),(0.5,math.sqrt(3)/2)]
def wheel7():
    pts=[(0.0,0.0)]
    for k in range(6):
        a=k*math.pi/3
        pts.append((math.cos(a),math.sin(a)))
    return pts
def mink(A,B,th):
    Br=[rot(q,th) for q in B]
    return [(a[0]+b[0],a[1]+b[1]) for a in A for b in Br]

print("="*74)
print("K3 [] W7  (proven u(21)=57 extremal): scan relative angle th for U>57")
print("="*74)
A=tri(); B=wheel7()
N=3600  # 0.05-deg resolution over full 360 (period 60 expected but scan all)
maxU=0; argmax=None; over=0; over_samples=[]
collisions=0
for i in range(N):
    th=2*math.pi*i/N
    P=mink(A,B,th)
    nd=ndistinct(P)
    if nd<21:
        collisions+=1
        continue
    U=count_unit(P)
    if U>maxU: maxU=U; argmax=math.degrees(th)
    if U>57:
        over+=1
        if len(over_samples)<10: over_samples.append((math.degrees(th),U,nd))
print(f"  scanned {N} angles: {collisions} had collisions (<21 distinct).")
print(f"  MAX U over 21-distinct configs = {maxU}  at th={argmax} deg")
print(f"  #configs with U>57 AND 21 distinct = {over}")
if maxU<=57 and over==0:
    print("  => ANGLE-RIGID at 57: no tuning beats the proven optimum (consistent")
    print("     with u(21)=57 PROVEN; same accidental<=>collision mechanism as THM-437).")
elif over>0:
    print("  !!! U>57 with 21 distinct points -- since u(21)=57 is PROVEN, this would")
    print("      indicate a COUNTING/COLLISION-TOL bug. Samples:", over_samples)

print("\n"+"="*74)
print("FINE near-collision sweep: do special angles add accidental edges BEFORE")
print("collapsing (i.e. is there a window U in (57, ...] with 21 distinct)?")
print("="*74)
# zoom near multiples of 60 deg and near other candidate accidental angles
hits=0
for base in range(0,360,1):
    for frac in range(-50,51):
        th=math.radians(base+frac*0.01)
        P=mink(A,B,th)
        if ndistinct(P)<21: continue
        U=count_unit(P)
        if U>57:
            hits+=1
            if hits<=5: print(f"   th={base+frac*0.01:.4f} U={U}")
print(f"  total (th, 21-distinct, U>57) hits in fine sweep: {hits}")
print("  => "+("none: confirms angle-rigidity at 57." if hits==0
              else "see samples (investigate)."))
print("\nDONE.")
