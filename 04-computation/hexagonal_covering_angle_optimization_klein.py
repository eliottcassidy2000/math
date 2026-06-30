#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE GEOMETRIC CORE: the rhombic-lattice covering density is MINIMIZED at 60deg (klein-S27).

Pursue 'the optimal LRC covering is the hexagonal/zeta_6 one'. A 2D LATTICE covering's fundamental cell is a
rhombus (two triangles); scan the rhombus angle theta. Covering density delta(theta) = pi R(theta)^2 / area,
R = covering radius (deep hole = Delaunay circumradius), area = sin theta. MIN at theta=60deg (rhombus = TWO
EQUILATERAL triangles = the hexagonal lattice = Kershner 2pi/sqrt27), beating theta=90deg (square = two
RIGHT triangles = pi/2). The convex-tiling classification (Reinhardt 1918 + Rao 2017) makes the search FINITE.
"""
import math

def covering_radius(theta, N=240):
    """lattice basis (1,0),(cos t, sin t); covering radius = max over fundamental cell of min dist to lattice."""
    ct,st=math.cos(theta),math.sin(theta)
    pts=[(i+ci*ct, ci*st) for i in range(-2,3) for ci in range(-2,3)]
    best=0.0
    for a in range(N+1):
        for b in range(N+1):
            x=a/N + (b/N)*ct; y=(b/N)*st
            d=min((x-px)**2+(y-py)**2 for px,py in pts)
            if d>best: best=d
    return math.sqrt(best)

print("="*80)
print(" Rhombic-lattice COVERING DENSITY vs rhombus angle theta (unit basis vectors)")
print("="*80)
print(f" {'theta':>7}{'cov.radius R':>13}{'area=sin':>10}{'density pi R^2/area':>20}{'  cell = 2 triangles'}")
results=[]
for deg in [50,55,60,65,70,75,80,85,90]:
    th=math.radians(deg); R=covering_radius(th); area=math.sin(th); dens=math.pi*R*R/area
    results.append((dens,deg))
    tri = "2 EQUILATERAL (optimal)" if deg==60 else ("2 right (45-45-90)" if deg==90 else "2 isoceles")
    print(f" {deg:>6}d{R:>13.5f}{area:>10.5f}{dens:>20.5f}   {tri}")
mn=min(results)
print(f"\n MIN density at theta = {mn[1]} deg, delta = {mn[0]:.5f}  (Kershner 2pi/sqrt27 = {2*math.pi/math.sqrt(27):.5f})")
print(f" square (90 deg) density = {[d for d,deg in results if deg==90][0]:.5f}  (= pi/2 = {math.pi/2:.5f}); hexagon WINS.")

print("\n"+"="*80)
print(" THE FINITE CLASSIFICATION (Reinhardt 1918 + Rao 2017) makes the covering search FINITE")
print("="*80)
print(" Convex polygons that monohedrally tile the plane: ALL triangles, ALL quadrilaterals,")
print(" exactly 3 hexagon types, exactly 15 pentagon types; NONE with >= 7 sides. (Rao 2017 closed it.)")
print(" For a LATTICE covering the Voronoi cell is CENTRALLY SYMMETRIC + convex with <=6 sides")
print("  => only the HEXAGON or the RECTANGLE (parallelogram). The angle scan above: hexagon (60deg,")
print("  rhombus=2 equilateral triangles) beats rectangle/square (90deg, 2 right triangles).")
print("  So among lattice coverings the optimum is the EQUILATERAL-TRIANGULAR (hexagonal) lattice -- a")
print("  FINITE check (hexagon vs rectangle), won by the hexagon (Kershner). The triangular lattice = the")
print("  A2 'god': the 3-fold (tridiagonal/A2) structure, p6m, the zeta_6-line's ambient lattice (HYP-3715).")

print("\n"+"="*80)
print(" APERIODIC (Socolar-Taylor) does NOT beat it")
print("="*80)
print(" Fejes Toth: the hexagonal lattice is the thinnest covering among ALL coverings (not just lattices).")
print(" The Socolar-Taylor aperiodic monotile is a DECORATED HEXAGON (forces aperiodicity) -- still in the")
print(" hexagonal family; aperiodicity gives no covering-density gain. So the optimum stays hexagonal/A2.")
