#!/usr/bin/env python3
"""
S631 Part 2 — unit distance as a TRIENEMENT (THM-389 structure): for each pair of points,
  dist < 1  -> arrow one way,   dist = 1 -> TIE (unit-distance edge),   dist > 1 -> arrow other way.
The unit-distance graph = the TIE subgraph; Erdos' problem = maximize ties.
Claim to test: in this framework the disproof's quantity is the TIE-DEGREE = #lattice vectors of
norm 1 = #modulus-1 elements; maximizing ties = maximizing the symmetry (the tie-graph is a CAYLEY
graph). So 'CM beats grid' = 'a lattice with more norm-1 vectors / bigger automorphism group'.
"""
import itertools, math, cmath
from collections import Counter

def tie_degree_lattice(gram_norm1_count):
    return gram_norm1_count  # placeholder; computed per lattice below

def trienement(pts, tol=1e-7):
    """counts (<1, =1, >1) and the tie-graph degree sequence."""
    n = len(pts); lt = eq = gt = 0
    deg = [0]*n
    for i, j in itertools.combinations(range(n), 2):
        d = abs(pts[i]-pts[j])
        if abs(d-1) < tol: eq += 1; deg[i]+=1; deg[j]+=1
        elif d < 1: lt += 1
        else: gt += 1
    return lt, eq, gt, deg

# triangular (Eisenstein) lattice patch
w = cmath.exp(1j*math.pi/3)
def disk(R): return [(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if i*i+i*j+j*j<=R*R]
patch = [i+j*w for (i,j) in disk(3)]
lt,eq,gt,deg = trienement(patch)
print("TRIANGULAR (Eisenstein CM Q(sqrt-3)) patch, %d points:" % len(patch))
print(f"   trienement: <1 arrows={lt}, =1 TIES(unit dist)={eq}, >1 arrows={gt}")
print(f"   tie-degree distribution: {dict(sorted(Counter(deg).items()))} (interior deg=6 = #units of Z[w])")
print("   => the tie-graph is the Cayley graph of the lattice; tie-degree = #norm-1 vectors = #roots of unity in Q(sqrt-3) = 6")

# square (Gaussian) lattice Z[i]: norm-1 vectors = +-1,+-i = 4 (units of Z[i])
print("\nSQUARE (Gaussian CM Q(i)) lattice: norm-1 vectors = units of Z[i] = {+-1,+-i} = 4 -> interior tie-degree 4")
print("   triangular (6) > square (4): MORE ties per point = denser unit-distance graph (why triangular wins at small n)")

# the disproof's quantity, restated: tie-degree growth via CM modulus-1 elements
print("\nDISPROOF restated in trienement terms:")
print("   grid: tie-degree bounded by lattice units (4 or 6) + a few extra rotations -> u(n)=n^{1+o(1)}")
print("   CM tower: tie-degree = #modulus-1 algebraic numbers of bounded height, GROWS with field degree")
print("   => 'CM beats grid' = 'the tie-graph (unit-distance graph) is a Cayley graph of a CM lattice")
print("      whose automorphism/unit supply is superlinear' — the symmetry (perspective/rigidity) IS the count.")

# concrete: number of representations of integers by x^2+xy+y^2 (Eisenstein norm) = tie-multiplicity at radius sqrt(m)
print("\nnorm-form representation counts r(m) for x^2+xy+y^2 (Eisenstein) — the tie-degree at distance sqrt(m):")
def r_eis(m):
    c=0
    B=int(m**0.5)+2
    for x in range(-B,B+1):
        for y in range(-B,B+1):
            if x*x+x*y+y*y==m: c+=1
    return c
vals=[(m,r_eis(m)) for m in [1,3,7,13,21,31,37,43,49,91]]
print("  ", vals)
print("   (high r(m) at m with many prime factors p=1 mod 3 = the unit-distance 'rich radii' = max ties;")
print("    the disproof pushes this multiplicity up via class field towers = many split primes.)")
