#!/usr/bin/env python3
"""
S634 — LRC and unit distance are ONE problem: chi of a unit-distance Cayley graph Cay(G,U).
  LRC:  G = Z/m, U = {±1,...,±band} (the resonance/tie set) -> the cycle/circulant; chi = sieve arity.
  Unit distance: G = R^2/lattice, U = unit sphere -> the plane unit-distance graph; chi = Hadwiger-Nelson.
Dimension of G is the knob. Shared mechanism: RIGIDITY (odd cycle in 1D, Moser spindle in 2D) forces
chi up. Shared key: the CM ROTATION (LRC multiplier / spindle angle arccos(5/6)) glues rotated copies.

This script: (1) the chi dimension-ladder; (2) the Moser spindle as a '2D multiplier' (rotation by
arccos(5/6)); (3) the transfer: a finite circulant/lattice subgraph lower-bounds chi.
"""
import itertools, math, cmath

def chromatic_number(adj, ub=8):
    N=len(adj)
    def colorable(k):
        col=[-1]*N
        def bt(v):
            if v==N: return True
            for c in range(k):
                if all(not(adj[v][u] and col[u]==c) for u in range(N)):
                    col[v]=c
                    if bt(v+1): return True
                    col[v]=-1
            return False
        return bt(0)
    for k in range(1,ub+1):
        if colorable(k): return k
    return ub+1

def chi_points(pts, tol=1e-7):
    N=len(pts); adj=[[0]*N for _ in range(N)]
    for i,j in itertools.combinations(range(N),2):
        if abs(abs(pts[i]-pts[j])-1)<tol: adj[i][j]=adj[j][i]=1
    return chromatic_number(adj), sum(adj[i][j] for i in range(N) for j in range(i+1,N))

def cycle_chi(m):
    adj=[[0]*m for _ in range(m)]
    for i in range(m): adj[i][(i+1)%m]=adj[(i+1)%m][i]=1
    return chromatic_number(adj)

print("(1) chi DIMENSION-LADDER of unit-distance Cayley graphs")
print("   1D  Cay(Z/m,{±1}) = C_m (LRC tie-graph):")
for m in (5,6,7,8,9): print(f"        m={m}: chi={cycle_chi(m)}  (2 even / 3 odd = the sieve arity, 2-adic seam)")
print("   2D  triangular (Eisenstein) unit-distance patches:")
w=cmath.exp(1j*math.pi/3)
def disk(R): return [(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if i*i+i*j+j*j<=R*R]
for R in (1,2,3):
    pts=[i+j*w for (i,j) in disk(R)]; chi,E=chi_points(pts)
    print(f"        triangular disk R={R} ({len(pts)} pts, {E} unit edges): chi={chi}")

print("\n(2) MOSER SPINDLE = a '2D multiplier': two unit rhombi glued by a rotation of arccos(5/6)")
A=0+0j; B=1+0j; C=cmath.exp(1j*math.pi/3); D=B+C            # rhombus far vertex, |D|=sqrt(3)
th=math.acos(5/6)                                            # the rotation gluing the far vertices at dist 1
r=cmath.exp(1j*th)
B2,C2,D2=r*B,r*C,r*D                                          # rotated rhombus
spindle=[A,B,C,D,B2,C2,D2]
chi,E=chi_points(spindle)
print(f"   spindle: |D-D'|={abs(D-D2):.4f} (=1 unit edge), {E} unit edges, chi={chi}  (=4: forces a 4th color)")
print(f"   the gluing rotation cos(th)=5/6 is the CM/multiplier 'perspective' in 2D — exactly the role")
print(f"   the LRC multiplier a in (Z/m)* plays on the shell (S625): rotated copies create the obstruction.")

print("\n(3) TRANSFER: rigidity raises chi in BOTH — odd cycle (1D) and spindle (2D) are the same move")
print(f"   1D: C_m odd  => not 2-colorable => chi=3  (single->multi sieve, S632)")
print(f"   2D: spindle  => not 3-colorable => chi=4  (Moser; the 2D 'odd structure')")
print("   => 'multi-sieve necessary' (LRC) and 'chi(plane)>=k' (Hadwiger-Nelson) are the SAME method:")
print("      exhibit a finite rigid subgraph whose chi is forced up. LRC circulants ARE such finite")
print("      unit-distance graphs; HN's spindle/rotation IS the LRC multiplier/perspective trick.")
