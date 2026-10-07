# FINITE-EXACT: unit-distance counts of triangular-lattice sections at the sqrt(7) layer.
# Points of Z[w] (w = e^{i pi/3}) as integer pairs (a,b) <-> a + b w ; norm N(a,b) = a^2 + a b + b^2.
from itertools import combinations
def mul(x, y):  # (a+bw)(c+dw) with w^2 = w - 1
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def nrm(x): a, b = x; return a*a + a*b + b*b
def edges(pts, D):
    return sum(1 for p, q in combinations(pts, 2) if nrm((p[0]-q[0], p[1]-q[1])) == D)
units = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
alpha, alphab = (2,1), (3,-1)        # 2+w and its conjugate 2+w^5 = 2 + (1-w) = 3 - w ; both norm 7
assert nrm(alpha) == 7 and nrm(alphab) == 7 and alpha != alphab
tri  = [(0,0),(1,0),(0,1)]                       # unit triangle
hex7 = [(0,0)] + units                           # hexagon + centre (7 pts, 12 edges)
print("check: e(tri) =", edges(tri, 1), " e(hex7) =", edges(hex7, 1))
def msum(P, Q):
    return sorted({add(mul(alpha, p), mul(alphab, q)) for p in P for q in Q})
for name, P, Q in (("K3 x K3 (N=9)", tri, tri), ("K3 x H7 (N=21)", tri, hex7), ("H7 x H7 (N=49)", hex7, hex7)):
    S = msum(P, Q)
    print(f"{name}: |S| = {len(S)}  sqrt7-unit distances = {edges(S, 7)}   (Cartesian-product prediction {len(P)*edges(Q,1)+len(Q)*edges(P,1)})")
# compare with Harborth penny numbers floor(3N - sqrt(12N-3)) and the AMP24 table u(N)
import math
uN = {9:18, 21:57}
for N in (9, 21, 49):
    print(f"N={N}: Harborth penny (D=1 sections) = {math.floor(3*N - math.sqrt(12*N-3))}", f"; u(N) = {uN[N]}" if N in uN else "")
# all pairwise distances realized in K3 x H7, as a sanity check on extra coincidences
S = msum(tri, hex7)
from collections import Counter
print("distance-norm census of K3 x H7:", sorted(Counter(nrm((p[0]-q[0], p[1]-q[1])) for p, q in combinations(S, 2)).items())[:12])
# Euclidean embedding check: coordinates
import cmath
w = cmath.exp(1j*math.pi/3)
Z = [complex(a + b*w.real, b*w.imag) for a, b in S]
cnt = sum(1 for z1, z2 in combinations(Z, 2) if abs(abs(z1 - z2) - math.sqrt(7)) < 1e-9)
print("float check: pairs at Euclidean distance sqrt(7):", cnt)
