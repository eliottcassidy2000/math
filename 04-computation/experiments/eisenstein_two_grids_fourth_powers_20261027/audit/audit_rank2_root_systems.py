"""(H) check: perpendicular roots in the rank-2 irreducible root systems, and the minimal vectors of their root lattices.
Coordinates in the A2 / Z^2 planes, exact integer Gram forms."""
import itertools, math
# A2 and G2 in the Eisenstein basis (1, w): Q(a,b) = a^2 - ab + b^2; G2 = short roots (norm 1) + long roots (norm 3)
Qw = lambda a, b: a*a - a*b + b*b
Bw = lambda x, y: 2*x[0]*y[0] - x[0]*y[1] - x[1]*y[0] + 2*x[1]*y[1]          # = 2<x,y>
box = [(a, b) for a in range(-3, 4) for b in range(-3, 4) if (a, b) != (0, 0)]
A2 = [p for p in box if Qw(*p) == 1]
G2 = A2 + [p for p in box if Qw(*p) == 3]
B2 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
perp = lambda R, B: any(B(x, y) == 0 for x, y in itertools.combinations(R, 2))
dot = lambda x, y: x[0]*y[0] + x[1]*y[1]
print("A2:", len(A2), "roots, perpendicular pair:", perp(A2, Bw))
print("B2:", len(B2), "roots, perpendicular pair:", perp(B2, dot), "; root lattice Z^2, minimal vectors +-e_i (square grid: free by (A1))")
print("G2:", len(G2), "roots, perpendicular pair:", perp(G2, Bw),
      "; root lattice = Z[w] (short roots generate), minimal vectors = the 6 short roots = A2 -> pays floor(Q/3) (THM-4820)")
print("   perpendicular pairs in G2 are short-long only:", all(Qw(*x) != Qw(*y) for x, y in itertools.combinations(G2, 2) if Bw(x, y) == 0))
