# FINITE-EXACT: Minkowski-product sections of the triangular lattice Z[w] at multi-class norms.
from itertools import combinations, product
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def nrm(x): a, b = x; return a*a + a*b + b*b
def conj(x): a, b = x; return (a + b, -b)       # conj(a + b w) = a + b w^5 = a + b(1 - w)
def edges(pts, D): return sum(1 for p, q in combinations(pts, 2) if nrm((p[0]-q[0], p[1]-q[1])) == D)
units = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
tri = [(0,0),(1,0),(0,1)]; hex7 = [(0,0)] + units; rhomb = [(0,0),(1,0),(0,1),(1,1)]
alpha = (2,1); beta = (3,1)                         # norms 7, 13
assert nrm(alpha) == 7 and nrm(beta) == 13
A, Ab, B, Bb = alpha, conj(alpha), beta, conj(beta)
dirs91 = [mul(A, B), mul(A, Bb), mul(Ab, B), mul(Ab, Bb)]
print("norm-91 directions:", dirs91, [nrm(d) for d in dirs91])
def msum(factors, dirs):
    pts = set()
    for choice in product(*factors):
        z = (0, 0)
        for d, p in zip(dirs, choice): z = add(z, mul(d, p))
        pts.add(z)
    return sorted(pts)
for name, factors, dirs, D in (
        ("K3^3 at D=91", [tri]*3, dirs91[:3], 91),
        ("K3^2 x W7 at D=91", [tri, tri, hex7], dirs91[:3], 91),
        ("K3^4 at D=91", [tri]*4, dirs91, 91),
        ("K3 x W7 at D=7", [tri, hex7], [A, Ab], 7)):
    S = msum(factors, dirs)
    n = 1
    for f in factors: n *= len(f)
    pred = sum(edges(f, 1) * (n // len(f)) for f in factors)
    print(f"{name}: |S| = {len(S)} (product {n}); unit distances = {edges(S, D)}; Cartesian prediction = {pred}; 3N = {3*len(S)}")
