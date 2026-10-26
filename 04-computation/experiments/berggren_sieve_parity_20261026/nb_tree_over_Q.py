# THM-4808 (G) over Q: the non-backtracking walk from the dart (1,0,1) -> (3,4,5) on the graph whose vertices are the
# Klein orbits of slopes and whose darts are slopes, reversal tau(t) = 2 - t, enumerates the Berggren tree.
# Checks, to depth DMAX: (i) the three NB continuations of a dart t are exactly U t, A t, D t; (ii) the head vertex of
# every dart t reached is the unsigned triangle Phi(t) (a primitive triple with V-orbit of 4 distinct slopes);
# (iii) no unsigned triangle is reached twice (tree); (iv) the degenerate vertex V.0 = {0, inf} carries the doubled
# edge to (3,4,5) (darts 2 twice) and a doubled loop (dart inf twice) in the multiset D(V.0) = [2, inf, 2, inf].
from fractions import Fraction as F
import sys
INF = 'inf'
def V(t):
    if t == INF or t == 0: return [0, INF, INF, 0] if False else [t, (0 if t == INF else INF), (INF if t == 0 else 0), t]
    return [t, -t, 1 / t, -1 / t]
def Dm(t): return INF if t == INF else t + 2
def U(t): return F(2) if t == INF else (INF if t == 0 else 2 - 1 / t)
def A(t): return F(2) if t == INF else (INF if t == 0 else 2 + 1 / t)
def tau(t): return INF if t == INF else 2 - t
def unsigned_triangle(t):
    t = F(t); m, n = t.numerator, t.denominator
    if m < 0: m = -m
    a, b, c = abs(m * m - n * n), 2 * abs(m * n), m * m + n * n
    from math import gcd
    g = gcd(gcd(a, b), c)
    return tuple(sorted((a // g, b // g))) + (c // g,)
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
ok_cont = ok_head = True; seen = {}; dup = 0; count = 0
frontier = [F(2)]                       # the dart from (1,0,1) into (3,4,5) is the slope t = 2
for depth in range(DMAX + 1):
    nxt = []
    for t in frontier:
        count += 1
        tri = unsigned_triangle(t)
        if tri in seen: dup += 1
        seen[tri] = t
        orbit = [t, -t, 1 / t, -1 / t]
        ok_head &= len(set(orbit)) == 4
        out = [Dm(x) for x in orbit]                       # darts leaving head(t)
        cont = [x for x in out if x != tau(t)]
        ok_cont &= sorted(cont) == sorted([U(t), A(t), Dm(t)])
        nxt += [U(t), A(t), Dm(t)]
    frontier = nxt
deg = sorted([Dm(0), Dm(INF), Dm(0), Dm(INF)], key=str)
print(f"depth <= {DMAX}: darts {count} (= (3^{DMAX+1}-1)/2: {count == (3**(DMAX+1)-1)//2}); (i) NB continuations = U,A,D: {ok_cont}; "
      f"(ii) generic head orbits: {ok_head}; (iii) repeated triangles: {dup}; (iv) degenerate out-multiset D(V.0) = {deg}, "
      f"tau(2) = {tau(F(2))} (the root's reverse dart points to (1,0,1))")
