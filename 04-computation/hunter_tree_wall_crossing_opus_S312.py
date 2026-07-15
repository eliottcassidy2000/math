# opus-2026-07-15-S312 -- THE HUNTER TREE BOUND AT THE SEVEN-COMB WALL.
# Hunter/Kounias (1976): mu(union A_i) <= sum mu(A_i) - sum_{(i,j) in T} mu(A_i cap A_j)
# for ANY spanning tree T on the index set. Hence
#   uncovered = mu(E) - mu(union (D_i cap E))
#            >= mu(E) - sum mu(D_i cap E) + max_T sum_T mu(D_i cap D_j cap E).
# With densities 2/13 + excess and pairwise ~ (4/169)muE:
#   tree sum ~ (m'-1)*4/169 vs needed (2m'-13)/13 = (2m'-13)*13/169:
#   coercive iff 4(m'-1) >= 13(2m'-13) iff m' <= 7.5.  THE WALL CROSSES AT 7.
# This script: exact Hunter bound on real radius-7 packets vs actual uncovered.
from fractions import Fraction
import random, itertools

DELTA = Fraction(1, 13)

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def mu(ivs): return sum(b - a for (a, b) in ivs)

def comb_teeth_in(x, a, b):
    import math
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def restrict(E, x):
    out = []
    for (a, b) in E: out.extend(comb_teeth_in(x, a, b))
    return out

def subtract_comb(ivs, x):
    out = []
    for (a, b) in ivs:
        cur = a
        for (lo, hi) in sorted(comb_teeth_in(x, a, b)):
            if lo > cur: out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return [iv for iv in out if iv[0] < iv[1]]

def max_tree_sum(W, m):
    # Prim max spanning tree on complete graph with weights W[(i,j)]
    intree = {0}; tot = Fraction(0)
    while len(intree) < m:
        best = None
        for i in intree:
            for j in range(m):
                if j not in intree:
                    w = W[(min(i,j), max(i,j))]
                    if best is None or w > best[0]: best = (w, j)
        tot += best[0]; intree.add(best[1])
    return tot

P = [1, 2, 3, 4, 5]
R = [6, 7, 8, 9, 10, 11, 12]
E = safe_set(P)
muE = mu(E)
print(f"prefix {P}: mu(E) = {muE} ~ {float(muE):.5f}; needed tree-sum > "
      f"sum(mass_i) - mu(E)")
print(f"asymptotic: tree 6*(4/169)muE = {float(24*muE/169):.5f} vs needed "
      f"(1/13)muE = {float(muE/13):.5f} -- margin x{24/13:.2f}\n")

random.seed(7)
packets = [[r + 13*random.randint(2, 40) for r in R] for _ in range(4)]
packets.append([499 + i for i in range(7)])            # adversarial consecutive
packets.append([6+13*2, 7+13*2, 8+13*2, 9+13*2, 10+13*2, 11+13*2, 12+13*2])  # small lifts h=2

for xs in packets:
    Ei = [restrict(E, x) for x in xs]
    masses = [mu(e) for e in Ei]
    W = {}
    for i, j in itertools.combinations(range(7), 2):
        # pair overlap: intersect Ei[i] intervals with comb j
        tot = Fraction(0)
        for (a, b) in Ei[i]:
            for (lo, hi) in comb_teeth_in(xs[j], a, b): tot += hi - lo
        W[(i, j)] = tot
    tree = max_tree_sum(W, 7)
    hunter = muE - sum(masses) + tree
    U = list(E)
    for x in xs: U = subtract_comb(U, x)
    actual = mu(U)
    print(f"xs={xs}")
    print(f"   sum mass = {float(sum(masses)):.5f}  max-tree = {float(tree):.5f}  "
          f"HUNTER >= {float(hunter):+.5f}  actual = {float(actual):.5f}  "
          f"{'COERCIVE' if hunter > 0 else 'not coercive'}")
