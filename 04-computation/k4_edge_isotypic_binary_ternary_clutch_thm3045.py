"""Exact referee for THM-3045.

The six coordinates are the edges 01,02,03,12,13,23 of K4.  We check the
orthogonal rational 1+[22]+[31] decomposition, its integral Smith defect,
the explicit (F2[matchings],F3) quotient, all S4 equivariance, and the exact
projector-integrality criterion.  Truth-bearing checks use explicit raises,
so ``python -O`` executes the same verification.
"""

from itertools import permutations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
PAIRS = ((0, 5), (1, 4), (2, 3))


def row(entries):
    return sp.Matrix(1, 6, entries)


c = row((1, 1, 1, 1, 1, 1))
p = []
n = []
for i, j in PAIRS:
    pv = [0] * 6
    nv = [0] * 6
    pv[i] = pv[j] = 1
    nv[i], nv[j] = 1, -1
    p.append(row(pv))
    n.append(row(nv))

aug = [p[0] - p[2], p[1] - p[2]]
basis_rows = [c, *aug, *n]
A = sp.Matrix.vstack(*basis_rows)

require(A.rank() == 6, "isotypic intersection lattice is not full rank")
require(abs(int(A.det())) == 24, "integral clutch index is not 24")
D = smith_normal_form(A, domain=ZZ)
smith = tuple(abs(int(D[i, i])) for i in range(6))
require(smith == (1, 1, 1, 2, 2, 6), "unexpected Smith form")

gram = A * A.T
require(gram[:1, 1:] == sp.zeros(1, 5), "constant line is not orthogonal")
require(gram[1:3, 3:] == sp.zeros(2, 3), "[22] and [31] are not orthogonal")
require(int(gram[:1, :1].det()) == 6, "constant Gram determinant changed")
require(int(gram[1:3, 1:3].det()) == 12, "A2 Gram determinant changed")
require(int(gram[3:, 3:].det()) == 8, "opposition-minus Gram determinant changed")

I6 = sp.eye(6)
O = sp.zeros(6)
for i, j in PAIRS:
    O[i, j] = O[j, i] = 1
Pminus = (I6 - O) / 2
Pplus = (I6 + O) / 2
Pzero = sp.ones(6, 6) / 6
Ptwo = Pplus - Pzero
projectors = (Pzero, Ptwo, Pminus)
for P in projectors:
    require(P * P == P, "projector is not idempotent")
require(sum(projectors, sp.zeros(6)) == I6, "projectors do not sum to identity")
for i in range(3):
    for j in range(i + 1, 3):
        require(projectors[i] * projectors[j] == sp.zeros(6), "projectors are not orthogonal")
require(tuple(P.rank() for P in projectors) == (1, 2, 3), "wrong isotypic ranks")


def phi(values):
    bits = tuple((values[i] + values[j]) % 2 for i, j in PAIRS)
    ternary = sum(values) % 3
    return bits, ternary


def projections_integral(values):
    v = sp.Matrix(values)
    return all(all(q.q == 1 for q in P * v) for P in projectors)


images = set()
kernel_count = 0
for values in product(range(6), repeat=6):
    image = phi(values)
    images.add(image)
    in_kernel = image == ((0, 0, 0), 0)
    require(projections_integral(values) == in_kernel, "projector/kernel criterion failed")
    if in_kernel:
        kernel_count += 1
require(len(images) == 24, "quotient map is not onto 24 classes")
require(kernel_count == 6**6 // 24, "wrong residue-kernel size")


def act_edge_vector(g, values):
    out = [0] * 6
    for i, edge in enumerate(EDGES):
        target = tuple(sorted((g[edge[0]], g[edge[1]])))
        out[EDGE_INDEX[target]] = values[i]
    return tuple(out)


def matching_action(g):
    ans = []
    pair_sets = [frozenset((EDGES[i], EDGES[j])) for i, j in PAIRS]
    for pair_set in pair_sets:
        moved = frozenset(tuple(sorted((g[a], g[b]))) for a, b in pair_set)
        ans.append(pair_sets.index(moved))
    return tuple(ans)


kernel_matching = []
for g in permutations(range(4)):
    mg = matching_action(g)
    if mg == (0, 1, 2):
        kernel_matching.append(g)
    for i in range(6):
        v = [0] * 6
        v[i] = 1
        before = phi(v)
        after = phi(act_edge_vector(g, v))
        expected_bits = [0, 0, 0]
        for old, new in enumerate(mg):
            expected_bits[new] = before[0][old]
        require(after == (tuple(expected_bits), before[1]), "S4 quotient equivariance failed")
require(len(kernel_matching) == 4, "matching action kernel is not V4")

# The C2 and C3 controls on the matching base.
transposition = (1, 0, 2, 3)
three_cycle = (1, 2, 0, 3)
mt = matching_action(transposition)
mc = matching_action(three_cycle)
require(sum(mt[i] == i for i in range(3)) == 1, "C2 does not fix one matching")
require(sum(mc[i] == i for i in range(3)) == 0, "C3 does not cycle the matchings")
require(tuple(mt[mt[i]] for i in range(3)) == (0, 1, 2), "C2 matching action is not order two")
require(tuple(mc[mc[mc[i]]] for i in range(3)) == (0, 1, 2), "C3 matching action is not order three")

print("THM3045 K4 EDGE ISOTYPIC INTEGRAL CLUTCH: PASS")
print("edge_order=01,02,03,12,13,23")
print("rational_blocks=1+[22]+[31] ranks=1,2,3")
print(f"orthogonal_gram_determinants={6},{12},{8}")
print(f"intersection_lattice_index={abs(int(A.det()))}")
print("smith=" + ",".join(map(str, smith)))
print("quotient=F2[three_matchings]+F3_trivial size=24")
print(f"residue_universe={6**6} kernel={kernel_count} images={len(images)}")
print("binary_map=kappa_m=x_e+x_opposite(e) mod2")
print("ternary_map=tau=sum_edges(x_e) mod3")
print("integral_all_projectors iff kappa=0 and tau=0")
print(f"S4_matching_kernel_size={len(kernel_matching)}")
print(f"C2_matching_action={mt} C3_matching_action={mc}")
print("prime_support=2(binary matching clutch),3(ternary scalar clutch)")
print("scope=no tree identification, no affine owner, no Keller/LRC exclusion")
