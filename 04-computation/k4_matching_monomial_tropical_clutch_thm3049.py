"""Exact referee for THM-3049.

The script checks the K4 matching-monomial cocharacter map, its exact
binary/ternary integral clutch, S4 equivariance, and sharp field-root and
matching-fibre cancellation boundaries.  Truth-bearing checks use explicit
raises, so ``python -O`` executes the same verification.
"""

from fractions import Fraction
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


# Matching-monomial exponent/cocharacter matrix.
B = sp.zeros(3, 6)
for m, (i, j) in enumerate(PAIRS):
    B[m, i] = B[m, j] = 1
require(B.rank() == 3, "matching-monomial cocharacter map is not onto")
require(B[:, (0, 1, 2)] == sp.eye(3), "chosen monomial section changed")

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

N = sp.Matrix.vstack(*n)
require(B * N.T == sp.zeros(3), "opposite-edge differences leave monomials")
require(N.rank() == 3, "monomial kernel rank changed")

aug = [p[0] - p[2], p[1] - p[2]]
lattice_basis = sp.Matrix.vstack(c, *aug, *n)
require(abs(int(lattice_basis.det())) == 24, "clutch index changed")
D = smith_normal_form(lattice_basis, domain=ZZ)
smith = tuple(abs(int(D[i, i])) for i in range(6))
require(smith == (1, 1, 1, 2, 2, 6), "unexpected clutch Smith form")


def matching_values(values):
    return tuple(int(z) for z in B * sp.Matrix(values))


def clutch(values):
    lam = matching_values(values)
    return tuple(z % 2 for z in lam), sum(lam) % 3


def projectors_integral(values):
    pair_sums = matching_values(values)
    return all(z % 2 == 0 for z in pair_sums) and sum(pair_sums) % 3 == 0


images = set()
kernel_count = 0
for values in product(range(6), repeat=6):
    image = clutch(values)
    images.add(image)
    in_kernel = image == ((0, 0, 0), 0)
    require(projectors_integral(values) == in_kernel, "root/projector criterion failed")
    kernel_count += int(in_kernel)
require(len(images) == 24, "clutch map is not onto")
require(kernel_count == 6**6 // 24, "wrong clutch kernel size")


def act_edge_vector(g, values):
    out = [0] * 6
    for i, edge in enumerate(EDGES):
        target = tuple(sorted((g[edge[0]], g[edge[1]])))
        out[EDGE_INDEX[target]] = values[i]
    return tuple(out)


def matching_action(g):
    pair_sets = [frozenset((EDGES[i], EDGES[j])) for i, j in PAIRS]
    ans = []
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
        before = clutch(v)
        after = clutch(act_edge_vector(g, v))
        expected = [0, 0, 0]
        for old, new in enumerate(mg):
            expected[new] = before[0][old]
        require(after == (tuple(expected), before[1]), "S4 clutch equivariance failed")
require(len(kernel_matching) == 4, "matching-action kernel is not V4")

# Independent binary and ternary value-lattice controls.
controls = {
    (2, 0, 0): ((0, 0, 0), 2),
    (1, 1, 1): ((1, 1, 1), 0),
}
for lam, expected in controls.items():
    values = (lam[0], lam[1], lam[2], 0, 0, 0)
    require(matching_values(values) == lam, "monomial section failed")
    require(clutch(values) == expected, "independence control changed")

# Zero valuation does not imply field roots: 3 is neither square nor cube mod 7.
nonzero_squares_mod7 = {pow(a, 2, 7) for a in range(1, 7)}
nonzero_cubes_mod7 = {pow(a, 3, 7) for a in range(1, 7)}
require(3 not in nonzero_squares_mod7, "Q7 square hostile failed")
require(3 not in nonzero_cubes_mod7, "Q7 cube hostile failed")

# The same zero valuation/clutch data can give a nonzero or zero fibre sum.
amplitude_live = (Fraction(1), Fraction(1), Fraction(1))
amplitude_cancel = (Fraction(1), Fraction(1), Fraction(-2))
require(sum(amplitude_live) == 3, "live amplitude control changed")
require(sum(amplitude_cancel) == 0, "cancellation hostile changed")
require(tuple(0 for _ in amplitude_live) == tuple(0 for _ in amplitude_cancel),
        "valuation-profile control changed")

print("THM3049 K4 MATCHING-MONOMIAL TROPICAL CLUTCH: PASS")
print("edge_order=01,02,03,12,13,23")
print("matching_map=lambda=(x01+x23,x02+x13,x03+x12)")
print("cocharacter_rank=3 kernel=opposite_edge_differences")
print(f"intersection_index={abs(int(lattice_basis.det()))}")
print("smith=" + ",".join(map(str, smith)))
print("quotient=F2[three_matchings]+F3_trivial size=24")
print(f"residue_universe={6**6} kernel={kernel_count} images={len(images)}")
print("square_value_obstruction=lambda mod2")
print("cube_total_value_obstruction=sum(lambda) mod3")
print(f"S4_matching_kernel_size={len(kernel_matching)}")
print("independence=lambda200->(000,2),lambda111->(111,0)")
print("Q7_unit_hostile=3 is neither square nor cube mod7")
print("fibre_sum_hostile=(1,1,1)->3 versus (1,1,-2)->0")
print("scope=uncontracted labelled monomials/cocharacters only")
print("no field-root, coloured-response, Keller, or LRC conclusion")
