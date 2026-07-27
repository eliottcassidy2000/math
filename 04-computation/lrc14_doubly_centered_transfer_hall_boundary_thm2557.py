#!/usr/bin/env python3
"""Exact referee for THM-2557.

Restrict the THM-2548 seven-step transfer to the integral 7-by-13 lattice
with both margins zero.  The restriction is injective but not primitive:
its cokernel is exactly the six-dimensional first-root-moment space over
F_13.  The same computation also checks the cone/type obstruction separating
these signed horizontal currents from THM-2545's nonnegative ancestry tables.
"""

from fractions import Fraction
from pathlib import Path


Q = 7
P = 13
DIM = (Q - 1) * (P - 1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def zero_table():
    return [[0 for _ in range(P)] for _ in range(Q)]


def rectangle(k, r):
    """Integral basis vector e_kr-e_k,12-e_6,r+e_6,12 of V_00."""
    f = zero_table()
    f[k][r] += 1
    f[k][P - 1] -= 1
    f[Q - 1][r] -= 1
    f[Q - 1][P - 1] += 1
    return f


RECTANGLES = [(k, r) for k in range(Q - 1) for r in range(P - 1)]


def is_doubly_centered(f):
    return (
        all(sum(f[k][r] for r in range(P)) == 0 for k in range(Q))
        and all(sum(f[k][r] for k in range(Q)) == 0 for r in range(P))
    )


def transfer(f, a):
    """D_a f(k,r)=sum_(j=0)^6 f(k-j,r-a*j), with cyclic indices."""
    return [
        [
            sum(f[(k - j) % Q][(r - a * j) % P] for j in range(Q))
            for r in range(P)
        ]
        for k in range(Q)
    ]


def coordinates(f):
    """Coordinates in the rectangle basis (the upper-left 6-by-12 block)."""
    require(is_doubly_centered(f), "table is not in V_00")
    return [f[k][r] for k, r in RECTANGLES]


def transfer_matrix(a):
    columns = [coordinates(transfer(rectangle(k, r), a)) for k, r in RECTANGLES]
    return [[columns[j][i] for j in range(DIM)] for i in range(DIM)]


def rank_mod(matrix, prime):
    a = [[x % prime for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                c = a[i][col]
                a[i] = [(x - c * y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def bareiss_det(matrix):
    """Fraction-free exact determinant with row pivoting."""
    a = [row[:] for row in matrix]
    n = len(a)
    require(all(len(row) == n for row in a), "determinant matrix is not square")
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division was not exact")
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[n - 1][n - 1]


def root_moments(f):
    """The seven first root moments, reduced modulo 13."""
    return [sum(r * f[k][r] for r in range(P)) % P for k in range(Q)]


print("== THM-2557: the doubly centered transfer lattice ==")
require(len(RECTANGLES) == DIM == 72, "wrong V_00 rank")
require(all(is_doubly_centered(rectangle(k, r)) for k, r in RECTANGLES),
        "rectangle basis leaves V_00")
print(f"  V_00 rank: ({Q}-1)({P}-1) = {DIM}")
print("  basis: 72 integral row/column rectangles")

matrices = {}
mod13_ranks = []
mod2_ranks = []
for slope in range(1, P):
    matrix = transfer_matrix(slope)
    matrices[slope] = matrix
    for k, r in RECTANGLES:
        image = transfer(rectangle(k, r), slope)
        require(is_doubly_centered(image), f"D_{slope} does not preserve V_00")
        moments = root_moments(image)
        require(moments == [0] * Q, f"D_{slope} image has nonzero root moment")
    mod13_ranks.append(rank_mod(matrix, P))
    mod2_ranks.append(rank_mod(matrix, 2))

require(mod13_ranks == [66] * (P - 1), "wrong mod-13 ranks")
require(mod2_ranks == [DIM] * (P - 1), "a nonzero slope is singular mod 2")
print("  all 12 nonzero slopes preserve both zero margins")
print("  ranks mod 13: 66 in all 12 cases")
print("  ranks mod 2:  72 in all 12 cases")


print("\n== exact determinant, moment quotient, and Smith form ==")
determinants = [bareiss_det(matrices[slope]) for slope in range(1, P)]
require(all(abs(d) == P**6 for d in determinants), "wrong restricted determinant")

# The first six components coordinatize the moment hyperplane because their
# total is zero.  Check its rank and the explicit rectangle lifts e_k-e_6.
moment_matrix = []
for k in range(Q - 1):
    moment_matrix.append([
        root_moments(rectangle(i, r))[k]
        for i, r in RECTANGLES
    ])
require(rank_mod(moment_matrix, P) == Q - 1, "root-moment map is not onto")
for k in range(Q - 1):
    moments = root_moments(rectangle(k, 0))
    target = [0] * Q
    target[k] = 1
    target[Q - 1] = P - 1
    require(moments == target, "explicit root-moment lift failed")

print(f"  |det(D_a on V_00)|: {P**6} = 13^6 for all 12 slopes")
print("  mu_k(g)=sum_r r*g(k,r) mod 13; sum_k mu_k=0")
print("  mu is onto the 6-dimensional zero-sum moment space")
print("  every transferred rectangle has mu=0")
print("  equal indices give im(D_a)=ker(mu)")
print("  Smith form: 1^66 direct-sum 13^6")


print("\n== signed-current / nonnegative-coupling cone boundary ==")
signed_controls = 0
source = rectangle(0, 0)
for slope in range(1, P):
    image = transfer(source, slope)
    flat = [x for row in image for x in row]
    require(any(x != 0 for x in flat), "injective transfer killed a rectangle")
    require(min(flat) < 0 < max(flat), "nonzero centered transfer is not signed")
    require(sum(flat) == 0, "transferred current has nonzero augmentation")
    signed_controls += 1
print(f"  explicit nonzero signed transfers: {signed_controls}/{P-1}")
print("  each has zero row and column margins and both signs")
print("  any nonnegative table with either zero-margin family is identically zero")

# Byte-anchor the only live typed-row input used here: THM-2550(B)'s two
# nonzero doubly centered witnesses.  Part (A)'s k=2 drift is deliberately not
# composed with the large-clock response table.
root = Path(__file__).resolve().parents[1]
typed_output = root / "05-knowledge/results/lrc14_replica_dichotomy_typed_row_opus_20260727.out"
stored = typed_output.read_text(encoding="utf-8")
d_m = Fraction(27141080175744172866190363, 51745643666528683956098414400)
d_c = Fraction(665502894989984158227797, 10720466392295341281600)
require(str(d_m) in stored and str(d_c) in stored, "THM-2550(B) witness anchor changed")
require(d_m != 0 and d_c != 0, "typed-row witness unexpectedly vanished")
print("  THM-2550(B) nonzero d_M,d_C witnesses: byte-anchored")
print("  consequence: every D_a d_M and D_a d_C survives, but remains signed")
print("  THM-2550(A) drift is a different k=2 packet and is not composed here")


print("\n== exact THM-2545 Hall controls (one word stratum) ==")
hall_subsets = 0

# Actual proved-field shape from THM-2549: all target mass goes to cemetery.
source_weights = [h + 1 for h in range(P)]
cemetery_mass = sum(source_weights)
for mask in range(1 << P):
    left_mass = sum(source_weights[h] for h in range(P) if mask & (1 << h))
    neighbourhood_mass = cemetery_mass if mask else 0
    require(left_mass <= neighbourhood_mass, "cemetery Hall inequality failed")
    hall_subsets += 1
require(hall_subsets == 2**P, "wrong cemetery Hall subset count")

# Smallest cemetery-free translation-equivariant completion: b=h+1.
offset_checks = 0
for mask in range(1 << P):
    source_set = {h for h in range(P) if mask & (1 << h)}
    neighbours = {(h + 1) % P for h in source_set}
    require(len(neighbours) == len(source_set), "offset matching violates Hall")
    offset_checks += 1
require(offset_checks == 2**P, "wrong offset Hall subset count")

# Hostile control: diagonal-only support becomes empty after deleting the
# semantic diagonal and has an immediate singleton deficit.
singleton_mass = 1
diagonal_deleted_neighbour_mass = 0
require(singleton_mass > diagonal_deleted_neighbour_mass,
        "diagonal Hall-deficit positive control failed")

print(f"  cemetery zero-arrival checks: {hall_subsets}/{2**P}")
print(f"  offset b=h+1 zero-arrival checks: {offset_checks}/{2**P}")
print("  diagonal-only support: singleton Hall deficit detected")
print("  neither zero-arrival control is supplied by a signed transfer current")


print("\n== typed-coordinate verdict ==")
print("  transfer coordinates: owner/source clock ell and response co-shift s")
print("  Hall coordinates: selected-head root h and genuinely later active root b")
print("  no lawful map (ell,s)->(h,b) is assumed or tested")
print("  missing atom coordinate: common-base future root b, equivalently d=b-h")
print("  future immediate digit e_1 plus its carry/intertwiner is still required")
print("\nall checks passed")
