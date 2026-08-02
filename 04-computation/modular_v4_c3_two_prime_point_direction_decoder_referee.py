#!/usr/bin/env python3
"""Exact referee for the V4-point/C3-direction two-prime decoder reflection."""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(matrix):
    a = [[Fraction(entry) for entry in row] for row in matrix]
    n = len(a)
    out = Fraction(1)
    for column in range(n):
        pivot = next((row for row in range(column, n) if a[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            out = -out
        value = a[column][column]
        out *= value
        for j in range(column, n):
            a[column][j] /= value
        for row in range(column + 1, n):
            scale = a[row][column]
            for j in range(column, n):
                a[row][j] -= scale * a[column][j]
    require(out.denominator == 1, "nonintegral determinant")
    return int(out)


def rank_mod(matrix, prime):
    a = [[entry % prime for entry in row] for row in matrix]
    rows = len(a)
    columns = len(a[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if a[row][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][column], -1, prime)
        a[rank] = [(inverse * value) % prime for value in a[rank]]
        for row in range(rows):
            if row != rank and a[row][column]:
                scale = a[row][column]
                a[row] = [(x - scale * y) % prime for x, y in zip(a[row], a[rank])]
        rank += 1
    return rank


def gcd_minors(matrix, size):
    value = 0
    for rows in combinations(range(len(matrix)), size):
        for columns in combinations(range(len(matrix[0])), size):
            minor = [[matrix[i][j] for j in columns] for i in rows]
            value = gcd(value, abs(determinant(minor)))
    return value


# Eisenstein integers a+b*w, w^2+w+1=0.
def eadd(x, y):
    return x[0] + y[0], x[1] + y[1]


def esum(values):
    total = (0, 0)
    for value in values:
        total = eadd(total, value)
    return total


def emul(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def eneg(x):
    return -x[0], -x[1]


def econj(x):
    a, b = x
    return a - b, -b


def enorm(x):
    a, b = x
    return a * a - a * b + b * b


def epow_w(power):
    return ((1, 0), (0, 1), (-1, -1))[power % 3]


def edeterminant(matrix):
    total = (0, 0)
    for perm in permutations(range(len(matrix))):
        inversions = sum(perm[i] > perm[j]
                         for i in range(len(perm)) for j in range(i + 1, len(perm)))
        term = (-1, 0) if inversions % 2 else (1, 0)
        for row, column in enumerate(perm):
            term = emul(term, matrix[row][column])
        total = eadd(total, term)
    return total


def add2(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def mat_vec2(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(2)) % 2 for i in range(2))


def det2(matrix):
    return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) % 2


points = tuple(product(range(2), repeat=2))
directions = tuple(point for point in points if point != (0, 0))

hadamard = [[(-1) ** ((chi[0] * x[0] + chi[1] * x[1]) % 2)
             for x in points] for chi in points]
require(determinant(hadamard) in (16, -16), "Walsh determinant changed")
minor_gcds = tuple(gcd_minors(hadamard, size) for size in range(1, 5))
require(minor_gcds == (1, 2, 4, 16), minor_gcds)
smith = (minor_gcds[0],
         minor_gcds[1] // minor_gcds[0],
         minor_gcds[2] // minor_gcds[1],
         minor_gcds[3] // minor_gcds[2])
require(smith == (1, 2, 2, 4), smith)
require(rank_mod(hadamard, 2) == 1, "characteristic-two point collapse changed")
require(all(rank_mod(hadamard, prime) == 4 for prime in (3, 5, 7)),
        "odd point decoder lost rank")

fourier3 = [[epow_w(row * column) for column in range(3)] for row in range(3)]
fourier_det = edeterminant(fourier3)
one_minus_w = (1, -1)
require(
    fourier_det == emul(emul(one_minus_w, one_minus_w), one_minus_w),
    (fourier_det, "expected (1-w)^3"),
)
require(enorm(fourier_det) == 27, (fourier_det, enorm(fourier_det)))
gram = [[esum(emul(econj(fourier3[k][i]), fourier3[k][j]) for k in range(3))
         for j in range(3)] for i in range(3)]
require(gram == [[(3 if i == j else 0, 0) for j in range(3)] for i in range(3)],
        "C3 Fourier Gram changed")
# Modulo the Eisenstein prime (1-w), w=1 in F3 and the matrix is all ones.
fourier_mod3 = [[1 for _column in range(3)] for _row in range(3)]
require(rank_mod(fourier_mod3, 3) == 1, "characteristic-three direction collapse changed")

joint_norm = abs(determinant(hadamard)) ** 6 * enorm(fourier_det) ** 4
require(joint_norm == 2**24 * 3**12, joint_norm)

# Every nonzero direction gives one K4 perfect matching x <-> x+d.
matchings = {}
for direction in directions:
    edges = frozenset(frozenset((x, add2(x, direction))) for x in points)
    require(len(edges) == 2 and set().union(*map(set, edges)) == set(points), direction)
    matchings[direction] = edges
require(len(set(matchings.values())) == 3, "direction/matching bijection changed")

flags = tuple((point, direction) for point in points for direction in directions)
flag_edges = tuple(frozenset((point, add2(point, direction)))
                   for point, direction in flags)
edge_multiplicities = {edge: flag_edges.count(edge) for edge in set(flag_edges)}
require(len(flags) == 12 and len(edge_multiplicities) == 6, "directed flag census changed")
require(set(edge_multiplicities.values()) == {2}, "flag-to-edge map is not two-to-one")

matrices = []
for entries in product(range(2), repeat=4):
    matrix = (entries[:2], entries[2:])
    if det2(matrix):
        matrices.append(matrix)
require(len(matrices) == 6, len(matrices))

affine_checks = 0
direction_permutations = set()
affine_permutations = set()
for matrix in matrices:
    direction_permutations.add(tuple(directions.index(mat_vec2(matrix, d)) for d in directions))
    for shift in points:
        image = {x: add2(mat_vec2(matrix, x), shift) for x in points}
        require(len(set(image.values())) == 4, "affine point action lost bijectivity")
        affine_permutations.add(tuple(points.index(image[x]) for x in points))
        for direction in directions:
            target_direction = mat_vec2(matrix, direction)
            edge_image = frozenset(frozenset((image[x] for x in edge))
                                   for edge in matchings[direction])
            require(edge_image == matchings[target_direction], "matching equivariance failed")
            affine_checks += 1
require(len(direction_permutations) == 6 and affine_checks == 72,
        (len(direction_permutations), affine_checks))
require(len(affine_permutations) == 24, "affine point action is not faithful")

# The marked C2*C3 quotient on the direction set.
S = ((0, 1), (1, 0))
T = ((0, 1), (1, 1))
def compose(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(2)) % 2
                       for j in range(2)) for i in range(2))
identity = ((1, 0), (0, 1))
require(compose(S, S) == identity, "C2 generator changed")
require(compose(compose(T, T), T) == identity and T != identity,
        "C3 generator changed")

generated = {identity}
while True:
    enlarged = generated | {
        compose(left, right)
        for left in generated
        for right in (S, T)
    }
    if enlarged == generated:
        break
    generated = enlarged
require(set(generated) == set(matrices), "marked generators do not generate GL(2,2)")

# Constant/augmentation splitting fails exactly at the cardinality primes.
point_constant_sum = len(points)
direction_constant_sum = len(directions)
require(tuple(p for p in (2, 3, 5, 7) if point_constant_sum % p == 0) == (2,),
        "point exceptional prime changed")
require(tuple(p for p in (2, 3, 5, 7) if direction_constant_sum % p == 0) == (3,),
        "direction exceptional prime changed")

print("MODULAR V4/C3 TWO-PRIME POINT-DIRECTION DECODER")
print(f"point_walsh_det={determinant(hadamard)} smith={smith} ranks_p2_p3={rank_mod(hadamard,2)},{rank_mod(hadamard,3)}")
print(f"direction_fourier_det={fourier_det} norm={enorm(fourier_det)} rank_at_1-w={rank_mod(fourier_mod3,3)}")
print(f"joint_decoder_norm=2^24*3^12={joint_norm}")
print(f"affine_group=V4_semidirect_S3;maps={len(matrices)*len(points)} matching_equivariance_checks={affine_checks}")
print("marked_free_factors=C2_order2,C3_order3;exceptional_primes=point:2,direction:3")
print("scope=abstract point-direction decoder only;no canonical physical PSL2Z,Keller,LRC,or quartic-origin map")
print("all_exact_checks=PASS")
