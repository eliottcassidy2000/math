#!/usr/bin/env python3
"""Exact companion for THM-3742.

The audit keeps three different projective objects separate:

* ordinary projectivization of a two-coordinate Pell state, which forgets
  the central sign and therefore has two seven-cycles;
* stereographic (half-angle) parametrization of the norm-one conic, which is
  a bijection with P^1(F_13) and retains all fourteen conic states; and
* a scalar linear observer, whose image has only seven or eight residues.

All finite-field universes are exhausted.  Integral identities are checked by
exact matrix algebra and on the declared first thirteen states of each Pell
coset.  No Python ``assert`` is used, so normal and optimized runs have the
same gates.
"""

from __future__ import annotations

from collections import Counter


P = 13
Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]
Projective = int | None  # None denotes infinity.


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def mat_mul(left: Matrix, right: Matrix, modulus: int | None = None) -> Matrix:
    entries = (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )
    if modulus is None:
        return entries
    return tuple(tuple(entry % modulus for entry in row) for row in entries)  # type: ignore[return-value]


def mat_pow(matrix: Matrix, exponent: int, modulus: int | None = None) -> Matrix:
    require(exponent >= 0, "nonnegative matrix exponent")
    result: Matrix = ((1, 0), (0, 1))
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = mat_mul(result, base, modulus)
        base = mat_mul(base, base, modulus)
        power //= 2
    return result


def transpose(matrix: Matrix) -> Matrix:
    return ((matrix[0][0], matrix[1][0]), (matrix[0][1], matrix[1][1]))


def mat_vec(matrix: Matrix, vector: Vector, modulus: int | None = None) -> Vector:
    result = (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )
    if modulus is None:
        return result
    return result[0] % modulus, result[1] % modulus


def determinant(matrix: Matrix, modulus: int | None = None) -> int:
    value = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return value if modulus is None else value % modulus


def inverse_mod(matrix: Matrix) -> Matrix:
    det_inverse = pow(determinant(matrix, P), -1, P)
    return (
        (matrix[1][1] * det_inverse % P, -matrix[0][1] * det_inverse % P),
        (-matrix[1][0] * det_inverse % P, matrix[0][0] * det_inverse % P),
    )


def scalar_mul(scalar: int, matrix: Matrix, modulus: int = P) -> Matrix:
    return tuple(tuple(scalar * entry % modulus for entry in row) for row in matrix)  # type: ignore[return-value]


def projectively_equal(left: Matrix, right: Matrix) -> bool:
    for i in range(2):
        for j in range(2):
            if right[i][j] % P:
                scale = left[i][j] * pow(right[i][j], -1, P) % P
                return left == scalar_mul(scale, right)
    return False


def mobius(matrix: Matrix, value: Projective) -> Projective:
    a, b = matrix[0]
    c, d = matrix[1]
    if value is None:
        if c % P == 0:
            return None
        return a * pow(c, -1, P) % P
    numerator = (a * value + b) % P
    denominator = (c * value + d) % P
    if denominator == 0:
        return None
    return numerator * pow(denominator, -1, P) % P


def fmt_projective(value: Projective) -> str:
    return "inf" if value is None else str(value)


def quadratic_character(value: int) -> int:
    value %= P
    if value == 0:
        return 0
    return 1 if pow(value, (P - 1) // 2, P) == 1 else -1


def q2(vector: Vector) -> int:
    x, s = vector
    return x * x - 2 * s * s


def q8_mod(vector: Vector) -> int:
    x, q = vector
    return (x * x - 8 * q * q) % P


def orbit(matrix: Matrix, seed: Vector) -> tuple[Vector, ...]:
    points: list[Vector] = []
    current = seed
    while current not in points:
        points.append(current)
        current = mat_vec(matrix, current, P)
    require(current == seed, "orbit must return to its seed")
    return tuple(points)


checks = 0


print("=== Integral signed Pell interlacing ===")
S: Matrix = ((1, 2), (1, 1))
A: Matrix = ((3, 4), (2, 3))
J: Matrix = ((1, 0), (0, -2))
minus_J = ((-1, 0), (0, 2))
require(mat_mul(transpose(S), mat_mul(J, S)) == minus_J, "S reverses the norm")
require(mat_mul(S, S) == A, "A=S^2")
require(mat_mul(transpose(A), mat_mul(J, A)) == J, "A preserves the norm")
checks += 3

positive = (1, 0)
negative = (1, 1)
positive_rows: list[tuple[int, int, int, int]] = []
negative_rows: list[tuple[int, int, int, int]] = []
for depth in range(13):
    x, s = positive
    require(q2(positive) == 1, "positive Pell norm")
    require(x % 2 == 1 and s % 2 == 0, "square-triangular parity")
    n, q = (x - 1) // 2, s // 2
    require(n * (n + 1) // 2 == q * q, "square-triangular identity")
    positive_rows.append((depth, n, q, x))

    xn, sn = negative
    require(q2(negative) == -1, "negative Pell norm")
    require(xn % 2 == 1 and sn % 2 == 1, "triangular-tower parity")
    nn, mn = (xn - 1) // 2, (sn - 1) // 2
    require(nn * (nn + 1) // 2 == mn * (mn + 1), "T_n=2T_m identity")
    require(negative == mat_vec(S, positive), "signed cosets interlace by S")
    negative_rows.append((depth, nn, mn, sn))
    checks += 7

    positive = mat_vec(A, positive)
    negative = mat_vec(A, negative)

print("S=[[1,2],[1,1]] reverses x^2-2s^2; A=S^2=[[3,4],[2,3]] preserves it.")
print("A^k(1,0) gives T_n=q^2; A^k(1,1)=S A^k(1,0) gives T_n=2T_m.")
print("square-triangular first q:", [row[2] for row in positive_rows[:6]])
print("triangular-tower first odd companions:", [row[3] for row in negative_rows[:6]])
print("Integral bounded replay: depths 0..12 on both norm cosets.")


print("\n=== Complete mod-13 norm-fibre audit ===")
M: Matrix = ((3, 8), (1, 3))
minus_identity: Matrix = ((P - 1, 0), (0, P - 1))
identity: Matrix = ((1, 0), (0, 1))
require(mat_pow(A, 7, P) == minus_identity, "A^7=-I")
require(mat_pow(A, 14, P) == identity, "A^14=I")
require(mat_pow(M, 7, P) == minus_identity, "M^7=-I")
require(mat_pow(M, 14, P) == identity, "M^14=I")
require(mat_pow(S, 7, P) == ((5, 0), (0, 5)), "S^7=5I")
require(mat_pow(S, 14, P) == minus_identity, "S^14=-I")
require(mat_pow(S, 28, P) == identity, "S^28=I")
checks += 7

squares = {value * value % P for value in range(P)}
require(2 not in squares and 8 not in squares and 5 not in squares, "2, 8, and 5 nonsquare")
require((-1) % P in squares, "-1 is square")
checks += 2

all_vectors = {(x, q) for x in range(P) for q in range(P)}
fiber_sizes: list[int] = []
projective_by_character: dict[int, set[Projective]] = {1: set(), -1: set()}
for norm in range(1, P):
    fiber = {(x, q) for x, q in all_vectors if q8_mod((x, q)) == norm}
    require(len(fiber) == 14, f"norm fibre {norm} has fourteen points")
    seed = min(fiber)
    fibre_orbit = set(orbit(M, seed))
    require(fibre_orbit == fiber, f"M sharply transitive on norm fibre {norm}")
    ordered = orbit(M, seed)
    require(all(ordered[(index + 7) % 14] == ((-ordered[index][0]) % P, (-ordered[index][1]) % P)
                for index in range(14)), "central sign occurs after seven steps")
    slopes = {None if x == 0 else q * pow(x, -1, P) % P for x, q in fiber}
    require(len(slopes) == 7, "ordinary projectivization is two-to-one")
    character = quadratic_character(norm)
    if not projective_by_character[character]:
        projective_by_character[character] = slopes
    require(projective_by_character[character] == slopes, "norm character determines projective orbit")
    fiber_sizes.append(len(fiber))
    checks += 5

require(projective_by_character[1].isdisjoint(projective_by_character[-1]), "two projective cycles disjoint")
require(projective_by_character[1] | projective_by_character[-1] == set(range(P)) | {None}, "cycles cover P1")
checks += 2
print("M^7=-I, M^14=I; all 12 nonzero norm fibres have 14 points and one M-orbit.")
print("Direct [x:q] quotient (reported as q/x) has two seven-cycles:")
print("  square norm:", [fmt_projective(v) for v in sorted(projective_by_character[1], key=lambda z: (z is None, z))])
print("  nonsquare norm:", [fmt_projective(v) for v in sorted(projective_by_character[-1], key=lambda z: (z is None, z))])
print("Norm +1 and norm -1 share the square-character cycle because -1 is square mod 13.")

active_offsets = {1, 2, -3 % P, -2 % P, -1 % P}
require(active_offsets == projective_by_character[-1] - {None, 3}, "THM-3713 five-offset near match")
checks += 1
print("THM-3713's {1,2,-3,-2,-1} equals the nonsquare cycle minus {inf,3}.")
print("That exact near-match lies on the wrong norm-character torsor.")


print("\n=== Stereographic half-angle P1 clock ===")
conic_plus = {(x, q) for x in range(P) for q in range(P) if q8_mod((x, q)) == 1}
require(len(conic_plus) == 14, "norm-one q-conic size")


def theta(point: Vector) -> Projective:
    x, q = point
    if (x + 1) % P == 0:
        require(q == 0, "unique theta point at infinity")
        return None
    return q * pow(x + 1, -1, P) % P


def theta_inverse(value: Projective) -> Vector:
    if value is None:
        return P - 1, 0
    denominator = (1 - 8 * value * value) % P
    require(denominator != 0, "half-angle denominator nonzero because 8 is nonsquare")
    x = (1 + 8 * value * value) * pow(denominator, -1, P) % P
    q = 2 * value * pow(denominator, -1, P) % P
    return x, q


parameters = {theta(point) for point in conic_plus}
require(parameters == set(range(P)) | {None}, "theta is onto P1")
for parameter in parameters:
    require(theta(theta_inverse(parameter)) == parameter, "theta inverse on P1")
    require(theta_inverse(parameter) in conic_plus, "theta inverse on conic")
    checks += 2

N: Matrix = ((4, 1), (8, 4))
for point in conic_plus:
    require(theta(mat_vec(M, point, P)) == mobius(N, theta(point)), "theta conjugates M to N")
    checks += 1

n_cycle: list[Projective] = []
current: Projective = 0
while current not in n_cycle:
    n_cycle.append(current)
    current = mobius(N, current)
require(current == 0 and len(n_cycle) == 14, "N is a single P1 fourteen-cycle")
require(set(n_cycle) == set(range(P)) | {None}, "N cycle covers P1")
antipode: Matrix = ((0, 1), (8, 0))
require(projectively_equal(mat_pow(N, 7, P), antipode), "N^7 is conic antipode in theta coordinates")
for point in conic_plus:
    negative_point = ((-point[0]) % P, (-point[1]) % P)
    require(theta(negative_point) == mobius(antipode, theta(point)), "theta retains central sign")
    checks += 1
checks += 3
print("theta(x,q)=q/(x+1), theta(-1,0)=inf, bijects x^2-8q^2=1 with P1(F13).")
print("theta M theta^-1 is N=[[4,1],[8,4]], det(N)=8 nonsquare.")
print("single 14-cycle:", [fmt_projective(value) for value in n_cycle])
print("N^7 is t -> 1/(8t), the antipodal state; theta is not ordinary vector projectivization.")


print("\n=== Existing mod-13 clock conjugacies ===")
g6: Matrix = ((0, 1), (-1 % P, 6))
conjugator: Matrix = ((11, 5), (5, 0))
require(determinant(conjugator, P) == 1, "g6 conjugator lies in SL2")
require(mat_mul(conjugator, g6, P) == mat_mul(M, conjugator, P), "P g6=M P")

singer: Matrix = ((1, 4), (2, 1))
p0: Matrix = ((0, 1), (2, 0))
singer_expression = scalar_mul(11, mat_mul(p0, mat_mul(mat_pow(singer, 5, P), inverse_mod(p0), P), P))
require(singer_expression == N, "half-angle clock is conjugate to Singer fifth power")
checks += 3
print("M is SL2-conjugate to THM-2619 g_6 by P=[[11,5],[5,0]].")
print("N=11 P0 Singer^5 P0^-1 for Singer=[[1,4],[2,1]], P0=[[0,1],[2,0]].")
print("Thus both clocks refine existing P1 mechanisms; the integral Pell typing is the new sidecar.")

reflection_j: Matrix = ((1, 0), (0, -1 % P))
reflection_r = mat_mul(inverse_mod(M), reflection_j, P)
require(reflection_r == ((3, 8), (12, 10)), "explicit Pell reflection")
require(mat_pow(reflection_r, 2, P) == identity, "Pell reflection is involutive")
require(mat_mul(reflection_r, mat_mul(M, reflection_r, P), P) == inverse_mod(M), "RMR=M^-1")
phase_orbit = orbit(M, (1, 0))
for phase, point in enumerate(phase_orbit):
    require(mat_vec(reflection_r, point, P) == phase_orbit[(-1 - phase) % 14], "dihedral phase law")
    checks += 1
checks += 3
print("R=M^-1 diag(1,-1)=[[3,8],[12,10]] obeys R^2=I and RMR=M^-1.")
print("On the norm-one orbit it sends phase r to -1-r mod 14 (and mod 7 after sign quotient).")


print("\n=== Complete scalar linear-observer classification ===")
observer_counts: dict[int, Counter[int]] = {}
for norm in range(1, P):
    fiber = {(x, q) for x in range(P) for q in range(P) if q8_mod((x, q)) == norm}
    counts: Counter[int] = Counter()
    for a in range(P):
        for b in range(P):
            if (a, b) == (0, 0):
                continue
            image = {(a * x + b * q) % P for x, q in fiber}
            require(len(image) in (7, 8), "linear observer image is seven or eight")
            criterion = quadratic_character(norm * (a * a - 5 * b * b)) == 1
            require((len(image) == 8) == criterion, "quadratic-character observer criterion")
            for value in range(P):
                fibre_multiplicity = sum((a * x + b * q) % P == value for x, q in fiber)
                predicted = 1 - quadratic_character(value * value - norm * (a * a - 5 * b * b))
                require(fibre_multiplicity == predicted, "exact observer fibre multiplicity")
                checks += 1
            counts[len(image)] += 1
            checks += 2
    require(counts == Counter({7: 84, 8: 84}), "balanced observer count on each fibre")
    observer_counts[norm] = counts
    checks += 1
print("For every e!=0 and L=a*x+b*q!=0 on x^2-8q^2=e:")
print("  #L^-1(y)=1-chi(y^2-e(a^2-5b^2));")
print("  |L(C_e)|=8 iff chi(e(a^2-5b^2))=+1; otherwise |L(C_e)|=7.")
print("Each fibre has 84 observers of each size (12*168 observer/fibre pairs exhausted).")


print("\n=== Square/triangular scalar residue hostile ===")
square_residues = {value * value % P for value in range(P)}
triangular_residues = {value * (value + 1) * pow(2, -1, P) % P for value in range(P)}
missing = set(range(P)) - square_residues - triangular_residues
require(square_residues == {0, 1, 3, 4, 9, 10, 12}, "square residues")
require(triangular_residues == {0, 1, 2, 3, 6, 8, 10}, "triangular residues")
require(missing == {5, 7, 11}, "square/triangular union misses three residues")
checks += 3
print("squares:", sorted(square_residues))
print("triangular numbers:", sorted(triangular_residues))
print("their union misses:", sorted(missing))


print("\n=== Constant-Cohn nonentry gate ===")
a, b = M[0]
c, d = M[1]
epsilon = 2 * a - d
require((a, b, c, d, epsilon, epsilon * epsilon - 1) == (3, 8, 1, 3, 3, 8), "Cohn specialization")
lower_constant_solutions = [h for h in range(-20, 21) if c * h == d - a and (4 * a - d) * h == b]
upper_constant_solutions = [h for h in range(-20, 21) if (a - d) * h == -c and b * h == 4 * a - d]
require(not lower_constant_solutions and not upper_constant_solutions, "no constant exposed closure")
require(b != 0 and c != 0, "THM-3736 nonconstant triangularity hypotheses fail")
checks += 3
print("R=M has (a,b,c,d)=(3,8,1,3), epsilon=3, epsilon^2-1=8.")
print("THM-3726 constant equations have no solution; b,c!=0 triggers THM-3736 nonconstant nonentry.")


print("\n=== Scope ===")
print("The Pell update supplies an exact typed C14 conic clock and a half-angle P1 address.")
print("It loses owner, root chart, temporal word, target positivity, and every JC polynomial target packet.")
print("No LRC(14), planar Jacobian, or constant-Cohn closure follows.")
print(f"CHECKS={checks}")
