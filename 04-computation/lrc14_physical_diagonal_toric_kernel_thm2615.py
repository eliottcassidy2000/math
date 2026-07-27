#!/usr/bin/env python3
"""Exact finite referee for THM-2615.

The analytic group-algebra kernel and positivity statements are proved in the
theorem.  This companion checks the complete finite F_13 Fourier/Radon
algebra, the two-permutation marginal hostile, and the explicit target
relations over an exact sample scalar word.  No third-party package is used.
"""


P = 13
MOD = 53
ZETA = 16  # primitive 13th root in F_53


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


require(pow(ZETA, P, MOD) == 1, "zeta^13")
require(all(pow(ZETA, j, MOD) != 1 for j in range(1, P)), "zeta primitive")


def zpow(exponent):
    return pow(ZETA, exponent % P, MOD)


def dft1(vector):
    return [
        sum(vector[t] * zpow(-q * t) for t in range(P)) * pow(P, -1, MOD) % MOD
        for q in range(P)
    ]


def dft2(table):
    inv = pow(P * P, -1, MOD)
    return [
        [
            sum(
                table[r][s] * zpow(-r * lam + s * nu)
                for r in range(P)
                for s in range(P)
            )
            * inv
            % MOD
            for nu in range(P)
        ]
        for lam in range(P)
    ]


# Delta-basis verification suffices for every joint series by linearity.
radon_basis_checks = 0
dephase_checks = 0
for left_residue in range(P):
    for right_residue in range(P):
        table = [
            [zpow(r * left_residue - s * right_residue) for s in range(P)]
            for r in range(P)
        ]
        transform = dft2(table)
        for lam in range(P):
            for nu in range(P):
                expected = int((lam, nu) == (left_residue, right_residue))
                require(transform[lam][nu] == expected, "2D residue selector")

        diagonal = [table[t][t] for t in range(P)]
        diagonal_transform = dft1(diagonal)
        for q in range(P):
            line_sum = sum(
                transform[lam][nu]
                for lam in range(P)
                for nu in range(P)
                if (lam - nu) % P == q
            ) % MOD
            require(diagonal_transform[q] == line_sum, "Radon line identity")
            radon_basis_checks += 1

        for deep_charge in range(P):
            dephased = [zpow(deep_charge * t) * diagonal[t] % MOD for t in range(P)]
            dephased_transform = dft1(dephased)
            for q in range(P):
                require(
                    dephased_transform[q]
                    == diagonal_transform[(q - deep_charge) % P],
                    "deep-charge dephasing",
                )
                dephase_checks += 1


# Two nonnegative permutation matrices with the same marginals and base axes.
def sigma(r):
    if r == 0:
        return 1
    if r == 1:
        return 12
    if r == 12:
        return 0
    return r


shift_matrix = [[int(s == (r + 1) % P) for s in range(P)] for r in range(P)]
sigma_matrix = [[int(s == sigma(r)) for s in range(P)] for r in range(P)]

for matrix in (shift_matrix, sigma_matrix):
    require([sum(row) for row in matrix] == [1] * P, "row marginals")
    require([sum(matrix[r][s] for r in range(P)) for s in range(P)] == [1] * P, "column marginals")

require(shift_matrix[0] == sigma_matrix[0], "base future row")
require(
    [shift_matrix[r][0] for r in range(P)] == [sigma_matrix[r][0] for r in range(P)],
    "base source column",
)

shift_diagonal = [shift_matrix[t][t] for t in range(P)]
sigma_diagonal = [sigma_matrix[t][t] for t in range(P)]
require(shift_diagonal == [0] * P, "shift diagonal")
require(sigma_diagonal == [int(2 <= t <= 11) for t in range(P)], "sigma diagonal")
sigma_colours = dft1(sigma_diagonal)
require(all(sigma_colours[q] != 0 for q in range(1, P)), "primitive diagonal colour")


# Exact THM-2309/2350 target-relation control on one typed scalar word.
# U={0,...,5}; B={j=6,a=7,b=8}; u0=0; ka=1; kb=2.
w = (1, 2, 3, 4, 5, 6, 13, 26, 39)
u0, ka, kb, a, b = 0, 1, 2, 7, 8


def unit_vector(index, coefficient=1):
    return tuple(coefficient if j == index else 0 for j in range(len(w)))


def add_vectors(*vectors):
    return tuple(sum(vector[j] for vector in vectors) for j in range(len(w)))


t_a = add_vectors(unit_vector(u0, w[a]), unit_vector(a, -w[u0]))
t_b = add_vectors(unit_vector(u0, w[b]), unit_vector(b, -w[u0]))


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def target_residue(vector):
    return ((vector[a] - vector[ka]) % P, (vector[b] - vector[kb]) % P)


require(dot(t_a, w) == 0 and dot(t_b, w) == 0, "exact target relations")
require(target_residue(t_a) == (12, 0), "a-axis residue")
require(target_residue(t_b) == (0, 12), "b-axis residue")

seen_residues = set()
for x in range(P):
    for y in range(P):
        relation = tuple(x * t_a[j] + y * t_b[j] for j in range(len(w)))
        require(dot(relation, w) == 0, "relation lattice combination")
        seen_residues.add(target_residue(relation))
require(len(seen_residues) == P * P, "target quotient not onto")


# The positive hostile has Laurent support {0,+relation,-relation}, all of
# whose physical exponents coincide because relation.w=0.
for relation in (t_a, t_b, add_vectors(t_a, t_b)):
    exponents = (0, dot(relation, w), -dot(relation, w))
    require(exponents == (0, 0, 0), "physical restriction did not collapse")
    residue = target_residue(relation)
    require(residue != (0, 0), "hostile residue accidentally zero")


print("THM-2615 exact finite referee")
print(f"field=F_{MOD}; zeta13={ZETA}")
print(f"radon_delta_basis_checks={radon_basis_checks}")
print(f"deep_charge_dephase_checks={dephase_checks}")
print("permutation_hostile: row_marginals=same column_marginals=same base_axes=same")
print("diagonals: shift=empty sigma={2,...,11}; primitive_colours=12/12")
print(f"sigma_primitive_values={','.join(str(sigma_colours[q]) for q in range(1, P))}")
print("target_relation_residues=169/169; physical_exponents=collapsed")
print("OK")
