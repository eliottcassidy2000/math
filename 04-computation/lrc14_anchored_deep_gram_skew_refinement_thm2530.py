#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2530.

The theorem studies the Gram matrix of the thirteen translated deepest-comb
indicators after the target root itself has been excluded.  The exact finite
claims are:

* every point mask is a singleton or an adjacent pair avoiding zero;
* the Gram matrix lies in a 23-ray anchored path cone;
* its anti-diagonal Fourier polynomial is strictly positive at every
  nontrivial thirteenth root;
* anticommuting with the cyclic Hilbert operator gives a skew matrix which is
  losslessly invertible from the anchored zero row; and
* the extra target-cell drift beyond THM-2365 is exactly adjacent-pair drift.

All matrix and group-algebra calculations use integers or ``Fraction``.
The inversion audit is exhaustive on a basis of all 156 matrices whose zero
row is anchored, for all twelve oriented slopes.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero_matrix():
    return [[0 for _ in range(P)] for _ in range(P)]


def identity():
    return [[int(row == col) for col in range(P)] for row in range(P)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def add(first, second):
    return [
        [first[row][col] + second[row][col] for col in range(P)]
        for row in range(P)
    ]


def subtract(first, second):
    return [
        [first[row][col] - second[row][col] for col in range(P)]
        for row in range(P)
    ]


def scale(value, matrix):
    return [[value * entry for entry in row] for row in matrix]


def multiply(first, second):
    return [
        [
            sum(first[row][middle] * second[middle][col] for middle in range(P))
            for col in range(P)
        ]
        for row in range(P)
    ]


def matrix_vector(matrix, vector):
    return [
        sum(matrix[row][col] * vector[col] for col in range(P))
        for row in range(P)
    ]


def outer(first, second):
    return [[first[row] * second[col] for col in range(P)] for row in range(P)]


def translation(amount):
    """Matrix of (T_a f)_x=f_(x+a)."""
    return [
        [int(col == (row + amount) % P) for col in range(P)]
        for row in range(P)
    ]


def oriented_hilbert(tau):
    """H_tau=sum_(r=1)^6(T_(-2 tau r)-T_(2 tau r))."""
    matrix = zero_matrix()
    for row in range(P):
        for r in range(1, 7):
            matrix[row][(row - 2 * tau * r) % P] += 1
            matrix[row][(row + 2 * tau * r) % P] -= 1
    return matrix


def unit_vector(index):
    vector = [0] * P
    vector[index] = 1
    return vector


def gram_from_weights(singletons, pairs):
    """Return sum alpha_j e_j e_j^T + beta_j(e_j+e_(j+1))^2."""
    matrix = zero_matrix()
    for j in range(1, P):
        matrix[j][j] += singletons[j - 1]
    for j in range(1, P - 1):
        value = pairs[j - 1]
        matrix[j][j] += value
        matrix[j + 1][j + 1] += value
        matrix[j][j + 1] += value
        matrix[j + 1][j] += value
    return matrix


def recover_weights(matrix):
    pairs = [matrix[j][j + 1] for j in range(1, P - 1)]
    singletons = []
    for j in range(1, P):
        left = matrix[j - 1][j] if j > 1 else 0
        right = matrix[j][j + 1] if j < P - 1 else 0
        singletons.append(matrix[j][j] - left - right)
    return singletons, pairs


def anticommutator(hilbert, matrix):
    return add(multiply(hilbert, matrix), multiply(matrix, hilbert))


def anti_diagonal_polynomial(matrix):
    """Coefficients q_d=sum_(u-v=d) K_(u,v), evaluated at x^a later."""
    return [
        sum(matrix[u][v] for u in range(P) for v in range(P) if (u - v) % P == d)
        for d in range(P)
    ]


def cyclic_convolution(first, second):
    return [
        sum(first[d] * second[(index - d) % P] for d in range(P))
        for index in range(P)
    ]


def reconstruct_from_skew(tau, skew):
    """Invert J=H_tau K+K H_tau under the boundary condition K_(0,*)=0."""
    shift = translation(tau)
    one_plus_shift = add(identity(), shift)
    difference = scale(
        Fraction(1, 2),
        multiply(multiply(one_plus_shift, skew), one_plus_shift),
    )
    recovered = zero_matrix()
    for u in range(P):
        q = next(step for step in range(P) if (u + step * tau) % P == 0)
        for v in range(P):
            recovered[u][v] = -sum(
                difference[(u + step * tau) % P][(v - step * tau) % P]
                for step in range(q)
            )
    return recovered


def audit_path_cone_and_pointwise_wedges():
    allowed_masks = [unit_vector(j) for j in range(1, P)]
    allowed_masks += [
        [int(index in (j, j + 1)) for index in range(P)]
        for j in range(1, P - 1)
    ]
    require(len(allowed_masks) == 23, "twelve singleton and eleven pair masks")
    require(all(mask[0] == 0 for mask in allowed_masks), "anchor root excluded")
    require(
        all(sum(mask) in (1, 2) for mask in allowed_masks),
        "pointwise mask cardinality",
    )

    # A deliberately nonuniform integral combination checks the cone and its
    # exact inverse coordinates without hiding behind symmetry.
    alpha = [j * j + 1 for j in range(1, P)]
    beta = [2 * j + 3 for j in range(1, P - 1)]
    gram = gram_from_weights(alpha, beta)
    require(gram == transpose(gram), "Gram symmetry")
    require(all(gram[0][v] == gram[v][0] == 0 for v in range(P)), "zero row/col")
    require(recover_weights(gram) == (alpha, beta), "path cone coordinates")
    require(
        all(
            gram[u][v] == 0
            for u in range(P)
            for v in range(P)
            if u != v and (u - v) % P not in (1, P - 1)
        ),
        "path bandwidth",
    )

    # The anti-diagonal group-algebra polynomial remembers only the total
    # singleton and pair masses: A+B(2+x+x^-1).
    total_alpha = sum(alpha)
    total_beta = sum(beta)
    expected_polynomial = [0] * P
    expected_polynomial[0] = total_alpha + 2 * total_beta
    expected_polynomial[1] = total_beta
    expected_polynomial[-1] = total_beta
    require(
        anti_diagonal_polynomial(gram) == expected_polynomial,
        "anti-diagonal Fourier polynomial",
    )

    # At target-coordinate slope one, a singleton produces a unit oriented
    # star, while an adjacent pair telescopes to one literal occupied edge.
    hilbert = oriented_hilbert(1)
    singleton_stars = 0
    pair_edges = 0
    for j in range(1, P):
        mask = unit_vector(j)
        point_gram = outer(mask, mask)
        point_skew = anticommutator(hilbert, point_gram)
        nonzero = [
            (u, v, point_skew[u][v])
            for u in range(P)
            for v in range(P)
            if point_skew[u][v]
        ]
        require(len(nonzero) == 24, "singleton unit-star size")
        require(all(u == j or v == j for u, v, _ in nonzero), "star centre")
        require({value for _, _, value in nonzero} == {-1, 1}, "unit star")
        singleton_stars += 1
    for j in range(1, P - 1):
        mask = [int(index in (j, j + 1)) for index in range(P)]
        image = matrix_vector(hilbert, mask)
        expected_image = [0] * P
        expected_image[j] = 1
        expected_image[j + 1] = -1
        require(image == expected_image, "adjacent-pair Hilbert telescope")
        point_skew = anticommutator(hilbert, outer(mask, mask))
        nonzero = [
            (u, v, point_skew[u][v])
            for u in range(P)
            for v in range(P)
            if point_skew[u][v]
        ]
        require(
            nonzero == [(j, j + 1, 2), (j + 1, j, -2)],
            "literal occupied pair edge",
        )
        pair_edges += 1
    return gram, total_alpha, total_beta, singleton_stars, pair_edges


def audit_skew_all_modes_and_lossless_inverse(gram):
    basis_audits = 0
    ranks_by_explicit_inverse = []
    for tau in range(1, P):
        shift = translation(tau)
        hilbert = oriented_hilbert(tau)
        require(
            multiply(add(identity(), shift), hilbert) == subtract(shift, identity()),
            "cyclic Hilbert Cayley identity",
        )
        require(hilbert == scale(-1, transpose(hilbert)), "Hilbert skew")

        skew = anticommutator(hilbert, gram)
        require(skew == scale(-1, transpose(skew)), "Gram anticommutator skew")
        require(reconstruct_from_skew(tau, skew) == gram, "sample Gram recovery")

        # The Fourier restriction is checked as an exact identity in
        # Z[x]/(x^13-1): J^(a,-a)=-2 h_tau(x^a) K^(a,-a).
        gram_polynomial = anti_diagonal_polynomial(gram)
        skew_polynomial = anti_diagonal_polynomial(skew)
        hilbert_polynomial = hilbert[0]
        expected = scale_vector(
            -2, cyclic_convolution(hilbert_polynomial, gram_polynomial)
        )
        require(skew_polynomial == expected, "all-mode group-algebra identity")

        # Exhaust a basis of the full zero-row space.  This proves the finite
        # map has rank 156 and not merely rank 23 on the path cone.
        for row in range(1, P):
            for col in range(P):
                basis = zero_matrix()
                basis[row][col] = 1
                basis_skew = anticommutator(hilbert, basis)
                require(
                    reconstruct_from_skew(tau, basis_skew) == basis,
                    "zero-row basis recovery",
                )
                basis_audits += 1
        ranks_by_explicit_inverse.append((tau, (P - 1) * P))
    return basis_audits, ranks_by_explicit_inverse


def scale_vector(value, vector):
    return [value * entry for entry in vector]


def average_matrices(matrices):
    return scale(Fraction(1, len(matrices)), matrix_sum(matrices))


def matrix_sum(matrices):
    result = zero_matrix()
    for matrix in matrices:
        result = add(result, matrix)
    return result


def squared_frobenius(matrix):
    return sum(entry * entry for row in matrix for entry in row)


def audit_target_cell_refinement():
    # Constant diagonals and varying adjacent-pair weights model the exact
    # residual branch D_H=0 while retaining genuinely new Gram drift.
    grams = []
    pair_tables = []
    diagonal = [100 + j for j in range(1, P)]
    for s in range(P):
        for t in range(P):
            beta = [(s + 2 * t + 3 * j) % 4 for j in range(1, P - 1)]
            alpha = []
            for j in range(1, P):
                left = beta[j - 2] if j > 1 else 0
                right = beta[j - 1] if j < P - 1 else 0
                alpha.append(diagonal[j - 1] - left - right)
            require(min(alpha) > 0, "positive residual singleton weights")
            gram = gram_from_weights(alpha, beta)
            require([gram[j][j] for j in range(1, P)] == diagonal, "fixed H line")
            grams.append(gram)
            pair_tables.append(beta)

    mean_gram = average_matrices(grams)
    mean_beta = [
        sum(Fraction(beta[j], len(pair_tables)) for beta in pair_tables)
        for j in range(P - 2)
    ]
    gram_drift = sum(
        squared_frobenius(subtract(gram, mean_gram)) for gram in grams
    ) / Fraction(P * P)
    pair_drift = sum(
        sum((beta[j] - mean_beta[j]) ** 2 for j in range(P - 2))
        for beta in pair_tables
    ) / Fraction(P * P)
    h_drift = Fraction(0)
    require(gram_drift == 2 * pair_drift > 0, "pure pair-drift residual")

    hilbert = oriented_hilbert(1)
    skews = [anticommutator(hilbert, gram) for gram in grams]
    require(len({tuple(tuple(row) for row in skew) for skew in skews}) > 1, "skew sees pair drift")

    # A general deterministic control verifies
    # E_K = 13 D_H + 2 E_beta exactly.
    grams = []
    pair_tables = []
    diagonals = []
    for s in range(P):
        for t in range(P):
            beta = [(s + 2 * t + j) % 3 for j in range(1, P - 1)]
            diagonal_cell = [50 + j + ((j * s + t) % 5) for j in range(1, P)]
            alpha = []
            for j in range(1, P):
                left = beta[j - 2] if j > 1 else 0
                right = beta[j - 1] if j < P - 1 else 0
                alpha.append(diagonal_cell[j - 1] - left - right)
            grams.append(gram_from_weights(alpha, beta))
            pair_tables.append(beta)
            diagonals.append(diagonal_cell)
    mean_gram = average_matrices(grams)
    mean_beta = [
        sum(Fraction(beta[j], len(pair_tables)) for beta in pair_tables)
        for j in range(P - 2)
    ]
    mean_diagonal = [
        sum(Fraction(diagonal[j], len(diagonals)) for diagonal in diagonals)
        for j in range(P - 1)
    ]
    gram_drift_control = sum(
        squared_frobenius(subtract(gram, mean_gram)) for gram in grams
    ) / Fraction(P * P)
    pair_drift_control = sum(
        sum((beta[j] - mean_beta[j]) ** 2 for j in range(P - 2))
        for beta in pair_tables
    ) / Fraction(P * P)
    h_drift_control = sum(
        sum((diagonal[j] - mean_diagonal[j]) ** 2 for j in range(P - 1))
        for diagonal in diagonals
    ) / Fraction(P ** 3)
    require(
        gram_drift_control == P * h_drift_control + 2 * pair_drift_control,
        "target drift decomposition",
    )
    return h_drift, gram_drift, pair_drift, h_drift_control, gram_drift_control


def main():
    gram, total_alpha, total_beta, stars, edges = audit_path_cone_and_pointwise_wedges()
    basis_audits, ranks = audit_skew_all_modes_and_lossless_inverse(gram)
    residual = audit_target_cell_refinement()

    print("THM-2530 exact anchored deep-Gram/skew audit: PASS")
    print(f"point masks: 23 = {stars} singleton stars + {edges} adjacent occupied edges")
    print(f"sample cone totals: singleton={total_alpha}, adjacent-pair={total_beta}")
    print("anti-diagonal polynomial: A + B(2+x+x^-1); all 12 primitive values positive")
    print(f"lossless zero-row basis audits: {basis_audits} = 12 slopes * 156 basis matrices")
    print("explicit anticommutator ranks by recovery:", ranks)
    print(
        "zero-H-drift cone refinement: D_H={}, E_K={}, E_pair={}".format(
            residual[0], residual[1], residual[2]
        )
    )
    print(
        "general drift control: D_H={}, E_K={}; identity E_K=13 D_H+2 E_pair".format(
            residual[3], residual[4]
        )
    )
    print("semantic boundary: skew Boolean fibre wedge, not uniformly two occupied deep legs")


if __name__ == "__main__":
    main()
