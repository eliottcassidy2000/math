#!/usr/bin/env python3
"""Exact referee for THM-2383.

All complex arithmetic is over Q(i), represented by pairs of Fractions.
No assertion is used, so the optimized run checks the same identities.
"""

from fractions import Fraction as Q


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(0))
I = (Q(0), Q(1))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cneg(a):
    return (-a[0], -a[1])


def csub(a, b):
    return cadd(a, cneg(b))


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def cdiv(a, b):
    den = b[0] * b[0] + b[1] * b[1]
    require(den != 0, "division by zero in Q(i)")
    return ((a[0] * b[0] + a[1] * b[1]) / den,
            (a[1] * b[0] - a[0] * b[1]) / den)


def cscale(a, scalar):
    return (a[0] * scalar, a[1] * scalar)


def vzero(dim):
    return tuple(ZERO for _ in range(dim))


def vadd(a, b):
    return tuple(cadd(x, y) for x, y in zip(a, b))


def vsub(a, b):
    return tuple(csub(x, y) for x, y in zip(a, b))


def vscale(a, scalar):
    return tuple(cmul(scalar, x) for x in a)


def inner(a, b):
    """Hermitian inner product linear in the first slot."""
    total = ZERO
    for x, y in zip(a, b):
        total = cadd(total, cmul(x, cconj(y)))
    return total


def chi(mask, x):
    return -1 if (mask & x).bit_count() % 2 else 1


def submasks(mask):
    result = []
    h = mask
    while True:
        result.append(h)
        if h == 0:
            return result
        h = (h - 1) & mask


def inverse_walsh(coefficients, rank, dim):
    size = 1 << rank
    values = []
    for x in range(size):
        value = vzero(dim)
        for support, coefficient in enumerate(coefficients):
            value = vadd(value, vscale(coefficient, (Q(chi(support, x)), Q(0))))
        values.append(value)
    return values


def walsh(values, rank, dim):
    size = 1 << rank
    coefficients = []
    for support in range(size):
        value = vzero(dim)
        for x in range(size):
            value = vadd(value, vscale(values[x], (Q(chi(support, x), size), Q(0))))
        coefficients.append(value)
    return coefficients


def cross_dirichlet(left, right, subset):
    size = len(left)
    shifts = submasks(subset)
    total = ZERO
    for x in range(size):
        for h in shifts:
            dl = vsub(left[x ^ h], left[x])
            dr = vsub(right[x ^ h], right[x])
            total = cadd(total, inner(dl, dr))
    return cscale(total, Q(1, size * len(shifts)))


def dirichlet(values, subset):
    answer = cross_dirichlet(values, values, subset)
    require(answer[1] == 0, "a squared Dirichlet energy was not real")
    return answer[0]


def maps_add(left, right):
    return [vadd(a, b) for a, b in zip(left, right)]


def maps_scale(values, scalar):
    return [vscale(value, scalar) for value in values]


def solve(matrix, vector):
    """Exact Gaussian elimination over Q(i)."""
    n = len(vector)
    augmented = [list(matrix[row]) + [vector[row]] for row in range(n)]
    for column in range(n):
        pivot = next((row for row in range(column, n)
                      if augmented[row][column] != ZERO), None)
        require(pivot is not None, "singular reference Gram matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [cdiv(value, pivot_value)
                             for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == ZERO:
                continue
            augmented[row] = [
                csub(augmented[row][j], cmul(factor, augmented[column][j]))
                for j in range(n + 1)
            ]
    return [augmented[row][-1] for row in range(n)]


RANK = 4
DIM = 3
SIZE = 1 << RANK
FULL = SIZE - 1


def gaussian_rational(seed, coordinate):
    real = Q(((seed + 2) * (coordinate + 3)) % 17 - 8, seed + coordinate + 5)
    imag = Q(((seed + 5) * (coordinate + 1)) % 19 - 9, seed + 2 * coordinate + 7)
    return (real, imag)


fhat = [
    tuple(gaussian_rational(support + 1, coordinate)
          for coordinate in range(DIM))
    for support in range(SIZE)
]

base_columns = (
    (ONE, ZERO, ONE),
    (ONE, ONE, ZERO),
    (ZERO, ONE, I),
)

rhat_bank = []
for reference in range(DIM):
    coefficients = []
    for support in range(SIZE):
        scalar = (Q((support + 1) * (reference + 2)), Q(reference, support + 2))
        coefficients.append(tuple(cmul(scalar, entry)
                                  for entry in base_columns[reference]))
    rhat_bank.append(coefficients)

f_values = inverse_walsh(fhat, RANK, DIM)
r_values_bank = [
    inverse_walsh(coefficients, RANK, DIM)
    for coefficients in rhat_bank
]

require(walsh(f_values, RANK, DIM) == fhat, "Walsh inversion failed for F")
for reference in range(DIM):
    require(walsh(r_values_bank[reference], RANK, DIM) == rhat_bank[reference],
            "Walsh inversion failed for a reference")

direct_checks = 0
polarization_checks = 0
mobius_cells = 0
reconstruction_cells = 0

cross_banks = []
for reference in range(DIM):
    bank = []
    for subset in range(SIZE):
        direct = cross_dirichlet(f_values, r_values_bank[reference], subset)
        spectral = ZERO
        for support in range(1, SIZE):
            if support & subset:
                spectral = cadd(
                    spectral,
                    cscale(inner(fhat[support], rhat_bank[reference][support]), Q(2)),
                )
        require(direct == spectral, "cross-Dirichlet 0/2 multiplier failed")
        bank.append(direct)
        direct_checks += 1

        d_f = dirichlet(f_values, subset)
        d_r = dirichlet(r_values_bank[reference], subset)
        d_plus = dirichlet(maps_add(f_values, r_values_bank[reference]), subset)
        d_iquad = dirichlet(
            maps_add(f_values, maps_scale(r_values_bank[reference], I)), subset
        )
        recovered = (Q(d_plus - d_f - d_r, 2),
                     Q(d_iquad - d_f - d_r, 2))
        require(recovered == direct, "two-quadrature polarization failed")
        polarization_checks += 2
    cross_banks.append(bank)

for reference in range(DIM):
    bank = cross_banks[reference]
    qsum = []
    for subset in range(SIZE):
        qsum.append(cscale(csub(bank[FULL], bank[FULL ^ subset]), Q(1, 2)))
    for support in range(1, SIZE):
        recovered = ZERO
        for subset in submasks(support):
            sign = -1 if (support.bit_count() - subset.bit_count()) % 2 else 1
            recovered = cadd(recovered, cscale(qsum[subset], Q(sign)))
        expected = inner(fhat[support], rhat_bank[reference][support])
        require(recovered == expected, "Möbius cross-Gram inversion failed")
        mobius_cells += 1

for support in range(1, SIZE):
    references = [rhat_bank[j][support] for j in range(DIM)]
    gram = [
        [inner(references[column], references[row]) for column in range(DIM)]
        for row in range(DIM)
    ]
    measured = [inner(fhat[support], references[row]) for row in range(DIM)]
    coefficients = solve(gram, measured)
    recovered = vzero(DIM)
    for coefficient, reference_vector in zip(coefficients, references):
        recovered = vadd(recovered, vscale(reference_vector, coefficient))
    require(recovered == fhat[support], "spanning-reference reconstruction failed")
    reconstruction_cells += 1

# Proper-span hostile: +/-e_2 have the same self-bank and zero e_1 cross-bank.
rank_h = 2
size_h = 1 << rank_h
support_h = 1
e1 = (ONE, ZERO)
e2 = (ZERO, ONE)
zero2 = (ZERO, ZERO)
plus_hat = [zero2 for _ in range(size_h)]
minus_hat = [zero2 for _ in range(size_h)]
ref_hat = [zero2 for _ in range(size_h)]
plus_hat[support_h] = e2
minus_hat[support_h] = vscale(e2, (Q(-1), Q(0)))
ref_hat[support_h] = e1
plus_values = inverse_walsh(plus_hat, rank_h, 2)
minus_values = inverse_walsh(minus_hat, rank_h, 2)
ref_values = inverse_walsh(ref_hat, rank_h, 2)
for subset in range(size_h):
    require(dirichlet(plus_values, subset) == dirichlet(minus_values, subset),
            "proper-span hostile changed the self-bank")
    require(cross_dirichlet(plus_values, ref_values, subset) == ZERO,
            "proper-span hostile was not reference-orthogonal")
    require(cross_dirichlet(minus_values, ref_values, subset) == ZERO,
            "proper-span hostile was not reference-orthogonal")

# One real quadrature cannot distinguish +i from -i against the same reference.
plus_i_hat = [tuple([ZERO]) for _ in range(size_h)]
minus_i_hat = [tuple([ZERO]) for _ in range(size_h)]
scalar_ref_hat = [tuple([ZERO]) for _ in range(size_h)]
plus_i_hat[support_h] = (I,)
minus_i_hat[support_h] = (cneg(I),)
scalar_ref_hat[support_h] = (ONE,)
plus_i_values = inverse_walsh(plus_i_hat, rank_h, 1)
minus_i_values = inverse_walsh(minus_i_hat, rank_h, 1)
scalar_ref_values = inverse_walsh(scalar_ref_hat, rank_h, 1)
for subset in range(size_h):
    require(dirichlet(plus_i_values, subset) == dirichlet(minus_i_values, subset),
            "one-quadrature hostile changed the self-bank")
    require(
        dirichlet(maps_add(plus_i_values, scalar_ref_values), subset)
        == dirichlet(maps_add(minus_i_values, scalar_ref_values), subset),
        "one real union distinguished opposite imaginary phases",
    )
require(
    dirichlet(maps_add(plus_i_values, maps_scale(scalar_ref_values, I)), support_h)
    != dirichlet(maps_add(minus_i_values, maps_scale(scalar_ref_values, I)), support_h),
    "the second quadrature failed to distinguish phase",
)

# Constants are invisible to every difference bank.
constant = ((Q(7, 5), Q(-3, 4)),) * DIM
shifted_f = [vadd(value, constant) for value in f_values]
for subset in range(SIZE):
    require(dirichlet(shifted_f, subset) == dirichlet(f_values, subset),
            "constant shift changed a self-bank")
    for reference in range(DIM):
        require(
            cross_dirichlet(shifted_f, r_values_bank[reference], subset)
            == cross_dirichlet(f_values, r_values_bank[reference], subset),
            "constant shift changed a cross-bank",
        )

# Forgetting coordinate labels identifies two different singleton supports.
unit_hat_0 = [tuple([ZERO]) for _ in range(SIZE)]
unit_hat_1 = [tuple([ZERO]) for _ in range(SIZE)]
unit_hat_0[1] = (ONE,)
unit_hat_1[2] = (ONE,)
unit_values_0 = inverse_walsh(unit_hat_0, RANK, 1)
unit_values_1 = inverse_walsh(unit_hat_1, RANK, 1)
for cardinality in range(RANK + 1):
    profile_0 = sorted(
        dirichlet(unit_values_0, subset)
        for subset in range(SIZE)
        if subset.bit_count() == cardinality
    )
    profile_1 = sorted(
        dirichlet(unit_values_1, subset)
        for subset in range(SIZE)
        if subset.bit_count() == cardinality
    )
    require(profile_0 == profile_1, "unlabelled support hostile failed")
require(dirichlet(unit_values_0, 1) != dirichlet(unit_values_1, 1),
        "labelled bank did not separate coordinates")

# THM-2370 clone boundary: global complement preserves every squared bank.
clone_rank = 5
clone_size = 1 << clone_rank
clone_f = []
clone_complement = []
for x in range(clone_size):
    weight = x.bit_count()
    clone_f.append(((Q(clone_rank - weight, clone_rank), Q(0)),))
    clone_complement.append(((Q(weight, clone_rank), Q(0)),))
clone_controls = 0
for subset in range(clone_size):
    require(dirichlet(clone_f, subset) == dirichlet(clone_complement, subset),
            "clone complementation changed the squared bank")
    clone_controls += 1
require(clone_f[-1] == (ZERO,) and clone_complement[-1] == (ONE,),
        "clone terminal orientation hostile failed")

print("THM-2383 polarized complete-subcube Gram tomography exact referee")
print(f"boolean_rank: {RANK}")
print(f"hilbert_dimension: {DIM}")
print(f"direct_cross_dirichlet_checks: {direct_checks}")
print(f"polarization_quadrature_checks: {polarization_checks}")
print(f"mobius_cross_gram_cells: {mobius_cells}")
print(f"spanning_reconstruction_cells: {reconstruction_cells}")
print("proper_span_hostile: PASS")
print("one_quadrature_hostile: PASS")
print("constant_mode_hostile: PASS")
print("unlabelled_support_hostile: PASS")
print(f"clone_terminal_orientation_controls: {clone_controls}")
print("VERDICT: labelled polarized subcube banks recover exactly the referenced nonconstant Walsh coefficients")
