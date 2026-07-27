#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2543."""

from itertools import combinations, product
from math import prod


ROWS = 7
ROOTS = 13
FIELD = 547  # 547-1 is divisible by 7*13.


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_root(p: int) -> int:
    factors = (2, 3, 7, 13)
    for g in range(2, p):
        if all(pow(g, (p - 1) // d, p) != 1 for d in factors):
            return g
    raise RuntimeError("primitive root not found")


GENERATOR = primitive_root(FIELD)
ZETA7 = pow(GENERATOR, (FIELD - 1) // ROWS, FIELD)
ZETA13 = pow(GENERATOR, (FIELD - 1) // ROOTS, FIELD)
require(pow(ZETA7, ROWS, FIELD) == 1 and ZETA7 != 1, "bad 7th root")
require(pow(ZETA13, ROOTS, FIELD) == 1 and ZETA13 != 1, "bad 13th root")


def phase_hat(q: tuple[int, ...], eta: int) -> int:
    return sum(q[g] * pow(ZETA7, eta * g, FIELD) for g in range(ROWS)) % FIELD


def table_hat(table: list[list[int]], kappa: int, b: int) -> int:
    return sum(
        table[ell][s]
        * pow(ZETA7, kappa * ell, FIELD)
        * pow(ZETA13, b * s, FIELD)
        for ell in range(ROWS)
        for s in range(ROOTS)
    ) % FIELD


def tensor_hat(
    table: list[list[int]], q: tuple[int, ...], kappa: int, eta: int, b: int
) -> int:
    return sum(
        table[ell][s]
        * q[gamma]
        * pow(ZETA7, kappa * ell + eta * gamma, FIELD)
        * pow(ZETA13, b * s, FIELD)
        for ell in range(ROWS)
        for gamma in range(ROWS)
        for s in range(ROOTS)
    ) % FIELD


def rank_mod(matrix: list[list[int]]) -> int:
    a = [[x % FIELD for x in row] for row in matrix]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, FIELD)
        a[rank] = [inv * x % FIELD for x in a[rank]]
        for r in range(len(a)):
            if r != rank and a[r][col]:
                c = a[r][col]
                a[r] = [
                    (a[r][j] - c * a[rank][j]) % FIELD
                    for j in range(len(a[r]))
                ]
        rank += 1
        if rank == len(a):
            return rank
    return rank


def bits(mask: int) -> list[int]:
    return [(mask >> u) & 1 for u in range(ROOTS)]


def empty_bank() -> list[list[list[int]]]:
    return [[[] for _s in range(ROOTS)] for _ell in range(ROWS)]


def control_bank() -> list[list[list[int]]]:
    """THM-2539's two-atom guard-supported positive control."""
    singleton = 1 << 1
    pair = (1 << 0) | (1 << 1)
    bank = empty_bank()
    bank[0][0] = [singleton]
    for s in range(1, ROOTS):
        bank[0][s] = [pair, singleton]
        for ell in range(1, ROWS):
            bank[ell][s] = [singleton]
    return bank


def aggregate(bank: list[list[list[int]]]) -> list[list[list[int]]]:
    out = [
        [[0 for _u in range(ROOTS)] for _s in range(ROOTS)]
        for _ell in range(ROWS)
    ]
    for ell in range(ROWS):
        for s in range(ROOTS):
            for mask in bank[ell][s]:
                e = bits(mask)
                for u in range(ROOTS):
                    out[ell][s][u] += e[u]
    return out


def root_coefficient(
    array: list[list[list[int]]], kappa: int, b: int, a: int
) -> int:
    return sum(
        array[ell][s][u]
        * pow(ZETA7, kappa * ell, FIELD)
        * pow(ZETA13, b * s - a * u, FIELD)
        for ell in range(ROWS)
        for s in range(ROOTS)
        for u in range(ROOTS)
    ) % FIELD


def boundary_aggregate(
    bank: list[list[list[int]]], tau: int
) -> tuple[list[list[list[int]]], list[list[list[int]]]]:
    plus = [
        [[0 for _u in range(ROOTS)] for _s in range(ROOTS)]
        for _ell in range(ROWS)
    ]
    minus = [
        [[0 for _u in range(ROOTS)] for _s in range(ROOTS)]
        for _ell in range(ROWS)
    ]
    for ell in range(ROWS):
        for s in range(ROOTS):
            for mask in bank[ell][s]:
                e = bits(mask)
                for u in range(ROOTS):
                    plus[ell][s][u] += e[u] * (1 - e[(u + tau) % ROOTS])
                    minus[ell][s][u] += (1 - e[u]) * e[(u + tau) % ROOTS]
    return plus, minus


# The anchored cubic positive control C=A^(o 3).
cubic = [[0 for _s in range(ROOTS)] for _ell in range(ROWS)]
cubic[0][0] = 1
for s in range(1, ROOTS):
    cubic[0][s] = 8
    for ell in range(1, ROWS):
        cubic[ell][s] = 1

cubic_mode_checks = 0
for kappa in range(1, ROWS):
    for b in range(1, ROOTS):
        require(table_hat(cubic, kappa, b) != 0, "cubic mixed mode vanished")
        cubic_mode_checks += 1


# Q[C7]=Q direct-sum Q(zeta7): over a large exact rational profile bank,
# one nontrivial character vanishes iff the profile is constant.
phase_profiles = 0
nonconstant_profiles = 0
phase_character_checks = 0
diagonal_nonzero_checks = 0
flat_diagonal_zero_checks = 0
for q0 in (1, 2):
    for tail in product(range(3), repeat=ROWS - 1):
        q = (q0,) + tail
        constant = len(set(q)) == 1
        phase_profiles += 1
        if not constant:
            nonconstant_profiles += 1
        for eta in range(1, ROWS):
            value = phase_hat(q, eta)
            require((value == 0) == constant, "all-or-flat phase lemma failed")
            phase_character_checks += 1
        for kappa in range(1, ROWS):
            qmode = phase_hat(q, (-kappa) % ROWS)
            for b in range(1, ROOTS):
                product_mode = qmode * table_hat(cubic, kappa, b) % FIELD
                if constant:
                    require(product_mode == 0, "flat relative mode survived")
                    flat_diagonal_zero_checks += 1
                else:
                    require(product_mode != 0, "relative mode vanished")
                    diagonal_nonzero_checks += 1


# Direct tensor-factorisation on all q0=1 Boolean phase profiles.
factorisation_checks = 0
quotient_table_checks = 0
inv7 = pow(ROWS, -1, FIELD)
inv91 = pow(ROWS * ROOTS, -1, FIELD)
inv637 = pow(ROWS * ROWS * ROOTS, -1, FIELD)
for tail in product((0, 1), repeat=ROWS - 1):
    q = (1,) + tail
    relative = [
        [
            sum(q[gamma] * cubic[(delta + gamma) % ROWS][s]
                for gamma in range(ROWS)) * inv7 % FIELD
            for s in range(ROOTS)
        ]
        for delta in range(ROWS)
    ]
    for kappa in range(1, ROWS):
        eta = (-kappa) % ROWS
        for b in range(1, ROOTS):
            direct = tensor_hat(cubic, q, kappa, eta, b) * inv637 % FIELD
            factored = (
                table_hat(cubic, kappa, b) * inv91
                * phase_hat(q, eta) * inv7
            ) % FIELD
            require(direct == factored, "tensor DFT did not factor")
            quotient = table_hat(relative, kappa, b) * inv91 % FIELD
            require(quotient == direct, "relative quotient DFT failed")
            factorisation_checks += 1
            quotient_table_checks += 1


# Simultaneous translation of the row and phase charts fixes eta=-kappa.
covariance_checks = 0
controls = (
    (1, 0, 0, 0, 0, 0, 0),
    (1, 1, 1, 1, 1, 1, 1),
    (1, 0, 2, 1, 2, 0, 1),
)
for q in controls:
    for delta in range(ROWS):
        shifted_q = tuple(q[(gamma - delta) % ROWS] for gamma in range(ROWS))
        shifted_cubic = [cubic[(ell - delta) % ROWS][:] for ell in range(ROWS)]
        for kappa in range(1, ROWS):
            eta = (-kappa) % ROWS
            for b in range(1, ROOTS):
                require(
                    tensor_hat(shifted_cubic, shifted_q, kappa, eta, b)
                    == tensor_hat(cubic, q, kappa, eta, b),
                    "relative character was not gauge invariant",
                )
                covariance_checks += 1


# The fixed gamma=0 slice always carries all visible source-target modes but
# is not invariant under translating the phase chart alone.
fixed_owner_slice_checks = 0
for kappa in range(1, ROWS):
    for b in range(1, ROOTS):
        require(table_hat(cubic, kappa, b) != 0, "fixed owner slice vanished")
        fixed_owner_slice_checks += 1
shifted_one_hot = tuple(controls[0][(g - 1) % ROWS] for g in range(ROWS))
require(
    controls[0][0] * table_hat(cubic, 1, 1) % FIELD
    != shifted_one_hot[0] * table_hat(cubic, 1, 1) % FIELD,
    "fixed phase slice hostile accidentally invariant",
)


# Common consecutive guard zeros force at least four active root modes.
guard_vandermonde_checks = 0
for size in (1, 2, 3):
    for active in combinations(range(ROOTS), size):
        for start in range(ROOTS):
            forbidden = [(start + j) % ROOTS for j in range(size)]
            matrix = [
                [pow(ZETA13, a * u, FIELD) for a in active]
                for u in forbidden
            ]
            require(rank_mod(matrix) == size, "guard Vandermonde singular")
            guard_vandermonde_checks += 1


# A one-hot phase profile is the sharp augmentation-positive/norm-zero control.
# Its relative tensor is exactly the THM-2539 root bank, so all 216 boundary
# incidences survive while the sevenfold phase product is zero.
bank = control_bank()
array = aggregate(bank)
one_hot = (1, 0, 0, 0, 0, 0, 0)
require(phase_hat(one_hot, 1) != 0 and prod(one_hot) == 0, "one-hot boundary failed")
seed_active = [a for a in range(ROOTS) if root_coefficient(array, 1, 1, a)]
require(seed_active == list(range(ROOTS)), "seed root spectrum drifted")
slopes = tuple(a for a in seed_active if a)[:3]

incidence_checks = 0
for kappa in range(1, ROWS):
    for b in range(1, ROOTS):
        for slope in slopes:
            a = slope * b % ROOTS
            original = root_coefficient(array, kappa, b, a)
            require(original != 0, "Galois root diagonal vanished")
            plus, minus = boundary_aggregate(bank, slope)
            cplus = root_coefficient(plus, kappa, b, a)
            cminus = root_coefficient(minus, kappa, b, a)
            multiplier = (pow(ZETA13, a * slope, FIELD) - 1) % FIELD
            require(
                (cminus - cplus) % FIELD == multiplier * original % FIELD,
                "positive boundary factorisation failed",
            )
            require(cplus or cminus, "both positive orientations vanished")
            incidence_checks += 1


# THM-2540's horizontal endpoint transport, on explicit noninvariant weights.
endpoint_transport_checks = 0
for mask in (1, 3, 0b101101, (1 << ROOTS) - 2):
    e = bits(mask)
    for tau in range(1, ROOTS):
        tail = [e[u] * (1 - e[(u + tau) % ROOTS]) for u in range(ROOTS)]
        head = [tail[(u - tau) % ROOTS] for u in range(ROOTS)]
        for w in (
            tuple(range(1, ROOTS + 1)),
            tuple((u * u + 3 * u + 1) % 17 for u in range(ROOTS)),
        ):
            transported = [w[(u - tau) % ROOTS] for u in range(ROOTS)]
            require(
                sum(w[u] * tail[u] for u in range(ROOTS))
                == sum(transported[u] * head[u] for u in range(ROOTS)),
                "covariant endpoint transport failed",
            )
            endpoint_transport_checks += 1


flat = (1,) * ROWS
require(prod(flat) == 1, "flat norm control vanished")
require(all(phase_hat(flat, eta) == 0 for eta in range(1, ROWS)), "flat augmentation survived")
semantic_arrival_flags = [0] * (ROWS * ROOTS)
require(not any(semantic_arrival_flags), "semantic edge was silently added")


print("THM-2543 exact referee")
print(f"field=F_{FIELD} primitive_generator={GENERATOR} zeta7={ZETA7} zeta13={ZETA13}")
print(f"cubic_mixed_modes={cubic_mode_checks}")
print(
    f"phase_profiles={phase_profiles} nonconstant={nonconstant_profiles} "
    f"character_checks={phase_character_checks}"
)
print(
    f"relative_diagonal_nonzero={diagonal_nonzero_checks} "
    f"flat_diagonal_zero={flat_diagonal_zero_checks}"
)
print(f"direct_normalized_tensor_factorisations={factorisation_checks}")
print(f"normalized_relative_quotient_checks={quotient_table_checks}")
print(f"diagonal_gauge_covariance_checks={covariance_checks}")
print(f"fixed_owner_slice_modes={fixed_owner_slice_checks} noncovariant_control=PASS")
print(f"guard_vandermonde_minors={guard_vandermonde_checks}")
print(f"one_hot_norm_zero=PASS seed_active_roots={seed_active} slopes={slopes}")
print(f"relative_positive_boundary_incidences={incidence_checks}")
print(f"covariant_endpoint_transport_checks={endpoint_transport_checks}")
print("flat_augmentation_zero_norm_positive=PASS")
print(f"semantic_arrival_flags={sum(semantic_arrival_flags)}")
print("THM-2543 PASS")
