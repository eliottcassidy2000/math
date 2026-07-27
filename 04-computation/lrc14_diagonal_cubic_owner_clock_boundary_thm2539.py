#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2539.

The cyclotomic checks are evaluated in F_547, whose multiplicative group
contains primitive seventh and thirteenth roots.  All other checks use exact
integers or Fraction.  The finite bank is a positive control for the theorem's
index algebra and, simultaneously, the sharp canonical-selector hostile.
"""

from fractions import Fraction
from itertools import combinations


ROWS = 7
Q = 13
PRIME = 547
ROOTS = range(Q)
MIXED_SOURCE = range(1, ROWS)
MIXED_TARGET = range(1, Q)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_root(p: int) -> int:
    factors = (2, 3, 7, 13)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in factors):
            return g
    raise RuntimeError("primitive root not found")


GENERATOR = primitive_root(PRIME)
ZETA7 = pow(GENERATOR, (PRIME - 1) // ROWS, PRIME)
ZETA13 = pow(GENERATOR, (PRIME - 1) // Q, PRIME)
require(pow(ZETA7, ROWS, PRIME) == 1 and ZETA7 != 1, "bad seventh root")
require(pow(ZETA13, Q, PRIME) == 1 and ZETA13 != 1, "bad thirteenth root")


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def bits(mask: int) -> list[int]:
    return [(mask >> u) & 1 for u in ROOTS]


def shift_vector(v: list[int], tau: int) -> list[int]:
    return [v[(u + tau) % Q] for u in ROOTS]


def cayley_vector(v: list[int], tau: int) -> list[int]:
    return [
        sum((1 if j % 2 else -1) * v[(u + j * tau) % Q] for j in range(1, Q))
        for u in ROOTS
    ]


def boundary(mask: int, tau: int) -> tuple[list[int], list[int]]:
    e = bits(mask)
    plus = [e[u] * (1 - e[(u + tau) % Q]) for u in ROOTS]
    minus = [(1 - e[u]) * e[(u + tau) % Q] for u in ROOTS]
    return plus, minus


def word(e: list[int], tau: int, anchor: int) -> tuple[int, ...]:
    return tuple(e[(anchor + j * tau) % Q] for j in ROOTS)


def selector(mask: int, tau: int) -> tuple[int, int]:
    e = bits(mask)
    anchor = max(ROOTS, key=lambda u: word(e, tau, u))
    run = next(j for j in range(1, Q) if not e[(anchor + j * tau) % Q])
    tail = (anchor + (run - 1) * tau) % Q
    head = (tail + tau) % Q
    require(e[tail] == 1 and e[head] == 0, "selector did not find a wall")
    return tail, head


def rank_mod(matrix: list[list[int]], p: int) -> int:
    if not matrix:
        return 0
    a = [[x % p for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, p)
        a[rank] = [(inv * x) % p for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                factor = a[r][col]
                a[r] = [
                    (a[r][c] - factor * a[rank][c]) % p for c in range(cols)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def empty_bank() -> list[list[list[int]]]:
    return [[[] for _s in range(Q)] for _ell in range(ROWS)]


def hostile_bank(singleton_root: int) -> list[list[list[int]]]:
    """Anchored two-atom bank; pair is {0,1}."""
    single = 1 << singleton_root
    pair = (1 << 0) | (1 << 1)
    bank = empty_bank()
    bank[0][0] = [single]
    for s in range(1, Q):
        bank[0][s] = [pair, single]
        for ell in range(1, ROWS):
            bank[ell][s] = [single]
    return bank


def aggregate(bank: list[list[list[int]]]) -> list[list[list[int]]]:
    answer = [[[0 for _u in ROOTS] for _s in range(Q)] for _ell in range(ROWS)]
    for ell in range(ROWS):
        for s in range(Q):
            for mask in bank[ell][s]:
                e = bits(mask)
                for u in ROOTS:
                    answer[ell][s][u] += e[u]
    return answer


def boundary_aggregate(
    bank: list[list[list[int]]], tau: int
) -> tuple[list[list[list[int]]], list[list[list[int]]]]:
    plus = [[[0 for _u in ROOTS] for _s in range(Q)] for _ell in range(ROWS)]
    minus = [[[0 for _u in ROOTS] for _s in range(Q)] for _ell in range(ROWS)]
    for ell in range(ROWS):
        for s in range(Q):
            for mask in bank[ell][s]:
                bp, bm = boundary(mask, tau)
                for u in ROOTS:
                    plus[ell][s][u] += bp[u]
                    minus[ell][s][u] += bm[u]
    return plus, minus


def selected_aggregate(
    bank: list[list[list[int]]], tau: int
) -> list[list[list[int]]]:
    selected = [[[0 for _u in ROOTS] for _s in range(Q)] for _ell in range(ROWS)]
    for ell in range(ROWS):
        for s in range(Q):
            for mask in bank[ell][s]:
                tail, _head = selector(mask, tau)
                selected[ell][s][tail] += 1
    return selected


def coefficient(
    array: list[list[list[int]]], kappa: int, b: int, a: int
) -> int:
    total = 0
    for ell in range(ROWS):
        x7 = pow(ZETA7, kappa * ell, PRIME)
        for s in range(Q):
            x713 = x7 * pow(ZETA13, b * s, PRIME) % PRIME
            for u in ROOTS:
                total += (
                    array[ell][s][u]
                    * x713
                    * pow(ZETA13, (-a * u) % Q, PRIME)
                )
    return total % PRIME


def mass_table(array: list[list[list[int]]]) -> list[list[int]]:
    return [[sum(array[ell][s]) for s in range(Q)] for ell in range(ROWS)]


def table_coefficient(table: list[list[int]], kappa: int, b: int) -> int:
    total = 0
    for ell in range(ROWS):
        for s in range(Q):
            total += (
                table[ell][s]
                * pow(ZETA7, kappa * ell, PRIME)
                * pow(ZETA13, b * s, PRIME)
            )
    return total % PRIME


# Exact root invariance of every late factor.
root_invariance_checks = 0
for length in (2, 4, 6):
    exponents = {length, 2 * length}
    exponents.update((3 + d) * length for d in range(ROWS))
    for denominator in (17, 19):
        for numerator in range(denominator):
            z = Fraction(numerator, denominator)
            for u in ROOTS:
                x = (z + u) / Q
                for exponent in exponents:
                    require(
                        frac((Q**exponent) * x)
                        == frac((Q ** (exponent - 1)) * z),
                        "late factor was not root invariant",
                    )
                    root_invariance_checks += 1


# Latin owner slot and the exact source/clock character on d=-ell.
assignments = [[(center + d) % ROWS for d in range(ROWS)] for center in range(ROWS)]
owner_slots = [row.index(0) for row in assignments]
require(owner_slots == [(-ell) % ROWS for ell in range(ROWS)], "owner slot drifted")
clock_character_checks = 0
for ell in range(ROWS):
    d = (-ell) % ROWS
    for kappa in MIXED_SOURCE:
        require(
            pow(ZETA7, kappa * ell, PRIME)
            == pow(ZETA7, (-kappa * d) % ROWS, PRIME),
            "source/clock character failed",
        )
        clock_character_checks += 1


# Every cell uses all phases.  One zero phase annihilates every product.
zero_phase_profiles = 0
for profile in range(1 << ROWS):
    q = [(profile >> gamma) & 1 for gamma in range(ROWS)]
    cell_products = []
    for center in range(ROWS):
        product = 1
        for d in range(ROWS):
            product *= q[(center + d) % ROWS]
        cell_products.append(product)
    expected = int(all(q))
    require(cell_products == [expected] * ROWS, "phase product depends on cell")
    if not expected:
        zero_phase_profiles += 1
require(zero_phase_profiles == 127, "zero-phase census drifted")


# Consecutive guard roots see an invertible Fourier minor whenever the active
# set has size at most three.
guard_vandermonde_checks = 0
for size in (1, 2, 3):
    for active in combinations(ROOTS, size):
        for start in ROOTS:
            forbidden = [(start + j) % Q for j in range(size)]
            matrix = [
                [pow(ZETA13, a * u, PRIME) for a in active] for u in forbidden
            ]
            require(rank_mod(matrix, PRIME) == size, "guard Vandermonde singular")
            guard_vandermonde_checks += 1
require(guard_vandermonde_checks == 4_901, "guard minor census drifted")


# The positive control and canonical-selector hostile share one exact bank.
# Every mask is supported on a chosen adjacent pair in a common guard-safe run.
bank = hostile_bank(singleton_root=1)
array = aggregate(bank)
require(
    all(array[ell][s][u] == 0 for ell in range(ROWS) for s in range(Q) for u in range(2, Q)),
    "control escaped common guard-safe pair",
)
require(max(array[ell][s][u] for ell in range(ROWS) for s in range(Q) for u in ROOTS) == 2,
        "aggregate unexpectedly Boolean")

table = mass_table(array)
require(table[0][0] == 1 and all(table[ell][0] == 0 for ell in range(1, ROWS)),
        "hostile lost delta anchor")
require(
    all(table[0][s] == 3 for s in range(1, Q))
    and all(table[ell][s] == 1 for ell in range(1, ROWS) for s in range(1, Q)),
    "hostile occupancy table drifted",
)

mixed_mean_checks = 0
for kappa in MIXED_SOURCE:
    for b in MIXED_TARGET:
        require(table_coefficient(table, kappa, b) != 0, "mixed mean vanished")
        require(coefficient(array, kappa, b, 0) != 0, "root-zero channel vanished")
        mixed_mean_checks += 1
require(mixed_mean_checks == 72, "mixed mean census drifted")


# One seed has all root colours in this control.  Select three; the exact
# Galois index diagonals a=lambda*b then survive on all 72 channels.
seed_active = [a for a in ROOTS if coefficient(array, 1, 1, a) != 0]
require(seed_active == list(ROOTS), "positive control root spectrum drifted")
slopes = tuple(a for a in seed_active if a != 0)[:3]
require(slopes == (1, 2, 3), "slope selection drifted")

incidence_checks = 0
orientation_pair_checks = 0
for kappa in MIXED_SOURCE:
    for b in MIXED_TARGET:
        for slope in slopes:
            a = slope * b % Q
            original = coefficient(array, kappa, b, a)
            require(original != 0, "Galois-uniform root diagonal vanished")
            plus, minus = boundary_aggregate(bank, slope)
            cplus = coefficient(plus, kappa, b, a)
            cminus = coefficient(minus, kappa, b, a)
            multiplier = (pow(ZETA13, a * slope, PRIME) - 1) % PRIME
            require(
                (cminus - cplus) % PRIME == multiplier * original % PRIME,
                "boundary incidence multiplier failed",
            )
            require(cplus != 0 or cminus != 0, "both positive orientations vanished")
            incidence_checks += 1
            orientation_pair_checks += 1
require(incidence_checks == orientation_pair_checks == 216, "incidence census drifted")


# Atomwise Cayley/signless boundary factorisation on every atom and slope.
atomwise_boundary_checks = 0
for ell in range(ROWS):
    for s in range(Q):
        for mask in bank[ell][s]:
            e = bits(mask)
            for tau in MIXED_TARGET:
                plus, minus = boundary(mask, tau)
                gradient = [e[(u + tau) % Q] - e[u] for u in ROOTS]
                require([minus[u] - plus[u] for u in ROOTS] == gradient,
                        "positive boundary difference failed")
                ce = cayley_vector(e, tau)
                signless = [ce[u] + ce[(u + tau) % Q] for u in ROOTS]
                require(signless == gradient, "Cayley signless incidence failed")
                tail, head = selector(mask, tau)
                require(ce[tail] + ce[head] == -1, "selected scalar failed")
                atomwise_boundary_checks += 1


# In the first hostile, the canonical occupied tail is root 1 for both the
# singleton and pair.  Its complete mixed source-target-root transform is
# zero, even though the occupancy table has all 72 mixed means.
selected = selected_aggregate(bank, 1)
plus, minus = boundary_aggregate(bank, 1)
require(selected == plus, "canonical selector differs from plus hostile")
selected_zero_checks = 0
for kappa in MIXED_SOURCE:
    for b in MIXED_TARGET:
        for a in ROOTS:
            require(coefficient(selected, kappa, b, a) == 0,
                    "canonical selected packet retained a mixed mode")
            require(coefficient(plus, kappa, b, a) == 0,
                    "fixed plus orientation retained a mixed mode")
            selected_zero_checks += 1
require(selected_zero_checks == 936, "selected hostile census drifted")

selected_table = mass_table(selected)
require(
    selected_table[0][0] == 1
    and all(selected_table[ell][0] == 0 for ell in range(1, ROWS))
    and all(selected_table[0][s] == 2 for s in range(1, Q))
    and all(selected_table[ell][s] == 1 for ell in range(1, ROWS) for s in range(1, Q)),
    "selected delta-plus-replica table drifted",
)
require(
    all(
        selected[ell][s][u] == 0
        for ell in range(ROWS)
        for s in range(Q)
        for u in ROOTS
        if u != 1
    ),
    "selected hostile lost exact single-root support",
)


# A reflected choice kills the other fixed orientation identically.
reflected_bank = hostile_bank(singleton_root=0)
_reflected_plus, reflected_minus = boundary_aggregate(reflected_bank, 1)
reflected_table = mass_table(reflected_minus)
require(reflected_table == selected_table, "reflected hostile table drifted")
require(
    all(
        reflected_minus[ell][s][u] == 0
        for ell in range(ROWS)
        for s in range(Q)
        for u in ROOTS
        if u != Q - 1
    ),
    "reflected hostile lost exact single-root support",
)
minus_zero_checks = 0
for kappa in MIXED_SOURCE:
    for b in MIXED_TARGET:
        for a in ROOTS:
            require(coefficient(reflected_minus, kappa, b, a) == 0,
                    "fixed minus orientation retained a mixed mode")
            minus_zero_checks += 1
require(minus_zero_checks == 936, "minus hostile census drifted")


print("THM-2539 exact referee")
print(f"field=F_{PRIME} primitive_generator={GENERATOR} zeta7={ZETA7} zeta13={ZETA13}")
print(f"root_invariance_samples={root_invariance_checks} PASS; equation11=universal")
print(f"owner_slots={owner_slots} source_clock_character_checks={clock_character_checks}")
print(f"zero_phase_profiles_annihilated={zero_phase_profiles} full_support_profiles=1")
print(f"guard_vandermonde_minors={guard_vandermonde_checks} sizes=1,2,3")
print("control_anchor=(1,0x6) nonanchor_values=(owner=3,replica=1)")
print(f"mixed_means_nonzero={mixed_mean_checks} simultaneous=PASS")
print(f"seed_active_roots={seed_active} selected_slopes={slopes}")
print(f"galois_uniform_boundary_incidences={incidence_checks} expected=216")
print(f"atomwise_positive_boundary_checks={atomwise_boundary_checks}")
print("aggregate_boolean=False aggregate_max=2 labelled_atoms_boolean=True")
print(f"canonical_plus_zero_mixed_root_checks={selected_zero_checks}")
print(f"reflected_minus_zero_mixed_root_checks={minus_zero_checks}")
print("selector_mass_table=delta_plus_replicas mixed_means=0")
print("hostile_zero_basis=exact_integer_table_plus_single_root_support")
print("scope=full_support_O9c_only; THM-2540_proved_separate_scope")
print("THM-2539 PASS")
