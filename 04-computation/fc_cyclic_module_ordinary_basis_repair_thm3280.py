#!/usr/bin/env python3
"""Exact regression checks for THM-3280's cyclic-module repair."""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations

from sympy import Matrix, diag, linear_eq_to_matrix, symbols


def require(condition, message):
    if not condition:
        raise AssertionError(message)


ledger = []
commutant_profiles = []
trace_pairing_profiles = []
cycle_subset_checks = 0
zero_degree_profiles = []
degree_bounds = []

for d in range(3, 11):
    r = d - 1

    # The general weighted nonzero cycle is diagonally similar over Qbar to a
    # scalar multiple of this cycle.  D is unchanged by the diagonal gauge.
    cycle = Matrix.zeros(r)
    for j in range(r - 1):
        cycle[j, j + 1] = 1
    cycle[r - 1, 0] = 1
    residue = diag(*[Fraction(-j, d) for j in range(1, r + 1)])

    entries = symbols(f"x0:{r*r}")
    unknown = Matrix(r, r, entries)
    equations = list(cycle * unknown - unknown * cycle)
    equations += list(residue * unknown - unknown * residue)
    coefficient_matrix, _ = linear_eq_to_matrix(equations, entries)
    nullity = r * r - coefficient_matrix.rank()
    require(nullity == 1, f"simultaneous commutant not scalar at d={d}")
    commutant_profiles.append((d, nullity))

    centralizer_basis = [cycle**k for k in range(r)]
    gram = Matrix(
        r,
        r,
        lambda i, j: (centralizer_basis[i] * centralizer_basis[j]).trace(),
    )
    determinant = gram.det()
    require(determinant != 0, f"centralizer trace pairing degenerate at d={d}")
    trace_pairing_profiles.append((d, abs(int(determinant))))

    successor = {0: r - 1, **{j: j - 1 for j in range(1, r)}}
    for mask in range(1, (1 << r) - 1):
        subset = {j for j in range(r) if mask & (1 << j)}
        require(
            any(successor[j] not in subset for j in subset),
            f"proper coordinate subset invariant at d={d}, subset={subset}",
        )
        cycle_subset_checks += 1

    if d <= 9:
        bound = Fraction((d - 1) ** 2, 8 * d)
        require(bound < 1, f"Pluecker bound not coercive at d={d}")
        degree_bounds.append((d, str(bound)))
        zero_count = 0
        for m in range(1, r):
            for subset in combinations(range(1, d), m):
                degree = Fraction(sum(subset), d) - Fraction(m, 2)
                if degree >= 0 and degree.denominator == 1:
                    require(degree == 0, f"positive degree profile before d=10: {d,m,subset}")
                    zero_count += 1
        zero_degree_profiles.append((d, zero_count))

d10_hostiles = []
d = 10
for m in range(1, d - 1):
    for subset in combinations(range(1, d), m):
        degree = Fraction(sum(subset), d) - Fraction(m, 2)
        if degree > 0 and degree.denominator == 1:
            d10_hostiles.append((m, subset, int(degree)))

require(
    d10_hostiles
    == [
        (4, (6, 7, 8, 9), 1),
        (5, (5, 6, 7, 8, 9), 1),
    ],
    f"unexpected d=10 hostile list: {d10_hostiles}",
)

ledger.append(f"C1:simultaneous_commutant_nullities={commutant_profiles}")
ledger.append(f"C2:centralizer_trace_gram_abs_determinants={trace_pairing_profiles}")
ledger.append(f"C3:proper_cycle_subset_checks={cycle_subset_checks};all_noninvariant=YES")
ledger.append(f"C4:integer_degree_zero_profiles={zero_degree_profiles}")
ledger.append(f"C5:pluecker_upper_bounds={degree_bounds};all_strictly_below_one=YES")
ledger.append(f"C6:d10_degree_one_bound_hostiles={d10_hostiles};submodule_claim=NO")

digest = sha256("\n".join(ledger).encode("utf-8")).hexdigest()

print("THM-3280 FC CYCLIC-MODULE ORDINARY-BASIS REPAIR AUDIT")
for row in ledger:
    print(row)
print(f"semantic_sha256={digest}")
print("CONCLUSION=BINOMIAL_VALUE_SPECIALIZATION_REPAIRED_DEGREES_3_TO_9;D10_PLUCKER_BOUND_OPEN")
