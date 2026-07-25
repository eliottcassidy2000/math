#!/usr/bin/env python3
"""Exact guard/danger hidden-state Bellman referee for THM-2233.

The state retains the guard bit and the three unit-core danger bits at every
multiplication-by-13 time.  Conditional on the four future bits, the actual
thirteen-root law has the exact one-coordinate marginals

    guard:  (10 - future_guard) / 13,
    danger: ( 2 - future_danger) / 13.

Every cross-coordinate coupling is allowed.  Each resulting four-Bernoulli
polytope is enumerated exactly once.  Bellman objectives are maximized over
its exact vertices, and an exact complementary-slackness dual certificate is
constructed for every distinct LP call.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product


P = 13
TARGET = Fraction(961, 6930)
BIT_VECTORS = tuple(product((0, 1), repeat=4))
COLUMNS = tuple(
    (Fraction(1),) + tuple(Fraction(bit) for bit in bits)
    for bits in BIT_VECTORS
)
RANK = 5
TERMINAL_KEY = len(BIT_VECTORS)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def invert_square(matrix):
    """Return the exact inverse of a square matrix, or None if singular."""
    n = len(matrix)
    aug = [
        list(row)
        + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for col in range(n):
        pivot = next((row for row in range(col, n) if aug[row][col]), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [entry / scale for entry in aug[col]]
        for row in range(n):
            if row == col:
                continue
            scale = aug[row][col]
            if scale:
                aug[row] = [
                    aug[row][j] - scale * aug[col][j]
                    for j in range(2 * n)
                ]
    return tuple(tuple(aug[i][n:]) for i in range(n))


def matrix_vector(matrix, vector):
    return tuple(
        sum((matrix[i][j] * vector[j] for j in range(len(vector))), Fraction())
        for i in range(len(matrix))
    )


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), Fraction())


def basis_matrix(basis):
    """Feature matrix whose columns are the selected Bernoulli atoms."""
    return tuple(
        tuple(COLUMNS[basis[col]][row] for col in range(RANK))
        for row in range(RANK)
    )


INVERTIBLE_BASES = []
for candidate_basis in combinations(range(len(BIT_VECTORS)), RANK):
    candidate_inverse = invert_square(basis_matrix(candidate_basis))
    if candidate_inverse is not None:
        INVERTIBLE_BASES.append((candidate_basis, candidate_inverse))
INVERTIBLE_BASES = tuple(INVERTIBLE_BASES)
INVERSE_BY_BASIS = dict(INVERTIBLE_BASES)


def marginal_rhs(key):
    if key == TERMINAL_KEY:
        return (
            Fraction(1),
            Fraction(5, 7),
            Fraction(1, 7),
            Fraction(1, 7),
            Fraction(1, 7),
        )
    bits = BIT_VECTORS[key]
    return (
        Fraction(1),
        Fraction(10 - bits[0], P),
        Fraction(2 - bits[1], P),
        Fraction(2 - bits[2], P),
        Fraction(2 - bits[3], P),
    )


@lru_cache(maxsize=None)
def coupling_vertices(key):
    """All distinct exact vertices, with every basis producing each vertex."""
    rhs = marginal_rhs(key)
    by_distribution = {}
    for basis, inverse in INVERTIBLE_BASES:
        weights = matrix_vector(inverse, rhs)
        if any(weight < 0 for weight in weights):
            continue
        distribution = [Fraction() for _ in BIT_VECTORS]
        for column, weight in zip(basis, weights):
            distribution[column] = weight
        distribution = tuple(distribution)
        by_distribution.setdefault(distribution, []).append(basis)
    require(by_distribution, f"no coupling vertices for key {key}")
    return tuple(
        (
            distribution,
            tuple(bases),
            tuple(
                (index, weight)
                for index, weight in enumerate(distribution)
                if weight
            ),
        )
        for distribution, bases in by_distribution.items()
    )


LP_CALLS = 0
DUAL_CERTIFICATES = 0


@lru_cache(maxsize=None)
def maximize_joint_bits_uncounted(values, key):
    """Exact primal vertex maximum plus exact basic dual certificate."""
    global DUAL_CERTIFICATES
    rhs = marginal_rhs(key)
    vertices = coupling_vertices(key)
    primal = None
    optimizers = []
    for distribution, bases, support in vertices:
        score = sum(
            (weight * values[index] for index, weight in support),
            Fraction(),
        )
        if primal is None or score > primal:
            primal = score
            optimizers = [(distribution, bases)]
        elif score == primal:
            optimizers.append((distribution, bases))

    require(primal is not None, "coupling polytope has no scored vertex")
    for distribution, bases in optimizers:
        for basis in bases:
            inverse = INVERSE_BY_BASIS[basis]
            basic_values = tuple(values[index] for index in basis)
            # If A_B is the feature matrix, A_B^T y = c_B, so
            # y=(A_B^{-1})^T c_B.
            dual = tuple(
                sum(
                    (inverse[row][col] * basic_values[row]
                     for row in range(RANK)),
                    Fraction(),
                )
                for col in range(RANK)
            )
            if any(dot(column, dual) < value
                   for column, value in zip(COLUMNS, values)):
                continue
            dual_value = dot(rhs, dual)
            require(dual_value == primal, "exact primal/dual value drift")
            require(
                dot(distribution, values) == primal,
                "selected primal distribution drift",
            )
            DUAL_CERTIFICATES += 1
            return primal
    raise RuntimeError("no exact dual-feasible optimal basis")


def maximize_joint_bits(values, key):
    global LP_CALLS
    LP_CALLS += 1
    return maximize_joint_bits_uncounted(tuple(values), key)


def active_incidences(profile):
    """Map unit times to the guard/checkpoint clause incidences."""
    first = profile[0]
    checkpoints = tuple(range(0, first + 1, 2))
    by_time = {0: (("guard", 0),)}
    mutable = {time: list(rows) for time, rows in by_time.items()}
    for clause, checkpoint in enumerate(checkpoints, start=1):
        for blocker, valuation in enumerate(profile):
            time = valuation - checkpoint
            require(time >= 0, "negative unit-time offset")
            mutable.setdefault(time, []).append((blocker, clause))
    return (
        1 + len(checkpoints),
        {time: tuple(rows) for time, rows in mutable.items()},
    )


@lru_cache(maxsize=None)
def profile_bound(profile):
    """Exact guard/checkpoint Bellman upper bound for one scalar profile."""
    require(1 <= profile[0] <= 5, "profile outside THM-2233 census")
    clause_count, by_time = active_incidences(profile)
    all_clauses = (1 << clause_count) - 1
    value = tuple(
        tuple(Fraction(int(satisfied == all_clauses)) for _ in BIT_VECTORS)
        for satisfied in range(1 << clause_count)
    )

    for time in range(max(by_time) + 1):
        incidences = by_time.get(time, ())
        hit_masks = []
        for bits in BIT_VECTORS:
            newly = 0
            for source, clause in incidences:
                if source == "guard":
                    active = bits[0]
                else:
                    active = bits[source + 1]
                if active:
                    newly |= 1 << clause
            hit_masks.append(newly)

        updated = []
        for later in range(1 << clause_count):
            row = []
            for future_key in range(len(BIT_VECTORS)):
                costs = tuple(
                    value[later | hit_masks[current]][current]
                    for current in range(len(BIT_VECTORS))
                )
                row.append(maximize_joint_bits(costs, future_key))
            updated.append(tuple(row))
        value = tuple(updated)

    return maximize_joint_bits(value[0], TERMINAL_KEY)


def all_profiles():
    return tuple(
        (first, middle, deepest)
        for deepest in range(2, 20)
        for middle in range(1, deepest)
        for first in range(1, middle + 1)
    )


def profiles_with_first(first):
    return tuple(
        (first, middle, deepest)
        for middle in range(first, 19)
        for deepest in range(middle + 1, 20)
    )


def digest_profiles(profiles):
    payload = ";".join(",".join(map(str, profile)) for profile in profiles)
    return sha256(payload.encode("ascii")).hexdigest()


def thm2229_raw_closed():
    closed = set()
    for profile in profiles_with_first(2):
        if (
            (profile[1] == 3 and profile[2] >= 6)
            or (profile[1] >= 5 and profile[2] != profile[1] + 2)
        ):
            closed.add(profile)
    for profile in profiles_with_first(3):
        if (
            (profile[1] == 3 and (profile[2] == 4 or profile[2] >= 6))
            or (profile[1] == 4 and profile[2] >= 7)
            or (profile[1] >= 6 and profile[2] != profile[1] + 2)
        ):
            closed.add(profile)
    closed.update(set(profiles_with_first(4)) - {(4, 6, 8)})
    closed.update(set(profiles_with_first(5)) - {(5, 7, 9)})
    return closed


def expected_guard_closed():
    """Closed-profile classification frozen independently of the recursion."""
    closed = set()
    for profile in profiles_with_first(2):
        if (
            (profile[1] == 2 and profile[2] >= 5)
            or (profile[1] == 3 and profile[2] >= 6)
            or (profile[1] >= 5 and profile[2] != profile[1] + 2)
        ):
            closed.add(profile)
    for profile in profiles_with_first(3):
        if (
            (profile[1] == 3 and (profile[2] == 4 or profile[2] >= 6))
            or (profile[1] == 4 and profile[2] >= 7)
            or (profile[1] >= 6 and profile[2] != profile[1] + 2)
        ):
            closed.add(profile)
    closed.update(set(profiles_with_first(4)) - {(4, 6, 8)})
    closed.update(profiles_with_first(5))
    return closed


def main():
    require(len(BIT_VECTORS) == 16, "four-bit state count drift")
    require(len(INVERTIBLE_BASES) > 0, "invertible basis census empty")
    for key in range(TERMINAL_KEY + 1):
        coupling_vertices(key)

    universe = set(all_profiles())
    require(len(universe) == 1140, "legal scalar-profile universe drift")
    scoped = tuple(
        profile for profile in all_profiles() if profile[0] <= 5
    )
    require(len(scoped) == 685, "first-depth-at-most-five census drift")

    bounds = {profile: profile_bound(profile) for profile in scoped}
    raw_closed = {
        profile for profile, bound in bounds.items() if bound < TARGET
    }
    expected_raw_closed = expected_guard_closed()
    require(
        raw_closed == expected_raw_closed,
        "guard-hidden raw closure classification drift",
    )

    low_exact = {
        (1, 1, 2),
        (1, 1, 3),
        (1, 2, 3),
        (2, 2, 3),
        (1, 1, 4),
        (1, 2, 4),
        (1, 3, 4),
        (2, 2, 4),
        (2, 3, 4),
        (3, 3, 4),
    }
    high_closed = {profile for profile in universe if profile[0] >= 6}
    hidden_parity_closed = {
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    }
    positive_set_closed = thm2229_raw_closed()
    current_closed = (
        low_exact | high_closed | hidden_parity_closed | positive_set_closed
    )
    require(len(current_closed) == 900, "pre-THM-2233 ledger drift")
    require(len(universe - current_closed) == 240, "240-profile residue drift")

    novel = raw_closed - current_closed
    combined_closed = current_closed | raw_closed
    remaining = tuple(sorted(universe - combined_closed))
    expected_novel = {
        (2, 2, deepest) for deepest in range(5, 20)
    } | {(5, 7, 9)}
    require(novel == expected_novel, "novel-profile classification drift")
    expected_remaining = (
        {
            profile
            for profile in universe
            if profile[0] == 1 and profile[2] >= 5
        }
        | {(2, 3, 5)}
        | {(2, 4, deepest) for deepest in range(5, 20)}
        | {(2, middle, middle + 2) for middle in range(5, 18)}
        | {(3, 3, 5), (3, 4, 5), (3, 4, 6)}
        | {(3, 5, deepest) for deepest in range(6, 20)}
        | {(3, middle, middle + 2) for middle in range(6, 18)}
        | {(4, 6, 8)}
    )
    require(
        set(remaining) == expected_remaining,
        "combined remaining-profile classification drift",
    )

    geometric = ((4, 6, 8), (5, 7, 9))
    expected_geometric_bounds = {
        (4, 6, 8): Fraction(800_688_435, 5_710_115_047),
        (5, 7, 9): Fraction(10_149_510_390, 74_231_495_611),
    }
    for profile, expected in expected_geometric_bounds.items():
        require(bounds[profile] == expected, f"geometric bound drift: {profile}")

    by_first = {
        first: {
            profile for profile in raw_closed if profile[0] == first
        }
        for first in range(1, 6)
    }
    novel_by_first = {
        first: {profile for profile in novel if profile[0] == first}
        for first in range(1, 6)
    }
    require(
        tuple(len(by_first[first]) for first in range(1, 6))
        == (0, 121, 107, 119, 105),
        "per-first-depth raw closure counts drift",
    )
    require(len(raw_closed) == 452, "raw closure total drift")
    require(len(raw_closed & low_exact) == 1, "low-exact overlap drift")
    require(
        len(raw_closed & positive_set_closed) == 436,
        "THM-2229 overlap drift",
    )
    require(
        raw_closed & hidden_parity_closed == hidden_parity_closed,
        "THM-2227 overlap drift",
    )
    require(len(remaining) == 224, "combined remaining ledger drift")
    raw_digest = digest_profiles(tuple(sorted(raw_closed)))
    novel_digest = digest_profiles(tuple(sorted(novel)))
    remaining_digest = digest_profiles(remaining)
    require(
        raw_digest
        == "1534c393cd40c008e958ac6376defeb6ab050d6a2fe82f0b6bc74756f7b3c19f",
        "raw-closed digest drift",
    )
    require(
        novel_digest
        == "1a038c930394ed7be3e7410d875c5cc09409543b1f3b3f5405ba2aeaff376079",
        "novel-profile digest drift",
    )
    require(
        remaining_digest
        == "d7617205a31504563f5f3d14eebc7d6dfd6148349f43a468949a8ecd850f8505",
        "remaining-profile digest drift",
    )

    print("THM-2233 GUARD-DANGER HIDDEN-STATE BELLMAN")
    print(f"target={TARGET}")
    print(f"bit_states={len(BIT_VECTORS)}")
    print(f"invertible_bases={len(INVERTIBLE_BASES)}")
    print(
        "coupling_vertex_counts="
        + repr(tuple(len(coupling_vertices(key))
                     for key in range(TERMINAL_KEY + 1)))
    )
    for profile in geometric:
        bound = bounds[profile]
        relation = "<" if bound < TARGET else ">="
        margin = TARGET - bound if bound < TARGET else bound - TARGET
        print(
            f"profile={profile} bound={bound} decimal={float(bound):.15f} "
            f"relation={relation}target absolute_margin={margin}"
        )
    for first in range(1, 6):
        print(
            f"lambda1={first} raw_closed={len(by_first[first])}/"
            f"{len(profiles_with_first(first))} "
            f"novel_vs_240_ledger={len(novel_by_first[first])}"
        )
    print(f"raw_closed_total={len(raw_closed)}")
    print(f"overlap_low_exact={len(raw_closed & low_exact)}")
    print(f"overlap_thm2227={len(raw_closed & hidden_parity_closed)}")
    print(f"overlap_thm2229={len(raw_closed & positive_set_closed)}")
    print(f"novel_vs_240_ledger={len(novel)}")
    print(f"novel_profiles={tuple(sorted(novel))}")
    print(f"combined_remaining={len(remaining)}")
    print("remaining_lambda1_1=all 165 profiles with lambda3>=5")
    print(
        "remaining_lambda1_2=(2,3,5); (2,4,c) for 5<=c<=19; "
        "(2,b,b+2) for 5<=b<=17"
    )
    print(
        "remaining_lambda1_3=(3,3,5),(3,4,5),(3,4,6); "
        "(3,5,c) for 6<=c<=19; (3,b,b+2) for 6<=b<=17"
    )
    print("remaining_high_first=((4,6,8),)")
    print(f"raw_closed_digest={raw_digest}")
    print(f"novel_digest={novel_digest}")
    print(f"remaining_digest={remaining_digest}")
    print(f"lp_calls={LP_CALLS}")
    print(f"lp_cache_hits={maximize_joint_bits_uncounted.cache_info().hits}")
    print(f"distinct_exact_lp_calls={maximize_joint_bits_uncounted.cache_info().misses}")
    print(f"exact_dual_certificates={DUAL_CERTIFICATES}")
    require(
        DUAL_CERTIFICATES == maximize_joint_bits_uncounted.cache_info().misses,
        "not every distinct LP received an exact dual certificate",
    )
    print("exact_primal_dual_vertex_parity=PASS")
    print("status=THM2233_GUARD_DANGER_HIDDEN_STATE_BELLMAN")


if __name__ == "__main__":
    main()
