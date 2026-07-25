#!/usr/bin/env python3
"""Exact deep-pair signed-transfer Bellman referee for THM-2266.

In the remaining scalar first-depth-one profiles ``(1,b,c)``, put

    X_(j,t) = 1_{D_(u_j)}(13^t x),       rho = -1/13.

The ordinary centered atoms

    Z_(j,t) = X_(j,t) - rho^t X_(j,0) + max(rho^t, 0)

are response-null, nonnegative, and dominate ``X_(j,t)``.  The time-zero
literal for the first blocker is retained directly.  For the two deep
negative literals, this script additionally centers their product:

    F  = X_(2,b-1) X_(3,c-1),
    s0 = (b-1) mod 2,
    d  = b-1-s0,
    F0 = X_(2,s0) X_(3,c-b+s0),
    W  = F - 13^(-d) F0.

Here ``d`` is even and ``F=F0 o T^d``, so ``integral R W=0``.  Therefore

    A = X_(1,0) + Z_(2,b-1) + Z_(3,c-1) - W

has the same signed response as ``X_(1,0)`` while satisfying

    A >= x1+x2+x3-x2*x3 >= 1

on the negative carrier.  The exact Bellman computation below maximizes
``E[1_P A]`` over every cross-core coupling with the exact one-coordinate
root marginals.  Every small rational coupling LP receives matching exact
primal and dual certificates.

The companion freezes the complete 150-profile table, the exact 120/30
subtarget split, the weakest subtarget margin, the reduced-product cutoff
758, and the exact boundary complement.  A subtarget Bellman value does not
by itself exclude a profile: the theorem concludes a finite relation-atlas
alternative through the remaining response ``integral R X_(1,0)``.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from math import gcd


P = 13
RHO = Fraction(-1, P)
TARGET = Fraction(961, 6930)
BIT_VECTORS = tuple(product((0, 1), repeat=3))
COLUMNS = tuple(
    (Fraction(1),) + tuple(Fraction(bit) for bit in bits)
    for bits in BIT_VECTORS
)
RANK = 4
TERMINAL_KEY = len(BIT_VECTORS)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def invert_square(matrix):
    """Return the exact inverse of a square matrix, or None if singular."""
    n = len(matrix)
    augmented = [
        list(row)
        + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            (
                row
                for row in range(column, n)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [
            entry / scale for entry in augmented[column]
        ]
        for row in range(n):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    augmented[row][j] - scale * augmented[column][j]
                    for j in range(2 * n)
                ]
    return tuple(tuple(augmented[i][n:]) for i in range(n))


def matrix_vector(matrix, vector):
    return tuple(
        sum(
            (matrix[i][j] * vector[j] for j in range(len(vector))),
            Fraction(),
        )
        for i in range(len(matrix))
    )


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), Fraction())


def basis_matrix(basis):
    return tuple(
        tuple(COLUMNS[basis[column]][row] for column in range(RANK))
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
            Fraction(1, 7),
            Fraction(1, 7),
            Fraction(1, 7),
        )
    future = BIT_VECTORS[key]
    return (
        Fraction(1),
        Fraction(2 - future[0], P),
        Fraction(2 - future[1], P),
        Fraction(2 - future[2], P),
    )


@lru_cache(maxsize=None)
def coupling_vertices(key):
    """All exact vertices, together with every basis producing each one."""
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
    """Exact primal vertex maximum and an independent exact basic dual."""
    global DUAL_CERTIFICATES
    rhs = marginal_rhs(key)
    primal = None
    optimizers = []
    for distribution, bases, support in coupling_vertices(key):
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
            # A_B^T y=c_B, hence y=(A_B^-1)^T c_B.
            dual = tuple(
                sum(
                    (
                        inverse[row][column] * basic_values[row]
                        for row in range(RANK)
                    ),
                    Fraction(),
                )
                for column in range(RANK)
            )
            if any(
                dot(column, dual) < value
                for column, value in zip(COLUMNS, values)
            ):
                continue
            require(dot(rhs, dual) == primal, "primal/dual value drift")
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


def strict_profiles():
    return {
        (1, middle, deepest)
        for deepest in range(5, 20)
        for middle in range(2, deepest)
    }


def interior_profiles():
    return {
        profile
        for profile in strict_profiles()
        if 3 <= profile[1] <= profile[2] - 2
    }


def boundary_profiles():
    return strict_profiles() - interior_profiles()


def pair_packet(profile):
    """Return all times and coefficients in the centered pair packet."""
    one, middle, deepest = profile
    require(one == 1, "profile is not first-depth-one")
    require(2 <= middle < deepest, "profile is not strict")
    require(5 <= deepest <= 19, "profile is outside the scalar ledger")

    current_two = middle - 1
    current_three = deepest - 1
    parity_base = current_two % 2
    shift_time = current_two - parity_base
    reference_two = parity_base
    reference_three = deepest - middle + parity_base
    epsilon = Fraction(1, P**shift_time)

    require(shift_time >= 0, "negative transfer shift")
    require(shift_time % 2 == 0, "pair transfer shift is not even")
    require(RHO**shift_time == epsilon, "pair transfer phase drift")
    require(
        reference_two + shift_time == current_two,
        "second-core pair shift drift",
    )
    require(
        reference_three + shift_time == current_three,
        "third-core pair shift drift",
    )

    return {
        "current_two": current_two,
        "current_three": current_three,
        "reference_two": reference_two,
        "reference_three": reference_three,
        "shift_time": shift_time,
        "epsilon": epsilon,
        "rho_two": RHO**current_two,
        "rho_three": RHO**current_three,
    }


def audit_pointwise_packet(profile):
    """Hostile bit audit for nonnegativity and carrier domination."""
    packet = pair_packet(profile)
    rho_two = packet["rho_two"]
    rho_three = packet["rho_three"]
    epsilon = packet["epsilon"]
    for (
        x_one,
        x_two_zero,
        x_two_current,
        x_three_zero,
        x_three_current,
        x_two_reference,
        x_three_reference,
    ) in product((0, 1), repeat=7):
        z_two = (
            x_two_current
            - rho_two * x_two_zero
            + max(rho_two, Fraction())
        )
        z_three = (
            x_three_current
            - rho_three * x_three_zero
            + max(rho_three, Fraction())
        )
        current_pair = x_two_current * x_three_current
        reference_pair = x_two_reference * x_three_reference
        score = (
            x_one
            + z_two
            + z_three
            - current_pair
            + epsilon * reference_pair
        )
        lower = (
            x_one
            + x_two_current
            + x_three_current
            - current_pair
        )
        require(z_two >= x_two_current, "second centered atom lost domination")
        require(
            z_three >= x_three_current,
            "third centered atom lost domination",
        )
        require(score >= lower, "deep-pair score lost carrier lower bound")
        require(score >= 0, "deep-pair score became negative")
        if x_one or x_two_current or x_three_current:
            require(score >= 1, "negative carrier score fell below one")


@lru_cache(maxsize=None)
def deep_pair_capacity(profile):
    """Robust capacity of the depth-one centered deep-pair signed dual."""
    packet = pair_packet(profile)
    _, middle, deepest = profile
    audit_pointwise_packet(profile)

    positive_by_time = {
        1: (0,),
        middle: (1,),
        deepest: (2,),
    }
    # Slot order: current core 2, current core 3, reference core 2,
    # reference core 3.
    slot_specs = (
        (1, packet["current_two"]),
        (2, packet["current_three"]),
        (1, packet["reference_two"]),
        (2, packet["reference_three"]),
    )
    slots_by_time = {}
    for slot, (core, time) in enumerate(slot_specs):
        slots_by_time.setdefault(time, []).append((slot, core))
    slots_by_time = {
        time: tuple(entries) for time, entries in slots_by_time.items()
    }

    @lru_cache(maxsize=None)
    def recurse(time, positive_hit, slot_mask, next_key):
        if time < 0:
            if not positive_hit:
                return Fraction()
            x_zero = BIT_VECTORS[next_key]
            current_two = (slot_mask >> 0) & 1
            current_three = (slot_mask >> 1) & 1
            reference_two = (slot_mask >> 2) & 1
            reference_three = (slot_mask >> 3) & 1
            value = (
                Fraction(x_zero[0] + current_two + current_three)
                - packet["rho_two"] * x_zero[1]
                - packet["rho_three"] * x_zero[2]
                + max(packet["rho_two"], Fraction())
                + max(packet["rho_three"], Fraction())
                - current_two * current_three
                + packet["epsilon"] * reference_two * reference_three
            )
            require(value >= 0, "terminal score became negative")
            return value

        costs = []
        for current_key, current_bits in enumerate(BIT_VECTORS):
            next_hit = positive_hit or any(
                current_bits[core]
                for core in positive_by_time.get(time, ())
            )
            next_mask = slot_mask
            for slot, core in slots_by_time.get(time, ()):
                if current_bits[core]:
                    next_mask |= 1 << slot
            costs.append(
                recurse(time - 1, next_hit, next_mask, current_key)
            )
        return maximize_joint_bits(tuple(costs), next_key)

    terminal_values = tuple(
        recurse(deepest, False, 0, terminal_key)
        for terminal_key in range(len(BIT_VECTORS))
    )
    return maximize_joint_bits(terminal_values, TERMINAL_KEY)


def digest_profile_set(profiles):
    payload = ";".join(
        ",".join(map(str, profile)) for profile in sorted(profiles)
    )
    return sha256(payload.encode("ascii")).hexdigest()


def digest_bounds(bounds):
    payload = ";".join(
        f"{','.join(map(str, profile))}:{bound.numerator}/{bound.denominator}"
        for profile, bound in sorted(bounds.items())
    )
    return sha256(payload.encode("ascii")).hexdigest()


def strict_relation_cutoff(gap):
    """Least N with 6/(4N)<gap."""
    candidate = 1
    while Fraction(3, 2 * candidate) >= gap:
        candidate += 1
    require(Fraction(3, 2 * candidate) < gap, "relation cutoff failed")
    if candidate > 1:
        require(
            Fraction(3, 2 * (candidate - 1)) >= gap,
            "relation cutoff is not minimal",
        )
    return candidate


def main():
    strict = strict_profiles()
    interior = interior_profiles()
    boundary = boundary_profiles()
    require(len(strict) == 150, "strict-profile census drift")
    require(len(interior) == 120, "interior-profile census drift")
    require(len(boundary) == 30, "boundary-profile census drift")

    expected_boundary = (
        {(1, 2, deepest) for deepest in range(5, 20)}
        | {(1, middle, middle + 1) for middle in range(4, 19)}
    )
    require(boundary == expected_boundary, "boundary decomposition drift")

    bounds = {
        profile: deep_pair_capacity(profile)
        for profile in sorted(strict)
    }
    subtarget = {
        profile for profile, bound in bounds.items() if bound < TARGET
    }
    nonsubtarget = strict - subtarget
    require(subtarget == interior, "subtarget profile set drift")
    require(nonsubtarget == boundary, "nonsubtarget profile set drift")

    weakest = max(subtarget, key=lambda profile: bounds[profile])
    require(weakest == (1, 3, 5), "weakest subtarget profile drift")
    expected_weakest_bound = Fraction(10_146_800_138, 74_231_495_611)
    expected_weakest_gap = Fraction(
        145_591_760_833,
        73_489_180_654_890,
    )
    require(
        bounds[weakest] == expected_weakest_bound,
        "weakest subtarget bound drift",
    )
    require(
        TARGET - bounds[weakest] == expected_weakest_gap,
        "weakest subtarget gap drift",
    )

    frozen_spot_checks = {
        (1, 4, 7): Fraction(
            283_293_422_086_229,
            2_120_125_746_145_771,
        ),
        (1, 4, 8): Fraction(
            47_679_834_772_180_172,
            358_301_251_098_635_299,
        ),
        (1, 10, 19): Fraction(
            152_389_911_763_114_035_508_249_271_108_093_757_321_845,
            1_150_805_888_298_461_593_768_523_506_367_492_533_368_331,
        ),
    }
    for profile, expected in frozen_spot_checks.items():
        require(bounds[profile] == expected, f"spot-check drift at {profile}")

    cutoffs = {
        profile: strict_relation_cutoff(TARGET - bounds[profile])
        for profile in sorted(subtarget)
    }
    cutoff_histogram = tuple(sorted(Counter(cutoffs.values()).items()))
    require(
        cutoff_histogram
        == (
            (240, 28),
            (241, 38),
            (264, 11),
            (265, 1),
            (268, 12),
            (297, 1),
            (298, 12),
            (299, 1),
            (342, 1),
            (467, 11),
            (468, 2),
            (584, 1),
            (758, 1),
        ),
        "profile-cutoff histogram drift",
    )
    uniform_cutoff = max(cutoffs.values())
    require(uniform_cutoff == 758, "uniform relation cutoff drift")
    require(
        cutoffs[weakest] == uniform_cutoff,
        "weakest row does not control the uniform cutoff",
    )
    relation_tail = Fraction(3, 2 * uniform_cutoff)
    relation_margin = expected_weakest_gap - relation_tail
    require(
        relation_margin
        == Fraction(
            124_783_729_079,
            55_704_798_936_406_620,
        ),
        "relation-atlas margin drift",
    )
    require(relation_margin > 0, "relation-atlas margin is not positive")
    require(
        expected_weakest_gap - Fraction(3, 2 * (uniform_cutoff - 1))
        == Fraction(
            -10_404_015_877,
            27_815_654_877_875_865,
        ),
        "minimal-cutoff hostile boundary drift",
    )
    atlas_bound = uniform_cutoff - 1
    oriented_atlas = tuple(
        (left, right)
        for left in range(1, atlas_bound + 1)
        for right in range(1, atlas_bound // left + 1)
        if gcd(left, right) == 1
    )
    unordered_atlas = tuple(
        pair for pair in oriented_atlas if pair[0] <= pair[1]
    )
    require(len(oriented_atlas) == 3_643, "oriented atlas census drift")
    require(len(unordered_atlas) == 1_822, "unordered atlas census drift")

    bound_digest = digest_bounds(bounds)
    subtarget_digest = digest_profile_set(subtarget)
    boundary_digest = digest_profile_set(boundary)
    cutoff_digest = digest_bounds(
        {
            profile: Fraction(cutoff)
            for profile, cutoff in cutoffs.items()
        }
    )
    require(
        bound_digest
        == "0f5e660eecf284406777d66781d04a2528fafac714b7ff99ca49b9462f890847",
        "complete bound-table digest drift",
    )
    require(
        subtarget_digest
        == "b5259a91c2a9fca337cce681e21455040a3c944fe8c63e86c38d2b803d899189",
        "subtarget-profile digest drift",
    )
    require(
        boundary_digest
        == "c5e141aba821f6f285f1c9c8c44211bb2843e7d0c31687c2ea59ca07efdd077f",
        "boundary-profile digest drift",
    )
    require(
        cutoff_digest
        == "3874ce2ba7b2c26e873c4a028883f94673d58403f3f0bf401f5b161282472d6f",
        "profile-cutoff digest drift",
    )

    cache = maximize_joint_bits_uncounted.cache_info()
    require(
        DUAL_CERTIFICATES == cache.misses,
        "not every distinct coupling LP received a dual certificate",
    )
    require(LP_CALLS == 118_970, "coupling LP call count drift")
    require(cache.misses == 77_646, "distinct coupling LP count drift")

    print(f"target={TARGET}")
    print(f"strict_profiles={len(strict)}")
    print(f"subtarget_profiles={len(subtarget)}")
    print(f"nonsubtarget_profiles={len(nonsubtarget)}")
    print(f"invertible_bases={len(INVERTIBLE_BASES)}")
    print(
        "root_vertex_counts="
        + str(tuple(len(coupling_vertices(key)) for key in range(8)))
    )
    print(f"terminal_vertex_count={len(coupling_vertices(TERMINAL_KEY))}")
    print(
        f"weakest_subtarget={weakest} bound={bounds[weakest]} "
        f"target_minus={expected_weakest_gap}"
    )
    print(
        f"uniform_reduced_product_cutoff={uniform_cutoff} "
        f"s_floor=-{relation_tail} "
        f"final_margin={relation_margin}"
    )
    print(
        f"small_relation_product_max={atlas_bound} "
        f"oriented_coprime_shapes={len(oriented_atlas)} "
        f"unordered_coprime_shapes={len(unordered_atlas)}"
    )
    print(f"profile_cutoff_histogram={cutoff_histogram}")
    print(
        "cutoff_757_hostile_difference="
        + str(
            expected_weakest_gap
            - Fraction(3, 2 * (uniform_cutoff - 1))
        )
    )
    print(
        "boundary_structure="
        "15 rows (1,2,c), 5<=c<=19; "
        "15 rows (1,b,b+1), 4<=b<=18"
    )
    print(f"bound_digest={bound_digest}")
    print(f"subtarget_digest={subtarget_digest}")
    print(f"boundary_digest={boundary_digest}")
    print(f"cutoff_digest={cutoff_digest}")
    print(
        f"lp_calls={LP_CALLS} distinct_lps={cache.misses} "
        f"dual_certificates={DUAL_CERTIFICATES}"
    )
    print("exact_primal_dual_vertex_parity=PASS")
    print("complete_profile_table_begin")
    for profile in sorted(strict):
        bound = bounds[profile]
        difference = TARGET - bound
        status = "SUBTARGET" if profile in subtarget else "NONSUBTARGET"
        cutoff = cutoffs.get(profile, "-")
        print(
            f"profile={profile} bound={bound} "
            f"target_minus={difference} status={status} "
            f"relation_cutoff={cutoff}"
        )
    print("complete_profile_table_end")


if __name__ == "__main__":
    main()
