#!/usr/bin/env python3
"""Exact multicore signed-transfer Bellman referee for THM-2239.

For a unit core u and X_t=1_{D_u}(13^t x), put rho=-1/13 and

    Y_t = X_t - rho^t X_0,
    Z_t = Y_t + max(rho^t, 0).

The signed residual R satisfies integral R=0 and
integral R X_t=rho^t integral R X_0, hence integral R Z_t=0.
Pointwise, Z_t is nonnegative and is at least one on {X_t=1}.

The script uses separately centered copies for the three blocker cores.  It
keeps every even checkpoint clause, the exact three future danger bits, and
the charge count for selected odd clauses.  Cross-core current incidence is
relaxed to every coupling with the exact one-coordinate root marginals.
Every rational coupling LP is solved by complete primal vertex enumeration
and an independent exact basic dual certificate.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product


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
            Fraction(1, 7),
            Fraction(1, 7),
            Fraction(1, 7),
        )
    bits = BIT_VECTORS[key]
    return (
        Fraction(1),
        Fraction(2 - bits[0], P),
        Fraction(2 - bits[1], P),
        Fraction(2 - bits[2], P),
    )


@lru_cache(maxsize=None)
def coupling_vertices(key):
    """All distinct exact vertices, together with every producing basis."""
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
    """Exact primal vertex maximum and an independent exact dual."""
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
            # A_B^T y=c_B, so y=(A_B^{-1})^T c_B.
            dual = tuple(
                sum(
                    (
                        inverse[row][col] * basic_values[row]
                        for row in range(RANK)
                    ),
                    Fraction(),
                )
                for col in range(RANK)
            )
            if any(
                dot(column, dual) < value
                for column, value in zip(COLUMNS, values)
            ):
                continue
            require(dot(rhs, dual) == primal, "exact primal/dual value drift")
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


def charge_packet(profile, odd_checkpoints):
    """Return selected centered charges and their expanded coefficients."""
    require(odd_checkpoints, "at least one odd checkpoint is required")
    require(
        all(checkpoint % 2 == 1 for checkpoint in odd_checkpoints),
        "charge checkpoints must be odd",
    )
    charges_by_time = {}
    correction = [Fraction() for _ in range(3)]
    shift = Fraction()
    for blocker, valuation in enumerate(profile):
        for checkpoint in odd_checkpoints:
            time = valuation - checkpoint
            require(time >= 0, "negative charge time")
            rho_power = RHO**time
            charges_by_time.setdefault(time, []).append(blocker)
            correction[blocker] += rho_power
            shift += max(rho_power, Fraction())
            # Z_t=X_t-rho^t X_0+(rho^t)_+ is nonnegative and >=1
            # on X_t=1.  Check all four pointwise bit pairs exactly.
            for x_zero in (0, 1):
                for x_time in (0, 1):
                    charge = (
                        x_time
                        - rho_power * x_zero
                        + max(rho_power, Fraction())
                    )
                    require(charge >= 0, "centered charge became negative")
                    if x_time:
                        require(charge >= 1, "active centered charge below one")
    return (
        {time: tuple(rows) for time, rows in charges_by_time.items()},
        tuple(correction),
        shift,
    )


@lru_cache(maxsize=None)
def signed_capacity(profile, odd_checkpoints):
    """Robust capacity of the separately centered signed dual."""
    first = profile[0]
    require(2 <= first <= 5, "profile outside signed-dual scope")
    even_checkpoints = tuple(range(0, first + 1, 2))
    clause_count = len(even_checkpoints)
    all_clauses = (1 << clause_count) - 1
    m = len(odd_checkpoints)
    z_max = 3 * m

    even_by_time = {}
    for blocker, valuation in enumerate(profile):
        for clause, checkpoint in enumerate(even_checkpoints):
            time = valuation - checkpoint
            require(time >= 0, "negative even-checkpoint time")
            even_by_time.setdefault(time, []).append((blocker, clause))
    even_by_time = {
        time: tuple(rows) for time, rows in even_by_time.items()
    }
    charge_by_time, correction, shift = charge_packet(
        profile, odd_checkpoints
    )
    horizon = max(profile)

    def new_clause_mask(time, bits):
        mask = 0
        for blocker, clause in even_by_time.get(time, ()):
            if bits[blocker]:
                mask |= 1 << clause
        return mask

    def charge_increment(time, bits):
        return sum(bits[blocker] for blocker in charge_by_time.get(time, ()))

    @lru_cache(maxsize=None)
    def recurse(time, mask, z_count, next_key):
        if time < 0:
            require(0 <= z_count <= z_max, "terminal charge count drift")
            if mask != all_clauses:
                return Fraction()
            x_zero = BIT_VECTORS[next_key]
            one_minus_q = (
                Fraction(1 - m + z_count)
                - sum(
                    correction[j] * x_zero[j]
                    for j in range(3)
                )
                + shift
            )
            return max(one_minus_q, Fraction())

        costs = []
        for current_key, current_bits in enumerate(BIT_VECTORS):
            next_mask = mask | new_clause_mask(time, current_bits)
            next_count = z_count + charge_increment(time, current_bits)
            require(next_count <= z_max, "charge count overflow")
            costs.append(
                recurse(time - 1, next_mask, next_count, current_key)
            )
        return maximize_joint_bits(tuple(costs), next_key)

    terminal_values = tuple(
        recurse(horizon, 0, 0, terminal_key)
        for terminal_key in range(len(BIT_VECTORS))
    )
    return maximize_joint_bits(terminal_values, TERMINAL_KEY)


def first_one_residue():
    return {
        (1, middle, deepest)
        for deepest in range(5, 20)
        for middle in range(1, deepest)
    }


def first_two_residue():
    return (
        {(2, 3, 5)}
        | {(2, 4, deepest) for deepest in range(5, 20)}
        | {(2, middle, middle + 2) for middle in range(5, 18)}
    )


def first_three_residue():
    return (
        {(3, 3, 5), (3, 4, 5), (3, 4, 6)}
        | {(3, 5, deepest) for deepest in range(6, 20)}
        | {(3, middle, middle + 2) for middle in range(6, 18)}
    )


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


def common_core_control():
    """Independent exact Markov evaluation for the (4,6,8) dual."""
    transition = (
        (Fraction(11, 13), Fraction(2, 13)),
        (Fraction(12, 13), Fraction(1, 13)),
    )
    stationary = (Fraction(6, 7), Fraction(1, 7))
    profile = (4, 6, 8)
    even = (0, 2, 4)
    odd = (1, 3)
    capacity = Fraction()
    unsigned = Fraction()
    for word in product((0, 1), repeat=9):
        probability = stationary[word[0]]
        for left, right in zip(word, word[1:]):
            probability *= transition[left][right]
        p_event = all(
            any(word[valuation - checkpoint] for valuation in profile)
            for checkpoint in even
        )
        if not p_event:
            continue
        unsigned += probability
        z_count = sum(
            word[valuation - checkpoint]
            for valuation in profile
            for checkpoint in odd
        )
        epsilon = tuple(
            -sum(RHO ** (valuation - checkpoint) for checkpoint in odd)
            for valuation in profile
        )
        one_minus_q = (
            Fraction(z_count - 1)
            + sum(epsilon[j] * word[0] for j in range(3))
        )
        capacity += probability * max(one_minus_q, Fraction())
    return unsigned, capacity


def main():
    global LP_CALLS, DUAL_CERTIFICATES
    require(len(first_one_residue()) == 165, "first-depth-one census drift")
    low_two = first_two_residue()
    low_three = first_three_residue()
    require(len(low_two) == 29, "first-depth-two census drift")
    require(len(low_three) == 29, "first-depth-three census drift")

    post_thm2233 = (
        first_one_residue()
        | low_two
        | low_three
        | {(4, 6, 8)}
    )
    require(len(post_thm2233) == 224, "post-THM2233 ledger drift")

    low_bounds = {
        profile: signed_capacity(profile, (1,))
        for profile in sorted(low_two | low_three)
    }
    passing_two = {
        profile for profile in low_two if low_bounds[profile] < TARGET
    }
    passing_three = {
        profile for profile in low_three if low_bounds[profile] < TARGET
    }
    require(passing_two == low_two, "not every first-depth-two row passed")
    require(
        passing_three == low_three - {(3, 4, 5)},
        f"first-depth-three pass set drift: {sorted(low_three-passing_three)}",
    )

    worst_two = max(passing_two, key=lambda profile: low_bounds[profile])
    worst_three = max(passing_three, key=lambda profile: low_bounds[profile])
    require(worst_two == (2, 4, 5), "first-depth-two worst row drift")
    require(worst_three == (3, 5, 6), "first-depth-three worst row drift")
    require(
        low_bounds[worst_two] == Fraction(8_945_166_533, 74_231_495_611),
        "first-depth-two worst bound drift",
    )
    require(
        low_bounds[worst_three]
        == Fraction(1_471_733_046_268, 12_545_122_758_259),
        "first-depth-three worst passing bound drift",
    )
    hostile = low_bounds[(3, 4, 5)]
    require(
        hostile == Fraction(17_878_637_620, 74_231_495_611),
        "hostile (3,4,5) bound drift",
    )
    require(hostile > TARGET, "hostile (3,4,5) unexpectedly passed")

    high = signed_capacity((4, 6, 8), (1, 3))
    expected_high = Fraction(
        17_322_925_655_936_326,
        358_301_251_098_635_299,
    )
    require(high == expected_high, "high-row signed capacity drift")
    require(high < TARGET, "high-row signed capacity missed target")
    high_gap = TARGET - high
    require(
        high_gap
        == Fraction(
            32_039_946_787_164_254_737,
            354_718_238_587_648_946_010,
        ),
        "high-row target gap drift",
    )

    # The six high-row charge times are odd.  If z=0, the positive-part
    # terminal raw value is -1+sum epsilon_j; this must remain negative.
    high_epsilon = tuple(
        -sum(RHO ** (valuation - checkpoint) for checkpoint in (1, 3))
        for valuation in (4, 6, 8)
    )
    require(sum(high_epsilon) < 1, "high-row zero-charge boundary drift")
    require(
        high_epsilon
        == (
            Fraction(170, 2197),
            Fraction(170, 371293),
            Fraction(170, 62748517),
        ),
        "high-row epsilon packet drift",
    )

    common_unsigned, common_capacity = common_core_control()
    require(
        common_unsigned == Fraction(916_159, 4_826_809),
        "common-core unsigned control drift",
    )
    require(
        common_capacity
        == Fraction(2_303_649_491_556_761, 51_185_893_014_090_757),
        "common-core signed control drift",
    )
    require(common_capacity <= high, "robust bound missed common-core control")

    newly_closed = passing_two | passing_three | {(4, 6, 8)}
    require(len(newly_closed) == 58, "new closure count drift")
    remaining = post_thm2233 - newly_closed
    expected_remaining = first_one_residue() | {(3, 4, 5)}
    require(remaining == expected_remaining, "remaining ledger drift")
    require(len(remaining) == 166, "remaining ledger count drift")

    # Hostile repair control: merely scaling the even-time centered atom
    # Y_2 by its active minimum leaves an inactive negative value.
    naive_scaled_inactive = -RHO**2 / (1 - RHO**2)
    require(
        naive_scaled_inactive == Fraction(-1, 168),
        "naive even-time scaling witness drift",
    )
    require(naive_scaled_inactive < 0, "naive scaling witness lost its sign")
    require(RHO**0 == 1, "time-zero transfer phase drift")
    require(
        Fraction(0) - RHO**0 * Fraction(1) + max(RHO**0, Fraction()) == 0,
        "time-zero constant-charge cancellation drift",
    )

    expected_bound_digest = (
        "e2c4650604218f8409bbc384c6144a96e010895c9ca6947b72922233dec5c24f"
    )
    expected_closed_digest = (
        "d056889baf99ba2972cd39329807810a1ae9e291b63fcdf960b8bf9713580249"
    )
    expected_remaining_digest = (
        "fa7c60161dd2a1a0de77bd3a1133d030ab65d4633f881a0833f67cac0cc25d09"
    )
    bound_digest = digest_bounds(low_bounds | {(4, 6, 8): high})
    closed_digest = digest_profile_set(newly_closed)
    remaining_digest = digest_profile_set(remaining)
    require(bound_digest == expected_bound_digest, "bound digest drift")
    require(closed_digest == expected_closed_digest, "closed digest drift")
    require(
        remaining_digest == expected_remaining_digest,
        "remaining digest drift",
    )

    cache = maximize_joint_bits_uncounted.cache_info()
    require(
        DUAL_CERTIFICATES == cache.misses,
        "not every distinct coupling LP received one dual certificate",
    )
    require(LP_CALLS == 45_617, "coupling LP call count drift")
    require(cache.misses == 39_295, "distinct coupling LP count drift")

    print(f"target={TARGET}")
    print(f"invertible_bases={len(INVERTIBLE_BASES)}")
    print(
        "root_vertex_counts="
        + str(tuple(len(coupling_vertices(key)) for key in range(8)))
    )
    print(f"terminal_vertex_count={len(coupling_vertices(TERMINAL_KEY))}")
    print(
        f"lambda1=2 closed={len(passing_two)}/{len(low_two)} "
        f"worst={worst_two} bound={low_bounds[worst_two]} "
        f"target_minus={TARGET-low_bounds[worst_two]}"
    )
    print(
        f"lambda1=3 closed={len(passing_three)}/{len(low_three)} "
        f"worst_passing={worst_three} bound={low_bounds[worst_three]} "
        f"target_minus={TARGET-low_bounds[worst_three]}"
    )
    print(
        f"hostile_profile=(3,4,5) bound={hostile} "
        f"bound_minus_target={hostile-TARGET}"
    )
    print(
        f"high_profile=(4,6,8) bound={high} target_minus={high_gap} "
        f"epsilon={high_epsilon}"
    )
    print(
        f"common_core_unsigned={common_unsigned} "
        f"common_core_signed={common_capacity}"
    )
    print(
        f"post_thm2233={len(post_thm2233)} newly_closed={len(newly_closed)} "
        f"remaining={len(remaining)}"
    )
    print(
        "remaining_structure=165 lambda1=1 profiles with lambda3>=5; "
        "(3,4,5)"
    )
    print(f"bound_digest={bound_digest}")
    print(f"closed_digest={closed_digest}")
    print(f"remaining_digest={remaining_digest}")
    print(f"naive_scaled_inactive_witness={naive_scaled_inactive}")
    print(
        f"lp_calls={LP_CALLS} distinct_lps={cache.misses} "
        f"dual_certificates={DUAL_CERTIFICATES}"
    )
    print("exact_primal_dual_vertex_parity=PASS")


if __name__ == "__main__":
    main()
