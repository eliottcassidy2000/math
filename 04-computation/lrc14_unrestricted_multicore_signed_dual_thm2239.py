#!/usr/bin/env python3
"""Exact signed-score Bellman referee for THM-2239.

The companion imports the already audited four-bit coupling engine of
THM-2233.  It adds the signed residual score

    q = number of odd checkpoints
        - sum_(odd k,j) (X_(j,lambda_j-k)
                         - (-1/13)^(lambda_j-k) X_(j,0))

and evaluates the pointwise majorant

    1_P (1-q)_+ + 5 1_N q_+.

Here P is the guard plus every even checkpoint carrier and N is every odd
checkpoint carrier.  All arithmetic after the finite coupling vertices have
been constructed is exact Fraction arithmetic.  Every distinct local
coupling maximization receives the exact primal/dual certificate supplied by
the THM-2233 engine.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256

import lrc14_guard_danger_hidden_state_bellman_thm2233 as base


RHO = Fraction(-1, 13)
ZERO = Fraction()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest_profiles(profiles):
    payload = ";".join(
        ",".join(map(str, profile)) for profile in sorted(profiles)
    )
    return sha256(payload.encode("ascii")).hexdigest()


def digest_bounds(rows):
    payload = ";".join(
        f"{','.join(map(str, profile))}={bound}"
        for profile, bound in sorted(rows.items())
    )
    return sha256(payload.encode("ascii")).hexdigest()


def current_residue():
    """The exact 224-profile residue produced by THM-2233."""
    universe = set(base.all_profiles())
    return (
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


def signed_closed_profiles():
    """Profiles claimed closed by the canonical odd-checkpoint dual."""
    return (
        {(2, 3, 5)}
        | {(2, 4, deepest) for deepest in range(5, 20)}
        | {(2, middle, middle + 2) for middle in range(5, 18)}
        | {(4, 6, 8)}
    )


def build_profile_data(profile):
    """Return clause masks and exact q-score increments for one profile."""
    first = profile[0]
    even_checkpoints = tuple(range(0, first + 1, 2))
    odd_checkpoints = tuple(range(1, first + 1, 2))
    positive_clause_count = 1 + len(even_checkpoints)
    odd_clause_count = len(odd_checkpoints)
    horizon = profile[2]

    hit_masks = [
        [0 for _ in base.BIT_VECTORS] for _ in range(horizon + 1)
    ]
    score_increments = [
        [ZERO for _ in base.BIT_VECTORS] for _ in range(horizon + 1)
    ]

    # Positive clause zero is the guard event C_H at time zero.
    for state, bits in enumerate(base.BIT_VECTORS):
        if bits[0]:
            hit_masks[0][state] |= 1

    # The remaining positive clauses are the even checkpoint unions.
    for clause, checkpoint in enumerate(even_checkpoints, start=1):
        for blocker, valuation in enumerate(profile):
            time = valuation - checkpoint
            require(time >= 0, "negative even-checkpoint time")
            for state, bits in enumerate(base.BIT_VECTORS):
                if bits[blocker + 1]:
                    hit_masks[time][state] |= 1 << clause

    # Negative clauses are the odd checkpoint unions.  At the same time,
    # construct q by subtracting each centered odd-checkpoint atom.
    for offset, checkpoint in enumerate(odd_checkpoints):
        clause = positive_clause_count + offset
        for blocker, valuation in enumerate(profile):
            time = valuation - checkpoint
            require(time >= 0, "negative odd-checkpoint time")
            correction = RHO**time
            for state, bits in enumerate(base.BIT_VECTORS):
                bit = bits[blocker + 1]
                if bit:
                    hit_masks[time][state] |= 1 << clause
                    score_increments[time][state] -= 1
                    score_increments[0][state] += correction

    positive_mask = (1 << positive_clause_count) - 1
    negative_mask = (
        ((1 << odd_clause_count) - 1) << positive_clause_count
    )
    return (
        horizon,
        odd_clause_count,
        positive_mask,
        negative_mask,
        tuple(tuple(row) for row in hit_masks),
        tuple(tuple(row) for row in score_increments),
    )


def profile_bound(profile):
    """Exact robust signed Bellman capacity for one valuation profile."""
    (
        horizon,
        odd_clause_count,
        positive_mask,
        negative_mask,
        hit_masks,
        score_increments,
    ) = build_profile_data(profile)

    @lru_cache(maxsize=None)
    def value(time, clauses, score, future_state):
        if time < 0:
            q_value = Fraction(odd_clause_count) + score
            positive = ZERO
            negative = ZERO
            if clauses & positive_mask == positive_mask:
                positive = max(Fraction(1) - q_value, ZERO)
            if clauses & negative_mask == negative_mask:
                negative = 5 * max(q_value, ZERO)
            return positive + negative

        costs = tuple(
            value(
                time - 1,
                clauses | hit_masks[time][current_state],
                score + score_increments[time][current_state],
                current_state,
            )
            for current_state in range(len(base.BIT_VECTORS))
        )
        return base.maximize_joint_bits(costs, future_state)

    terminal_costs = tuple(
        value(horizon, 0, ZERO, future_state)
        for future_state in range(len(base.BIT_VECTORS))
    )
    answer = base.maximize_joint_bits(terminal_costs, base.TERMINAL_KEY)
    return answer, value.cache_info().currsize


def main():
    residue = current_residue()
    claimed = signed_closed_profiles()
    require(len(residue) == 224, "THM-2233 residue count drift")
    require(len(claimed) == 30, "signed closure count drift")
    require(claimed <= residue, "claimed profile outside THM-2233 residue")
    require(
        digest_profiles(residue)
        == "d7617205a31504563f5f3d14eebc7d6dfd6148349f43a468949a8ecd850f8505",
        "THM-2233 residue digest drift",
    )
    require(
        digest_profiles(claimed)
        == "d7330636404abfd41aea31c08e72e0edbdf96b25f46d2909394d842d79d99859",
        "signed closure digest drift",
    )

    bounds = {}
    total_states = 0
    for profile in sorted(claimed):
        bound, state_count = profile_bound(profile)
        bounds[profile] = bound
        total_states += state_count
        require(
            bound < base.TARGET,
            f"claimed signed profile does not close: {profile}",
        )

    worst_profile = max(bounds, key=lambda profile: bounds[profile])
    worst_bound = bounds[worst_profile]
    require(worst_profile == (2, 3, 5), "worst closed profile drift")
    require(
        worst_bound == Fraction(7_265_183_507, 74_231_495_611),
        "worst closed bound drift",
    )

    high_profile = (4, 6, 8)
    high_bound = bounds[high_profile]
    require(
        high_bound
        == Fraction(13_063_069_772_240_011, 358_301_251_098_635_299),
        "high-first signed capacity drift",
    )

    # In the high-first row all six odd-checkpoint times are odd.  Hence
    # every centering correction is nonpositive and q<=0 on the full odd
    # carrier.  This freezes the pointwise sign mechanism independently of
    # the numerical Bellman value.
    high_odd_times = tuple(
        valuation - checkpoint
        for checkpoint in (1, 3)
        for valuation in high_profile
    )
    require(
        high_odd_times == (3, 5, 7, 1, 3, 5),
        "high-first odd schedule drift",
    )
    require(
        all(time % 2 == 1 for time in high_odd_times),
        "high-first correction lost its nonpositive sign",
    )

    remaining = residue - claimed
    require(len(remaining) == 194, "post-signed residue count drift")
    require(
        digest_profiles(remaining)
        == "4b939be77a2e6aa0fff0ac367269677afb693bd2707a381a3ef82b9e284a16dc",
        "post-signed residue digest drift",
    )
    expected_remaining = (
        {
            profile
            for profile in set(base.all_profiles())
            if profile[0] == 1 and profile[2] >= 5
        }
        | {(3, 3, 5), (3, 4, 5), (3, 4, 6)}
        | {(3, 5, deepest) for deepest in range(6, 20)}
        | {(3, middle, middle + 2) for middle in range(6, 18)}
    )
    require(remaining == expected_remaining, "remaining profile shape drift")

    # A surviving control records that the canonical all-ones dual is not a
    # universal scalar closure.  No optimality among failed rows is claimed.
    hostile_profile = (3, 17, 19)
    hostile_bound, hostile_states = profile_bound(hostile_profile)
    total_states += hostile_states
    require(
        hostile_bound
        == Fraction(
            229_985_997_674_091_853_361_506_127_343_787_104_565_593,
            1_150_805_888_298_461_593_768_523_506_367_492_533_368_331,
        ),
        "hostile surviving bound drift",
    )
    require(
        hostile_bound > base.TARGET,
        "hostile canonical-dual control unexpectedly closed",
    )

    bound_digest = digest_bounds(bounds)
    require(
        bound_digest
        == "c535642cce0db21aaab42a3eb6a733065bbb797ff10d531f9a372decc98c7621",
        "closed-bound digest drift",
    )
    require(total_states == 61_748, "Bellman state count drift")
    require(base.LP_CALLS == 54_483, "local LP call count drift")
    require(
        base.maximize_joint_bits_uncounted.cache_info().misses == 44_631,
        "distinct local LP count drift",
    )

    print("THM-2239 UNRESTRICTED MULTICORE SIGNED-DUAL PROFILE EXCLUSION")
    print(f"target={base.TARGET}")
    print(f"input_residue={len(residue)}")
    print(f"signed_closed={len(claimed)}")
    print(
        "signed_closed_shape=(2,3,5); "
        "(2,4,c),5<=c<=19; "
        "(2,b,b+2),5<=b<=17; "
        "(4,6,8)"
    )
    print(f"worst_closed_profile={worst_profile}")
    print(
        f"worst_closed_bound={worst_bound} "
        f"decimal={float(worst_bound):.15f}"
    )
    print(f"worst_closed_margin={base.TARGET - worst_bound}")
    print(
        f"high_first_profile={high_profile} "
        f"bound={high_bound} decimal={float(high_bound):.15f}"
    )
    print(f"high_first_margin={base.TARGET - high_bound}")
    print("high_first_negative_term=POINTWISE_ZERO")
    print(f"remaining_profiles={len(remaining)}")
    print("remaining_shape=lambda1=1:165; lambda1=3:29")
    print(f"input_residue_digest={digest_profiles(residue)}")
    print(f"signed_closed_digest={digest_profiles(claimed)}")
    print(f"remaining_digest={digest_profiles(remaining)}")
    print(f"closed_bound_digest={bound_digest}")
    print(
        f"hostile_profile={hostile_profile} "
        f"bound={hostile_bound} decimal={float(hostile_bound):.15f}"
    )
    print(f"hostile_excess={hostile_bound - base.TARGET}")
    print(f"bellman_states={total_states}")
    print(f"lp_calls={base.LP_CALLS}")
    print(
        "distinct_exact_lp_calls="
        f"{base.maximize_joint_bits_uncounted.cache_info().misses}"
    )
    print(f"exact_dual_certificates={base.DUAL_CERTIFICATES}")
    require(
        base.DUAL_CERTIFICATES
        == base.maximize_joint_bits_uncounted.cache_info().misses,
        "not every distinct local LP received an exact dual certificate",
    )
    print("exact_primal_dual_vertex_parity=PASS")
    print("status=THM2239_PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
