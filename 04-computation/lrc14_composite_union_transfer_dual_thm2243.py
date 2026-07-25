#!/usr/bin/env python3
"""Exact composite-union transfer-dual referee for THM-2243.

For a depth-three valuation profile, all four checkpoint unions are shifts
of one Boolean function U.  With rho=-1/13 and r=rho^2=1/169, the score

    q = 1 - (U_1-r U_3)/(1-r)

is R-orthogonal.  It vanishes on the negative carrier U_1=U_3=1, while

    (1-q)_+ = U_1 (169-U_3)/168.

The script prices the resulting positive-carrier payoff using the audited
four-bit arbitrary-coupling engine from THM-2233.  All arithmetic after the
finite coupling vertices have been constructed is exact Fraction arithmetic.
Every distinct local maximization receives that engine's exact primal/dual
basic certificate.
"""

from fractions import Fraction
from hashlib import sha256

import lrc14_guard_danger_hidden_state_bellman_thm2233 as base


TARGET = Fraction(961, 6930)
RHO = Fraction(-1, 13)
R = RHO**2
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


def post_thm2233_residue():
    """Reconstruct the exact 224-profile residue of THM-2233."""
    return (
        first_one_residue()
        | first_two_residue()
        | first_three_residue()
        | {(4, 6, 8)}
    )


def build_hit_masks(profile):
    """Return exact guard/U_0/U_1/U_2/U_3 incidence masks."""
    require(profile[0] == 3, "THM-2243 only prices depth three")
    horizon = profile[2]
    hit_masks = [
        [0 for _ in base.BIT_VECTORS] for _ in range(horizon + 1)
    ]

    # Clause zero is the guard at time zero.
    for state, bits in enumerate(base.BIT_VECTORS):
        if bits[0]:
            hit_masks[0][state] |= 1

    # Clause k+1 is U_k.  The incidence time lambda_j-k is exactly the
    # orbit identity U_k=U o S^(3-k), resolved into the three core bits.
    for checkpoint in range(4):
        clause = checkpoint + 1
        for blocker, valuation in enumerate(profile):
            time = valuation - checkpoint
            require(time >= 0, "negative checkpoint time")
            for state, bits in enumerate(base.BIT_VECTORS):
                if bits[blocker + 1]:
                    hit_masks[time][state] |= 1 << clause

    return horizon, tuple(tuple(row) for row in hit_masks)


def terminal_payoff(clauses):
    """The exact sign-safe payoff 1_P (1-q)_+."""
    guard = bool(clauses & 1)
    union_bits = tuple(
        int(bool(clauses & (1 << (checkpoint + 1))))
        for checkpoint in range(4)
    )
    u0, u1, u2, u3 = union_bits
    q_value = Fraction(1) - (Fraction(u1) - R * Fraction(u3)) / (1 - R)

    # Pointwise truth-table audits for the nonlinear centered score.
    if u1 and u3:
        require(q_value == 0, "q does not vanish on the negative carrier")
    expected_hinge = Fraction(u1 * (169 - u3), 168)
    require(
        max(Fraction(1) - q_value, ZERO) == expected_hinge,
        "composite-union hinge identity drift",
    )

    return expected_hinge if guard and u0 and u2 else ZERO


def profile_bound(profile):
    """Exact robust Bellman capacity for one depth-three profile."""
    horizon, hit_masks = build_hit_masks(profile)

    # The row indexed by a four-bit state is initially constant because the
    # terminal payoff depends only on the accumulated clause mask.
    value = tuple(
        tuple(terminal_payoff(clauses) for _ in base.BIT_VECTORS)
        for clauses in range(1 << 5)
    )

    for time in range(horizon + 1):
        updated = []
        for later_clauses in range(1 << 5):
            row = []
            for future_state in range(len(base.BIT_VECTORS)):
                costs = tuple(
                    value[
                        later_clauses | hit_masks[time][current_state]
                    ][current_state]
                    for current_state in range(len(base.BIT_VECTORS))
                )
                row.append(base.maximize_joint_bits(costs, future_state))
            updated.append(tuple(row))
        value = tuple(updated)

    answer = base.maximize_joint_bits(value[0], base.TERMINAL_KEY)
    state_count = (horizon + 1) * (1 << 5) * len(base.BIT_VECTORS)
    return answer, state_count


def main():
    residue = post_thm2233_residue()
    depth_three = first_three_residue()
    hostile_profile = (3, 4, 5)
    claimed = depth_three - {hostile_profile}

    require(len(residue) == 224, "post-THM-2233 residue count drift")
    require(len(depth_three) == 29, "depth-three residue count drift")
    require(len(claimed) == 28, "claimed closure count drift")
    require(
        digest_profiles(residue)
        == "d7617205a31504563f5f3d14eebc7d6dfd6148349f43a468949a8ecd850f8505",
        "post-THM-2233 residue digest drift",
    )
    require(
        digest_profiles(depth_three)
        == "d1cd2ef956e3b929735390fc3f92c262a1fab945c65662eb0594c666fcfdd1be",
        "depth-three residue digest drift",
    )
    require(
        digest_profiles(claimed)
        == "057e60270d63c8e192adfb14c96f70c50b4f20fc68f216a41e8472be5fe1f9ea",
        "claimed closure digest drift",
    )

    # The identity U_k=U o S^(3-k) is reflected by these exact shifts.
    for profile in depth_three:
        for checkpoint in range(4):
            for valuation in profile:
                require(
                    valuation - checkpoint
                    == (valuation - 3) + (3 - checkpoint),
                    "composite orbit time identity drift",
                )

    # Audit the complete U_0,...,U_3 truth table before any Bellman call.
    for word in range(1 << 4):
        clauses = sum(
            ((word >> checkpoint) & 1) << (checkpoint + 1)
            for checkpoint in range(4)
        )
        terminal_payoff(clauses)
        terminal_payoff(clauses | 1)

    bounds = {}
    state_count = 0
    for profile in sorted(depth_three):
        bound, states = profile_bound(profile)
        bounds[profile] = bound
        state_count += states

    for profile in claimed:
        require(
            bounds[profile] < TARGET,
            f"claimed composite-union profile does not close: {profile}",
        )
    require(
        bounds[hostile_profile] > TARGET,
        "hostile profile unexpectedly closes",
    )

    worst_closed_profile = max(claimed, key=lambda profile: bounds[profile])
    worst_closed_bound = bounds[worst_closed_profile]
    require(
        worst_closed_profile == (3, 4, 6),
        "worst closed profile drift",
    )
    require(
        worst_closed_bound == Fraction(18956451, 270301304),
        "worst closed bound drift",
    )
    hostile_bound = bounds[hostile_profile]
    require(
        hostile_bound == Fraction(8907541, 62377224),
        "hostile bound drift",
    )
    require(
        digest_bounds(bounds)
        == "d0d86f8b91dcfb6fead5f83da4c81f416cec6d189c9e119065dfe7d60fe3bb6e",
        "depth-three bound digest drift",
    )

    remaining = residue - claimed
    expected_remaining = (
        first_one_residue()
        | first_two_residue()
        | {hostile_profile, (4, 6, 8)}
    )
    require(len(remaining) == 196, "standalone THM-2243 residue count drift")
    require(remaining == expected_remaining, "remaining profile shape drift")
    require(
        digest_profiles(remaining)
        == "d33487ee4439a501a2492a4bfe477cf7a7112eeaf0598e58f8271294347525a1",
        "remaining profile digest drift",
    )
    canonical_remaining = first_one_residue() | {hostile_profile}
    require(len(canonical_remaining) == 166, "canonical residue count drift")
    require(
        digest_profiles(canonical_remaining)
        == "fa7c60161dd2a1a0de77bd3a1133d030ab65d4633f881a0833f67cac0cc25d09",
        "canonical residue digest drift",
    )

    require(
        base.DUAL_CERTIFICATES
        == base.maximize_joint_bits_uncounted.cache_info().misses,
        "not every distinct LP received an exact dual certificate",
    )
    require(
        base.LP_CALLS == state_count + len(depth_three),
        "Bellman LP-call count drift",
    )

    print("THM-2243 COMPOSITE-UNION TRANSFER-DUAL DEPTH-THREE EXCLUSION")
    print(f"rho={RHO} r=rho^2={R}")
    print("orbit_identity=U_k=U_o_S^(3-k):PASS")
    print("negative_carrier_q=POINTWISE_ZERO")
    print("hinge_identity=(1-q)_+=U1*(169-U3)/168:PASS")
    print(f"target={TARGET}")
    print(f"post_thm2233_residue={len(residue)}")
    print(f"depth_three_profiles={len(depth_three)}")
    print(f"composite_union_closed={len(claimed)}")
    print(f"hostile_profile={hostile_profile}")
    print(
        f"hostile_bound={hostile_bound} "
        f"decimal={float(hostile_bound):.15f}"
    )
    print(f"hostile_excess={hostile_bound - TARGET}")
    print(f"worst_closed_profile={worst_closed_profile}")
    print(
        f"worst_closed_bound={worst_closed_bound} "
        f"decimal={float(worst_closed_bound):.15f}"
    )
    print(f"worst_closed_margin={TARGET - worst_closed_bound}")
    print(f"standalone_remaining_profiles={len(remaining)}")
    print(
        "standalone_remaining_shape=lambda1=1:165; lambda1=2:29; "
        "exceptional=(3,4,5),(4,6,8)"
    )
    print(f"canonical_remaining_with_thm2239={len(canonical_remaining)}")
    print("canonical_remaining_shape=lambda1=1:165; exceptional=(3,4,5):1")
    print(f"post_thm2233_residue_digest={digest_profiles(residue)}")
    print(f"depth_three_digest={digest_profiles(depth_three)}")
    print(f"closed_digest={digest_profiles(claimed)}")
    print(f"standalone_remaining_digest={digest_profiles(remaining)}")
    print(f"canonical_remaining_digest={digest_profiles(canonical_remaining)}")
    print(f"bound_digest={digest_bounds(bounds)}")
    print(f"bellman_states={state_count}")
    print(f"lp_calls={base.LP_CALLS}")
    print(
        "distinct_exact_lp_calls="
        f"{base.maximize_joint_bits_uncounted.cache_info().misses}"
    )
    print(f"exact_dual_certificates={base.DUAL_CERTIFICATES}")
    print("exact_primal_dual_vertex_parity=PASS")
    print("status=THM2243_PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
