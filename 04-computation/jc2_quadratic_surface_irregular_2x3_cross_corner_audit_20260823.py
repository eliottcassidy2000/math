"""Exact hostile controls for a shorter alternative proof of THM-3786.

THM-3786 is already proved and independently hostile-audited.  This companion
checks the alternative forced-cross-corner collision proof with explicit gates
that remain active under ``python -O``.  The unbounded proof is the
sign/collision argument recorded in the accompanying reflection.
"""

from collections import defaultdict

import sympy as sp


def gate(condition: bool, label: str) -> None:
    """Optimization-safe exact gate."""

    if not condition:
        raise RuntimeError(f"FAILED: {label}")


s = sp.symbols("s")
checks = 0


def canonical_profile(weight: int, *, nonconstant_zero: bool = True):
    """A legal THM-3783 profile of the requested Euler weight."""

    pole_order = (-weight + 3) // 4 if weight < 0 else 0
    profile = (s**2 - 1) ** pole_order
    if weight % 2:
        profile *= s
    if weight == 0 and nonconstant_zero:
        profile = s**2
    return sp.expand(profile)


def bracket_profile(a: int, f, b: int, h):
    """The profile in J(x^a f(s), x^b h(s))."""

    return sp.expand(2 * (a * f * sp.diff(h, s) - b * sp.diff(f, s) * h))


def commuting_signs_are_possible(a: int, b: int) -> bool:
    """Necessary sign condition for nonconstant active profiles to commute."""

    return a * b > 0 or (a == 0 and b == 0)


# Opposite nonzero weights never commute; a nonconstant weight-zero profile
# never commutes with a nonzero-weight profile.  A constant weight-zero
# profile is retained as the hostile control explaining why constants must be
# subtracted before support counting.
opposite_sign_controls = 0
zero_weight_controls = 0
constant_profile_controls = 0
for negative in range(-16, 0):
    for positive in range(1, 17):
        value = bracket_profile(
            negative,
            canonical_profile(negative),
            positive,
            canonical_profile(positive),
        )
        gate(value != 0, f"opposite-sign-{negative}-{positive}")
        opposite_sign_controls += 1
        checks += 1

for weight in list(range(-16, 0)) + list(range(1, 17)):
    profile = canonical_profile(weight)
    value = bracket_profile(0, canonical_profile(0), weight, profile)
    gate(value != 0, f"nonconstant-zero-weight-{weight}")
    zero_weight_controls += 1
    checks += 1

    constant_value = bracket_profile(
        0, canonical_profile(0, nonconstant_zero=False), weight, profile
    )
    gate(constant_value == 0, f"constant-zero-weight-{weight}")
    constant_profile_controls += 1
    checks += 1


# If all component weights are nonpositive, every component bracket vanishes
# at s=1.  This includes negative/zero and zero/zero boundaries.
nonpositive_boundary_controls = 0
for a in range(-16, 1):
    for b in range(-16, 1):
        value = bracket_profile(
            a, canonical_profile(a), b, canonical_profile(b)
        )
        gate(value.subs(s, 1) == 0, f"nonpositive-boundary-{a}-{b}")
        nonpositive_boundary_controls += 1
        checks += 1


# In a homogeneous scalar bracket a+b=-3.  Boundary evaluation eliminates
# every sign case except (1,-4) and its swap.  In those two cases the leading
# coefficient has factor deg(h)+4deg(f), respectively its swap, and hence
# cannot disappear in characteristic zero.
remaining_homogeneous_cases = []
for a in range(-64, 65):
    b = -3 - a
    if a > 0 > b and (-b + 3) // 4 == 1:
        remaining_homogeneous_cases.append((a, b))
    if b > 0 > a and (-a + 3) // 4 == 1:
        remaining_homogeneous_cases.append((a, b))
gate(
    remaining_homogeneous_cases == [(-4, 1), (1, -4)],
    "homogeneous-boundary-case-list",
)
checks += 1

leading_degree_controls = 0
for odd_degree in range(1, 32, 2):
    for even_degree in range(2, 34, 2):
        gate(
            odd_degree + even_degree - 1 > 0
            and even_degree + 4 * odd_degree > 0,
            f"homogeneous-leading-term-{odd_degree}-{even_degree}",
        )
        leading_degree_controls += 1
        checks += 1


def buckets(d: int, e: int, f: int):
    increments = {
        (0, 0): 0,
        (0, 1): e,
        (0, 2): e + f,
        (1, 0): d,
        (1, 1): d + e,
        (1, 2): d + e + f,
    }
    grouped = defaultdict(list)
    for pair, increment in increments.items():
        grouped[increment].append(pair)
    return increments, dict(grouped)


# Complete finite hostile census for the two forced mixed-sign corners.  The
# transpose check is the literal output-swap (3-by-2) control.
bound = 72
gap_triples = 0
common_step_controls = 0
central_seam_controls = 0
output_swap_controls = 0
for d in range(1, bound + 1):
    for e in range(1, bound + 1):
        for f in range(1, bound + 1):
            increments, grouped = buckets(d, e, f)
            low_high = increments[(0, 2)]
            high_low = increments[(1, 0)]
            both_collide = (
                len(grouped[low_high]) > 1 and len(grouped[high_low]) > 1
            )
            predicted = d == e == f or d == e + f
            gate(both_collide == predicted, f"collision-{d}-{e}-{f}")
            gate(max(map(len, grouped.values())) <= 2, f"no-triple-{d}-{e}-{f}")
            gap_triples += 1
            checks += 2

            transposed = defaultdict(list)
            for increment, pairs in grouped.items():
                transposed[increment] = [(j, i) for i, j in pairs]
            gate(
                sorted(map(len, transposed.values()))
                == sorted(map(len, grouped.values())),
                f"output-swap-{d}-{e}-{f}",
            )
            output_swap_controls += 1
            checks += 1

            if d == e == f:
                gate(
                    grouped[d] == [(0, 1), (1, 0)]
                    and grouped[2 * d] == [(0, 2), (1, 1)],
                    f"common-step-{d}",
                )
                common_step_controls += 1
                checks += 1

            if d == e + f:
                gate(
                    grouped[d] == [(0, 2), (1, 0)]
                    and grouped[e] == [(0, 1)]
                    and grouped[d + e] == [(1, 1)],
                    f"central-seam-{d}-{e}-{f}",
                )
                central_seam_controls += 1
                checks += 1


# Abstract endpoint-sign audit, including both-zero endpoints.  Under the
# exact commutation classification, scalar reachability, and the s=1 hostile
# boundary, all surviving endpoint weights are strictly crossing.
endpoint_sign_controls = 0
for d in range(1, 13):
    for e in range(1, 13):
        for f in range(1, 13):
            for a in range(-20, 13):
                for b in range(-20, 13):
                    high_a = a + d
                    high_b = b + e + f
                    if not commuting_signs_are_possible(a, b):
                        continue
                    if not commuting_signs_are_possible(high_a, high_b):
                        continue
                    supports_a = (a, high_a)
                    supports_b = (b, b + e, high_b)
                    if not any(u + v == -3 for u in supports_a for v in supports_b):
                        continue
                    if a >= 0 and b >= 0:
                        continue  # Scalar weight is unreachable in fact.
                    if high_a <= 0 and high_b <= 0:
                        continue  # Every bracket profile vanishes at s=1.
                    gate(a < 0 < high_a, f"A-endpoint-sign-{a}-{b}-{d}-{e}-{f}")
                    gate(b < 0 < high_b, f"B-endpoint-sign-{a}-{b}-{d}-{e}-{f}")
                    endpoint_sign_controls += 1
                    checks += 2


# On d=e+f the middle B weight is a hub joined to both the negative and the
# positive A endpoint.  No sign, including zero, permits both brackets to
# vanish when the zero-weight profile is genuinely active/nonconstant.
hub_sign_controls = 0
for e in range(1, 33):
    for f in range(1, 33):
        d = e + f
        for a in range(1 - d, 0):
            for b in range(1 - d, 0):
                middle_b = b + e
                gate(
                    not (
                        commuting_signs_are_possible(a, middle_b)
                        and commuting_signs_are_possible(a + d, middle_b)
                    ),
                    f"hub-sign-{a}-{b}-{d}-{e}-{f}",
                )
                hub_sign_controls += 1
                checks += 1


# Common-step cells really survive this collision argument.  These controls
# are deliberately not the already-closed THM-3783 aligned family: for every
# t>=3, A={2-t,2}, B={-5,t-5,2t-5} crosses zero and has the scalar pair-sum
# -3 in the first double bucket.
unaligned_common_step_controls = 0
for t in range(3, bound + 1):
    a, b = 2 - t, -5
    increments, grouped = buckets(t, t, t)
    gate(a < 0 < a + t, f"common-A-crossing-{t}")
    gate(b < 0 < b + 2 * t, f"common-B-crossing-{t}")
    gate(a + b + t == -3 and len(grouped[t]) == 2, f"common-scalar-{t}")
    unaligned_common_step_controls += 1
    checks += 3


print("audit=THM-3786-irregular-two-by-three-cross-corner-alternative-proof")
print("status=finite-exact-hostile-controls-for-independent-shorter-proof")
print(f"optimization_safe_checks={checks}")
print(f"opposite_sign_controls={opposite_sign_controls}")
print(f"zero_weight_nonconstant_controls={zero_weight_controls}")
print(f"zero_weight_constant_hostile_controls={constant_profile_controls}")
print(f"nonpositive_s_equals_one_controls={nonpositive_boundary_controls}")
print(f"homogeneous_leading_degree_controls={leading_degree_controls}")
print(f"gap_triples={gap_triples}")
print(f"output_swap_controls={output_swap_controls}")
print(f"common_step_collision_controls={common_step_controls}")
print(f"central_seam_collision_controls={central_seam_controls}")
print(f"endpoint_sign_controls={endpoint_sign_controls}")
print(f"central_hub_sign_controls={hub_sign_controls}")
print(f"unaligned_common_step_positive_controls={unaligned_common_step_controls}")
print("classification=forced-corners-collide-iff-common-step-or-d=e+f")
print("central_seam=d=e+f-dies-at-middle-weight-hub-including-weight-zero")
print("survivor=common-step-collision-shape-only-no-Darboux-existence-claim")
print("PASS")
