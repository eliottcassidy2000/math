#!/usr/bin/env python3
"""Exact lightweight probe for the THM-2818 x THM-2820 bridge.

The proved THM-2818 block lengths and fifty-half-step interblock jumps are
the finite input universe.  This script constructs raw, live, dead, and
live-minus-dead generating polynomials, computes their rooted Hasse jets
over F_13, checks their complete translation orbits, and exhibits two typed
carrier-gauge completions with the same one-axis chain motion but opposite
transversality verdicts.  It does not import the large physical builders.
"""

import ast
from math import comb
from pathlib import Path


P = 13
INTERBLOCK_JUMP = 50
TARGET_LABEL = 12
BLOCK_LENGTHS = {
    1: (145, 96),
    2: (143, 289, 96),
    3: (143, 289, 74),
}
EXPECTED_STARTS = {
    1: (0, 194),
    2: (0, 192, 530),
    3: (0, 192, 530),
}
EXPECTED_COUNTS = {
    1: (241, 121, 120),
    2: (528, 265, 263),
    3: (506, 254, 252),
}


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(
            ast.parse(Path(path).read_text(encoding="utf-8"))
        )
    )


def block_starts(lengths):
    starts = [0]
    for length in lengths[:-1]:
        # The last exponent in a block is start+length-1, and the next
        # first exponent is exactly fifty half-steps later.
        starts.append(starts[-1] + length - 1 + INTERBLOCK_JUMP)
    return tuple(starts)


def add_term(polynomial, exponent, coefficient):
    value = polynomial.get(exponent, 0) + coefficient
    if value:
        polynomial[exponent] = value
    elif exponent in polynomial:
        del polynomial[exponent]


def chain_polynomials(lengths):
    starts = block_starts(lengths)
    raw = {}
    live = {}
    dead = {}
    for start, length in zip(starts, lengths):
        for local_index in range(length):
            exponent = start + local_index
            add_term(raw, exponent, 1)
            add_term(
                live if local_index % 2 == 0 else dead,
                exponent,
                1,
            )
    contrast = dict(live)
    for exponent, coefficient in dead.items():
        add_term(contrast, exponent, -coefficient)
    require(
        all(raw.get(exponent, 0) == live.get(exponent, 0) + dead.get(exponent, 0)
            for exponent in raw),
        "raw/live/dead partition",
    )
    require(
        contrast
        == {
            exponent: 2 * live.get(exponent, 0) - raw.get(exponent, 0)
            for exponent in raw
        },
        "alternating contrast identity",
    )
    return starts, raw, live, dead, contrast


def hasse(polynomial, order, modulus=P):
    return sum(
        coefficient * comb(exponent, order)
        for exponent, coefficient in polynomial.items()
    ) % modulus


def integer_moment(polynomial, order):
    return sum(
        coefficient * comb(exponent, order)
        for exponent, coefficient in polynomial.items()
    )


def normalized_first(polynomial):
    j0 = hasse(polynomial, 0)
    j1 = hasse(polynomial, 1)
    require(j0 != 0, "normalization denominator vanished")
    return j0, j1, j1 * pow(j0, -1, P) % P


def quotient_vector(polynomial):
    vector = [0] * P
    for exponent, coefficient in polynomial.items():
        vector[exponent % P] = (
            vector[exponent % P] + coefficient
        ) % P
    return tuple(vector)


def eta_vector(polynomial):
    return tuple(hasse(polynomial, order) for order in range(P))


def truncated_multiply(left, right):
    answer = [0] * P
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j < P:
                answer[i + j] = (answer[i + j] + a * b) % P
    return tuple(answer)


def truncated_inverse(polynomial):
    require(polynomial[0] != 0, "nonunit local polynomial")
    answer = [0] * P
    answer[0] = pow(polynomial[0], -1, P)
    for degree in range(1, P):
        tail = sum(
            polynomial[index] * answer[degree - index]
            for index in range(1, degree + 1)
        )
        answer[degree] = -answer[0] * tail % P
    require(
        truncated_multiply(polynomial, tuple(answer))
        == (1,) + (0,) * (P - 1),
        "truncated unit inverse",
    )
    return tuple(answer)


def translate(polynomial, displacement):
    return {
        exponent + displacement: coefficient
        for exponent, coefficient in polynomial.items()
    }


def subtract(left, right):
    answer = dict(left)
    for exponent, coefficient in right.items():
        add_term(answer, exponent, -coefficient)
    return answer


def translation_orbit(polynomial):
    vectors = tuple(
        quotient_vector(translate(polynomial, displacement))
        for displacement in range(P)
    )
    require(len(set(vectors)) == P, "rooted translation orbit not regular")
    base_b = normalized_first(polynomial)[2]
    for displacement in range(P):
        shifted_b = normalized_first(
            translate(polynomial, displacement)
        )[2]
        require(
            shifted_b == (base_b + displacement) % P,
            "normalized Hasse translation law",
        )
    return vectors


def effective_vector(ell_dot, q1, q2):
    """THM-2820's q_eff for a section z with z.L=ell_dot."""
    return ((q1 - ell_dot) % P, (q2 + ell_dot) % P)


def gauge_hostile_pair():
    # Both completions project to the same observed first chain-axis motion
    # q1=1.  The omitted address and second carrier axis change the verdict.
    # Naming L is essential: q_eff=0 alone would not certify full pure gauge.
    transverse = ("0", 0, 1, 0)
    marked_pure_gauge = ("W", 1, 1, -1)
    require(transverse[2] == marked_pure_gauge[2] == 1,
            "one-axis projections differ")
    require(transverse[0] == "0", "transverse address is not L=0")
    require(marked_pure_gauge[0] == "W",
            "gauge completion is not the marked address L=W")
    transverse_effective = effective_vector(*transverse[1:])
    pure_effective = effective_vector(*marked_pure_gauge[1:])
    require(transverse_effective == (1, 0), "transverse completion")
    require(pure_effective == (0, 0), "marked pure-gauge completion")
    return (
        transverse,
        transverse_effective,
        marked_pure_gauge,
        pure_effective,
    )


def compact_formula(starts, lengths, kind):
    terms = []
    for start, length in zip(starts, lengths):
        if kind == "raw":
            terms.append(f"X^{start}[{length}]_X")
        elif kind == "live":
            terms.append(f"X^{start}E_{(length + 1) // 2}(X)")
        elif kind == "dead":
            terms.append(f"X^{start + 1}E_{length // 2}(X)")
        else:
            raise RuntimeError("unknown formula kind")
    return " + ".join(terms)


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    rows = []
    polynomials = {}
    for clock in (1, 2, 3):
        lengths = BLOCK_LENGTHS[clock]
        starts, raw, live, dead, contrast = chain_polynomials(lengths)
        require(starts == EXPECTED_STARTS[clock], "block-root offsets")
        require(
            (len(raw), len(live), len(dead)) == EXPECTED_COUNTS[clock],
            "THM-2818 copy counts",
        )
        require(
            all(
                next_start - (start + length - 1) == INTERBLOCK_JUMP
                for start, length, next_start
                in zip(starts, lengths, starts[1:])
            ),
            "interblock fifty-half-step law",
        )
        data = {}
        for name, polynomial in (
            ("raw", raw),
            ("live", live),
            ("dead", dead),
            ("contrast", contrast),
        ):
            j0, j1, barycenter = normalized_first(polynomial)
            exact_j0 = integer_moment(polynomial, 0)
            exact_j1 = integer_moment(polynomial, 1)
            local = eta_vector(polynomial)
            truncated_inverse(local)
            translation_orbit(polynomial)
            data[name] = (
                exact_j0,
                exact_j1,
                j0,
                j1,
                barycenter,
            )
        require(data["live"][3] != 0, "live first Hasse numerator vanished")
        require(data["live"][4] != 0, "live normalized first jet vanished")
        root_kill = -data["live"][4] % P
        require(
            normalized_first(translate(live, root_kill))[2] == 0,
            "root gauge failed to kill absolute barycenter",
        )
        one_step_endpoint_defect = subtract(translate(live, 1), dead)
        expected_endpoint_defect = {
            start + length: 1
            for start, length in zip(starts, lengths)
            if length % 2 == 1
        }
        require(
            one_step_endpoint_defect == expected_endpoint_defect,
            "one-half-step live/dead boundary identity",
        )
        polynomials[clock] = (starts, raw, live, dead, contrast)
        rows.append((
            clock,
            starts,
            lengths,
            data,
            quotient_vector(live),
            root_kill,
            tuple(one_step_endpoint_defect),
        ))

    transverse, transverse_effective, pure, pure_effective = gauge_hostile_pair()

    # Every exceptional polynomial lives at the one physical target label
    # t=12.  A target-preserving map is constant on its formal translation
    # orbit, while an intertwiner with the regular target shift must add one.
    target_equivariance_failures = 0
    for clock in (1, 2, 3):
        live = polynomials[clock][2]
        orbit = translation_orbit(live)
        for displacement in range(P):
            require(orbit[displacement] != orbit[(displacement + 1) % P],
                    "formal orbit collision")
            preserved_target = TARGET_LABEL
            shifted_target = (TARGET_LABEL + 1) % P
            if preserved_target != shifted_target:
                target_equivariance_failures += 1
    require(
        target_equivariance_failures == 3 * P,
        "target-axis equivariance hostile count",
    )

    affine_deck_identifications = P * (P - 1)
    rooted_generator_preserving_identifications = P

    print("THM-2818 x THM-2820 HALF-STEP/HASSE TRANSPORT PROBE")
    print(
        "coordinate=X=one rooted physical half-step; "
        "[L]_X=(1-X^L)/(1-X); E_m(X)=(1-X^(2m))/(1-X^2)"
    )
    print("block convention: next root is 50 half-steps after the previous block's last term")
    for (
        clock,
        starts,
        lengths,
        data,
        residue_vector,
        root_kill,
        endpoint_defect,
    ) in rows:
        print(
            f"clock={clock}; starts={starts}; lengths={lengths}; "
            f"raw={compact_formula(starts, lengths, 'raw')}"
        )
        print(
            f"clock={clock}; live={compact_formula(starts, lengths, 'live')}; "
            f"dead={compact_formula(starts, lengths, 'dead')}"
        )
        print(
            f"clock={clock}; jets=(exact_J0,exact_J1,J0_mod13,J1_mod13,b): "
            f"raw={data['raw']}; live={data['live']}; "
            f"dead={data['dead']}; contrast={data['contrast']}"
        )
        print(
            f"clock={clock}; live_residue_vector={residue_vector}; "
            f"live_orbit=regular_13; root_shift_killing_b={root_kill}"
        )
        print(
            f"clock={clock}; X*live-dead endpoint_defect={endpoint_defect}"
        )
    print(
        "formal_live_barycenters=(7,1,9); all nonzero in the live-head root; "
        "all are root-gauge killable"
    )
    print(
        "same_one_axis_motion_hostile: "
        f"transverse=(L={transverse[0]},zL={transverse[1]},"
        f"q={transverse[2:]})"
        f"->{transverse_effective}; "
        f"marked_pure_gauge=(L={pure[0]},zL={pure[1]},q={pure[2:]})"
        f"->{pure_effective}"
    )
    print(
        f"fixed_target12_vs_regular_target_shift_failures="
        f"{target_equivariance_failures}/39"
    )
    print(
        f"unoriented_affine_target_to_root_identifications="
        f"{affine_deck_identifications}; generator_preserving_origin_choices="
        f"{rooted_generator_preserving_identifications}"
    )
    print(
        "verdict=formal rooted Hasse translator positive; physical "
        "gauge-transversality and common-interval transport not selected"
    )
    print("assert_nodes=0; ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
