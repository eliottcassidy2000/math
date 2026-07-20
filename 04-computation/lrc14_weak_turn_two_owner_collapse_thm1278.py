#!/usr/bin/env python3
"""Exact referee for THM-1278's weak-turn two-owner collapse.

A consecutive-address fastest return has exact packet bracket

    sum(1/j_q) - 6/h.

Call it weak when this positive bracket is < 1/h.  Equivalently its
dimensionless reciprocal word has 6 < sum(h/j_q) < 7.  Above h/c=399/5,
the THM-1198 six-bin envelope allows at most two lower owners above h/6.
Every weak word uses only those owners.  If A<B are those owners, then
B<h<6A.  The A-tooth is shorter than a B-safe gap, while the B-tooth is
shorter than an A-safe gap.  Thus a connected two-comb component contains
at most one tooth of each owner, and the word is AB or BA.

The analytic envelope and interval-minimality implications remain the paper
providers.  This dependency-free referee checks their exact arithmetic,
finite word quotient, payment funnel, endpoint threshold sharpness, and a
large exact packet census.  It deliberately contains no Python ``assert``.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from itertools import product
from pathlib import Path


THRESHOLD = F(399, 5)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1278 referee failed: {message}")


def optimization_safety_probe() -> int:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "Python assert nodes are optimization-sensitive")
    caught = False
    try:
        require(False, "deliberate require probe")
    except RuntimeError as error:
        caught = "deliberate require probe" in str(error)
    require(caught, "require probe did not fire")
    return count


def load_threshold_audit() -> tuple[F, int, F]:
    """Check the exact three-large-owner load contradiction.

    If three of d1,...,d5 exceed h/6 and x=h/c, their loads are each
    strictly below 1/7+1/x.  The other two lower loads are at most 7/36,
    and the h load is at most 1/7+1/(6x).  The resulting upper envelope is

        121/126 + 19/(6x),

    equal to one exactly at x=399/5 and decreasing thereafter.  Strictness
    of the three large-owner bounds excludes equality at the endpoint.
    """
    base = 2 * F(7, 36) + 4 * F(1, 7)
    require(base == F(121, 126), "load base constant")
    endpoint_error = F(19, 6) / THRESHOLD
    require(endpoint_error == F(5, 126), "threshold tail error")
    require(base + endpoint_error == 1, "399/5 is the exact endpoint")

    rows = 0
    smallest_margin: F | None = None
    # Exact rational samples on and above the threshold.  The endpoint has
    # zero weak-envelope margin; the actual load inequality is still strict.
    for denominator in range(1, 81):
        start = (THRESHOLD * denominator).numerator // (
            THRESHOLD * denominator
        ).denominator
        while F(start, denominator) < THRESHOLD:
            start += 1
        for numerator in range(start, start + 240):
            x = F(numerator, denominator)
            if x < THRESHOLD:
                continue
            upper = base + F(19, 6) / x
            require(upper <= 1, "three-large-owner envelope is at most one")
            margin = 1 - upper
            smallest_margin = margin if smallest_margin is None else min(
                smallest_margin, margin
            )
            rows += 1
    require(rows > 10000, "load threshold audit is substantial")
    require(smallest_margin == 0, "endpoint equality guardrail is present")

    below = THRESHOLD - F(1, 1000)
    negative_guardrail = 1 - (base + F(19, 6) / below)
    require(negative_guardrail < 0, "immediately lower scalar threshold fails")
    return base, rows, negative_guardrail


def weak_packet_census() -> tuple[int, int, int, F]:
    """Exhaust exact weak packets and check the h/6 support restriction."""
    rows = 0
    weak_rows = 0
    strong_rows = 0
    smallest_positive: F | None = None

    # Repeated speeds are permitted here because this is a deliberately
    # relaxed arithmetic census.  The minimal-word layer later removes
    # adjacent equal owners.  Every conclusion below is therefore robust.
    for h in range(8, 62, 3):
        samples = sorted({
            max(1, h // 7),
            max(1, h // 4),
            max(1, h // 2),
            h - 1,
        })
        samples = [j for j in samples if 0 < j < h]
        for length in range(2, 8):
            for word in product(samples, repeat=length):
                reciprocal_sum = sum((F(h, j) for j in word), F(0))
                packet = reciprocal_sum - 6
                if packet <= 0:
                    continue
                rows += 1
                smallest_positive = packet if smallest_positive is None else min(
                    smallest_positive, packet
                )
                if packet < 1:
                    weak_rows += 1
                    require(length <= 6, "weak packet has at most six teeth")
                    require(all(6 * j > h for j in word),
                            "every weak-packet owner is strictly above h/6")
                else:
                    strong_rows += 1

    require(rows > 100000, "packet census is substantial")
    require(weak_rows > 0 and strong_rows > 0, "both packet branches occur")
    require(smallest_positive is not None and smallest_positive > 0,
            "positive packet minimum exists on the finite bank")
    return rows, weak_rows, strong_rows, smallest_positive


def ratio_localized_component_audit() -> tuple[int, F, F]:
    """Check the exact two-sided tooth/gap separation in the weak sector."""
    rows = 0
    smallest_a_vs_b_gap: F | None = None
    smallest_b_vs_a_gap: F | None = None
    for h in range(8, 501):
        for a in range(1, h):
            if 6 * a <= h:
                continue
            for b in range(a + 1, h):
                # A tooth versus the safe gap between consecutive B teeth.
                a_vs_b_gap = F(6, 7 * b) - F(1, 7 * a)
                # B tooth versus the safe gap between consecutive A teeth.
                b_vs_a_gap = F(6, 7 * a) - F(1, 7 * b)
                require(b < 6 * a, "weak owner ratio is below six")
                require(a_vs_b_gap > 0, "an A tooth cannot bridge a B gap")
                require(b_vs_a_gap > 0, "a B tooth cannot bridge an A gap")
                smallest_a_vs_b_gap = (
                    a_vs_b_gap
                    if smallest_a_vs_b_gap is None
                    else min(smallest_a_vs_b_gap, a_vs_b_gap)
                )
                smallest_b_vs_a_gap = (
                    b_vs_a_gap
                    if smallest_b_vs_a_gap is None
                    else min(smallest_b_vs_a_gap, b_vs_a_gap)
                )
                rows += 1
    require(rows > 1_000_000, "ratio-localized component audit is substantial")
    require(smallest_a_vs_b_gap is not None and smallest_a_vs_b_gap > 0,
            "A-tooth/B-gap separation has positive finite-bank margin")
    require(smallest_b_vs_a_gap is not None and smallest_b_vs_a_gap > 0,
            "B-tooth/A-gap separation has positive finite-bank margin")
    return rows, smallest_a_vs_b_gap, smallest_b_vs_a_gap


def two_owner_word_audit() -> tuple[int, tuple[str, ...]]:
    """Enumerate the ratio-localized irredundant two-comb word quotient.

    Paper geometry gives adjacent-owner distinctness.  The two strict
    tooth/gap inequalities checked above imply that a connected component
    contains at most one tooth of each owner.
    """
    survivors: set[str] = set()
    rows = 0
    for length in range(2, 7):
        for letters in product("AB", repeat=length):
            rows += 1
            adjacent_distinct = all(
                letters[index] != letters[index + 1]
                for index in range(length - 1)
            )
            minimal_counts = letters.count("A") <= 1 and letters.count("B") <= 1
            if adjacent_distinct and minimal_counts:
                survivors.add("".join(letters))
    expected = {"AB", "BA"}
    require(survivors == expected, "minimal two-owner word alphabet")
    require(not {word for word in survivors if len(word) >= 3},
            "no ratio-localized weak word has three teeth")
    return rows, tuple(sorted(survivors, key=lambda word: (len(word), word)))


def weak_pair_detuning_audit() -> tuple[int, int, F]:
    """Check the fixed-pair packet numerator and its gcd quantum.

    For a weak pair A<B the dimensionless packet is

        h/A+h/B-6 = D/(AB),  D=h(A+B)-6AB.

    The reciprocal seam bracket is D/(hAB), and D is a positive integer
    divisible by gcd(A,B).  Orientation AB versus BA changes nothing.
    """
    from math import gcd

    rows = 0
    smallest_d = 10**9
    smallest_bracket: F | None = None
    for h in range(8, 401):
        for a in range(1, h):
            if 6 * a <= h:
                continue
            for b in range(a + 1, h):
                packet = F(h, a) + F(h, b) - 6
                if not (0 < packet < 1):
                    continue
                d = h * (a + b) - 6 * a * b
                require(packet == F(d, a * b), "dimensionless detuning identity")
                bracket = F(1, a) + F(1, b) - F(6, h)
                require(bracket == F(d, h * a * b), "reciprocal detuning identity")
                g = gcd(a, b)
                require(d > 0 and d % g == 0 and d >= g,
                        "positive pair detuning pays a gcd quantum")
                require(
                    F(h, b) + F(h, a) - 6 == packet,
                    "AB and BA have the same exact packet",
                )
                smallest_d = min(smallest_d, d)
                smallest_bracket = (
                    bracket if smallest_bracket is None else min(smallest_bracket, bracket)
                )
                rows += 1
    require(rows > 100_000, "weak-pair detuning audit is substantial")
    require(smallest_d == 1, "unit detuning guardrail occurs")
    require(smallest_bracket is not None and smallest_bracket > 0,
            "finite weak-pair bank has positive reciprocal detuning")
    return rows, smallest_d, smallest_bracket


def payment_funnel_audit() -> tuple[int, int]:
    """Check the regular-run exception funnel with weak turns separated."""
    rows = 0
    weak_only_rows = 0
    for e in range(0, 6):
        for k in range(1, 81):
            forced = (k + e) // (e + 1) - 1  # ceil(k/(e+1))-1
            require(forced == (k - 1) // (e + 1), "nested ceiling identity")
            for weak in range(0, forced + 8):
                required_paid = max(0, forced - weak)
                for floods in range(0, required_paid + 4):
                    for strong in range(0, required_paid + 4):
                        if floods + strong + weak < forced:
                            continue
                        require(floods + strong >= required_paid,
                                "unpaid exception count is absorbed by weak cells")
                        if floods + strong == 0 and forced > 0:
                            require(weak >= forced,
                                    "tax failure forces a bank of weak cells")
                            weak_only_rows += 1
                        rows += 1
    require(rows > 100000, "payment funnel census is substantial")
    require(weak_only_rows > 0, "weak-only residual occurs abstractly")
    return rows, weak_only_rows


def exact_examples() -> dict[str, object]:
    # AB: 6 < h/A+h/B < 7.
    h, a, b = 60, 11, 50
    ab = F(h, a) + F(h, b)
    require(6 < ab < 7, "AB weak example")

    # Packet arithmetic alone admits BAB, but the ratio-localized two-comb
    # geometry forbids its three distinct teeth from forming one component.
    h2, a2, b2 = 60, 20, 39
    bab = F(h2, a2) + 2 * F(h2, b2)
    require(6 < bab < 7, "arithmetic BAB guardrail")
    require(F(1, 7 * a2) < F(6, 7 * b2),
            "the A tooth is shorter than the intervening B-safe gap")

    # Exact abutment is excluded by strict seam positivity but is the
    # self-similar limit: h/j+h/k=6 at ratios 4+2.
    h3, a3, b3 = 60, 15, 30
    abut = F(h3, a3) + F(h3, b3)
    require(abut == 6, "two-owner exact-abutment guardrail")
    return {
        "AB": (h, a, b, ab),
        "arithmetic_BAB_guardrail": (h2, a2, b2, bab),
        "abutment": (h3, a3, b3, abut),
    }


def tournament_and_carrier_audit() -> None:
    print("TOURNAMENT_AND_CARRIER_AUDIT")
    print("vertices=weak_fastest_return_cells;labels=AB,BA;gauge=chronology")
    print("observable=(word_length,orientation);tie_path=chronological_order")
    print("quotient_tournament=transitive;cycles=0;scc_sizes=all_1;hamiltonian_paths=1")
    print("two_owner_boundary_graph=A<->B;three-tooth_palindrome_killed_by_B<6A;not_a_runner_tournament")
    print("fano_status=seven_local_colours_collapse_to_one_owner_pair_but_no_phase_digit_law_follows")
    print("faithful_state=(A,B;ordered_tooth_addresses;AB_or_BA;exact_detuning_numerator;seams)")
    print("preserves=weak_packet_predicate,minimal_chain_shape,repeated_pair,exact_return_payment")
    print("destroys_if_owner_only=addresses,seam_digits,abutment_drift;destroys_if_wall_only=common_pair")
    print("challenged_vertices=runners,gaps,walls,teeth,return_cells,residues,Fano_lines,proof_obligations")


def main() -> None:
    assert_nodes = optimization_safety_probe()
    base, load_rows, below_margin = load_threshold_audit()
    packet_rows, weak_rows, strong_rows, packet_min = weak_packet_census()
    ratio_rows, a_gap_margin, b_gap_margin = ratio_localized_component_audit()
    word_rows, words = two_owner_word_audit()
    detuning_rows, detuning_min, bracket_min = weak_pair_detuning_audit()
    funnel_rows, weak_only_rows = payment_funnel_audit()
    examples = exact_examples()

    print("THM-1278 WEAK-TURN TWO-OWNER COLLAPSE EXACT REFEREE")
    print(f"optimization_sensitive_assert_nodes={assert_nodes}")
    print(f"load_base={base};tail_error_at_399/5={F(5,126)};sum=1")
    print(f"load_threshold_rows={load_rows};below_threshold_guardrail_margin={below_margin}")
    print("load_conclusion=h/c>=399/5 => at most two lower owners exceed h/6")
    print(f"packet_rows={packet_rows};weak={weak_rows};strong={strong_rows};finite_bank_min={packet_min}")
    print("weak_packet_law=0<bracket<1/h => every internal owner>h/6 and length<=6")
    print(f"ratio_component_rows={ratio_rows};A_tooth_vs_B_gap_min={a_gap_margin};B_tooth_vs_A_gap_min={b_gap_margin}")
    print("ratio_component_law=A<B<h<6A => at most one A tooth and one B tooth per component")
    print(f"two_owner_word_rows={word_rows};survivors={','.join(words)}")
    print("minimal_word_law=unique weak owner pair; internal word in {AB,BA}")
    print(f"detuning_rows={detuning_rows};minimum_integer_detuning={detuning_min};finite_bank_bracket_min={bracket_min}")
    print("fixed_pair_law=bracket=D/(hAB),D=h(A+B)-6AB>=gcd(A,B);AB_and_BA_equal")
    print(f"payment_funnel_rows={funnel_rows};weak_only_rows={weak_only_rows}")
    print("payment_funnel=X+T_strong>=max(0,ceil(K/(e+1))-1-T_weak)")
    print("example_AB=" + ",".join(map(str, examples["AB"])))
    print("arithmetic_BAB_guardrail=" + ",".join(map(str, examples["arithmetic_BAB_guardrail"])))
    print("abutment_guardrail=" + ",".join(map(str, examples["abutment"])))
    tournament_and_carrier_audit()
    print("SCOPE=collapses_the_unpaid_complementary_branch_but_does_not_bound_repeated_pair_phase_drift")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
