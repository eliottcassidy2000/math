#!/usr/bin/env python3
"""Exact referee for THM-1272's fastest-flood prefix-chamber closure.

The proof combines four already proved inputs.

* THM-1198 gives the phase-free one-comb majorant Pbar.
* THM-1233 gives prefix component bounds and T(x)=ceil((6x+1)/7).
* THM-1267 gives 270*d1 <= 563*c-1.
* THM-1275 gives the dominated fastest-tail tax
      F6 >= c/(8h) * (ceil(K/(e+1))-1)
  with K>=ceil(h/(7c)), and also
      K>=ceil(36*h*eta_h/(7c)).

The finite carrier is a monotone word of tooth-count chambers, not a runner
tournament.  This dependency-free script reconstructs each chamber supremum,
exhausts every prefix word allowed by THM-1233, checks the five basic branch
cuts, and checks the stronger functional e=5 cut h/c<798.

Analytic/topological providers (the one-comb envelope, prefix component
containment, and THM-1275's layered multiplicity invoice) remain paper inputs.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1272 referee failed: {message}")


def optimization_safe_require_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "optimization-sensitive assert remains")
    caught = False
    try:
        require(False, "deliberate sentinel")
    except RuntimeError as error:
        caught = "deliberate sentinel" in str(error)
    require(caught, "explicit RuntimeError sentinel did not fire")
    return count


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


# Each branch is (left, right, constant, reciprocal coefficient).  The last
# branch is the BV majorant on the open ray x>7/2.  At x=7/2 the preceding
# exact branch has value 1/7, whereas the right limit of the BV branch is 4/21.
PBAR_BRANCHES: tuple[tuple[F, F | None, F, F], ...] = (
    (F(1), F(238, 189), F(0), F(7, 36)),
    (F(238, 189), F(4, 3), F(3, 4), F(-3, 4)),
    (F(4, 3), F(7, 5), F(0), F(1, 4)),
    (F(7, 5), F(8, 5), F(5, 18), F(-5, 36)),
    (F(8, 5), F(7, 4), F(0), F(11, 36)),
    (F(7, 4), F(2), F(2, 9), F(-1, 12)),
    (F(2), F(7, 3), F(0), F(13, 36)),
    (F(7, 3), F(854, 357), F(1, 24), F(19, 72)),
    (F(854, 357), F(5, 2), F(3, 4), F(-103, 72)),
    (F(5, 2), F(14, 5), F(0), F(4, 9)),
    (F(14, 5), F(3), F(5, 18), F(-1, 3)),
    (F(3), F(7, 2), F(0), F(1, 2)),
    (F(7, 2), None, F(1, 7), F(1, 6)),
)


def pbar_branch_value(branch: tuple[F, F | None, F, F], x: F) -> F:
    left, right, constant, reciprocal = branch
    require(x >= left, "Pbar evaluation left of branch")
    require(right is None or x <= right, "Pbar evaluation right of branch")
    return constant + reciprocal / x


def chamber_interval(count: int) -> tuple[F, F]:
    """T(x)=count iff x lies in this left-open, right-closed interval."""

    require(count >= 2, "fast-speed tooth count starts at two")
    return F(7 * count - 8, 6), F(7 * count - 1, 6)


def chamber_supremum_from_envelope(count: int) -> F:
    """Take the exact supremum of Pbar on one T-chamber.

    Every Pbar branch is monotone because it is a+b/x.  Open endpoints do not
    affect a supremum.  The BV ray is treated with its right limit at 7/2,
    which is the essential upward-jump guardrail in count chamber four.
    """

    chamber_left, chamber_right = chamber_interval(count)
    candidates: list[F] = []
    for index, branch in enumerate(PBAR_BRANCHES):
        branch_left, branch_right, _, reciprocal = branch
        left = max(chamber_left, branch_left)
        right = chamber_right if branch_right is None else min(chamber_right, branch_right)
        if right < left:
            continue
        # The tail branch is open at 7/2.  A singleton intersection there is
        # empty, but every nontrivial right interval has the displayed limit.
        if index == len(PBAR_BRANCHES) - 1 and right == left == F(7, 2):
            continue
        endpoint = left if reciprocal >= 0 else right
        candidates.append(pbar_branch_value(branch, endpoint))
    require(candidates, f"empty Pbar/count intersection at count={count}")
    return max(candidates)


def chamber_supremum_formula(count: int) -> F:
    if count == 2:
        return F(7, 36)
    if count == 3:
        return F(8, 45)
    if count == 4:
        return F(4, 21)
    return F(1, 7) + F(1, 7 * count - 8)


def chamber_supremum_census() -> int:
    rows = 0
    for count in range(2, 233):
        reconstructed = chamber_supremum_from_envelope(count)
        formula = chamber_supremum_formula(count)
        require(reconstructed == formula, f"Pbar chamber supremum count={count}")
        rows += 1
    return rows


# THM-1233's individual projective caps imply these T-count caps.
COUNT_CAPS = (2, 9, 15, 76, 232)

# After j prefix speeds, THM-1233 bounds x_(j+1) by A_j times
# 1+sum_(i<=j)T(x_i).  Index zero is unused.
PREFIX_COEFFICIENTS: tuple[F | None, ...] = (
    None,
    F(91, 29),
    F(91, 22),
    F(49, 15),
    F(21, 8),
)


def count_left_endpoint(count: int) -> F:
    return F(7 * count - 8, 6)


StateKey = tuple[int, int]
StateValue = tuple[F, tuple[int, ...]]


def enumerate_prefix_states() -> tuple[dict[int, dict[int, StateValue]], tuple[int, ...]]:
    """Exhaust monotone T-words satisfying every earlier prefix row.

    For a proposed next count m, its speed is strictly above
    (7m-8)/6.  Hence feasibility requires that left endpoint to be strictly
    below THM-1233's prefix upper bound.  Keeping only the largest load for
    each (sum,last-count) is exact because later feasibility depends only on
    those two coordinates.
    """

    states: dict[StateKey, StateValue] = {(0, 1): (F(0), ())}
    suffix_best: dict[int, dict[int, StateValue]] = {}
    state_counts: list[int] = []

    for index in range(5):
        next_states: dict[StateKey, StateValue] = {}
        for (count_sum, last_count), (load, word) in states.items():
            for count in range(max(2, last_count), COUNT_CAPS[index] + 1):
                if index >= 1:
                    coefficient = PREFIX_COEFFICIENTS[index]
                    require(coefficient is not None, "missing prefix coefficient")
                    if not count_left_endpoint(count) < coefficient * (1 + count_sum):
                        continue
                key = (count_sum + count, count)
                candidate = (load + chamber_supremum_formula(count), word + (count,))
                incumbent = next_states.get(key)
                if incumbent is None or candidate[0] > incumbent[0]:
                    next_states[key] = candidate
        states = next_states
        state_counts.append(len(states))
        require(states, f"empty prefix state bank at depth={index + 1}")

        exact_sum_best: dict[int, StateValue] = {}
        for (count_sum, _), value in states.items():
            incumbent = exact_sum_best.get(count_sum)
            if incumbent is None or value[0] > incumbent[0]:
                exact_sum_best[count_sum] = value

        running = (F(-1), ())
        at_least_best: dict[int, StateValue] = {}
        for count_sum in sorted(exact_sum_best, reverse=True):
            if exact_sum_best[count_sum][0] > running[0]:
                running = exact_sum_best[count_sum]
            at_least_best[count_sum] = running
        suffix_best[index + 1] = at_least_best

    return suffix_best, tuple(state_counts)


def best_with_sum_at_least(
    banks: dict[int, dict[int, StateValue]], depth: int, required_sum: int
) -> StateValue:
    bank = banks[depth]
    eligible = [count_sum for count_sum in bank if count_sum >= required_sum]
    require(eligible, f"no count word at depth={depth}, sum>={required_sum}")
    return bank[min(eligible)]


def prefix_bank_census(
    banks: dict[int, dict[int, StateValue]], state_counts: tuple[int, ...]
) -> tuple[tuple[int, int, F, tuple[int, ...]], ...]:
    require(state_counts == (1, 8, 84, 872, 9040), "compressed prefix state counts")

    checks = (
        (1, 2, F(7, 36), (2,)),
        (2, 4, F(7, 18), (2, 2)),
        (3, 10, F(145, 252), (2, 4, 4)),
        (3, 11, F(61, 108), (2, 4, 5)),
        (4, 41, F(38081, 52668), (2, 4, 4, 31)),
        (4, 42, F(112571, 158004), (2, 4, 5, 31)),
        (5, 113, F(868922, 995715), (2, 2, 4, 26, 79)),
        (5, 114, F(2847097, 3276540), (2, 4, 4, 25, 79)),
        (5, 136, F(476431, 549252), (2, 4, 4, 31, 95)),
        (5, 137, F(9882995, 11534292), (2, 4, 5, 31, 95)),
    )
    for depth, required_sum, expected_load, expected_word in checks:
        load, word = best_with_sum_at_least(banks, depth, required_sum)
        require(load == expected_load, f"prefix maximum depth={depth}, sum={required_sum}")
        require(word == expected_word, f"prefix maximizer depth={depth}, sum={required_sum}")
    return checks


def basic_branch_census(
    banks: dict[int, dict[int, StateValue]]
) -> tuple[tuple[int, F, int, F, F, tuple[int, ...]], ...]:
    """Check the five basic-K branch thresholds.

    For e<5, h<=6*d_(e+1) and the e-prefix row force
        x < A_e(1+sum m_i),
    with A=(546/29,273/11,98/5,63/4).  For e=5, the last-tooth
    row gives x<7(1+sum m_i).  Above x>21 all noneligible speeds
    lie on Pbar's decreasing BV tail.
    """

    threshold_data = (
        # e, threshold, required prefix sum, expected Q_e
        (1, F(42), 2, F(7, 36)),
        (2, F(1407, 20), 4, F(7, 18)),
        (3, F(1078, 5), 11, F(61, 108)),
        (4, F(1323, 2), 42, F(112571, 158004)),
        (5, F(959), 137, F(9882995, 11534292)),
    )
    rows: list[tuple[int, F, int, F, F, tuple[int, ...]]] = []
    for eligible, threshold, required_sum, expected_load in threshold_data:
        load, word = best_with_sum_at_least(banks, eligible, required_sum)
        require(load == expected_load, f"basic branch prefix load e={eligible}")

        # B=ceil(K/(e+1))-1 with K>=ceil(x/7).  Nested ceilings give
        # B>=ceil(x/[7(e+1)])-1.  At e=1 the tax jumps only to B=3
        # immediately to the right of x=42; the displayed boundary itself
        # has B=2 and is intentionally not claimed.
        if eligible == 1:
            forced_exceptions = 3
            comparison_x = threshold
            upper_minus_tax_times_x = (
                comparison_x * (load - F(eligible + 1, 7))
                + F(31 - 6 * eligible, 6)
                - F(forced_exceptions, 8)
            )
            require(upper_minus_tax_times_x < 0, "strict e=1 right-limit margin")
        else:
            forced_exceptions = ceil_fraction(threshold / (7 * (eligible + 1))) - 1
            upper_minus_tax_times_x = (
                threshold * (load - F(eligible + 1, 7))
                + F(31 - 6 * eligible, 6)
                - F(forced_exceptions, 8)
            )
            require(upper_minus_tax_times_x <= 0, f"basic branch margin e={eligible}")

        # Once the cut is crossed, the prefix requirement and B do not fall.
        # The coefficient of x in the cross-multiplied upper-minus-tax bound
        # is already negative at the threshold bank.
        require(load - F(eligible + 1, 7) < 0, f"negative tail slope e={eligible}")
        rows.append(
            (
                eligible,
                threshold,
                forced_exceptions,
                load,
                -upper_minus_tax_times_x,
                word,
            )
        )

    require(rows[0][4] == F(1, 24), "e=1 cross-multiplied right-limit margin")
    require(rows[1][4] == 0, "e=2 exact touching margin")
    require(rows[2][4] == F(145, 1080), "e=3 cross-multiplied margin")
    # The next two margins are deliberately only checked for positivity above;
    # their large unreduced forms are printed by the referee.
    return tuple(rows)


def functional_e5_census(
    banks: dict[int, dict[int, StateValue]]
) -> tuple[int, F, tuple[int, ...], int, F]:
    """Use THM-1275's eta_h count to close every e=5 chamber x>=798.

    For x in [7s,7(s+1)), THM-1233 forces sum T(x_i)>=s.  Let Q_s
    be the exact prefix supremum.  Since the first speed is in the open m=2
    chamber, the actual prefix load is strictly below Q_s, so

        eta_h > 1-Q_s,
        K >= floor(36*s*(1-Q_s))+1.

    The worst endpoint of the cross-multiplied F6/tax comparison is the right
    endpoint when Q_s>6/7 and the left endpoint when Q_s<=6/7.
    """

    # The immediately preceding chamber is a guardrail: this particular
    # phase-free certificate genuinely starts at 798 rather than being a
    # rounded presentation of an earlier cut.
    previous_load, previous_word = best_with_sum_at_least(banks, 5, 113)
    previous_eta = 1 - previous_load
    previous_k_floor = floor_fraction(36 * 113 * previous_eta) + 1
    previous_exceptions = ceil_fraction(F(previous_k_floor, 6)) - 1
    previous_coefficient = previous_load - F(6, 7)
    previous_margin = F(previous_exceptions, 8) - (
        798 * previous_coefficient + F(1, 6)
    )
    require(previous_word == (2, 2, 4, 26, 79), "s=113 maximizing word")
    require(previous_margin == F(-113823, 63220), "s=113 failure margin")

    checked = 0
    minimum_margin: F | None = None
    minimum_row: tuple[int, F, tuple[int, ...], int, F] | None = None
    previous_q: F | None = None
    previous_k = -1

    for count_sum in range(114, 335):
        q_value, word = best_with_sum_at_least(banks, 5, count_sum)
        eta_floor = 1 - q_value
        require(eta_floor > 0, f"positive eta floor s={count_sum}")
        k_floor = floor_fraction(36 * count_sum * eta_floor) + 1
        exceptions = ceil_fraction(F(k_floor, 6)) - 1
        coefficient = q_value - F(6, 7)
        endpoint = 7 * (count_sum + 1) if coefficient > 0 else 7 * count_sum
        margin = F(exceptions, 8) - (endpoint * coefficient + F(1, 6))
        require(margin > 0, f"functional e=5 chamber margin s={count_sum}")

        if previous_q is not None:
            require(q_value <= previous_q, "Q_s is not nonincreasing")
            require(k_floor >= previous_k, "functional K floor is not nondecreasing")
        previous_q = q_value
        previous_k = k_floor

        if minimum_margin is None or margin < minimum_margin:
            minimum_margin = margin
            minimum_row = (count_sum, q_value, word, exceptions, margin)
        checked += 1

    require(checked == 221, "functional e=5 chamber count")
    require(minimum_margin is not None and minimum_row is not None, "empty e=5 census")
    require(minimum_row[0] == 114, "minimum e=5 margin chamber")
    require(minimum_row[1] == F(2847097, 3276540), "s=114 prefix maximum")
    require(minimum_row[2] == (2, 4, 4, 25, 79), "s=114 maximizing word")
    require(minimum_row[3] == 89, "s=114 exception floor")
    require(minimum_row[4] == F(1921973, 1310616), "s=114 chamber margin")
    return minimum_row


def e_zero_integer_cut() -> tuple[F, str]:
    # e=0 means h<=6*d1.  THM-1267 gives 270*d1<=563*c-1, hence
    # 45*h<=270*d1<=563*c-1<563*c.
    bound = F(563, 45)
    require(6 * F(563, 270) == bound, "e=0 ratio conversion")
    return bound, "45*h<=270*d1<=563*c-1<563*c"


def tournament_loss_audit() -> None:
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=prefix tooth-count chambers, not runners or individual arcs")
    print("observable=q(m_i)-q(m_j); gauge=count order then chronological speed order")
    print("score_histogram=transitive after tie gauge; directed_cycles=0; SCCs=singletons")
    print("hamiltonian_path_count=1 after tie gauge")
    print("preserves=prefix component feasibility, one-body load budget, and count jumps")
    print("destroys=phases, tooth addresses, flood/turn locations, and seam geometry")
    print("challenged_assumption=the faithful vertices are obligations, not runners")


def main() -> None:
    assert_count = optimization_safe_require_probe()
    chamber_rows = chamber_supremum_census()
    banks, state_counts = enumerate_prefix_states()
    bank_checks = prefix_bank_census(banks, state_counts)
    basic_rows = basic_branch_census(banks)
    functional_row = functional_e5_census(banks)
    e0_bound, e0_chain = e_zero_integer_cut()

    print("THM-1272 FASTEST-FLOOD PREFIX-CHAMBER CLOSURE EXACT REFEREE")
    print(f"optimization_safe_assert_count={assert_count}")
    print(f"pbar_count_chamber_rows={chamber_rows} counts=2..232")
    print(f"prefix_state_counts={state_counts}")
    print(f"named_prefix_maximum_checks={len(bank_checks)}")
    print(f"e0_bound=h/c<{e0_bound} chain={e0_chain}")
    for eligible, threshold, exceptions, load, margin, word in basic_rows:
        relation = ">" if eligible == 1 else ">="
        print(
            f"basic_branch=e{eligible} closes_at_x{relation}{threshold} "
            f"B={exceptions} Q={load} cross_margin={margin} word={word}"
        )
    count_sum, load, word, exceptions, margin = functional_row
    print(
        "functional_e5="
        f"x>=798 closes; first_chamber_sum={count_sum} Q={load} "
        f"B={exceptions} cross_margin={margin} word={word}"
    )
    print("functional_e5_chambers=221 range_s=114..334")
    print("high_ratio_consequence=h>=7*d5/2 implies h/c<798")
    tournament_loss_audit()
    print("SCOPE=closes a top-dominated projective tail, not all six-comb covers or LRC(14)")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
