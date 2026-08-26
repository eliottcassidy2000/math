#!/usr/bin/env python3
"""Exact primary audit for THM-4213.

The load-bearing theorem is algebraic: contextual prefix defects telescope
under ordinal sum, and certified Hamilton-normalized floors therefore obey a
weighted semigroup law.  This audit checks that symbolic law, replays it on
exact transferred tournament data, checks the tail-five word floors on a
labelled hostile/control bank, and emits the canonical good-suffix hostile.
"""

from __future__ import annotations

import hashlib
import itertools

import sympy as sp

import tournament_ordinal_cocycle_parity_thm4184 as base


TAIL_FLOOR = 10764


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def labelled_through_three() -> tuple[base.TournamentData, ...]:
    answer: list[base.TournamentData] = []
    for order in range(1, 4):
        width = order * (order - 1) // 2
        for integer in range(1 << width):
            bits = "" if width == 0 else format(integer, f"0{width}b")
            answer.append(base.tournament_data(base.parse(bits, order)))
    return tuple(answer)


def ordinal(*factors: base.TournamentData) -> base.TournamentData:
    need(bool(factors), "ordinal product needs a factor")
    value = factors[0]
    for factor in factors[1:]:
        value = base.ordinal_data(value, factor)
    return value


def defect(
    prefix: base.TournamentData,
    middle: base.TournamentData,
    right: base.TournamentData,
) -> int:
    return base.remainder(base.ordinal_data(prefix, middle), right) - base.remainder(
        middle, right
    )


def word_floor(cycle_blocks: int) -> int:
    need(cycle_blocks >= 1, "positive number of tail-five blocks")
    return TAIL_FLOOR * (9**cycle_blocks - 1) // 8


def main() -> None:
    # The exact telescoping identity is formal and uses no inequality.
    r_adb_c, r_db_c, r_b_c = sp.symbols("r_adb_c r_db_c r_b_c")
    same_identity = sp.expand(
        (r_adb_c - r_b_c)
        - ((r_adb_c - r_db_c) + (r_db_c - r_b_c))
    )
    need(same_identity == 0, "symbolic contextual-defect telescope")

    alpha, beta, h_d, h_b, h_c = sp.symbols(
        "alpha beta h_d h_b h_c", nonnegative=True
    )
    composed_floor = sp.expand(
        alpha * (h_d * h_b) ** 2 * h_c**2
        + beta * h_b**2 * h_c**2
    )
    expected_floor = sp.expand(
        (alpha * h_d**2 + beta) * h_b**2 * h_c**2
    )
    need(composed_floor == expected_floor, "weighted certified-floor law")

    one = base.tournament_data(base.transitive(1))
    p2 = base.tournament_data(base.transitive(2))
    p3 = base.tournament_data(base.transitive(3))
    cycle = base.tournament_data(base.parse("101", 3))
    q = {
        tail: ordinal(cycle, base.tournament_data(base.transitive(tail)))
        for tail in range(5, 8)
    }
    contexts = labelled_through_three()
    need(len(contexts) == 11, "labelled context count through order three")

    # Exact transferred-data controls for the formal telescope.
    prefix_bank = (one, p2, cycle, q[5])
    telescope_checks = 0
    for left, second in itertools.product(prefix_bank, repeat=2):
        composite = base.ordinal_data(left, second)
        for middle, right in itertools.product(contexts, repeat=2):
            lhs = defect(composite, middle, right)
            rhs = defect(left, base.ordinal_data(second, middle), right) + defect(
                second, middle, right
            )
            need(lhs == rhs, "exact contextual-defect telescope control")
            telescope_checks += 1

    # Direct finite controls of the inherited THM-4212 floor.
    tail_floor_checks = 0
    for tail in range(5, 8):
        for middle, right in itertools.product(contexts, repeat=2):
            value = defect(q[tail], middle, right)
            floor = TAIL_FLOOR * middle.hamilton**2 * right.hamilton**2
            need(value >= floor, "tail-five inherited floor control")
            if value == floor:
                need(
                    tail == 5
                    and len(middle.out) == len(right.out) == 1,
                    "tail-floor equality boundary control",
                )
            tail_floor_checks += 1

    # Quantitative multi-cycle controls.  The first word is source-free and
    # has two nontrivial strong components, so it is outside C3 ▷ P_n.
    two_cycle = ordinal(q[5], q[5])
    two_cycle_floor = word_floor(2)
    two_cycle_checks = 0
    for middle, right in itertools.product(contexts, repeat=2):
        value = defect(two_cycle, middle, right)
        floor = two_cycle_floor * middle.hamilton**2 * right.hamilton**2
        need(value > floor, "two-cycle strict geometric floor control")
        two_cycle_checks += 1

    padded_word = ordinal(p2, q[5], p3, q[6], one)
    padded_floor = word_floor(2)
    padded_checks = 0
    for middle, right in itertools.product(contexts, repeat=2):
        value = defect(padded_word, middle, right)
        floor = padded_floor * middle.hamilton**2 * right.hamilton**2
        need(value > floor, "transitively padded word floor control")
        padded_checks += 1

    # The repaired boundary: a good terminal Q5 is not enough if an arbitrary
    # (non-neutral) factor is multiplied on the left.
    hostile = ordinal(cycle, q[5])
    hostile_value = defect(hostile, one, one)
    need(len(hostile.out) == 11, "hostile order")
    need(hostile.hamilton == 9, "hostile Hamilton count")
    need(hostile_value == -338580, "good-suffix hostile value")

    # Small exact OS+ controls for the proved corollary X ▷ P_a.
    osplus_word = ordinal(two_cycle, one)
    no_sink_contexts = tuple(data for data in contexts if not base.has_sink(data.out))
    need(len(no_sink_contexts) == 2, "labelled no-sink order-three controls")
    osplus_checks = 0
    for right in no_sink_contexts:
        need(base.remainder(osplus_word, right) > 0, "OS+ corollary control")
        osplus_checks += 1

    coefficients = tuple(word_floor(m) for m in range(1, 7))
    need(
        coefficients
        == (10764, 107640, 979524, 8826480, 79449084, 715052520),
        "geometric floor coefficients",
    )

    lines = [
        "theorem=THM-4213",
        "symbolic_contextual_telescope=PASS",
        "symbolic_weighted_floor_law=alpha*H(D)^2+beta",
        f"labelled_contexts_through_order_three={len(contexts)}",
        f"exact_telescope_checks={telescope_checks}",
        f"tail_five_floor_checks={tail_floor_checks}",
        f"two_cycle_floor_checks={two_cycle_checks}",
        f"transitively_padded_word_checks={padded_checks}",
        "geometric_floor_formula=10764*(9^m-1)/8",
        "geometric_floor_coefficients_m1_to_m6="
        + ",".join(str(value) for value in coefficients),
        f"two_cycle_order={len(two_cycle.out)}",
        f"two_cycle_hamilton={two_cycle.hamilton}",
        f"two_cycle_singleton_defect={defect(two_cycle, one, one)}",
        f"good_suffix_hostile_order={len(hostile.out)}",
        f"good_suffix_hostile_hamilton={hostile.hamilton}",
        f"good_suffix_hostile_singleton_defect={hostile_value}",
        f"osplus_no_sink_controls={osplus_checks}",
        "equality_boundary=m=1,n1=5,no_transitive_padding,B=C=P1",
    ]
    payload = "\n".join(lines) + "\n"
    print(payload, end="")
    print(f"semantic_sha256={hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
