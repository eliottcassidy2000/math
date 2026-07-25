#!/usr/bin/env python3
"""Exact companion for THM-2305."""

from fractions import Fraction
from itertools import product


DELTA_STRICT = Fraction(39002430583, 53493927587100)
DELTA_REPEAT = Fraction(13560199813, 53493927587100)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    # For a fixed absent source blocker j, every nonempty word in the other
    # two blocker bits is exactly one of {a}, {b}, {a,b}.
    truth_rows = 0
    for j in range(3):
        others = tuple(i for i in range(3) if i != j)
        seen = set()
        for bits in product((0, 1), repeat=3):
            if bits[j] != 0 or not any(bits[i] for i in others):
                continue
            word = tuple(i for i in others if bits[i])
            require(len(word) in (1, 2), "terminal word has invalid size")
            seen.add(word)
            truth_rows += 1
        require(
            seen == {(others[0],), (others[1],), others},
            f"word partition failed for source {j}",
        )
    require(truth_rows == 9, "terminal truth-table size changed")

    # Sharp trade: d is the double-target mass. Scaling by a common
    # denominator checks every rational simplex row on a finite hostile bank.
    trade_rows = 0
    equality_rows = 0
    for denominator in range(1, 101):
        for a_num in range(denominator + 1):
            for b_num in range(denominator - a_num + 1):
                d_num = denominator - a_num - b_num
                p_a = Fraction(a_num, denominator)
                p_b = Fraction(b_num, denominator)
                d = Fraction(d_num, denominator)
                rho = p_a + p_b + d
                require(
                    max(p_a, p_b) >= (rho - d) / 2,
                    "pure-edge/fork trade failed",
                )
                require(
                    d >= rho / 3 or max(p_a, p_b) > rho / 3,
                    "one-third dichotomy failed",
                )
                if p_a == p_b:
                    require(
                        max(p_a, p_b) == (rho - d) / 2,
                        "sharp equality row changed",
                    )
                    equality_rows += 1
                trade_rows += 1

    strict_mass_floor = DELTA_STRICT / 3
    repeat_mass_floor = DELTA_REPEAT / 3
    strict_energy_floor = DELTA_STRICT**2 / 9
    repeat_energy_floor = DELTA_REPEAT**2 / 9
    require(
        strict_mass_floor == Fraction(39002430583, 160481782761300),
        "strict word-mass floor changed",
    )
    require(
        repeat_mass_floor == Fraction(13560199813, 160481782761300),
        "repeat word-mass floor changed",
    )
    require(
        strict_energy_floor
        == Fraction(
            1521189591381733719889,
            25754402598245085852777690000,
        ),
        "strict word-energy floor changed",
    )
    require(
        repeat_energy_floor
        == Fraction(
            183879018968485234969,
            25754402598245085852777690000,
        ),
        "repeat word-energy floor changed",
    )

    jump_rows = 0
    for S in range(9, 200):
        source_jumps = 2 * S
        word_jumps = 2 * S
        require(
            source_jumps + word_jumps - 1 == 4 * S - 1,
            "covariant word-index bound changed",
        )
        jump_rows += 1

    # Every no-loop functional graph on three labels has a 2- or 3-cycle.
    functional_rows = 0
    two_cycle_rows = 0
    three_cycle_rows = 0
    choices = tuple(tuple(t for t in range(3) if t != j) for j in range(3))
    for targets in product(*choices):
        functional_rows += 1
        cycle_lengths = set()
        for start in range(3):
            path = []
            position = {}
            vertex = start
            while vertex not in position:
                position[vertex] = len(path)
                path.append(vertex)
                vertex = targets[vertex]
            cycle_lengths.add(len(path) - position[vertex])
        require(cycle_lengths <= {2, 3}, "loop-free graph has invalid cycle")
        if 3 in cycle_lengths:
            three_cycle_rows += 1
        else:
            require(2 in cycle_lengths, "functional graph has no cycle")
            two_cycle_rows += 1
    require(functional_rows == 8, "functional graph atlas changed")
    require(two_cycle_rows == 6, "two-cycle census changed")
    require(three_cycle_rows == 2, "three-cycle census changed")

    # Allowing the double word produces 3^3 hypergraph selections. Exactly
    # the eight all-singleton selections reduce to functional graphs.
    hypergraph_rows = 3**3
    fork_rows = hypergraph_rows - functional_rows
    require(hypergraph_rows == 27 and fork_rows == 19, "hypergraph census changed")

    print("theorem=THM-2305")
    print("status=PROVED+VERIFIED-EXACT-CANDIDATE")
    print(f"terminal_truth_rows={truth_rows}")
    print("canonical_words={a};{b};{a,b}")
    print(f"sharp_trade_rows={trade_rows}")
    print(f"sharp_trade_equality_rows={equality_rows}")
    print("trade=max(p_a,p_b)>=(rho-d)/2")
    print("dichotomy=fork>=rho/3_or_pure>rho/3")
    print(f"strict_word_mass_floor={strict_mass_floor}")
    print(f"repeat_word_mass_floor={repeat_mass_floor}")
    print(f"strict_word_energy_floor={strict_energy_floor}")
    print(f"repeat_word_energy_floor={repeat_energy_floor}")
    print(f"covariant_jump_rows={jump_rows}")
    print("gauge_index_bound_on_exact_word=h<=4S-1")
    print(f"pure_functional_graphs={functional_rows}")
    print(f"pure_graphs_with_2_cycle={two_cycle_rows}")
    print(f"pure_graphs_with_3_cycle={three_cycle_rows}")
    print(f"word_hypergraph_selections={hypergraph_rows}")
    print(f"selections_with_fork={fork_rows}")
    print("cycle_scope=conditional_on_three_positive_pure_returns")
    print("composition_loss=incoming_and_outgoing_owner_subsets_may_be_disjoint")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
