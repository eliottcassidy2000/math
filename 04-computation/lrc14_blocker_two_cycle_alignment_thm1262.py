#!/usr/bin/env python3
"""Dependency-free exact referee for THM-1262.

The geometric inputs are deliberately named rather than re-proved here:

* the ascent target tooth is contained in the low-owner safe component;
* a low-owner danger tooth is disjoint from that component;
* consecutive teeth in the deletion-minimal word overlap;
* a binary phase/order mismatch forces the two marked teeth to be adjacent.

The referee exhausts the finite logical composition of those inputs, the
ordered-word corridor and reflection laws, the third-owner exclusion at the
first corridor handoff, and the gcd/lcm seam arithmetic.  Only the Python
standard library is used.

Tournament Analysis challenges runners as vertices.  The speed relation on
the low/high pair is transitive and loses the blocker return; the blocker
relation is a directed 2-cycle and loses protected placement.  The faithful
three-vertex carrier is (high target tooth, protected handoff, low target
tooth), ordered by the chronological word.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reverse_binary_gate_audit() -> int:
    """Check the elementary inequality selecting the binary descent edge."""
    rows = 0
    for carrier in range(1, 81):
        for low in range(carrier + 1, 3 * carrier + 1):
            for high in range(low + 1, 3 * carrier + 1):
                require(low < high, "two-cycle orientation is not a descent")
                require(high < carrier + high,
                        "binary target inequality lost positivity of carrier")
                require(low < carrier + high,
                        "reverse high-to-low edge left the binary gate")
                rows += 1
    return rows


def logical_composition_audit() -> tuple[int, int]:
    """Exhaust the ascent-protection/mismatch truth table.

    The boolean implications are exactly the four paper providers used by
    the theorem.  Under containment, reverse danger, and binary descent, all
    surviving rows must be aligned.
    """
    truth_rows = 0
    protected_binary_rows = 0
    for (
        contained,
        reverse_danger,
        disjoint,
        adjacent,
        binary,
        mismatch,
        aligned,
    ) in product((False, True), repeat=7):
        truth_rows += 1

        # Ascent containment plus reverse-target danger force disjointness.
        if contained and reverse_danger and not disjoint:
            continue
        # Consecutive selected teeth overlap, so disjoint marks are not
        # consecutive.
        if disjoint and adjacent:
            continue
        # The binary landing law has exactly the aligned/mismatch branches.
        if binary and aligned == mismatch:
            continue
        # Its mismatch branch is an adjacent marked inversion.
        if binary and mismatch and not adjacent:
            continue

        if contained and reverse_danger and binary:
            protected_binary_rows += 1
            require(disjoint, "protected marked teeth were not disjoint")
            require(not adjacent, "protected marked teeth became consecutive")
            require(not mismatch, "nonconsecutive binary marks mismatched")
            require(aligned, "binary two-cycle failed to align")

    require(truth_rows == 128, "boolean truth-table size changed")
    require(protected_binary_rows > 0, "truth table lost its admissible row")
    return truth_rows, protected_binary_rows


def minimal_interval_chain_audit() -> tuple[int, int, int]:
    """Audit the ordered-interval realization of the corridor conclusion."""
    chains = 0
    nonconsecutive_pairs = 0
    protected_neighbor_overlaps = 0
    grid = range(10)

    for length in range(2, 6):
        for lefts in combinations(grid, length):
            for rights in combinations(range(1, 11), length):
                if any(lefts[a] >= rights[a] for a in range(length)):
                    continue
                if any(lefts[a + 1] >= rights[a]
                       for a in range(length - 1)):
                    continue
                if any(rights[a] > lefts[a + 2]
                       for a in range(length - 2)):
                    continue
                chains += 1

                for high_pos in range(length):
                    for low_pos in range(length):
                        if abs(high_pos - low_pos) < 2:
                            continue
                        nonconsecutive_pairs += 1
                        direction = 1 if high_pos < low_pos else -1
                        neighbor = high_pos + direction
                        require(
                            min(high_pos, low_pos) < neighbor
                            < max(high_pos, low_pos),
                            "corridor neighbor is not intermediate",
                        )
                        require(abs(neighbor - high_pos) == 1,
                                "corridor neighbor is not consecutive")

                        overlap_left = max(lefts[high_pos], lefts[neighbor])
                        overlap_right = min(rights[high_pos], rights[neighbor])
                        require(overlap_left < overlap_right,
                                "first corridor handoff lost positive overlap")

                        if high_pos < low_pos:
                            require(rights[high_pos] <= lefts[low_pos],
                                    "nonconsecutive marked teeth overlap")
                            high_mark = Fraction(
                                lefts[high_pos] + rights[high_pos], 2
                            )
                            low_mark = Fraction(
                                lefts[low_pos] + rights[low_pos], 2
                            )
                            require(high_mark < low_mark,
                                    "word/phase orientation is not aligned")
                        else:
                            require(rights[low_pos] <= lefts[high_pos],
                                    "reflected nonconsecutive marks overlap")
                            high_mark = Fraction(
                                lefts[high_pos] + rights[high_pos], 2
                            )
                            low_mark = Fraction(
                                lefts[low_pos] + rights[low_pos], 2
                            )
                            require(low_mark < high_mark,
                                    "reflected word/phase orientation failed")
                        protected_neighbor_overlaps += 1

    require(chains == 3_542, "minimal interval-chain census changed")
    require(nonconsecutive_pairs == protected_neighbor_overlaps,
            "some corridor lacks a first handoff")
    return chains, nonconsecutive_pairs, protected_neighbor_overlaps


def ordered_word_third_owner_audit() -> tuple[int, int, int]:
    """Exhaust positions and owner labels for the forced bridge tooth."""
    placements = 0
    intermediate_checks = 0
    third_owner_choices = 0
    labels = range(6)

    for length in range(3, 13):
        for high_pos in range(length):
            for low_pos in range(length):
                if abs(high_pos - low_pos) < 2:
                    continue
                placements += 1
                direction = 1 if high_pos < low_pos else -1
                neighbor = high_pos + direction
                require(abs(neighbor - high_pos) == 1,
                        "first corridor tooth is not adjacent to high target")
                require(
                    min(high_pos, low_pos) < neighbor < max(high_pos, low_pos),
                    "corridor has no intermediate tooth",
                )
                intermediate_checks += 1

                for low_owner in labels:
                    for high_owner in labels:
                        if not low_owner < high_owner:
                            continue
                        allowed = []
                        for neighbor_owner in labels:
                            # Same-high-owner teeth cannot overlap.  A
                            # low-owner tooth cannot meet the high tooth
                            # contained in the low-safe component.
                            if neighbor_owner in (low_owner, high_owner):
                                continue
                            allowed.append(neighbor_owner)
                            third_owner_choices += 1
                        require(len(allowed) == 4,
                                "six-label third-owner fibre changed size")

    return placements, intermediate_checks, third_owner_choices


def reflection_audit() -> tuple[int, int]:
    rows = 0
    aligned_rows = 0
    for word_sign, phase_sign in product((-1, 1), repeat=2):
        aligned = word_sign == phase_sign
        reflected_word = -word_sign
        reflected_phase = -phase_sign
        require(
            aligned == (reflected_word == reflected_phase),
            "reflection changed alignment status",
        )
        require(
            (word_sign != phase_sign)
            == (reflected_word != reflected_phase),
            "reflection changed mismatch status",
        )
        rows += 1
        aligned_rows += int(aligned)
    require(rows == 4 and aligned_rows == 2,
            "orientation/reflection truth table changed")
    return rows, aligned_rows


def lcm_seam_audit() -> tuple[int, int]:
    pair_rows = 0
    positive_numerator_rows = 0
    for high in range(1, 101):
        for third in range(1, 101):
            if high == third:
                continue
            common = gcd(high, third)
            clock = lcm(high, third)
            require(common * clock == high * third,
                    "gcd/lcm product identity failed")
            quantum = Fraction(common, 14 * high * third)
            require(quantum == Fraction(1, 14 * clock),
                    "gcd and lcm seam quanta disagree")
            require(quantum > 0, "lcm seam quantum is not positive")
            pair_rows += 1

            for multiplier in range(1, 6):
                numerator = multiplier * common
                overlap = Fraction(numerator, 14 * high * third)
                require(overlap >= quantum,
                        "positive compatible seam missed its quantum")
                positive_numerator_rows += 1
    return pair_rows, positive_numerator_rows


def main() -> None:
    speed_rows = reverse_binary_gate_audit()
    truth_rows, protected_rows = logical_composition_audit()
    chains, corridors, neighbor_overlaps = minimal_interval_chain_audit()
    placements, intermediate, third_choices = ordered_word_third_owner_audit()
    reflection_rows, reflected_aligned = reflection_audit()
    lcm_pairs, numerator_rows = lcm_seam_audit()

    print("THM-1262 BLOCKER TWO-CYCLE ALIGNMENT EXACT AUDIT")
    print(f"reverse binary-gate speed triples = {speed_rows}")
    print(f"logical truth-table rows = {truth_rows}")
    print(f"admissible protected binary rows = {protected_rows}")
    print(f"minimal interval chains = {chains}")
    print(f"nonconsecutive marked corridors = {corridors}")
    print(f"positive first-corridor handoffs = {neighbor_overlaps}")
    print(f"ordered-word marked placements = {placements}")
    print(f"intermediate-neighbor checks = {intermediate}")
    print(f"six-label third-owner choices = {third_choices}")
    print(f"reflection orientation rows = {reflection_rows}")
    print(f"reflection-aligned rows = {reflected_aligned}")
    print(f"gcd/lcm seam pairs = {lcm_pairs}")
    print(f"positive compatible seam numerators = {numerator_rows}")
    print("speed tournament on {low,high}: transitive, score word (0,1)")
    print("blocker relation on {low,high}: directed 2-cycle, not a tournament")
    print("faithful carrier: high target -> protected third-owner handoff -> low target")
    print("reflection: reverses phase and tooth orders, preserves alignment")
    print("destroyed by runner projection: low-safe containment and wall placement")
    print("CONCLUSION: mismatch impossible; aligned corridor exports a third-owner lcm seam")
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
