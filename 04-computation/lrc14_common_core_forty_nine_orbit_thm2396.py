#!/usr/bin/env python3
"""Exact companion for THM-2396.

The script exhausts a relaxation of every generic 49-root orbit in the
THM-2393 common-core no-clean branch.  All words are finite subsets of
Z/49Z.  Translates which share one physical orbit are deliberately allowed
to vary independently, so emptiness of this larger universe is a valid
obstruction to the physical branch.
"""

from collections import Counter, defaultdict
from itertools import combinations
from math import gcd


MODULUS = 49
ALL_SEVEN = (1 << 7) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ordinary_word(unit: int, start: int) -> tuple[int, ...]:
    """Vertical graph of a seven-term unit comb word on Z/49Z.

    The returned value at residue r is the unique s such that
    j=r+7s lies in

        {j : unit*j-start in {0,...,6} mod 49}.
    """

    require(gcd(unit, 7) == 1, "ordinary word needs a seven-unit")
    inverse = pow(unit, -1, MODULUS)
    values = [-1] * 7
    for offset in range(7):
        index = inverse * (start + offset) % MODULUS
        residue = index % 7
        require(values[residue] == -1, "ordinary word repeated one bin")
        values[residue] = index // 7
    require(all(value >= 0 for value in values), "ordinary word missed a bin")
    return tuple(values)


def guard_word(unit: int, start: int) -> tuple[tuple[int, int], ...]:
    """Two-sheet graph of a fourteen-term guard word on Z/49Z."""

    require(gcd(unit, 7) == 1, "guard word needs a seven-unit")
    inverse = pow(unit, -1, MODULUS)
    values: list[list[int]] = [[] for _ in range(7)]
    for offset in range(14):
        index = inverse * (start + offset) % MODULUS
        values[index % 7].append(index // 7)
    require(all(len(value) == 2 for value in values), "guard sheet count")
    return tuple(tuple(sorted(value)) for value in values)


ordinary_bank = tuple(
    sorted(
        {
            ordinary_word(unit, start)
            for unit in range(1, MODULUS)
            if unit % 7
            for start in range(MODULUS)
        }
    )
)
guard_bank = tuple(
    sorted(
        {
            guard_word(unit, start)
            for unit in range(1, MODULUS)
            if unit % 7
            for start in range(MODULUS)
        }
    )
)

# The raw (unit,start) parameter banks each have 42*49 elements.  Reversal
# of a consecutive block identifies pairs, leaving 1,029 distinct words.
require(len(ordinary_bank) == 1029, "ordinary word-bank size")
require(len(guard_bank) == 1029, "guard word-bank size")


# Normalize the quotient C1 word to the horizontal graph A(r)=0.  The
# correlated common-core words then have speeds 13 and 169=22 mod 49.
A_WORD = ordinary_word(1, 0)
B_BANK = tuple(sorted({ordinary_word(13, start) for start in range(49)}))
C_BANK = tuple(sorted({ordinary_word(22, start) for start in range(49)}))
require(A_WORD == (0,) * 7, "A normalization")
require(len(B_BANK) == len(C_BANK) == 49, "fixed-speed translate bank")


candidate_count_distribution: Counter[int] = Counter()
packing_count_distribution: Counter[int] = Counter()
required_c_rows_distribution: Counter[int] = Counter()

relevant_guard_top_pairs = 0
q_packings = 0
b_packing_tests = 0
nonzero_hole_rejections = 0
c_capacity_boundaries = 0
c_word_tests = 0
survivors = 0


def four_q_packings(
    candidates: tuple[tuple[int, ...], ...],
    guard: tuple[tuple[int, int], ...],
    safe_rows: tuple[int, ...],
) -> tuple[
    tuple[tuple[tuple[int, ...], ...], tuple[int, ...]], ...
]:
    """Choose four q graphs disjoint from the guard and one another.

    On every safe row the guard contributes two addresses.  A legal
    four-q choice contributes four further distinct addresses, leaving
    the unique hole returned with the packing.
    """

    base_occupancy = [
        (1 << guard[row][0]) | (1 << guard[row][1])
        for row in safe_rows
    ]
    out = []

    def recurse(
        first_index: int,
        chosen: tuple[tuple[int, ...], ...],
        occupancy: tuple[int, ...],
    ) -> None:
        if len(chosen) == 4:
            holes = []
            for occupied in occupancy:
                missing = ALL_SEVEN ^ occupied
                require(missing.bit_count() == 1, "safe row is not transversal")
                holes.append(missing.bit_length() - 1)
            out.append((chosen, tuple(holes)))
            return

        for index in range(first_index, len(candidates)):
            word = candidates[index]
            updated = list(occupancy)
            legal = True
            for row_index, row in enumerate(safe_rows):
                bit = 1 << word[row]
                if updated[row_index] & bit:
                    legal = False
                    break
                updated[row_index] |= bit
            if legal:
                recurse(index + 1, chosen + (word,), tuple(updated))

    recurse(0, (), tuple(base_occupancy))
    return tuple(out)


for top_row in range(7):
    safe_rows = tuple(row for row in range(7) if row != top_row)

    for guard in guard_bank:
        # A and B must partition the two guard addresses on the top q*
        # bin.  Since A(top)=0, the guard must contain zero there.
        if 0 not in guard[top_row]:
            continue

        other_top_address = next(
            iter(set(guard[top_row]) - {0})
        )
        b_words = tuple(
            word for word in B_BANK if word[top_row] == other_top_address
        )
        require(len(b_words) == 7, "wrong prescribed-top B bank")
        relevant_guard_top_pairs += 1

        # THM-2391 puts every lower q point inside the guard on the top
        # bin.  THM-2394 makes the six remaining bins transversals, so a
        # lower q point must avoid both guard addresses there.
        candidates = tuple(
            word
            for word in ordinary_bank
            if word[top_row] in guard[top_row]
            and all(word[row] not in guard[row] for row in safe_rows)
        )
        candidate_count_distribution[len(candidates)] += 1

        packings = four_q_packings(candidates, guard, safe_rows)
        packing_count_distribution[len(packings)] += 1
        q_packings += len(packings)

        for _, holes in packings:
            for b_word in b_words:
                b_packing_tests += 1

                # In normalized coordinates THM-2394's hole trichotomy is
                #
                #   hole=B, or hole=A=C=0.
                #
                # A nonzero hole different from B is immediately illegal.
                nonzero_mismatches = tuple(
                    row
                    for row, hole in zip(safe_rows, holes)
                    if hole != b_word[row] and hole != 0
                )
                if nonzero_mismatches:
                    nonzero_hole_rejections += 1
                    continue

                required_c_rows = tuple(
                    row
                    for row, hole in zip(safe_rows, holes)
                    if hole == 0 and hole != b_word[row]
                )
                required_c_rows_distribution[len(required_c_rows)] += 1
                c_capacity_boundaries += 1

                # The final relaxed check allows C every independent
                # speed-22 translate.  It must agree with A=0 on every
                # required row.
                local_survivors = 0
                for c_word in C_BANK:
                    c_word_tests += 1
                    if all(c_word[row] == 0 for row in required_c_rows):
                        local_survivors += 1
                survivors += local_survivors


require(relevant_guard_top_pairs == 2058, "relevant guard/top count")
require(
    candidate_count_distribution
    == Counter({22: 588, 25: 588, 26: 294, 33: 588}),
    "q-candidate distribution",
)
require(
    packing_count_distribution == Counter({0: 1176, 1: 294, 3: 588}),
    "four-q packing distribution",
)
require(q_packings == 2058, "four-q packing total")
require(b_packing_tests == 14406, "B/packing test total")
require(nonzero_hole_rejections == 14386, "first rejection count")
require(c_capacity_boundaries == 20, "C-capacity boundary count")
require(
    required_c_rows_distribution == Counter({4: 18, 5: 2}),
    "required C-row distribution",
)
require(c_word_tests == 980, "C-word test count")


# Independent capacity audit for the last twenty cases.  A C translate can
# agree with A on at most two rows, while every residual case requires four
# or five.  The detailed mask bank is retained rather than inferred from the
# final zero count.
c_zero_masks: defaultdict[tuple[int, ...], int] = defaultdict(int)
for c_word in C_BANK:
    zero_rows = tuple(row for row in range(7) if c_word[row] == 0)
    c_zero_masks[zero_rows] += 1

require(len(c_zero_masks) == 13, "C zero-mask bank size")
require(max(map(len, c_zero_masks)) == 2, "C/A agreement exceeded two rows")
require(
    Counter(
        {
            size: sum(
                multiplicity
                for rows, multiplicity in c_zero_masks.items()
                if len(rows) == size
            )
            for size in (0, 1, 2)
        }
    )
    == Counter({0: 10, 1: 29, 2: 10}),
    "C zero-mask multiplicity distribution",
)
require(survivors == 0, "relaxed 49-orbit survivor")


print("theorem=THM-2396")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print("modulus=49; bins=7; addresses_per_bin=7")
print(
    f"ordinary_words={len(ordinary_bank)};"
    f" guard_words={len(guard_bank)};"
    f" B_translates={len(B_BANK)}; C_translates={len(C_BANK)}"
)
print(f"relevant_guard_top_pairs={relevant_guard_top_pairs}")
print(
    "q_candidate_distribution="
    + ",".join(
        f"{count}:{candidate_count_distribution[count]}"
        for count in sorted(candidate_count_distribution)
    )
)
print(
    "q_packing_distribution="
    + ",".join(
        f"{count}:{packing_count_distribution[count]}"
        for count in sorted(packing_count_distribution)
    )
)
print(f"q_packings={q_packings}; B_packing_tests={b_packing_tests}")
print(
    f"nonzero_hole_rejections={nonzero_hole_rejections};"
    f" C_capacity_boundaries={c_capacity_boundaries}"
)
print(
    "required_C_rows="
    + ",".join(
        f"{count}:{required_c_rows_distribution[count]}"
        for count in sorted(required_c_rows_distribution)
    )
)
print(
    f"C_zero_masks={len(c_zero_masks)};"
    f" max_A_C_agreement={max(map(len, c_zero_masks))};"
    f" C_word_tests={c_word_tests}"
)
print("relaxed_49_orbit_survivors=0")
print("consequence=delta_zero_common_core_branch_empty")
print("row_decrement=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
