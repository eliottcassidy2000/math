#!/usr/bin/env python3
"""Exact relaxed-orbit dual certificate for THM-2417.

The script uses only integer arithmetic and Python bit sets (stored as ints).
It reconstructs the full relaxed mask families rather than trusting the
offline optimization that found the displayed dual weights.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple


N = 343
WORD_LENGTH = 49
INV13 = pow(13, -1, N)
FULL = (1 << N) - 1
GUARD = (1 << 98) - 1
UNITS = tuple(d for d in range(1, N) if d % 7)
ORIENTED_STEPS = tuple(d for d in UNITS if d < N - d)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def progression_mask(start: int, step: int) -> int:
    mask = 0
    for j in range(WORD_LENGTH):
        mask |= 1 << ((start + j * step) % N)
    return mask


def points(mask: int) -> Iterable[int]:
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask -= bit


def parse_table(text: str) -> List[List[int]]:
    rows = [
        [int(entry) for entry in line.split()]
        for line in text.strip().splitlines()
    ]
    require(all(len(row) == 35 for row in rows), "bad certificate row")
    return rows


# Rows are indexed by the safe residues in increasing order; columns are
# j=14,...,48 in x=r+7j.  Reflections supply t=6,5,4 from t=0,1,2.
PARITY_T0 = parse_table(
    """
7 0 2 2 0 18 3 7 0 8 10 8 18 9 13 1 8 10 14 18 13 11 13 4 2 11 10 13 1 13 2 2 6 5 13
13 4 7 8 10 8 6 13 11 20 22 13 14 11 14 4 22 22 13 12 10 14 4 22 16 13 12 10 12 3 15 13 0 10 8
9 7 0 1 5 7 2 12 10 12 13 5 7 9 11 2 8 7 3 6 6 12 16 13 10 11 13 10 12 10 3 5 7 2 10
9 10 3 11 0 0 3 9 16 13 16 4 12 3 12 5 9 14 14 11 3 13 13 13 20 16 12 6 14 3 5 0 0 5 0
7 10 12 20 16 17 6 0 10 13 17 18 14 12 2 10 13 17 24 15 19 12 13 9 11 22 10 18 7 13 3 8 10 4 14
7 11 6 7 0 4 5 13 11 15 13 2 14 0 4 21 18 14 12 12 7 3 21 18 8 12 8 7 3 14 11 0 12 1 1
"""
)

PARITY_T1_ROW = [
    35, 36, 35, 30, 24, 35, 35, 50, 57, 62, 64, 52, 54, 49, 51, 56, 65,
    72, 62, 55, 48, 50, 56, 59, 70, 62, 54, 45, 36, 38, 29, 34, 32, 38, 33,
]
PARITY_T1 = [PARITY_T1_ROW[:] for _ in range(6)]

PARITY_T2 = parse_table(
    """
7 7 7 7 6 8 8 7 7 9 10 9 8 9 9 9 12 13 13 10 10 8 9 9 9 10 10 10 7 8 4 7 7 7 8
9 8 6 5 6 9 8 8 8 9 10 7 8 8 10 10 12 13 13 10 9 7 9 8 10 9 8 8 8 9 6 6 6 7 8
8 8 6 7 6 9 8 8 8 9 10 7 8 8 11 10 12 13 11 10 9 8 8 7 9 7 8 8 8 9 6 5 7 7 8
7 7 7 8 6 8 6 10 10 9 11 9 8 9 9 9 13 13 12 9 9 8 7 9 10 11 8 7 8 8 7 7 6 7 7
7 7 5 6 5 8 6 9 9 9 11 8 8 8 8 9 10 13 11 9 8 8 8 8 11 10 7 9 8 9 6 7 6 7 9
8 8 7 7 5 9 8 9 9 10 10 9 8 8 10 8 12 11 12 10 8 9 8 9 10 10 8 6 8 8 6 5 6 7 9
"""
)

PARITY_T3 = parse_table(
    """
8 8 8 7 6 9 9 10 9 11 11 8 9 9 9 10 12 14 12 10 10 8 9 9 11 10 9 10 7 9 6 6 6 7 8
7 8 8 8 7 9 9 7 8 11 10 9 8 9 8 10 14 14 14 10 10 10 10 9 11 10 11 11 7 9 7 8 5 7 7
9 8 7 7 5 11 9 9 10 8 10 8 11 10 10 12 12 14 13 11 11 8 9 8 10 9 8 8 8 9 6 7 6 8 9
9 8 7 7 6 9 9 8 9 9 10 8 9 8 11 12 13 14 12 11 11 9 9 8 11 9 9 9 9 11 6 6 6 8 9
6 7 6 8 7 9 7 11 11 10 11 10 10 10 11 10 14 14 14 10 8 9 8 9 10 11 8 7 9 9 7 8 8 7 7
8 7 6 7 6 9 7 10 10 10 12 9 9 8 10 10 12 14 12 11 9 9 9 8 11 11 9 10 9 9 6 7 8 8 9
"""
)

BASE_PARITY_TABLES = {
    0: PARITY_T0,
    1: PARITY_T1,
    2: PARITY_T2,
    3: PARITY_T3,
}

CONTIGUOUS_T2 = parse_table(
    """
2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 3 4 4 4 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 2
3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 3 4 4 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 2
2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3
2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 4 4 4 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3
2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 4 4 4 3 4 4 4 4 4 4 4 4 3 3 3 3 3 4 4
2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 4 4 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3
"""
)

CONTIGUOUS_T3 = parse_table(
    """
3 2 2 2 2 2 2 3 3 3 3 2 2 3 3 2 3 3 2 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 1
2 2 2 2 2 2 2 2 2 3 3 3 3 3 2 2 3 3 3 2 3 3 3 3 3 3 3 3 2 2 2 2 2 3 1
2 2 2 2 2 2 2 3 3 3 3 2 3 3 3 2 2 3 2 2 3 3 3 2 3 3 3 3 2 2 2 2 2 2 1
1 2 2 2 2 2 2 2 3 3 3 2 3 3 3 2 3 3 2 2 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2
1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 2 2 3 2 2 3 3 2 3 3 3 3 3 2 2 2 2 2 2 2
1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 2 3 3 3 2 3 3 3 2 3 3 3 3 2 2 2 2 2 2 3
"""
)


def top_mask(t: int) -> int:
    return sum(1 << (t + 7 * j) for j in range(49))


def support_mask(t: int, partition: str, side: int) -> int:
    if partition == "contiguous":
        js = range(0, 7) if side == 0 else range(7, 14)
    elif partition == "parity":
        js = range(0, 14, 2) if side == 0 else range(1, 14, 2)
    else:
        raise RuntimeError("unknown partition")
    return sum(1 << (t + 7 * j) for j in js)


def table_weights(t: int, table: Sequence[Sequence[int]]) -> Dict[int, int]:
    safe_rows = [r for r in range(7) if r != t]
    require(len(table) == len(safe_rows), "bad number of safe rows")
    return {
        r + 7 * j: table[row_index][j - 14]
        for row_index, r in enumerate(safe_rows)
        for j in range(14, 49)
    }


def reflect_weights(weights: Dict[int, int]) -> Dict[int, int]:
    # x -> 97-x preserves the guard and sends Q_t to Q_(6-t).
    return {(97 - x) % N: value for x, value in weights.items()}


def parity_weights(t: int) -> Dict[int, int]:
    if t <= 3:
        return table_weights(t, BASE_PARITY_TABLES[t])
    source = 6 - t
    return reflect_weights(table_weights(source, BASE_PARITY_TABLES[source]))


def contiguous_weights(t: int) -> Dict[int, int]:
    if t == 2:
        return table_weights(t, CONTIGUOUS_T2)
    if t == 3:
        return table_weights(t, CONTIGUOUS_T3)
    if t == 4:
        return reflect_weights(table_weights(2, CONTIGUOUS_T2))
    qmask = top_mask(t)
    return {
        x: 1
        for x in range(98, N)
        if not ((qmask >> x) & 1) and x // 7 < 48
    }


def mask_weight(mask: int, weights: Dict[int, int]) -> int:
    return sum(weights.get(x, 0) for x in points(mask))


def valuation(n: int, p: int) -> int:
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def strict_sign(v: int, width: int, s: int) -> int:
    denominator = 49 * 1009
    residue = (v * (1001 + 1009 * s)) % denominator
    distance = min(residue, denominator - residue)
    return 14 * distance - width * denominator


def main() -> None:
    all_progressions = {
        progression_mask(start, step)
        for step in ORIENTED_STEPS
        for start in range(N)
    }
    require(len(ORIENTED_STEPS) == 147, "oriented step count")
    require(len(all_progressions) == 50421, "progression mask count")
    require(all(mask.bit_count() == 49 for mask in all_progressions), "mask size")

    # Translation of the support changes starts, not steps, so t=0 suffices
    # to reconstruct the four quotient-step families.
    quotient_steps: Dict[Tuple[str, int], Tuple[int, ...]] = {}
    blocker_masks: Dict[Tuple[str, int], Tuple[int, ...]] = {}
    q0 = top_mask(0)
    for partition in ("contiguous", "parity"):
        for side in (0, 1):
            support = support_mask(0, partition, side)
            steps = tuple(
                d
                for d in UNITS
                if any(
                    progression_mask(start, d) & q0 == support
                    for start in range(N)
                )
            )
            physical = {
                progression_mask(start, (INV13 * d) % N)
                for d in steps
                for start in range(N)
            }
            require(len(steps) == 14, "quotient blocker step count")
            require(len(physical) == 2401, "physical blocker family count")
            quotient_steps[(partition, side)] = steps
            blocker_masks[(partition, side)] = tuple(physical)

    require(
        {d % 49 for d in quotient_steps[("contiguous", 0)]} == {1, 48},
        "contiguous quotient slopes",
    )
    require(
        {d % 49 for d in quotient_steps[("parity", 0)]} == {2, 47},
        "parity quotient slopes",
    )
    require(
        quotient_steps[("contiguous", 0)]
        == quotient_steps[("contiguous", 1)],
        "contiguous side step equality",
    )
    require(
        quotient_steps[("parity", 0)] == quotient_steps[("parity", 1)],
        "parity side step equality",
    )

    print("THM-2417 middle-depth-two relaxed orbit dual")
    print(f"N={N} progression_masks={len(all_progressions)}")
    print(
        "quotient_steps=14 per side; physical_blocker_masks=2401 per side"
    )

    expected_parity = {
        0: (1998, 339, 320, 320),
        1: (9978, 1702, 1584, 1584),
        2: (1765, 301, 280, 280),
        3: (1911, 326, 303, 303),
        4: (1765, 301, 280, 280),
        5: (9978, 1702, 1584, 1584),
        6: (1998, 339, 320, 320),
    }
    expected_contiguous = {
        0: (204, 35, 30, 30),
        1: (204, 35, 30, 30),
        2: (735, 127, 113, 113),
        3: (513, 89, 78, 78),
        4: (735, 127, 113, 113),
        5: (204, 35, 30, 30),
        6: (204, 35, 30, 30),
    }

    for t in range(7):
        qmask = top_mask(t)
        universe = FULL & ~GUARD & ~qmask
        require(universe.bit_count() == 210, "universe size")
        require((GUARD & qmask).bit_count() == 14, "top guard size")
        q_family = tuple(
            mask
            for mask in all_progressions
            if not (mask & qmask & ~GUARD)
        )
        require(len(q_family) == 490, "lower q family count")
        q_steps = {
            d
            for d in UNITS
            if any(
                not (progression_mask(start, d) & qmask & ~GUARD)
                for start in range(N)
            )
        }
        require(len(q_steps) == 28, "lower q step count")
        require({d % 49 for d in q_steps} == {1, 2, 47, 48}, "q slopes")

        for partition in ("contiguous", "parity"):
            weights = (
                contiguous_weights(t)
                if partition == "contiguous"
                else parity_weights(t)
            )
            require(set(weights) <= set(points(universe)), "weight domain")
            require(all(value >= 0 for value in weights.values()), "negative weight")
            total = sum(weights.values())
            q_cap = max(mask_weight(mask, weights) for mask in q_family)
            b_caps = tuple(
                max(
                    mask_weight(mask, weights)
                    for mask in blocker_masks[(partition, side)]
                )
                for side in (0, 1)
            )
            budget = 4 * q_cap + b_caps[0] + b_caps[1]
            require(
                total > budget,
                f"dual inequality did not close: {partition=} {t=} "
                f"{total=} {q_cap=} {b_caps=} {budget=}",
            )
            require(total <= 5 * q_cap + b_caps[0] + b_caps[1], "five-q hostile")
            if partition == "contiguous":
                require(
                    (total, q_cap, *b_caps) == expected_contiguous[t],
                    "contiguous certificate",
                )
            else:
                require(
                    (total, q_cap, *b_caps) == expected_parity[t],
                    "parity certificate",
                )
            print(
                f"{partition:10s} t={t} U=210 q_masks=490 "
                f"b_masks=2401,2401 total={total} "
                f"caps={q_cap},{b_caps[0]},{b_caps[1]} "
                f"budget={budget} margin={total-budget}"
            )

    # THM-2414's live W=8 atlas is an M=1 local packet, so it is outside
    # this M=2 universe.  Its displayed M=2 W=7 hostile fails the exact
    # full-bin premises at s=6 and s=13.
    require(valuation(7, 7) == 1, "W8 top depth")
    require(valuation(8281, 7) == 2, "W8 high depth")
    require(valuation(49, 7) == 2, "W7 top depth")
    require(valuation(57967, 7) == 3, "W7 high depth")
    require(strict_sign(49, 1, 6) < 0, "W7 s=6 top")
    require(strict_sign(8, 2, 6) < 0, "W7 s=6 guard")
    require(strict_sign(4, 1, 6) > 0, "W7 s=6 C1")
    require(strict_sign(8, 1, 6) > 0, "W7 s=6 C2")
    require(strict_sign(4459, 1, 6) > 0, "W7 s=6 C3")
    require(strict_sign(49, 1, 13) < 0, "W7 s=13 top")
    require(strict_sign(56, 1, 13) < 0, "W7 s=13 lower q")
    require(strict_sign(8, 2, 13) > 0, "W7 s=13 guard")
    require(strict_sign(4, 1, 13) > 0, "W7 s=13 C1")
    require(strict_sign(8, 1, 13) > 0, "W7 s=13 C2")
    require(strict_sign(4459, 1, 13) > 0, "W7 s=13 C3")
    print("control: THM-2414 W8 has M=1 and is outside the M=2 universe")
    print("control: THM-2414 W7 has M=2 but fails both full-bin premises")
    print("control: every stored certificate is specific to four lower q masks")
    print("VERIFIED")


if __name__ == "__main__":
    main()
