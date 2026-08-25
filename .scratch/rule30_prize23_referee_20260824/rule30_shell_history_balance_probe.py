#!/usr/bin/env python3
"""Finite-exact Rule 30 Prize-2/3 shell and suffix-history audit.

Scratch only.  The universal statements audited here are elementary finite-word
identities; every Rule 30 numerical statement stops at ``2**18`` center bits.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib


HORIZON = 1 << 18
EXPECTED_CENTER_SHA256 = (
    "75be8407c265798fa046baa3ba0f9336e4cbe27bff0be21aeba3e87a7681bea4"
)
HISTORY_DEPTHS = (0, 1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 30, 32, 33, 34)
EXPECTED_RIGHT_SPECIAL_CONTEXT = "011110100110000100110111001000101"
EXPECTED_RIGHT_SPECIAL_FOLLOWERS = ((45559, 1), (83154, 0))


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def center_prefix(stop: int) -> bytearray:
    row = 1
    answer = bytearray(stop)
    for time in range(stop):
        answer[time] = (row >> time) & 1
        row ^= (row << 1) | (row << 2)
    return answer


def history_profile(word: bytes | bytearray, history: int):
    length = len(word) - history
    require(length > 0, ("history exceeds word", history, len(word)))
    counts: dict[int, list[int]] = {}
    if history == 0:
        ones = sum(word)
        counts[0] = [len(word) - ones, ones]
    else:
        mask = (1 << history) - 1
        context = 0
        for value in word[:history]:
            context = ((context << 1) | value) & mask
        for index in range(history, len(word)):
            pair = counts.get(context)
            if pair is None:
                pair = [0, 0]
                counts[context] = pair
            pair[word[index]] += 1
            context = ((context << 1) | word[index]) & mask
    error = sum(min(pair) for pair in counts.values())
    mismatches = sum(pair[0] > 0 and pair[1] > 0 for pair in counts.values())
    variation = sum(abs(pair[1] - pair[0]) for pair in counts.values())
    maximum_bias = max(abs(pair[1] - pair[0]) for pair in counts.values())
    tail_discrepancy = sum(1 if value else -1 for value in word[history:])
    require(variation == length - 2 * error, ("mass identity", history))
    require(abs(tail_discrepancy) <= variation, ("discrepancy bound", history))
    return (
        history,
        length,
        error,
        mismatches,
        len(counts),
        variation,
        maximum_bias,
        tail_discrepancy,
    )


def right_special_witness(word: bytes | bytearray, history: int):
    mask = (1 << history) - 1
    context = 0
    for value in word[:history]:
        context = ((context << 1) | value) & mask
    first: dict[int, tuple[int, int]] = {}
    for index in range(history, len(word)):
        value = word[index]
        old = first.get(context)
        if old is not None and old[0] != value:
            return format(context, f"0{history}b"), (old[1], old[0]), (index, value)
        first.setdefault(context, (value, index))
        context = ((context << 1) | value) & mask
    return None


def shell_row(center: bytearray, scale: int):
    start = 1 << scale
    stop = 1 << (scale + 1)
    running = 0
    minimum = maximum = 0
    minimum_at = maximum_at = 0
    for offset, value in enumerate(center[start:stop], 1):
        running += 1 if value else -1
        if running < minimum:
            minimum, minimum_at = running, offset
        if running > maximum:
            maximum, maximum_at = running, offset
    amplitude = max(-minimum, maximum)
    return (
        scale,
        stop - start,
        running,
        minimum,
        minimum_at,
        maximum,
        maximum_at,
        amplitude,
    )


def sparse_flip_control(power: int):
    """Balanced/easy word with a linearly long right-special context."""

    stop = 2 * power - 1
    word = bytearray(
        (index & 1) ^ int(index > 0 and (index & (index - 1)) == 0)
        for index in range(stop)
    )
    history = (power - 2) // 2
    left_follower = power
    right_follower = 2 * power - 2
    require(
        word[left_follower - history : left_follower]
        == word[right_follower - history : right_follower],
        ("sparse control contexts", power),
    )
    require(
        word[left_follower] != word[right_follower],
        ("sparse control followers", power),
    )
    discrepancy = sum(1 if value else -1 for value in word)
    return power, stop, discrepancy, history, left_follower, right_follower


def main() -> None:
    center = center_prefix(HORIZON)
    center_hash = hashlib.sha256(center).hexdigest()
    require(center_hash == EXPECTED_CENTER_SHA256, "frozen center hash")

    discrepancy = 2 * sum(center) - len(center)
    shell_rows = tuple(shell_row(center, scale) for scale in range(18))
    profiles = tuple(history_profile(center, history) for history in HISTORY_DEPTHS)
    witness = right_special_witness(center, 33)
    require(witness is not None, "length-33 right-special witness")
    require(witness[0] == EXPECTED_RIGHT_SPECIAL_CONTEXT, "right-special context")
    require(
        (witness[1], witness[2]) == EXPECTED_RIGHT_SPECIAL_FOLLOWERS,
        "right-special followers",
    )
    require(history_profile(center, 34)[2] == 0, "depth-34 exact finite lookup")
    require(
        all(
            row[3] == (1 << row[0]) and row[4] == (1 << row[0])
            for row in profiles if row[0] <= 12
        ),
        "every short context is observed and ambiguous",
    )

    alternating = bytearray(index & 1 for index in range(1 << 16))
    thue_morse = bytearray(index.bit_count() & 1 for index in range(1 << 16))
    controls = (
        ("alternating", 2 * sum(alternating) - len(alternating),
         tuple(history_profile(alternating, h)[2:5] for h in (0, 1, 8, 32))),
        ("thue_morse", 2 * sum(thue_morse) - len(thue_morse),
         tuple(history_profile(thue_morse, h)[2:5] for h in (0, 1, 8, 32))),
    )
    sparse = tuple(sparse_flip_control(1 << exponent) for exponent in range(5, 13))

    print("RULE30_PRIZE23_SHELL_HISTORY_SCRATCH")
    print(f"universe=center_t=0..{HORIZON-1}")
    print(f"center_sha256={center_hash}")
    print(f"center_ones={sum(center)};prefix_discrepancy={discrepancy}")
    print(f"shell_rows={shell_rows}")
    print(
        "history_columns=(h,length,optimal_error,mismatch_fibres,observed_contexts,"
        "conditional_L1,max_fibre_bias,tail_discrepancy)"
    )
    print(f"history_profiles={profiles}")
    print(f"right_special_h33={witness}")
    print(f"controls={controls}")
    print(f"sparse_flip_controls={sparse}")
    print(
        "scope=FINITE-EXACT through 2^18; universal shell/predictor identities are "
        "proved separately; no balance limit or complexity lower bound"
    )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
