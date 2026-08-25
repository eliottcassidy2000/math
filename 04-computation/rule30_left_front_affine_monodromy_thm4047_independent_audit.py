#!/usr/bin/env python3
"""No-import packed-orbit audit for THM-4047's finite physical bank."""

from __future__ import annotations

from collections import Counter
import hashlib


MAX_R = 100_000
CHECK_TIME = 2_000_000
GLOBAL_PERIOD = 32
EXPECTED_ZEROS = (2, 7, 28, 399, 53207, 58286, 87866)
EXPECTED_DOUBLINGS = (3, 8, 29, 400, 87867)
EXPECTED_IDENTITIES = (53208, 58287)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def packed_step(state: int, mask: int) -> int:
    shifted = state << 1
    return ((shifted << 1) ^ shifted ^ state ^ (shifted & state)) & mask


def least_period(values: tuple[int, ...]) -> int:
    for period in (1, 2, 4, 8, 16, 32):
        if all(values[i] == values[i % period] for i in range(len(values))):
            return period
    raise RuntimeError("period exceeds audit window")


def main() -> None:
    # Slow scalar positive control for the packed orientation and boundary.
    control_width = 96
    scalar = [1] + [0] * (control_width - 1)
    packed = 1
    control_mask = (1 << control_width) - 1
    for time in range(256):
        require(
            packed == sum(bit << r for r, bit in enumerate(scalar)),
            ("scalar/packed control", time),
        )
        scalar = [
            (scalar[r - 2] if r >= 2 else 0)
            ^ (scalar[r - 1] if r >= 1 else 0)
            ^ scalar[r]
            ^ ((scalar[r - 1] if r >= 1 else 0) & scalar[r])
            for r in range(control_width)
        ]
        packed = packed_step(packed, control_mask)

    # Completely direct physical evolution; no theorem script is imported and
    # no fitted per-column start or phase is supplied.
    mask = (1 << (MAX_R + 1)) - 1
    state = 1
    for _ in range(CHECK_TIME):
        state = packed_step(state, mask)
    states = [state]
    for _ in range(GLOBAL_PERIOD):
        state = packed_step(state, mask)
        states.append(state)
    require(states[0] == states[GLOBAL_PERIOD], "closed strip does not repeat")
    require(states[0] != states[GLOBAL_PERIOD // 2], "claimed global period is not least")

    words: list[tuple[int, ...]] = []
    periods: list[int] = []
    for r in range(MAX_R + 1):
        values = tuple((states[offset] >> r) & 1 for offset in range(GLOBAL_PERIOD))
        period = least_period(values)
        words.append(values[:period])
        periods.append(period)

    zeros = tuple(r for r, word in enumerate(words) if not any(word))
    doublings: list[int] = []
    identities: list[int] = []
    for r in range(2, MAX_R + 1):
        if not any(words[r - 1]):
            translation = sum(words[r - 2]) & 1
            (doublings if translation else identities).append(r)

    require(zeros == EXPECTED_ZEROS, ("zero census", zeros))
    require(tuple(doublings) == EXPECTED_DOUBLINGS, ("doubling census", doublings))
    require(tuple(identities) == EXPECTED_IDENTITIES, ("identity census", identities))

    byte_width = (MAX_R + 8) // 8
    payload = b"".join(item.to_bytes(byte_width, "little") for item in states[:-1])
    period_histogram = tuple(sorted(Counter(periods).items()))
    print("RULE30_LEFT_FRONT_AFFINE_MONODROMY_THM4047_INDEPENDENT_AUDIT")
    print(f"universe=direct_single_seed_left_front_strip_r=0..{MAX_R}")
    print(f"physical_repeat=(time={CHECK_TIME},period={GLOBAL_PERIOD},least_global=True)")
    print(f"zero_columns={zeros}")
    print(f"period_doubling_columns={tuple(doublings)}")
    print(f"identity_monodromy_columns={tuple(identities)}")
    print(f"period_histogram={period_histogram}")
    print(f"cycle_sha256={hashlib.sha256(payload).hexdigest()}")
    print("audit=PASS;direct packed orbit;scalar orientation control;no primary import")


if __name__ == "__main__":
    main()
