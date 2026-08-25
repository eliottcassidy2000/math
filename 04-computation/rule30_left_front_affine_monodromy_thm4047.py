#!/usr/bin/env python3
"""Exact physical left-front certificates for THM-4047.

For the single-seed Rule 30 orbit let ell_r(t)=a_t(-t+r).  The script uses
the affine recurrence

    ell_r(t+1) = ell_(r-2)(t) + ell_(r-1)(t)
                   + (1+ell_(r-1)(t))*ell_r(t)

over F_2.  Ultimately periodic lower columns give an exact one-block affine
monodromy for the next column.  Reset blocks forget the transient; the rare
multiplier-one blocks obtain their missing phase from the packed physical
orbit.  A final packed-state repeat verifies the entire finite strip at once.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib


MAX_R = 100_000


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


@dataclass(frozen=True)
class Steady:
    start: int
    word: tuple[int, ...]  # value at absolute time t is word[t mod p]
    monodromy: tuple[int, int]
    queried_phase: bool

    def value(self, time: int) -> int:
        require(time >= self.start, ("pre-steady query", time, self.start))
        return self.word[time % len(self.word)]


def one_step(a: int, b: int, x: int) -> int:
    return a ^ b ^ x ^ (b & x)


def absolute_word(start: int, values: tuple[int, ...]) -> tuple[int, ...]:
    """Return the least-period word indexed by absolute time residue."""
    raw_period = len(values)
    period = raw_period
    for candidate in range(1, raw_period + 1):
        if raw_period % candidate == 0 and all(
            values[i] == values[i % candidate] for i in range(raw_period)
        ):
            period = candidate
            break
    word: list[int | None] = [None] * period
    for offset, value in enumerate(values):
        residue = (start + offset) % period
        previous = word[residue]
        require(previous is None or previous == value, ("phase conflict", start, period))
        word[residue] = value
    require(all(value is not None for value in word), ("unfilled residue", start, period))
    return tuple(int(value) for value in word)


def main() -> None:
    # Base columns are exact: ell_0 is constantly one, while ell_1 is one
    # from time one onward.  Their displayed monodromies describe their own
    # one-step startup maps and are not used as lower-driver certificates.
    steady = [
        Steady(0, (1,), (1, 0), False),
        Steady(1, (1,), (0, 1), False),
    ]

    # Bit r of ``front`` is the physical ell_r(t).  The truncated vector is
    # closed because every update uses only r, r-1, and r-2.
    mask = (1 << (MAX_R + 1)) - 1
    front = 1
    simulation_time = 0

    def advance_to(target_time: int) -> None:
        nonlocal front, simulation_time
        require(target_time >= simulation_time, ("nonmonotone query", target_time, simulation_time))
        while simulation_time < target_time:
            shifted = front << 1
            front = ((shifted << 1) ^ shifted ^ front ^ (shifted & front)) & mask
            simulation_time += 1

    phase_queries: list[tuple[int, int, int, tuple[int, int]]] = []
    for r in range(2, MAX_R + 1):
        left = steady[r - 2]
        right = steady[r - 1]
        common_start = max(left.start, right.start)
        common_period = max(len(left.word), len(right.word))
        require(
            common_period % len(left.word) == 0
            and common_period % len(right.word) == 0,
            ("non-dyadic lower periods", r),
        )

        image_zero = 0
        image_one = 1
        for time in range(common_start, common_start + common_period):
            a = left.value(time)
            b = right.value(time)
            image_zero = one_step(a, b, image_zero)
            image_one = one_step(a, b, image_one)
        multiplier = image_zero ^ image_one
        constant = image_zero

        if multiplier == 0:
            # One driver block resets both possible incoming states to C.
            start = common_start + common_period
            length = common_period
            initial = constant
            queried = False
        else:
            # The lower periodic tails retain one physical phase bit.
            start = common_start
            length = common_period if constant == 0 else 2 * common_period
            advance_to(start)
            initial = (front >> r) & 1
            queried = True
            phase_queries.append((r, start, initial, (multiplier, constant)))

        values: list[int] = []
        x = initial
        for time in range(start, start + length):
            values.append(x)
            x = one_step(left.value(time), right.value(time), x)
        require(x == initial, ("period wrap", r, start, length))
        word = absolute_word(start, tuple(values))
        item = Steady(start, word, (multiplier, constant), queried)

        # The wrapped recurrence is an all-future certificate, not a tail fit.
        check_length = 2 * max(common_period, len(word))
        for time in range(start, start + check_length):
            require(
                item.value(time + 1)
                == one_step(left.value(time), right.value(time), item.value(time)),
                ("wrapped recurrence mismatch", r, time),
            )
        steady.append(item)

    zeros = tuple(r for r, item in enumerate(steady) if not any(item.word))
    doublings = tuple(r for r, item in enumerate(steady) if item.monodromy == (1, 1))
    identities = tuple(
        r for r, item in enumerate(steady) if r >= 2 and item.monodromy == (1, 0)
    )
    resets = tuple(r for r, item in enumerate(steady) if r >= 2 and item.monodromy[0] == 0)

    # Audit the exact clock law: multiplier one occurs exactly behind a zero
    # column, and its translation is the previous column's weight parity.
    for r in range(2, MAX_R + 1):
        previous_is_zero = not any(steady[r - 1].word)
        require(
            (steady[r].monodromy[0] == 1) == previous_is_zero,
            ("zero/multiplier clock", r),
        )
        if previous_is_zero:
            expected_constant = sum(steady[r - 2].word) & 1
            require(
                steady[r].monodromy[1] == expected_constant,
                ("weight parity clock", r, expected_constant),
            )

    # A single exact repeat of the closed physical strip proves its future
    # periodicity by determinism and independently checks every stored phase.
    global_period = max(len(item.word) for item in steady)
    global_start = max(item.start for item in steady)
    advance_to(global_start)
    state_at_global_start = front
    for r, item in enumerate(steady):
        require(
            ((front >> r) & 1) == item.value(global_start),
            ("physical/certificate mismatch", r, global_start),
        )
    advance_to(global_start + global_period)
    require(front == state_at_global_start, ("packed global repeat", global_start, global_period))

    rows = tuple(
        (r, item.start, item.word, item.monodromy, item.queried_phase)
        for r, item in enumerate(steady)
    )
    period_histogram = tuple(sorted(Counter(len(item.word) for item in steady).items()))
    diagonal_after_stored_certificate_start = tuple(
        r for r, item in enumerate(steady) if r >= item.start
    )
    print("RULE30_LEFT_FRONT_AFFINE_MONODROMY_THM4047")
    print(f"universe=single_seed_left_front_columns_r=0..{MAX_R}")
    print("certificate=all_times_after_each_start;wrapped_recurrence_plus_physical_phase")
    print(f"zero_columns={zeros}")
    print(f"zero_column_natural_ranks={tuple(enumerate(zeros, start=1))}")
    print(f"period_doubling_columns={doublings}")
    print(f"identity_monodromy_columns={identities}")
    print(f"reset_column_count={len(resets)}")
    print(f"physical_phase_queries={tuple(phase_queries)}")
    print(f"period_histogram={period_histogram}")
    print(
        "center_diagonal_at_or_after_stored_certificate_start_columns="
        f"{diagonal_after_stored_certificate_start}"
    )
    print(f"global_repeat=(start={global_start},period={global_period})")
    print(f"packed_simulation_steps={simulation_time}")
    print(f"certificate_sha256={hashlib.sha256(repr(rows).encode('ascii')).hexdigest()}")
    print("scope=PROVED universal affine criterion;FINITE-EXACT physical bank;Rule30 prizes OPEN")


if __name__ == "__main__":
    main()
