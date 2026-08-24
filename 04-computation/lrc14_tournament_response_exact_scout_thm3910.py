#!/usr/bin/env python3
"""Exact AP11 mixed-incidence scout for the THM-3878 scale-one survivors.

This is the exact THM-3910 hostile/positive-control companion.  It independently builds open danger combs with Fraction endpoints
and compares the one/two-way incidence quotient with the load-bearing
three-way response

    mu(G_AP11 \\ (D_(tp) union D_(tq))).

The scan deliberately excludes t=11 when p=1, because that duplicates the
largest AP11 speed.  Positive measure is only an AP11-body control and makes
no claim for an arbitrary eleven-speed component.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
import sys


sys.stdout.reconfigure(newline="\n")

Interval = tuple[Fraction, Fraction]

SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)


def danger(frequency: int) -> list[Interval]:
    """Merged open 1/14-danger comb clipped to [0,1]."""
    radius = Fraction(1, 14 * frequency)
    pieces: list[Interval] = []
    for k in range(frequency + 1):
        centre = Fraction(k, frequency)
        left = max(Fraction(0), centre - radius)
        right = min(Fraction(1), centre + radius)
        if left < right:
            pieces.append((left, right))
    return merge(pieces)


def merge(intervals: list[Interval]) -> list[Interval]:
    merged: list[list[Fraction]] = []
    for left, right in sorted(intervals):
        # Touching open intervals remain distinct topologically, but merging
        # them is measure-exact, which is all this scout uses.
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return [(left, right) for left, right in merged]


def complement(intervals: list[Interval]) -> list[Interval]:
    answer: list[Interval] = []
    cursor = Fraction(0)
    for left, right in intervals:
        if cursor < left:
            answer.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        answer.append((cursor, Fraction(1)))
    return answer


def intersect(lefts: list[Interval], rights: list[Interval]) -> list[Interval]:
    answer: list[Interval] = []
    i = j = 0
    while i < len(lefts) and j < len(rights):
        left = max(lefts[i][0], rights[j][0])
        right = min(lefts[i][1], rights[j][1])
        if left < right:
            answer.append((left, right))
        if lefts[i][1] < rights[j][1]:
            i += 1
        elif rights[j][1] < lefts[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return answer


def measure(intervals: list[Interval]) -> Fraction:
    return sum((right - left for left, right in intervals), Fraction(0))


def fmt(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def parity_hostile() -> dict[str, object]:
    """Same one/two moments and cut metric, different LRC target atom 100."""
    even = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
    odd = ((1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1))

    def moments(atoms: tuple[tuple[int, int, int], ...]):
        one = tuple(sum(atom[i] for atom in atoms) for i in range(3))
        two = tuple(
            sum(atom[i] * atom[j] for atom in atoms)
            for i, j in ((0, 1), (0, 2), (1, 2))
        )
        cut = tuple(one[i] + one[j] - 2 * pair for (i, j), pair in zip(
            ((0, 1), (0, 2), (1, 2)), two
        ))
        target = sum(atom == (1, 0, 0) for atom in atoms)
        triple = sum(atom == (1, 1, 1) for atom in atoms)
        return one, two, cut, target, triple

    me = moments(even)
    mo = moments(odd)
    if me[:3] != mo[:3] or (me[3], mo[3]) != (0, 1):
        raise RuntimeError("parity hostile failed")
    return {
        "one": me[0],
        "two": me[1],
        "cut": me[2],
        "target_even_quarters": me[3],
        "target_odd_quarters": mo[3],
        "triple_even_quarters": me[4],
        "triple_odd_quarters": mo[4],
    }


def main() -> None:
    if len(SCALE1) != 57 or len(set(SCALE1)) != 57:
        raise RuntimeError("survivor universe mismatch")

    core_danger = merge(sum((danger(w) for w in range(1, 12)), []))
    core_safe = complement(core_danger)
    m0 = measure(core_safe)
    if m0 != Fraction(10931, 194040):
        raise RuntimeError(f"AP11 safe mass mismatch: {m0}")

    # Cache all combs used in the bounded t scan.
    t_min, t_max = 12, 42
    frequencies = {
        t * x
        for t in range(t_min, t_max + 1)
        for p, q in SCALE1
        for x in (p, q)
    }
    comb = {frequency: danger(frequency) for frequency in frequencies}

    pairwise_fibres: dict[tuple[Fraction, Fraction, Fraction], list[dict[str, object]]] = defaultdict(list)
    records: list[dict[str, object]] = []
    for t in range(t_min, t_max + 1):
        for p, q in SCALE1:
            dp = comb[t * p]
            dq = comb[t * q]
            gp = intersect(core_safe, dp)
            gq = intersect(core_safe, dq)
            pq = intersect(dp, dq)
            gpq = intersect(gp, dq)
            mp = measure(gp)
            mq = measure(gq)
            cpq = measure(pq)
            mpq = measure(gpq)
            residual = m0 - mp - mq + mpq
            if residual < 0:
                raise RuntimeError("negative inclusion-exclusion residual")
            record = {
                "t": t,
                "p": p,
                "q": q,
                "mp": mp,
                "mq": mq,
                "cpq": cpq,
                "mpq": mpq,
                "residual": residual,
            }
            records.append(record)
            # This is the complete one/two-way quotient for indicators
            # (G,D_tp,D_tq), since all three one-masses are fixed except the
            # displayed G-p/G-q intersections and Cpq.
            pairwise_fibres[(mp, mq, cpq)].append(record)

    collisions = []
    for key, fibre in pairwise_fibres.items():
        triple_values = {record["mpq"] for record in fibre}
        if len(triple_values) > 1:
            ordered = sorted(fibre, key=lambda r: (r["t"], r["p"], r["q"]))
            first = ordered[0]
            second = next(r for r in ordered[1:] if r["mpq"] != first["mpq"])
            collisions.append((key, first, second, len(fibre), len(triple_values)))

    minimum = min(records, key=lambda r: (r["residual"], r["t"], r["p"], r["q"]))
    zero_count = sum(record["residual"] == 0 for record in records)
    positive_count = len(records) - zero_count

    # The 58th survivor is the scale-two pair (1,9).  Directly test the AP11
    # body at every legal odd t in the same bounded window.  Here the core is
    # 2*(1,...,11), not the unscaled AP11 core above.
    core2_danger = merge(sum((danger(2 * w) for w in range(1, 12)), []))
    core2_safe = complement(core2_danger)
    scale2_records: list[dict[str, object]] = []
    for t in range(11, t_max + 1, 2):
        dt = danger(t)
        d9t = danger(9 * t)
        gt = intersect(core2_safe, dt)
        g9t = intersect(core2_safe, d9t)
        gt9t = intersect(gt, d9t)
        residual = measure(core2_safe) - measure(gt) - measure(g9t) + measure(gt9t)
        scale2_records.append({"t": t, "residual": residual})
    scale2_minimum = min(scale2_records, key=lambda r: (r["residual"], r["t"]))
    scale2_zero_count = sum(record["residual"] == 0 for record in scale2_records)

    def serial(record: dict[str, object]) -> dict[str, object]:
        return {
            key: fmt(value) if isinstance(value, Fraction) else value
            for key, value in record.items()
        }

    witness = None
    if collisions:
        key, first, second, fibre_size, triple_count = collisions[0]
        witness = {
            "pairwise_key": tuple(fmt(x) for x in key),
            "first": serial(first),
            "second": serial(second),
            "fibre_size": fibre_size,
            "triple_values": triple_count,
        }

    parity = parity_hostile()
    semantic = {
        "scope": "AP11 body, 57 THM-3878 scale-one certificate survivors, 12<=t<=42",
        "rows": len(records),
        "positive": positive_count,
        "zero": zero_count,
        "collision_fibres": len(collisions),
        "minimum": serial(minimum),
        "scale2_rows": len(scale2_records),
        "scale2_zero": scale2_zero_count,
        "scale2_minimum": serial(scale2_minimum),
        "first_collision": witness,
        "parity_hostile": parity,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("LRC14_TOURNAMENT_RESPONSE_EXACT_SCOUT_20260823")
    print(f"scope=AP11;scale1_survivors=57;t={t_min}..{t_max};rows={len(records)}")
    print(f"ap11_safe_mass={fmt(m0)};positive_responses={positive_count};zero_responses={zero_count}")
    print("minimum_response=" + json.dumps(serial(minimum), sort_keys=True, separators=(",", ":")))
    print(f"scale2_AP11_odd_t_rows={len(scale2_records)};zero_responses={scale2_zero_count};minimum=" + json.dumps(serial(scale2_minimum), sort_keys=True, separators=(",", ":")))
    print(f"pairwise_fibres={len(pairwise_fibres)};triple_split_fibres={len(collisions)}")
    print("first_native_collision=" + json.dumps(witness, sort_keys=True, separators=(",", ":")))
    print("parity_hostile=" + json.dumps(parity, sort_keys=True, separators=(",", ":")))
    print(f"semantic_sha256={digest}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
