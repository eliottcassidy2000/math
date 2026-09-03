#!/usr/bin/env python3
"""Clean-room exact replay of the finite h=420 attachment-address hostile."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from math import ceil, floor
import sys


H = 420
COMMON = tuple(11 + 1680 * k for k in range(7)) + (525, 945, 1365, 1575)


def interval(w: int, n: int) -> tuple[Fraction, Fraction]:
    return Fraction(14 * n - 1, 14 * w), Fraction(14 * n + 1, 14 * w)


def anchor(k: int) -> tuple[Fraction, Fraction]:
    return Fraction(14 * k + 1, 28 * H), Fraction(14 * k + 13, 28 * H)


def trace(k: int, speeds: tuple[int, ...]):
    left, right = anchor(k)
    teeth = []
    for w in speeds:
        lo = floor(w * left - Fraction(1, 14)) - 1
        hi = ceil(w * right + Fraction(1, 14)) + 1
        for n in range(lo, hi + 1):
            a, b = interval(w, n)
            if a < right and left < b:
                teeth.append((w, n, a, b))
    cursor = left
    chain = []
    while True:
        live = [t for t in teeth if t[2] < cursor < t[3]]
        if not live:
            return "missing", tuple(chain), cursor
        winner = max(live, key=lambda t: (t[3], -t[2], -t[0], -t[1]))
        chain.append(winner)
        cursor = winner[3]
        if cursor > right:
            return ("span" if len(chain) == 1 else "renew"), tuple(chain), cursor


def labels(extra: int) -> dict[int, str]:
    ans = dict(
        zip(
            COMMON,
            [f"D{i}" for i in range(7)] + [f"C{i}" for i in range(4)],
            strict=True,
        )
    )
    ans[extra] = "P"
    return ans


def row(extra: int):
    speeds = COMMON + (extra,)
    names = labels(extra)
    status_map = []
    traces = []
    coarse = Counter()
    for k in range(2 * H):
        status, chain, cursor = trace(k, speeds)
        epsilon = int(k >= H)
        word = tuple(names[t[0]] for t in chain)
        coarse[(epsilon, status, word if status != "missing" else ())] += 1
        status_map.append(status)
        traces.append((status, tuple((names[t[0]], t[1]) for t in chain), cursor))
    return tuple(status_map), tuple(traces), tuple(sorted(coarse.items(), key=repr))


def dist(x: Fraction) -> Fraction:
    r = x % 1
    return min(r, 1 - r)


def clipped(w: int, shift: Fraction):
    lo = floor(w * shift - Fraction(1, 14)) - 1
    hi = ceil(Fraction(w, 2) + w * shift + Fraction(1, 14)) + 1
    for n in range(lo, hi + 1):
        a = Fraction(n, w) - shift - Fraction(1, 14 * w)
        b = Fraction(n, w) - shift + Fraction(1, 14 * w)
        a, b = max(Fraction(), a), min(Fraction(1, 2), b)
        if a < b:
            yield a, b


def core_energy(speeds: tuple[int, ...]) -> Fraction:
    events = defaultdict(lambda: [0, 0, 0])
    for w in speeds:
        for a, b in clipped(w, Fraction()):
            events[a][0] += 1
            events[b][0] -= 1
        for a, b in clipped(w, Fraction(1, 2)):
            events[a][1] += 1
            events[b][1] -= 1
    for a, b in clipped(2 * H, Fraction()):
        events[a][2] += 1
        events[b][2] -= 1
    events[Fraction()]
    events[Fraction(1, 2)]
    lower = upper = anchor_depth = 0
    energy = mass = Fraction()
    walls = sorted(events)
    for x, y in zip(walls, walls[1:]):
        dl, du, da = events[x]
        lower += dl
        upper += du
        anchor_depth += da
        if anchor_depth == 0:
            width = y - x
            mass += width
            energy += width * (lower - upper) ** 2
    assert mass == Fraction(3, 7)
    return energy


def signed_current(t: Fraction, speeds: tuple[int, ...]):
    a = sum(dist(w * t) < Fraction(1, 14) for w in speeds)
    b = sum(dist(w * (t + Fraction(1, 2))) < Fraction(1, 14) for w in speeds)
    return a, b, a - b


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    sa, ta, qa = row(1287)
    sb, tb, qb = row(9009)
    assert sa == sb
    assert qa == qb
    changed = tuple(k for k in range(2 * H) if ta[k] != tb[k])
    role_changed = tuple(
        k
        for k in changed
        if tuple(x[0] for x in ta[k][1]) != tuple(x[0] for x in tb[k][1])
    )
    assert len(changed) == 209 and len(role_changed) == 208
    assert all(ta[k][0] == "missing" for k in changed)
    assert ta[9][:2] == (
        "missing",
        (("C3", 17), ("P", 14), ("D5", 92), ("C2", 15)),
    )
    assert tb[9][:2] == ("missing", (("C3", 17), ("D4", 73)))
    assert ta[9][2] == Fraction(211, 19110)
    assert tb[9][2] == Fraction(1023, 94234)

    speeds_a = COMMON + (1287,)
    speeds_b = COMMON + (9009,)
    x = Fraction(239, 22050)
    left_test = (Fraction(1021, 94234) + x) / 2
    right_test = (x + Fraction(766, 70637)) / 2
    assert signed_current(left_test, speeds_a) == (3, 1, 2)
    assert signed_current(left_test, speeds_b) == (2, 1, 1)
    assert signed_current(right_test, speeds_a) == (2, 1, 1)
    assert signed_current(right_test, speeds_b) == (1, 1, 0)

    ea = core_energy(speeds_a)
    eb = core_energy(speeds_b)
    assert ea == Fraction(103565349065690759276041319, 79794403580513459456309691)
    assert eb == Fraction(58621224727881646861397, 45136410657303198493260)
    assert ea != eb

    # Independent exact check of the live h=420,u=3 resonant side split.
    h, u, L = 420, 3, 20
    wall = Fraction(1, 14 * u)
    delta = Fraction(1, 28 * h)
    inward, outward = wall - delta, wall + delta
    assert dist(2 * h * inward) == dist(2 * h * outward) == Fraction(1, 14)
    assert dist(u * inward) == Fraction(1, 14) - Fraction(1, 196 * L)
    assert dist(u * outward) == Fraction(1, 14) + Fraction(1, 196 * L)
    assert dist(13 * u * inward) == Fraction(1, 14) + Fraction(13, 196 * L)
    assert dist(13 * u * outward) == Fraction(1, 14) - Fraction(13, 196 * L)
    component_width = Fraction(3, 49 * L * u)
    inward_slack = Fraction(1, 7 * u) - delta - component_width
    outward_slack = Fraction(1, 91 * u) - delta - component_width
    assert inward_slack == Fraction(28 * L - 13, 196 * L * u) > 0
    assert outward_slack == Fraction(28 * L - 169, 2548 * L * u) > 0

    print("LRC14_PARTIAL_ATTACHMENT_ADDRESS_NONFACTORIZATION_INDEPENDENT_REFEREE=PASS")
    print(f"coarse_sha256={sha256(repr(qa).encode()).hexdigest()}")
    print(
        f"status_map_equal={sa == sb};changed={len(changed)};"
        f"role_changed={len(role_changed)};changed_nonmissing=0"
    )
    print(f"k9_A={ta[9]};k9_B={tb[9]}")
    print(f"core_energy_A={ea};core_energy_B={eb}")
    print(
        "resonant_h420_u3="
        f"wall:{wall};in:{inward};out:{outward};"
        f"in_slack:{inward_slack};out_slack:{outward_slack}"
    )
    print("CHECKS=1742_PLUS_EXACT_SWEEPS")


if __name__ == "__main__":
    main()
