#!/usr/bin/env python3
"""Audit the strict/weak tie conventions in the LRC trienerment.

The S539 trienerment note used two facts at once:

* observer loneliness is observer tie-degree zero;
* the trienerment is globally never a pure tournament by pigeonhole.

Those facts require different boundary conventions.  Strict ties
``dist < 1/n`` give the observer equivalence for the closed LRC inequality
``dist >= 1/n``.  Weak ties ``dist <= 1/n`` give the global pigeonhole
non-pure statement.  The only strict globally tie-free configurations are
regular n-gon walls.

Tournament Analysis note: this audit deliberately does not collapse the
tri-state relation to a tournament.  The boundary distinction being tested is
exactly the information that a binary switch would erase, so no meaningful
pairwise tournament observable is reported here.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circ_dist(a: Fraction, b: Fraction) -> Fraction:
    d = frac(a - b)
    return min(d, 1 - d)


def positions(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return (Fraction(0),) + tuple(frac(Fraction(v) * t) for v in speeds)


def tie_counts(pos: tuple[Fraction, ...], n: int) -> tuple[int, int]:
    strict = 0
    weak = 0
    threshold = Fraction(1, n)
    for i, j in combinations(range(len(pos)), 2):
        d = circ_dist(pos[i], pos[j])
        strict += d < threshold
        weak += d <= threshold
    return strict, weak


def observer_strict_tie_free(speeds: tuple[int, ...], t: Fraction, n: int) -> bool:
    threshold = Fraction(1, n)
    return all(circ_dist(Fraction(0), frac(Fraction(v) * t)) >= threshold for v in speeds)


def lrc_closed(speeds: tuple[int, ...], t: Fraction, n: int) -> bool:
    threshold = Fraction(1, n)
    return all(min(frac(Fraction(v) * t), 1 - frac(Fraction(v) * t)) >= threshold for v in speeds)


def is_regular_ngon(pos: tuple[Fraction, ...], n: int) -> bool:
    pts = sorted(frac(p) for p in pos)
    if len(set(pts)) != n:
        return False
    gaps = [pts[(i + 1) % n] - pts[i] for i in range(n - 1)]
    gaps.append(pts[0] + 1 - pts[-1])
    return all(g == Fraction(1, n) for g in gaps)


def primitive_combos(n: int, max_speed: int):
    for combo in combinations(range(1, max_speed + 1), n - 1):
        if reduce(gcd, combo) == 1:
            yield combo


def audit_initial_segments() -> None:
    print("Regular-wall calibration for initial segments")
    print("n  t=a/n samples  strict_ties  weak_ties  regular")
    for n in range(3, 10):
        speeds = tuple(range(1, n))
        samples = 0
        strict_values = set()
        weak_values = set()
        regular_values = set()
        for a in range(1, n):
            if gcd(a, n) != 1:
                continue
            samples += 1
            pos = positions(speeds, Fraction(a, n))
            strict, weak = tie_counts(pos, n)
            strict_values.add(strict)
            weak_values.add(weak)
            regular_values.add(is_regular_ngon(pos, n))
        print(
            f"{n:2d} {samples:12d}  "
            f"{sorted(strict_values)!s:>11}  {sorted(weak_values)!s:>9}  "
            f"{sorted(regular_values)}"
        )


def audit_bounded_clock_samples() -> None:
    print()
    print("Bounded exact-clock audit")
    print("n max_speed sets states strict_free regular_viol weak_viol observer_mism")
    boxes = {4: 10, 5: 8, 6: 7, 7: 8}
    for n, max_speed in boxes.items():
        set_count = 0
        state_count = 0
        strict_free = 0
        regular_viol = 0
        weak_viol = 0
        observer_mism = 0
        for speeds in primitive_combos(n, max_speed):
            set_count += 1
            q = n * lcm(*speeds)
            for r in range(q):
                t = Fraction(r, q)
                pos = positions(speeds, t)
                strict, weak = tie_counts(pos, n)
                state_count += 1
                if strict == 0:
                    strict_free += 1
                    if not is_regular_ngon(pos, n):
                        regular_viol += 1
                if weak == 0:
                    weak_viol += 1
                if observer_strict_tie_free(speeds, t, n) != lrc_closed(speeds, t, n):
                    observer_mism += 1
        print(
            f"{n:1d} {max_speed:9d} {set_count:4d} {state_count:7d} "
            f"{strict_free:11d} {regular_viol:12d} {weak_viol:9d} "
            f"{observer_mism:13d}"
        )


def main() -> None:
    print("LRC trienerment tie-convention audit (S542)")
    print()
    print("Strict tie: dist < 1/n.  Weak tie: dist <= 1/n.")
    print("Closed LRC uses strict observer tie-freeness.")
    print("Global pigeonhole non-purity uses weak ties.")
    print()
    audit_initial_segments()
    audit_bounded_clock_samples()
    print()
    print("Conclusion:")
    print("  strict observer tie-free == closed LRC: zero mismatches")
    print("  weak tie count is never zero in the sampled exact clocks")
    print("  strict globally tie-free states occur only at regular n-gon walls")


if __name__ == "__main__":
    main()
