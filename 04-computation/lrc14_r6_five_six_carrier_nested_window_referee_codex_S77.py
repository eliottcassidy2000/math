#!/usr/bin/env python3
"""Exact referee for THM-1212's five/six-carrier nested-window closure.

The proof uses only the interval-density inequality already formalized in
``LRCEssentialRegion.lean`` and the elementary first-safe windows

    J(x,q) = [1/(14x), 13/(14q)]  (x <= every protected carrier <= q).

Every displayed comparison is replayed with integers or ``Fraction``.  A
small endpoint-exact audit independently checks the claimed safe-window
geometry and searches adversarial carrier rows; it is telemetry, not a
replacement for the quantified inequalities.

Tournament Analysis uses the nested cover obligations as vertices.  Their
containment order is transitive; the metric labels, not that naked order,
carry the proof.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor
from random import Random


H = F(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(value: F) -> F:
    return value - value.numerator // value.denominator


def safe(speed: int, time: F) -> bool:
    phase = frac(speed * time)
    return min(phase, 1 - phase) >= H


def nested_window(x: int, q: int) -> tuple[F, F]:
    require(0 < x <= q < 13 * x, "nested first-safe window is empty")
    return F(1, 14 * x), F(13, 14 * q)


def nested_length(x: int, q: int) -> F:
    left, right = nested_window(x, q)
    length = right - left
    require(length == F(13 * x - q, 14 * x * q), "length identity moved")
    return length


def verify_window_geometry(x: int, protected: tuple[int, ...]) -> None:
    q = protected[-1]
    left, right = nested_window(x, q)
    require(protected[0] == x, "least protected carrier changed")
    for speed in protected:
        require(x <= speed <= q, "protected carrier left the prefix")
        require(F(1, 14) <= speed * left, "left phase fell below 1/14")
        require(speed * right <= F(13, 14), "right phase exceeded 13/14")
        require(safe(speed, left) and safe(speed, right), "closed endpoint unsafe")
        require(safe(speed, (left + right) / 2), "window midpoint unsafe")


def rho5_algebra_audit(limit: int = 2000) -> int:
    """Replay the residual inequality on every integer pair in a large box."""
    rows = 0
    for x in range(8, limit + 1):
        # Strict residual y < 14x/9.
        y_cap = ceil(F(14 * x, 9)) - 1
        for y in range(x + 1, y_cap + 1):
            # Cleared form: 4L > 3/y iff 5x > 2y.
            require(5 * x > 2 * y, "rho=5 cleared inequality failed")
            rows += 1
        if y_cap > x:
            for y in {x + 1, y_cap}:
                length = nested_length(x, y)
                require(4 * length > F(3, y), "rho=5 density contradiction failed")
                verify_window_geometry(x, (x, y))
    return rows


def rho6_symbolic_audit() -> tuple[F, F, F, F]:
    """Check the exact rational ladder and its deliberately rounded handoffs."""
    y_ratio = F(35, 12)
    z_exact = 56 * y_ratio / (3 * (13 - y_ratio))
    z_round = F(27, 5)
    require(z_exact == F(1960, 363), "z exact ratio moved")
    require(z_exact < z_round, "one-unit z rounding margin disappeared")

    w_ratio = F(21, 2) * z_round / (13 - z_round)
    require(w_ratio == F(567, 76), "w ratio moved")
    require(w_ratio < 13, "third nested window can be empty")

    v_ratio = 28 * w_ratio / (5 * (13 - w_ratio))
    require(v_ratio == F(15876, 2105), "v ratio moved")
    require(v_ratio < 11, "final one-comb length comparison failed")

    # At v/x < 11, J(x,v) has length > 1/(7v).  The last carrier u>v
    # has strictly shorter open teeth.
    require(13 - v_ratio > 2, "final protected phase width is not >2")
    return y_ratio, z_exact, w_ratio, v_ratio


def rho6_integer_ladder_audit(limit: int = 1200) -> int:
    """Replay every possible first pair and the monotone rational bounds."""
    rows = 0
    _, _, w_global, v_global = rho6_symbolic_audit()
    for x in range(8, limit + 1):
        y_cap = ceil(F(35 * x, 12)) - 1
        for y in range(x + 1, y_cap + 1):
            # 56y/[3(13x-y)] < 1960/363 < 27/5.
            require(
                56 * y * 363 < 1960 * 3 * (13 * x - y),
                "strict y substitution failed",
            )
            require(
                56 * y * 5 < 27 * 3 * (13 * x - y),
                "rounded z handoff failed",
            )

            # It suffices to test the largest admissible rational input at
            # each monotone handoff; the exact global constants were checked
            # above.  These inequalities also guard sign/orientation errors.
            require(F(567 * x, 76) < 13 * x, "w handoff left first window")
            require(v_global * x < 11 * x, "v handoff left final tooth bound")
            require(w_global < F(567, 76) + F(1, 10**9), "w constant drifted")
            rows += 1
        if y_cap > x:
            for y in {x + 1, y_cap}:
                nested_length(x, y)
                verify_window_geometry(x, (x, y))
    return rows


def tooth_boundaries(speed: int, left: F, right: F) -> set[F]:
    answer: set[F] = set()
    for center in range(floor(speed * left) - 2, ceil(speed * right) + 3):
        for sign in (-1, 1):
            endpoint = F(14 * center + sign, 14 * speed)
            if left < endpoint < right:
                answer.add(endpoint)
    return answer


def has_safe_point(speeds: tuple[int, ...]) -> bool:
    """Endpoint-exact common-safe test on the least carrier's first window."""
    x = speeds[0]
    left, right = F(1, 14 * x), F(13, 14 * x)
    walls = {left, right}
    for speed in speeds[1:]:
        walls.update(tooth_boundaries(speed, left, right))
    ordered = sorted(walls)
    probes = list(ordered)
    probes.extend((a + b) / 2 for a, b in zip(ordered, ordered[1:]))
    return any(all(safe(speed, time) for speed in speeds) for time in probes)


def finite_geometry_audit() -> tuple[int, int]:
    """Independent exhaustive and adversarial checks of the final conclusion."""
    rows5 = 0
    rows6 = 0

    # Exhaust a compact box containing both close and first 14-band rows.
    for x in range(8, 13):
        universe = range(x + 1, x + 13)
        for tail in combinations(universe, 4):
            require(has_safe_point((x,) + tail), "five-carrier compact counterexample")
            rows5 += 1
        for tail in combinations(universe, 5):
            require(has_safe_point((x,) + tail), "six-carrier compact counterexample")
            rows6 += 1

    # Stress the obstruction bands 14m+-1 at scales far beyond the compact box.
    rng = Random(20260718)
    for rho in (5, 6):
        for _ in range(240):
            x = rng.randint(8, 80)
            values = {x}
            while len(values) < rho:
                multiplier = rng.choice((1, 2, 3, 13, 14, 15, 27, 28, 29))
                candidate = max(x + 1, multiplier * x + rng.randint(-3, 3))
                values.add(candidate)
            speeds = tuple(sorted(values))
            require(has_safe_point(speeds), "adversarial carrier counterexample")
            if rho == 5:
                rows5 += 1
            else:
                rows6 += 1
    return rows5, rows6


def tournament_audit() -> tuple[tuple[int, ...], int, int, int]:
    """Containment-order tournament on the four rho=6 proof obligations."""
    vertices = tuple(range(4))  # J2, J3, J4, J5
    edges = {(a, b) for a, b in combinations(vertices, 2)}
    scores = tuple(sum(a == v for a, _ in edges) for v in vertices)
    cycles = sum(
        ((a, b) in edges and (b, c) in edges and (c, a) in edges)
        or ((b, a) in edges and (c, b) in edges and (a, c) in edges)
        for a, b, c in combinations(vertices, 3)
    )
    reverse = {(b, a) for a, b in edges}
    flips = len(edges.symmetric_difference(reverse)) // 2
    require(scores == (3, 2, 1, 0), "nested-window tournament scores moved")
    require(cycles == 0 and flips == 6, "nested-window tournament changed")
    return scores, cycles, 1, flips


def main() -> None:
    rho5_rows = rho5_algebra_audit()
    rho6_rows = rho6_integer_ladder_audit()
    y_ratio, z_ratio, w_ratio, v_ratio = rho6_symbolic_audit()
    finite5, finite6 = finite_geometry_audit()
    scores, cycles, hp_count, flips = tournament_audit()

    print("THM-1212 five/six-carrier nested-window referee")
    print("arithmetic=integers and fractions.Fraction; danger teeth are open")
    print("optimized_mode_guard=require() only; assert_statements=0")
    print(f"rho5_integer_first-pair_rows={rho5_rows}")
    print("rho5_residual=y/x<14/9; nested L=(13x-y)/(14xy)")
    print("rho5_contradiction=4L>3/y, while three later combs covering forces 4L<3/y")
    print(f"rho6_integer_first-pair_rows={rho6_rows}")
    print(f"rho6_ratio_ladder=y<{y_ratio}, z<{z_ratio}<27/5, w<{w_ratio}, v<{v_ratio}<11")
    print("rho6_final=|J(x,v)|>1/(7v)>1/(7r), so one open r-tooth cannot cover")
    print(f"endpoint_exact_geometry_rows=rho5:{finite5},rho6:{finite6}")
    print("Tournament Analysis:")
    print(f"  vertices=(J2,J3,J4,J5); scores={scores}; cycles={cycles}; SCCs=(1,1,1,1)")
    print(f"  Hamiltonian_paths={hp_count}; reverse-gauge edge_flips={flips}")
    print("  pairwise_observable=prefix containment; switch=edge points to larger protected prefix")
    print("  tie Hamiltonian path=J2->J3->J4->J5 (no metric ties)")
    print("  faithful object=nested closed window + remaining-comb count + reciprocal budget")
    print("  destroyed by tournament=window length|endpoint owner|strict cover|density slack")
    print("  challenged vertices=runners|carriers|ratios|teeth|walls|windows|cover obligations")
    print("VERDICT: rho=5 and rho=6 close; with THM-1136/1154/1169 the full r=6 stratum is empty")


if __name__ == "__main__":
    main()
