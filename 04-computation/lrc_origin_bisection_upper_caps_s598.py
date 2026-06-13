#!/usr/bin/env python3
"""S598: origin-bisection laws as upper/lower cap certificates for C'.

For S = S' union {v=nw}, each v-danger cap is bisected by its center
j/(nw).  Lemma F used this locally: if a small owner pins one endpoint to the
center, the component must fit in one half-cap of length 1/(n^2 w).

This script records the aggregate version.  In each primary n-clock cell there
are w upper half-caps and w lower half-caps, each half of length 1/(n^2 w), so
the total upper capacity and total lower capacity are both 1/n^2.  If the
small-owner congruence law forces more than 1/n^2 safe mass into the upper (or
lower) halves of one primary cell, then D_{nw} cannot cover G(S') and S is
loose.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import floor, gcd
import random


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


@dataclass(frozen=True)
class Component:
    a: F
    b: F
    length: F
    left_owners: tuple[tuple[int, int, int], ...]
    right_owners: tuple[tuple[int, int, int], ...]

    def left_safe_owners(self) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.left_owners if o[2] == 1)

    def right_safe_owners(self) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.right_owners if o[2] == -1)

    def small_left(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.left_safe_owners() if o[0] < n)

    def small_right(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.right_safe_owners() if o[0] < n)


def components(sp: tuple[int, ...], n: int) -> list[Component]:
    """Exact components of G(S')={t: ||ut||>1/n for all u in sp}."""
    threshold = F(1, n)
    pts: dict[F, list[tuple[int, int, int]]] = {}
    for u in sp:
        for k in range(u + 1):
            for eps in (1, -1):
                pts.setdefault(F(k * n + eps, n * u) % 1, []).append((u, k, eps))
    order = sorted(pts)
    out: list[Component] = []
    for i, a in enumerate(order):
        b = order[(i + 1) % len(order)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(dist(u * mid) > threshold for u in sp):
            out.append(Component(a, b, length, tuple(pts[a]), tuple(pts[b])))
    return out


def primitive(vs: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in vs))


def endpoint_pinned(owner: tuple[int, int, int], w: int, n: int) -> bool:
    """Small-owner endpoint lies on the v=nw center lattice iff this holds."""
    u, k, eps = owner
    return (w * (k * n + eps)) % u == 0


def right_cell_of(x: F, n: int) -> int:
    """Primary n-cell containing a tiny move to the right of x."""
    y = (x % 1) * n
    if y.denominator == 1:
        return y.numerator % n
    return floor(y) % n


def left_cell_of(x: F, n: int) -> int:
    """Primary n-cell containing a tiny move to the left of x."""
    y = (x % 1) * n
    if y.denominator == 1:
        return (y.numerator - 1) % n
    return floor(y) % n


def split_interval(a: F, length: F) -> list[tuple[F, F]]:
    b = a + length
    if b <= 1:
        return [(a, b)]
    return [(a, F(1)), (F(0), b - 1)]


def overlap(lo1: F, hi1: F, lo2: F, hi2: F) -> F:
    lo = max(lo1, lo2)
    hi = min(hi1, hi2)
    return hi - lo if hi > lo else F(0)


def n_cell_loads(comps: list[Component], n: int) -> tuple[F, ...]:
    loads = [F(0) for _ in range(n)]
    for comp in comps:
        for lo, hi in split_interval(comp.a, comp.length):
            for r in range(n):
                loads[r] += overlap(lo, hi, F(r, n), F(r + 1, n))
    return tuple(loads)


def dual_cell_cap_exit(comps: list[Component], n: int) -> tuple[bool, int, F]:
    loads = n_cell_loads(comps, n)
    cap = F(2, n * n)
    best_r = max(range(n), key=lambda r: (loads[r] - cap, loads[r], -r))
    return loads[best_r] > cap, best_r, loads[best_r] - cap


def local_cover_exit(component: Component, n: int, w: int) -> str | None:
    full_radius = F(2, n * n * w)
    half_radius = F(1, n * n * w)
    if component.length > full_radius:
        return "Bprime_long_component"

    small_left = component.small_left(n)
    small_right = component.small_right(n)
    if small_left and small_right:
        return "LemmaC_two_small_owners"
    for owner in small_left + small_right:
        if not endpoint_pinned(owner, w, n):
            return "LemmaE_small_pin_off_lattice"
    if (small_left or small_right) and component.length > half_radius:
        return "LemmaF_half_cap"
    return None


def origin_bisection_loads(
    comps: list[Component], n: int, w: int
) -> tuple[tuple[F, ...], tuple[F, ...], Counter[str]]:
    """Return aggregate upper/lower forced loads per primary cell.

    A pinned small left owner forces the component into the upper half-cap based
    at its left endpoint.  A pinned small right owner forces it into the lower
    half-cap ending at its right endpoint.  Off-lattice small owners are not
    counted here because they are already Lemma E local exits.
    """
    upper = [F(0) for _ in range(n)]
    lower = [F(0) for _ in range(n)]
    shapes: Counter[str] = Counter()
    for comp in comps:
        pinned_left = [o for o in comp.small_left(n) if endpoint_pinned(o, w, n)]
        pinned_right = [o for o in comp.small_right(n) if endpoint_pinned(o, w, n)]
        if pinned_left:
            upper[right_cell_of(comp.a, n)] += comp.length
            shapes["upper_forced_components"] += 1
        if pinned_right:
            lower[left_cell_of(comp.b, n)] += comp.length
            shapes["lower_forced_components"] += 1
        if pinned_left and pinned_right:
            shapes["two_sided_forced_components"] += 1
    return tuple(upper), tuple(lower), shapes


def origin_bisection_exit(
    comps: list[Component], n: int, w: int
) -> tuple[bool, str, int, F, tuple[F, ...], tuple[F, ...], Counter[str]]:
    upper, lower, shapes = origin_bisection_loads(comps, n, w)
    cap = F(1, n * n)
    options = []
    for r in range(n):
        options.append(("upper", r, upper[r] - cap, upper[r]))
        options.append(("lower", r, lower[r] - cap, lower[r]))
    side, cell, surplus, _load = max(options, key=lambda x: (x[2], x[3], x[0], -x[1]))
    return surplus > 0, side, cell, surplus, upper, lower, shapes


def route_row(vs: tuple[int, ...], n: int) -> dict[str, object]:
    mults = [v for v in vs if v % n == 0]
    if not mults:
        return {"route": "no_multiple"}
    v = mults[0]
    w = v // n
    sp = tuple(x for x in vs if x != v)
    comps = components(sp, n)

    total_hit, total_cell, total_surplus = dual_cell_cap_exit(comps, n)
    bis_hit, bis_side, bis_cell, bis_surplus, upper, lower, shapes = origin_bisection_exit(
        comps, n, w
    )

    local_reasons = Counter()
    first_local = None
    for comp in comps:
        reason = local_cover_exit(comp, n, w)
        if reason:
            local_reasons[reason] += 1
            if first_local is None:
                first_local = reason

    if total_hit:
        route = "dual_total_cell_cap"
    elif bis_hit:
        route = f"origin_bisection_{bis_side}_cap"
    elif first_local is not None:
        route = first_local
    else:
        route = "residual"

    return {
        "route": route,
        "w": w,
        "components": len(comps),
        "total_hit": total_hit,
        "total_cell": total_cell,
        "total_surplus": total_surplus,
        "bis_hit": bis_hit,
        "bis_side": bis_side,
        "bis_cell": bis_cell,
        "bis_surplus": bis_surplus,
        "upper": upper,
        "lower": lower,
        "shapes": shapes,
        "local_reasons": local_reasons,
        "first_local": first_local,
    }


def sample_configs(n: int, trials: int, seed: int) -> list[tuple[int, ...]]:
    rng = random.Random(seed)
    m = n - 1
    rows: list[tuple[int, ...]] = []
    for _ in range(trials):
        others = set(rng.sample([x for x in range(1, n + 10) if x % n != 0], m - 1))
        w = rng.randint(1, 4)
        v = n * w
        if v in others:
            continue
        row = primitive(tuple(sorted(others | {v})))
        if len(row) == m and any(x % n == 0 for x in row):
            rows.append(row)
    return rows


def fmt_frac(x: F) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def tournament_fingerprint(route_totals: Counter[str]) -> str:
    vertices = [
        ("dual_total_cell_cap", 5, route_totals["dual_total_cell_cap"]),
        (
            "origin_bisection_upper_cap",
            4,
            route_totals["origin_bisection_upper_cap"],
        ),
        (
            "origin_bisection_lower_cap",
            4,
            route_totals["origin_bisection_lower_cap"],
        ),
        ("LemmaF_half_cap", 3, route_totals["LemmaF_half_cap"]),
        ("LemmaE_small_pin_off_lattice", 2, route_totals["LemmaE_small_pin_off_lattice"]),
        ("LemmaC_two_small_owners", 2, route_totals["LemmaC_two_small_owners"]),
        ("Bprime_long_component", 1, route_totals["Bprime_long_component"]),
        ("residual", 0, route_totals["residual"]),
    ]

    def beats(a: tuple[str, int, int], b: tuple[str, int, int]) -> bool:
        return (a[1], a[2], a[0]) > (b[1], b[2], b[0])

    scores = Counter({v[0]: 0 for v in vertices})
    c3 = 0
    for a, b in combinations(vertices, 2):
        scores[(a if beats(a, b) else b)[0]] += 1
    for a, b, c in combinations(vertices, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    order = sorted(vertices, key=lambda x: (x[1], x[2], x[0]), reverse=True)
    return (
        f"score_hist={dict(sorted(Counter(scores.values()).items()))}; "
        f"directed_3_cycles={c3}; SCCs={(1,) * len(vertices)}; "
        f"Hamiltonian_paths=1; tie_HP={' > '.join(v[0] for v in order)}"
    )


def print_named_rows() -> None:
    named = [
        ("unit_shift_AP_n14", 14, tuple(range(2, 15))),
        ("apex_AP_replace_13_n14", 14, tuple(range(1, 13)) + (14,)),
        ("near_AP_double_14_n14", 14, tuple(range(1, 12)) + (13, 28)),
    ]
    print("NAMED ROWS")
    for name, n, row in named:
        info = route_row(row, n)
        print(
            f"{name}: route={info['route']} total_surplus={fmt_frac(info['total_surplus'])} "
            f"bis={info['bis_side']}@{info['bis_cell']} surplus={fmt_frac(info['bis_surplus'])} "
            f"forced={dict(info['shapes'])}"
        )
    print()


def main() -> None:
    print("S598 origin-bisection upper/lower cap certificates for Cprime")
    print()
    print("PROVED CERTIFICATE")
    print("If a small left owner pins a component endpoint to j/(nw), any cover puts")
    print("that component in the upper half-cap.  Per primary n-cell, total upper")
    print("capacity is 1/n^2.  The analogous right-endpoint law gives lower caps.")
    print()
    print_named_rows()

    route_totals: Counter[str] = Counter()
    aggregate_shapes: Counter[str] = Counter()
    print("DETERMINISTIC SAMPLE")
    print(" n  rows  total_cap  bisect_route  bisect_any  local_after  residual")
    for n in [6, 8, 10, 12, 14]:
        rows = sample_configs(n, 5000, 598 + n)
        totals = Counter()
        bis_any = 0
        local_after = 0
        for row in rows:
            info = route_row(row, n)
            route = str(info["route"])
            totals[route] += 1
            route_totals[route] += 1
            if info["bis_hit"]:
                bis_any += 1
            if (
                not info["total_hit"]
                and not info["bis_hit"]
                and info["first_local"] is not None
            ):
                local_after += 1
            aggregate_shapes.update(info["shapes"])
        print(
            f"{n:2d} {len(rows):5d} {totals['dual_total_cell_cap']:10d} "
            f"{totals['origin_bisection_upper_cap'] + totals['origin_bisection_lower_cap']:12d} "
            f"{bis_any:10d} {local_after:11d} {totals['residual']:8d}"
        )
    print()
    print("ROUTE TOTALS")
    for key in sorted(route_totals):
        print(f"{key}: {route_totals[key]}")
    print()
    print("FORCED HALF-CAP SHAPES")
    for key in sorted(aggregate_shapes):
        print(f"{key}: {aggregate_shapes[key]}")
    print()
    print("TOURNAMENT ANALYSIS")
    print("vertices: cap certificates and local owner criteria, not runners or arcs")
    print("pair observable: (certificate_rank, sampled_route_count, name)")
    print("switch: stronger aggregate certificate beats later/local criterion")
    print(tournament_fingerprint(route_totals))
    print()
    print("ASSUMPTION CHALLENGE")
    print("Alternative vertices considered: runners, cap centers, components, owners,")
    print("primary n-cells, upper/lower half-caps, and proof obligations.")
    print("Chosen vertices are proof certificates.  The quotient preserves whether")
    print("owner bisection forces mass into one-sided cap capacity.  It destroys")
    print("exact phase order inside a cell; any mixed fiber must be lifted by")
    print("endpoint owner, side, and residue labels.")


if __name__ == "__main__":
    main()
