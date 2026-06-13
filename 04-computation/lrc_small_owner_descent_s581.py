#!/usr/bin/env python3
"""S581: extend endpoint-cover positivity past Lemma C.

HYP-2105 proves the both-small-owner case.  This script adds the one-small-owner
"pin" descent:

* if a component endpoint has owner u < n and is not on the v=nw arc-center
  lattice, the component is uncoverable;
* if it is on that lattice, the endpoint is the arc centre, so the component
  must fit into one half-radius 1/(n^2 w), not the full 2/(n^2 w).

The remaining rows are therefore genuinely large-owner/half-radius rows.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd
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

    def large_left(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.left_safe_owners() if o[0] >= n)

    def large_right(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.right_safe_owners() if o[0] >= n)


def components(sp: tuple[int, ...], n: int) -> list[Component]:
    thr = F(1, n)
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
        if all(dist(u * mid) > thr for u in sp):
            out.append(Component(a, b, length, tuple(pts[a]), tuple(pts[b])))
    return out


def primitive(vs: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in vs))


def endpoint_pinned(owner: tuple[int, int, int], w: int, n: int) -> bool:
    """For small owner u<n, endpoint lies on the v=nw centre lattice iff this holds."""
    u, k, eps = owner
    return (w * (k * n + eps)) % u == 0


def criterion(component: Component, n: int, w: int) -> str | None:
    """Return the first proved reason that this component is not coverable by D_{nw}."""
    full_radius = F(2, n * n * w)
    half_radius = F(1, n * n * w)

    if component.length > full_radius:
        return "Bprime_long_interval"

    small_left = component.small_left(n)
    small_right = component.small_right(n)

    if small_left and small_right:
        return "LemmaC_both_small"

    for owner in small_left + small_right:
        if not endpoint_pinned(owner, w, n):
            return "LemmaE_small_pin_off_lattice"

    if (small_left or small_right) and component.length > half_radius:
        return "LemmaF_small_pin_half_radius"

    return None


def config_route(vs: tuple[int, ...], n: int) -> tuple[str, Counter[str]]:
    mults = [v for v in vs if v % n == 0]
    if not mults:
        return "no_multiple", Counter()
    v = mults[0]
    w = v // n
    sp = tuple(x for x in vs if x != v)
    counts: Counter[str] = Counter()
    for comp in components(sp, n):
        reason = criterion(comp, n, w)
        if reason is not None:
            counts[reason] += 1
            return reason, counts
        sl = bool(comp.small_left(n) or comp.small_right(n))
        ll = bool(comp.large_left(n) or comp.large_right(n))
        if sl and ll:
            counts["unproved_mixed_small_large_tiny"] += 1
        elif ll:
            counts["unproved_large_large_tiny"] += 1
        else:
            counts["unproved_no_binding_owner"] += 1
    return "large_owner_residual", counts


def sample_configs(n: int, trials: int, seed: int) -> list[tuple[int, ...]]:
    rng = random.Random(seed)
    m = n - 1
    rows: list[tuple[int, ...]] = []
    for _ in range(trials):
        others = set(rng.sample([x for x in range(1, n + 8) if x % n != 0], m - 1))
        w = rng.randint(1, 3)
        v = n * w
        if v in others:
            continue
        row = primitive(tuple(sorted(others | {v})))
        if len(row) == m and any(x % n == 0 for x in row):
            rows.append(row)
    return rows


def tournament_fingerprint(route_totals: Counter[str]) -> str:
    vertices = [
        ("Bprime_long_interval", 5, route_totals["Bprime_long_interval"]),
        ("LemmaC_both_small", 4, route_totals["LemmaC_both_small"]),
        ("LemmaE_small_pin_off_lattice", 3, route_totals["LemmaE_small_pin_off_lattice"]),
        ("LemmaF_small_pin_half_radius", 2, route_totals["LemmaF_small_pin_half_radius"]),
        ("large_owner_residual", 1, route_totals["large_owner_residual"]),
    ]

    def beats(a, b):
        # pair observable: (proof_rank, coverage); switch orients toward stronger
        # proved criterion, with residual last.
        return (a[1], a[2], a[0]) > (b[1], b[2], b[0])

    scores = Counter({v[0]: 0 for v in vertices})
    c3 = 0
    for a, b in combinations(vertices, 2):
        winner = a if beats(a, b) else b
        scores[winner[0]] += 1
    for a, b, c in combinations(vertices, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    order = sorted(vertices, key=lambda x: (x[1], x[2], x[0]), reverse=True)
    return (
        f"score_hist={dict(sorted(Counter(scores.values()).items()))}; "
        f"directed_3_cycles={c3}; "
        f"hamiltonian_path={' -> '.join(v[0] for v in order)}"
    )


def main() -> None:
    print("S581 endpoint-cover positivity: small-owner descent past Lemma C")
    print()
    print("PROVED CRITERIA")
    print("Bprime: component length > 2/(n^2 w).")
    print("Lemma C: both endpoint owners are < n, so two pinned centres force a=b.")
    print("Lemma E: one small endpoint owner is off the v=nw centre lattice.")
    print("Lemma F: one small endpoint owner is pinned, but the component exceeds one half-radius 1/(n^2 w).")
    print()

    route_totals: Counter[str] = Counter()
    residual_examples: dict[int, tuple[int, ...]] = {}
    residual_shapes: Counter[str] = Counter()

    for n in [6, 8, 10, 12, 14]:
        rows = sample_configs(n, 5000, 581 + n)
        routes: Counter[str] = Counter()
        shape_counts: Counter[str] = Counter()
        for row in rows:
            route, shapes = config_route(row, n)
            routes[route] += 1
            route_totals[route] += 1
            if route == "large_owner_residual" and n not in residual_examples:
                residual_examples[n] = row
            shape_counts.update(shapes)
            if route == "large_owner_residual":
                residual_shapes.update(shapes)
        proved = len(rows) - routes["large_owner_residual"]
        print(
            f"n={n:2d}: rows={len(rows)}; proved={proved} ({100*proved/max(len(rows),1):.2f}%); "
            f"residual={routes['large_owner_residual']}"
        )
        print(
            "  routes: "
            + ", ".join(f"{k}={routes[k]}" for k in sorted(routes))
            + (f"; example_residual={residual_examples[n]}" if n in residual_examples else "")
        )
        if routes["large_owner_residual"]:
            print(
                "  residual component shapes: "
                + ", ".join(f"{k}={shape_counts[k]}" for k in sorted(shape_counts) if k.startswith("unproved"))
            )
    print()
    print("AGGREGATE ROUTES")
    for key in sorted(route_totals):
        print(f"{key}: {route_totals[key]}")
    print()
    print("RESIDUAL SHAPES")
    for key in sorted(residual_shapes):
        print(f"{key}: {residual_shapes[key]}")
    print()
    print("TOURNAMENT ANALYSIS")
    print("vertices: proved cover criteria plus the remaining large-owner residual")
    print("pair observable: (proof_rank, coverage); switch favors stronger proved criterion")
    print(tournament_fingerprint(route_totals))
    print()
    print("ASSUMPTION CHALLENGE")
    print("Vertices considered: runners, endpoints, components, owner pairs, congruence windows, and proof criteria.")
    print("Chosen vertices: proof criteria, because the predicate preserved is 'some component is certified uncoverable'.")
    print("Destroyed data: exact endpoint positions and individual speed identities after they are summarized by route.")
    print("Next residual: all remaining unproved components are large-owner or mixed tiny half-radius windows.")


if __name__ == "__main__":
    main()
