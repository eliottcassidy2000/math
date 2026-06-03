#!/usr/bin/env python3
"""S593: dual pigeonhole cap face for the C' multiple-of-n residual.

For S = S' union {v=nw}, the v-danger arcs have total length 2/n.  More
locally, in every n-clock cell [r/n,(r+1)/n) they have exactly w caps, each of
length 2/(n^2 w), hence total capacity 2/n^2.  Therefore S is loose whenever
G(S') has more than 2/n^2 safe mass in at least one n-clock cell.

This is the aggregate dual to B': B' finds one component too long for one cap;
the cell-cap test finds too many short components in one clock cell for all caps
in that cell.  It is a proved criterion, not an enumeration claim.
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
            out.append(Component(a, length, tuple(pts[a]), tuple(pts[b])))
    return out


def primitive(vs: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in vs))


def endpoint_pinned(owner: tuple[int, int, int], w: int, n: int) -> bool:
    u, k, eps = owner
    return (w * (k * n + eps)) % u == 0


def local_cover_exit(component: Component, n: int, w: int) -> str | None:
    """S581 local criteria proving this component cannot fit one v=nw cap."""
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
                cell_lo = F(r, n)
                cell_hi = F(r + 1, n)
                loads[r] += overlap(lo, hi, cell_lo, cell_hi)
    return tuple(loads)


def cap_face_exit(comps: list[Component], n: int) -> tuple[bool, F, int, tuple[F, ...]]:
    """Return whether an n-clock cell exceeds its v=nw danger capacity."""
    loads = n_cell_loads(comps, n)
    cap = F(2, n * n)
    surpluses = tuple(load - cap for load in loads)
    best_r = max(range(n), key=lambda r: (surpluses[r], loads[r], -r))
    best = surpluses[best_r]
    return best > 0, best, best_r, loads


def safe_measure_from_components(comps: list[Component]) -> F:
    return sum((comp.length for comp in comps), F(0))


def route_row(vs: tuple[int, ...], n: int) -> dict[str, object]:
    mults = [v for v in vs if v % n == 0]
    if not mults:
        return {"route": "no_multiple"}
    v = mults[0]
    w = v // n
    sp = tuple(x for x in vs if x != v)
    comps = components(sp, n)
    cap_hit, cap_surplus, cap_cell, loads = cap_face_exit(comps, n)
    global_surplus = safe_measure_from_components(comps) - F(2, n)

    local_reasons = Counter()
    first_local = None
    for comp in comps:
        reason = local_cover_exit(comp, n, w)
        if reason is not None:
            local_reasons[reason] += 1
            if first_local is None:
                first_local = reason

    if cap_hit:
        route = "dual_cell_cap"
    elif global_surplus > 0:
        route = "global_cap_mass"
    elif first_local is not None:
        route = first_local
    else:
        route = "cap_face_residual"

    return {
        "route": route,
        "v": v,
        "w": w,
        "components": len(comps),
        "safe_mass_sp": safe_measure_from_components(comps),
        "global_surplus": global_surplus,
        "cap_surplus": cap_surplus,
        "cap_cell": cap_cell,
        "loads": loads,
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


def exact_m(vs: tuple[int, ...]) -> tuple[F, F]:
    cands = {F(0), F(1, 2)}
    for a, b in combinations(vs, 2):
        for d in (a + b, abs(a - b)):
            if d:
                for j in range(d + 1):
                    cands.add(F(j, d))
    for v in vs:
        for j in range(2 * v + 1):
            cands.add(F(j, 2 * v))
    best = F(-1)
    best_t = F(0)
    for t in cands:
        m = min(dist(v * t) for v in vs)
        if m > best:
            best = m
            best_t = t
    return best, best_t


def cell_tournament(loads: tuple[F, ...]) -> str:
    n = len(loads)
    scores = Counter({r: 0 for r in range(n)})
    c3 = 0

    def beats(a: int, b: int) -> bool:
        return (loads[a], -a) > (loads[b], -b)

    for a, b in combinations(range(n), 2):
        winner = a if beats(a, b) else b
        scores[winner] += 1
    for a, b, c in combinations(range(n), 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    score_hist = Counter(scores.values())
    order = sorted(range(n), key=lambda r: (loads[r], -r), reverse=True)
    edge_flips_vs_clock_order = sum(
        1 for a, b in combinations(range(n), 2) if beats(b, a)
    )
    return (
        f"score_hist={dict(sorted(score_hist.items()))}; directed_3_cycles={c3}; "
        f"SCCs={(1,) * n}; Hamiltonian_paths=1; "
        f"edge_flips_vs_clock_order={edge_flips_vs_clock_order}; "
        f"tie_HP={' -> '.join(str(r) for r in order)}"
    )


def fmt_frac(x: F) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def print_named_rows() -> None:
    named = [
        ("apex_AP_replace_13_n14", 14, tuple(range(1, 13)) + (14,)),
        ("unit_shift_AP_n14", 14, tuple(range(2, 15))),
        ("near_AP_double_14_n14", 14, tuple(range(1, 12)) + (13, 28)),
    ]
    print("NAMED MULTIPLE-OF-14 ROWS")
    for name, n, row in named:
        info = route_row(row, n)
        m, t = exact_m(row)
        print(
            f"{name}: route={info['route']}; M={fmt_frac(m)} at t={fmt_frac(t)}; "
            f"safe_mass(G(S'))={fmt_frac(info['safe_mass_sp'])}; "
            f"cell_cap_surplus={fmt_frac(info['cap_surplus'])} in cell {info['cap_cell']}; "
            f"global_surplus={fmt_frac(info['global_surplus'])}"
        )


def main() -> None:
    print("S593 dual pigeonhole cap face for C' multiple-of-n rows")
    print()
    print("PROVED LEMMA")
    print("For v=nw, D_v has exactly w caps in each n-clock cell, each of length 2/(n^2 w).")
    print("Thus each n-clock cell has total cap capacity 2/n^2.")
    print("If mu(G(S') cap [r/n,(r+1)/n)) > 2/n^2 for some r, then G(S') is not contained in D_v, so S is loose.")
    print()

    aggregate = Counter()
    cap_added_after_local = Counter()
    best_n14: tuple[F, tuple[int, ...], dict[str, object]] | None = None

    for n in [6, 8, 10, 12, 14]:
        rows = sample_configs(n, 2500, 593 + n)
        routes = Counter()
        overlap = Counter()
        best_surplus = F(-1)
        best_row = None
        best_info = None
        for row in rows:
            info = route_row(row, n)
            route = str(info["route"])
            routes[route] += 1
            aggregate[route] += 1
            local = info.get("first_local")
            cap = route == "dual_cell_cap"
            if cap and local is None:
                cap_added_after_local[n] += 1
            if cap and local is not None:
                overlap["cap_and_local"] += 1
            elif cap:
                overlap["cap_only"] += 1
            elif local is not None:
                overlap["local_only"] += 1
            else:
                overlap["neither"] += 1
            surplus = info.get("cap_surplus", F(-1))
            if isinstance(surplus, F) and surplus > best_surplus:
                best_surplus = surplus
                best_row = row
                best_info = info
        proved = len(rows) - routes["cap_face_residual"]
        print(
            f"n={n:2d}: rows={len(rows)}; routed={proved} ({100*proved/max(len(rows),1):.2f}%); "
            f"dual_cell_cap={routes['dual_cell_cap']}; global_cap_mass={routes['global_cap_mass']}; "
            f"cap_face_residual={routes['cap_face_residual']}"
        )
        print(
            "  routes: "
            + ", ".join(f"{k}={routes[k]}" for k in sorted(routes))
        )
        print(
            "  cap/local overlap: "
            + ", ".join(f"{k}={overlap[k]}" for k in sorted(overlap))
        )
        if best_row is not None and best_info is not None:
            print(
                f"  max cell surplus row={best_row}; surplus={fmt_frac(best_surplus)}; "
                f"cell={best_info['cap_cell']}; components={best_info['components']}"
            )
            if n == 14:
                best_n14 = (best_surplus, best_row, best_info)
    print()
    print("AGGREGATE ROUTES")
    for key in sorted(aggregate):
        print(f"{key}: {aggregate[key]}")
    print()
    print("CAP ADDS BEFORE LOCAL OWNER DESCENT")
    for n in sorted(cap_added_after_local):
        print(f"n={n}: {cap_added_after_local[n]} rows had a cell-cap proof before any S581 local component exit")
    if not cap_added_after_local:
        print("none in this deterministic sample")
    print()
    print_named_rows()
    print()
    if best_n14 is not None:
        _, row, info = best_n14
        print("N=14 CELL-LOAD TOURNAMENT ANALYSIS")
        print(f"representative row={row}")
        print("vertices: the 14 n-clock cells, not runners or components")
        print("pair observable: (cell safe load, then clock order); switch points toward larger load/surplus")
        print("tie Hamiltonian path: cells sorted by decreasing load, clock order breaks ties")
        print(cell_tournament(info["loads"]))
        cap = F(2, 14 * 14)
        top = sorted(range(14), key=lambda r: (info["loads"][r], -r), reverse=True)[:5]
        print(
            "top cell loads: "
            + ", ".join(f"cell {r}: {fmt_frac(info['loads'][r])} (cap {fmt_frac(cap)})" for r in top)
        )
    print()
    print("ASSUMPTION CHALLENGE")
    print("Considered vertices: runners, gaps, components, cap centres, n-clock cells, owners, and proof obligations.")
    print("Chosen vertices: n-clock cells for the cap lemma, plus proof routes for coverage.")
    print("Preserved predicate: whether D_{nw} has enough cap capacity to cover G(S') in every primary clock cell.")
    print("Destroyed data: exact owner congruences and component adjacency inside a cell.")
    print("Challenged assumption: a cover proof must find one component too long; aggregate short components can overload one cell.")


if __name__ == "__main__":
    main()
