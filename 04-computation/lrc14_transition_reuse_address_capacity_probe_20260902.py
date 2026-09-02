#!/usr/bin/env python3
"""Exact address-capacity audit for THM-4335 minority renewals.

For the anchor {2h}, this program distinguishes an addressed proper crossing

    (u,n) -> (v,m)

from the quotient which remembers only the ordered speed pair u -> v.  It
checks the determinant parametrization, the exact anchor-component congruence,
the orbit-count formula, reflection/owner-bit symmetry, and every transition
used by the two frozen h=420 controls.  It also prints two hostile families:
an all-or-none gcd/CRT family and a coprime many-determinant family.

All decisions use integer or Fraction arithmetic.  This is not an LRC census.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, floor, gcd


@dataclass(frozen=True)
class Tooth:
    w: int
    n: int
    left: Fraction
    right: Fraction


@dataclass(frozen=True)
class Event:
    u: int
    n: int
    v: int
    m: int
    d: int
    determinant: int
    q: int
    component: int | None
    owner_bit: int | None
    entry_remainder: int
    exit_remainder: int


def frac_part(x: Fraction) -> Fraction:
    return x - floor(x)


def ceil_fraction(x: Fraction) -> int:
    return -floor(-x)


def tooth(w: int, n: int) -> Tooth:
    return Tooth(
        w,
        n,
        Fraction(14 * n - 1, 14 * w),
        Fraction(14 * n + 1, 14 * w),
    )


def anchor_component(h: int, k: int) -> tuple[Fraction, Fraction]:
    assert 0 <= k < 2 * h
    return Fraction(14 * k + 1, 28 * h), Fraction(14 * k + 13, 28 * h)


def determinant_indices(u: int, v: int) -> tuple[int, ...]:
    """Normalized d=D/g for proper crossings u -> v."""
    g = gcd(u, v)
    U, V = u // g, v // g
    lo = floor(abs(U - V) / 14) + 1
    hi = floor(Fraction(U + V - 1, 14))
    return tuple(range(lo, hi + 1))


def event_for_addresses(h: int, u: int, n: int, v: int, m: int) -> Event:
    assert u > 0 and v > 0 and u % 2 == v % 2 == 1 and u != v
    D = u * m - v * n
    assert abs(u - v) < 14 * D < u + v
    g = gcd(u, v)
    assert D % g == 0
    q = u + v - 14 * D
    assert q > 0 and q % (2 * g) == 0

    # A consecutive greedy use forces the whole handoff window
    # [left(T(v,m)), right(T(u,n))] into one closed anchor component: v was
    # not eligible at the previous frontier, and u does not reach past R_k.
    Y = 2 * h * (14 * m - 1)
    X = 2 * h * (14 * n + 1)
    k_entry, r_entry = divmod(Y, 14 * v)
    k_exit, r_exit = divmod(X, 14 * u)
    if (
        k_entry == k_exit
        and v <= r_entry < 13 * v
        and u < r_exit <= 13 * u
    ):
        assert 0 <= k_exit < 2 * h
        epsilon = int(k_exit >= h)
        N_u = 2 * n - epsilon * u
        N_v = 2 * m - epsilon * v
        assert N_u % 2 == N_v % 2 == epsilon
        assert u * N_v - v * N_u == 2 * D
        component: int | None = k_exit
        owner_bit: int | None = epsilon
    else:
        component = None
        owner_bit = None
    return Event(
        u,
        n,
        v,
        m,
        D // g,
        D,
        q,
        component,
        owner_bit,
        r_entry,
        r_exit,
    )


def addressed_events(h: int, u: int, v: int) -> tuple[Event, ...]:
    """All proper addressed crossings u -> v in one circle period."""
    g = gcd(u, v)
    U, V = u // g, v // g
    inverse = pow(V, -1, U) if U > 1 else 0
    events: list[Event] = []
    for d in determinant_indices(u, v):
        n0 = (-d * inverse) % U if U > 1 else 0
        m0 = (d + V * n0) // U
        assert U * m0 - V * n0 == d
        for s in range(g):
            n = n0 + U * s
            m = m0 + V * s
            assert 0 <= n < u and 0 <= m <= v
            events.append(event_for_addresses(h, u, n, v, m))
    assert len({(e.n, e.m) for e in events}) == len(events)
    return tuple(events)


def orbit_count_for_d(h: int, u: int, v: int, d: int) -> int:
    """O(1) count of addressed d-events whose exit wall is anchor-safe."""
    g = gcd(u, v)
    U, V = u // g, v // g
    inverse = pow(V, -1, U) if U > 1 else 0
    n0 = (-d * inverse) % U if U > 1 else 0
    m0 = (d + V * n0) // U
    delta = gcd(g, 2 * h)
    G = g // delta

    # Under t -> 2ht, the g handoff windows form delta copies of a shifted
    # 1/G-grid.  Their common scaled length is lambda.  Containment in an
    # anchor component is start phase in [1/14,13/14-lambda].
    rho = frac_part(Fraction(h * (14 * m0 - 1), 7 * delta * V))
    q = g * (U + V - 14 * d)
    lam = Fraction(h * q, 7 * u * v)
    if lam > Fraction(6, 7):
        return 0
    lower = Fraction(G, 14) - rho
    upper = G * (Fraction(13, 14) - lam) - rho
    base = max(0, floor(upper) - ceil_fraction(lower) + 1)
    assert 0 <= base <= G
    return delta * base


def brute_events(h: int, u: int, v: int) -> tuple[Event, ...]:
    ans: list[Event] = []
    for n in range(u):
        a = tooth(u, n)
        for m in range(v + 1):
            b = tooth(v, m)
            if a.left < b.left < a.right < b.right:
                ans.append(event_for_addresses(h, u, n, v, m))
    return tuple(ans)


def audit_parametrization() -> int:
    checks = 0
    for h in range(2, 21):
        for u in range(1, 40, 2):
            for v in range(1, 40, 2):
                if u == v:
                    continue
                formula = addressed_events(h, u, v)
                brute = brute_events(h, u, v)
                assert set(formula) == set(brute), (h, u, v, formula, brute)
                by_d: dict[int, int] = {}
                for e in formula:
                    if e.component is not None:
                        by_d[e.d] = by_d.get(e.d, 0) + 1
                for d in determinant_indices(u, v):
                    assert by_d.get(d, 0) == orbit_count_for_d(h, u, v, d)

                reverse = addressed_events(h, v, u)
                reflected = {
                    (e.d, v - e.m, u - e.n, 2 * h - 1 - e.component)
                    for e in formula
                    if e.component is not None
                }
                reverse_live = {
                    (e.d, e.n, e.m, e.component)
                    for e in reverse
                    if e.component is not None
                }
                assert reflected == reverse_live, (
                    h,
                    u,
                    v,
                    sorted(reflected - reverse_live),
                    sorted(reverse_live - reflected),
                )
                checks += 1
    return checks


def meeting_teeth(w: int, left: Fraction, right: Fraction) -> list[Tooth]:
    lo = floor(w * left - Fraction(1, 14))
    hi = ceil(w * right + Fraction(1, 14))
    return [
        tooth(w, n)
        for n in range(lo - 1, hi + 2)
        if tooth(w, n).left < right and left < tooth(w, n).right
    ]


def greedy_chain(h: int, k: int, speeds: tuple[int, ...]) -> tuple[Tooth, ...] | None:
    left, right = anchor_component(h, k)
    teeth = [t for w in speeds for t in meeting_teeth(w, left, right)]
    chain: list[Tooth] = []
    cursor = left
    while True:
        eligible = [t for t in teeth if t.left < cursor < t.right]
        if not eligible:
            return None
        nxt = max(eligible, key=lambda t: (t.right, -t.left, -t.w, -t.n))
        if chain and nxt.right <= chain[-1].right:
            return None
        chain.append(nxt)
        cursor = nxt.right
        if cursor > right:
            return tuple(chain)


def actual_transitions(
    h: int, speeds: tuple[int, ...]
) -> tuple[tuple[int, Tooth, Tooth, Event], ...]:
    uses: list[tuple[int, Tooth, Tooth, Event]] = []
    for k in range(2 * h):
        chain = greedy_chain(h, k, speeds)
        if chain is None:
            continue
        for one, two in zip(chain, chain[1:]):
            event = event_for_addresses(h, one.w, one.n, two.w, two.n)
            assert event.component == k
            uses.append((k, one, two, event))
    # This is the central no-reuse assertion after addresses are retained.
    assert len({(a.w, a.n, b.w, b.n) for _, a, b, _ in uses}) == len(uses)
    return tuple(uses)


def pair_capacity(h: int, speeds: tuple[int, ...]) -> tuple[int, int, int]:
    filtered = 0
    unfiltered = 0
    distinct_component_slots: set[tuple[int, int, int]] = set()
    for i, u in enumerate(speeds):
        for v in speeds[i + 1 :]:
            for a, b in ((u, v), (v, u)):
                events = addressed_events(h, a, b)
                unfiltered += len(events)
                for e in events:
                    if e.component is not None:
                        filtered += 1
                        distinct_component_slots.add((min(u, v), max(u, v), e.component))
    return filtered, unfiltered, len(distinct_component_slots)


def print_h10_control() -> None:
    h, speeds = 10, (3, 13)
    uses = actual_transitions(h, speeds)
    print("h10_control")
    print(" uses", tuple((k, a.w, a.n, b.w, b.n, e.d, e.q, e.owner_bit) for k, a, b, e in uses))
    print(" capacity", pair_capacity(h, speeds))


def print_h420_control(extra: int) -> None:
    h = 420
    common = tuple(11 + 1680 * k for k in range(7)) + (525, 945, 1365, 1575)
    speeds = common + (extra,)
    uses = actual_transitions(h, speeds)
    transition_components = sorted({k for k, _, _, _ in uses})
    used_pairs: dict[tuple[int, int], list[tuple[int, int, int, int]]] = {}
    for k, a, b, e in uses:
        key = (min(a.w, b.w), max(a.w, b.w))
        used_pairs.setdefault(key, []).append((k, a.w, b.w, e.q))

    pair_summaries = []
    for (u, v), rows in sorted(used_pairs.items()):
        events = addressed_events(h, u, v) + addressed_events(h, v, u)
        live = tuple(e for e in events if e.component is not None)
        pair_summaries.append(
            (
                (u, v),
                gcd(u, v),
                len(determinant_indices(u, v)),
                len(live),
                len(events),
                len({e.component for e in live}),
                len(rows),
            )
        )

    print(f"h420_P{extra}")
    print(" transition_components", tuple(transition_components))
    print(" uses", tuple((k, a.w, a.n, b.w, b.n, e.d, e.q, e.owner_bit) for k, a, b, e in uses))
    print(" used_pair_multiplicities", tuple((p, tuple(rows)) for p, rows in sorted(used_pairs.items())))
    print(" used_pair_gcd_d_filtered_unfiltered_component_actual", tuple(pair_summaries))
    print(" capacity_filtered_unfiltered_pair_component", pair_capacity(h, speeds))


def print_gcd_crt_family() -> None:
    print("gcd_crt_family")
    for a in (1, 3, 9, 27):
        u, v = 13 * a, 3 * a
        rows: list[tuple[int, int, tuple[int, ...], tuple[int, ...]]] = []
        for H in (8, 10):
            h = H * a
            forward = addressed_events(h, u, v)
            reverse = addressed_events(h, v, u)
            kf = tuple(e.component for e in forward if e.component is not None)
            kr = tuple(e.component for e in reverse if e.component is not None)
            rows.append((H, len(kf) + len(kr), kf[:4], kr[:4]))
            if H == 8:
                assert not kf and not kr
            else:
                assert kf == tuple(20 * s + 6 for s in range(a))
                assert kr == tuple(20 * s + 13 for s in range(a))

        # Common dilation gives two actual pair-renewal components per copy.
        uses_bare = actual_transitions(10 * a, (3 * a, 13 * a))
        direct_bare = [x for x in uses_bare if {x[1].w, x[2].w} == {3 * a, 13 * a}]
        assert len(direct_bare) == 2 * a

        # Adding speed 1 makes the row primitive.  On every component wholly
        # inside the speed-1 safe band, the same pair transition is unchanged.
        uses_primitive = actual_transitions(10 * a, (1, 3 * a, 13 * a))
        direct_primitive = [
            x for x in uses_primitive if {x[1].w, x[2].w} == {3 * a, 13 * a}
        ]
        safe_band_direct = []
        for row in direct_primitive:
            left, right = anchor_component(10 * a, row[0])
            if Fraction(1, 14) <= left and right <= Fraction(13, 14):
                safe_band_direct.append(row)
        print(
            " a", a,
            "gcd", a,
            "lcm", 39 * a,
            "phase_rows", tuple(rows),
            "bare_actual", len(direct_bare),
            "primitive_actual", len(direct_primitive),
            "primitive_safe_band", len(safe_band_direct),
        )


def print_coprime_actual_family() -> None:
    print("coprime_actual_family")
    for N in (1, 2, 5, 10, 25):
        u, v = 14 * N + 1, 14 * N + 3
        h = 42 * N
        forward = addressed_events(h, u, v)
        reverse = addressed_events(h, v, u)
        live_forward = tuple(e for e in forward if e.component is not None)
        live_reverse = tuple(e for e in reverse if e.component is not None)
        uses = actual_transitions(h, (u, v))
        actual_d = {(a.w, e.d) for _, a, _, e in uses}
        forced = set(range(N + 1, 2 * N + 1))
        assert gcd(u, v) == 1
        assert determinant_indices(u, v) == tuple(range(1, 2 * N + 1))
        assert len(forward) == len(reverse) == 2 * N
        assert {(u, d) for d in forced} <= actual_d
        assert {(v, d) for d in forced} <= actual_d
        print(
            " N", N,
            "u_v_h", (u, v, h),
            "gcd_lcm", (1, u * v),
            "candidate_events", len(live_forward) + len(live_reverse),
            "actual_events", len(uses),
            "proved_subfamily_events", 2 * N,
            "actual_k_head", tuple(k for k, _, _, _ in uses[:8]),
        )


def main() -> None:
    print("parametrization_orbit_reflection_checks", audit_parametrization())
    print_h10_control()
    print_h420_control(1287)
    print_h420_control(9009)
    print_gcd_crt_family()
    print_coprime_actual_family()
    print("PASS")


if __name__ == "__main__":
    main()
