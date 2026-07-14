#!/usr/bin/env python3
"""Exact audit for THM-771: the seven-exception sheet endpoint defect.

The bad condition is STRICT, ||w(t0+k)/c|| < 1/14.  This convention is
essential: an endpoint-transfer moment is itself a closed LRC witness, not a
place where an exiting owner can be silently handed to an entering owner.

Tournament Analysis
-------------------
Vertices are bad-set owners (exceptions), not runners by default and not the
sheet vertices.  For owners a,b the directed observable is the number of
cyclic sheet adjacencies whose endpoints are privately owned a then b.  The
switch is the sign of the a->b versus b->a difference; ties use the owner-order
Hamiltonian path 1->...->7.  We report scores, directed 3-cycles, SCCs, edge
flips, and Hamiltonian-path counts.  This quotient preserves oriented private
splices but destroys capacity slack, ramification surplus, and simultaneous
overlap.  Accordingly it is diagnostic only; the theorem-facing carrier is
the full owner-by-sheet incidence deck and its exact defect F=Q+Omega-sigma.

All arithmetic and all assertions are exact.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
import random


DELTA = F(1, 14)


def circ(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def bad_set(w: int, c: int, t0: F) -> frozenset[int]:
    return frozenset(
        k for k in range(c) if circ(F(w) * (t0 + k) / c) < DELTA
    )


def owner_data(w: int, c: int) -> tuple[int, int, int, int]:
    """Return g, effective grid size C, reduced winding u, maximum capacity."""
    g = gcd(w, c)
    C = c // g
    u = w // g
    capacity = g * ((C + 6) // 7)
    return g, C, u, capacity


def is_unramified(w: int, c: int) -> bool:
    return owner_data(w, c)[1] % 7 == 0


def is_event(w: int, c: int, t0: F) -> bool:
    """Exact integral-grid endpoint event in the 7 | C layer."""
    _g, C, u, _cap = owner_data(w, c)
    assert C % 7 == 0
    return (F(u) * t0 - F(C, 14)) % 1 == 0


def event_set(w: int, c: int) -> tuple[F, ...]:
    """The single endpoint-transfer coset, of mesh 1/u (not 1/w)."""
    _g, C, u, _cap = owner_data(w, c)
    assert C % 7 == 0
    return tuple(sorted({(F(C, 14) + j) / u % 1 for j in range(u)}))


def cyclic_gaps(points: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(
        (points[(i + 1) % len(points)] - points[i]) % 1
        for i in range(len(points))
    )


def deck_stats(W: tuple[int, ...], c: int, t0: F) -> dict[str, object]:
    sets = tuple(bad_set(w, c, t0) for w in W)
    capacities = tuple(owner_data(w, c)[3] for w in W)
    loads = tuple(sum(k in B for B in sets) for k in range(c))
    Q = sum(A - len(B) for A, B in zip(capacities, sets))
    overlap = sum(max(d - 1, 0) for d in loads)
    sigma = sum(capacities) - c
    free = sum(d == 0 for d in loads)
    assert free == Q + overlap - sigma
    return {
        "sets": sets,
        "capacities": capacities,
        "loads": loads,
        "Q": Q,
        "overlap": overlap,
        "sigma": sigma,
        "free": free,
    }


def chamber_radius(W: tuple[int, ...], c: int, t0: F) -> F:
    """A positive radius on which every strict bad/nonbad incidence is fixed."""
    return min(
        abs(circ(F(w) * (t0 + k) / c) - DELTA) * c / w
        for w in W
        for k in range(c)
    )


def mirror_capacity(a: int, W: tuple[int, ...]) -> int:
    """The withdrawn THM-767 KCL diagnostic, retained only to falsify it locally."""
    ans = 0
    for b in W:
        if b == a:
            continue
        d = gcd(a, b)
        if (a + b) % (14 * d) == 0:
            ans += d
    return ans


def tournament_fingerprint(stats: dict[str, object]) -> dict[str, object]:
    sets = stats["sets"]
    loads = stats["loads"]
    c = len(loads)
    n = len(sets)
    private = []
    for k in range(c):
        owners = [a for a, B in enumerate(sets) if k in B]
        private.append(owners[0] if len(owners) == 1 else None)

    obs: dict[tuple[int, int], int] = defaultdict(int)
    for k in range(c):
        a, b = private[k], private[(k + 1) % c]
        if a is not None and b is not None and a != b:
            obs[(a, b)] += 1

    adj = [[False] * n for _ in range(n)]
    flips = 0
    for a, b in combinations(range(n), 2):
        if obs[(a, b)] >= obs[(b, a)]:
            adj[a][b] = True
        else:
            adj[b][a] = True
            flips += 1

    scores = tuple(sum(adj[a]) for a in range(n))
    cycles3 = 0
    for a, b, d in combinations(range(n), 3):
        cycles3 += int(
            (adj[a][b] and adj[b][d] and adj[d][a])
            or (adj[a][d] and adj[d][b] and adj[b][a])
        )

    reach = [[a == b or adj[a][b] for b in range(n)] for a in range(n)]
    for k in range(n):
        for a in range(n):
            for b in range(n):
                reach[a][b] = reach[a][b] or (reach[a][k] and reach[k][b])
    unseen = set(range(n))
    scc = []
    while unseen:
        a = min(unseen)
        comp = {b for b in unseen if reach[a][b] and reach[b][a]}
        unseen -= comp
        scc.append(len(comp))

    # Exact directed Hamiltonian-path count by endpoint/subset DP.
    dp: dict[tuple[int, int], int] = {}
    for a in range(n):
        dp[(1 << a, a)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + count
    hp = sum(dp.get(((1 << n) - 1, last), 0) for last in range(n))
    return {
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": cycles3,
        "scc_sizes": tuple(sorted(scc, reverse=True)),
        "edge_flips": flips,
        "hamiltonian_paths": hp,
        "private_word": tuple(0 if a is None else a + 1 for a in private),
    }


def exact_formula_audit() -> None:
    checked = 0
    for c in (7, 14, 21, 28, 35, 42, 49, 56, 63, 70):
        for w in range(1, 3 * c + 1):
            if w % c == 0 or not is_unramified(w, c):
                continue
            g, _C, _u, _cap = owner_data(w, c)
            samples = {F(a, q) for q in range(1, 11) for a in range(q)}
            samples.update(event_set(w, c))
            for t0 in samples:
                expected = c // 7 - (g if is_event(w, c, t0) else 0)
                assert len(bad_set(w, c, t0)) == expected
                checked += 1
    print(f"unramified_count_formula_checks={checked} failures=0")


def defect_identity_audit() -> None:
    rng = random.Random(771)
    checked = 0
    for _ in range(12000):
        c = rng.randint(2, 60)
        universe = [w for w in range(1, 10 * c) if w % c]
        W = tuple(rng.sample(universe, 7))
        q = rng.randint(1, 80)
        t0 = F(rng.randrange(q), q)
        deck_stats(W, c, t0)
        checked += 1
    print(f"capacity_defect_identity_checks={checked} failures=0")


def print_examples() -> None:
    P = tuple(range(1, 7))

    c7 = 7
    W7 = (1, 4, 5, 6, 8, 9, 10)
    x7 = F(1, 7)
    s7 = deck_stats(W7, c7, x7)
    assert all(is_unramified(w, c7) for w in W7)
    assert (s7["Q"], s7["overlap"], s7["sigma"], s7["free"]) == (0, 0, 0, 0)
    assert chamber_radius(W7, c7, x7) == F(1, 140)
    assert mirror_capacity(10, W7) == 0 < 10
    fp7 = tournament_fingerprint(s7)

    event_x = F(1, 8)  # an event of owner w=4
    assert is_event(4, c7, event_x)
    assert min(circ(F(p) * event_x) for p in P) >= DELTA
    e7 = deck_stats(W7, c7, event_x)
    free_sheets = [k for k, d in enumerate(e7["loads"]) if d == 0]
    assert free_sheets == [5]
    witness = (event_x + free_sheets[0]) / c7
    V7 = tuple(sorted(c7 * p for p in P) + sorted(W7))
    clearance = min(circ(F(v) * witness) for v in V7)
    assert witness == F(41, 56) and clearance == DELTA

    print("C7 STATIC ZERO-DEFECT CHAMBER")
    print(f"  P={P} W={W7} t0={x7}")
    print(f"  bad_sets={[sorted(B) for B in s7['sets']]}")
    print(f"  capacities={s7['capacities']} loads={s7['loads']}")
    print(
        f"  Q={s7['Q']} Omega={s7['overlap']} sigma={s7['sigma']} "
        f"F={s7['free']} chamber_radius={chamber_radius(W7,c7,x7)}"
    )
    print(f"  withdrawn_KCL_capacity(owner=10)={mirror_capacity(10,W7)} < current=10")
    print(f"  tournament={fp7}")
    print("C7 EVENT PIERCE")
    print(
        f"  event_owner=4 t0={event_x} free_sheets={free_sheets} "
        f"witness={witness} exact_clearance={clearance}"
    )

    c21 = 21
    W21 = (1, 2, 3, 4, 7, 49, 56)
    x21 = F(1, 7)
    s21 = deck_stats(W21, c21, x21)
    V21 = tuple(sorted(c21 * p for p in P) + sorted(W21))
    assert reduce(gcd, V21) == 1
    assert tuple(len(B) for B in s21["sets"]) == (3, 3, 3, 3, 7, 7, 7)
    assert (s21["Q"], s21["overlap"], s21["sigma"], s21["free"]) == (0, 12, 12, 0)
    assert chamber_radius(W21, c21, x21) == F(1, 112)
    fp21 = tournament_fingerprint(s21)
    print("C21 RAMIFIED SURPLUS/OVERLAP OBSTRUCTION")
    print(f"  P={P} W={W21} t0={x21} primitive_gcd={reduce(gcd,V21)}")
    print(f"  bad_sizes={[len(B) for B in s21['sets']]} capacities={s21['capacities']}")
    print(f"  loads={s21['loads']}")
    print(
        f"  Q={s21['Q']} Omega={s21['overlap']} sigma={s21['sigma']} "
        f"F={s21['free']} chamber_radius={chamber_radius(W21,c21,x21)}"
    )
    print(f"  tournament={fp21}")

    # Scaling changes raw w and g together but leaves the event mesh controlled
    # by u=w/g.  This is the exact falsifier to THM-767's raw-w density sentence.
    c = 105
    w = 150
    g, C, u, _cap = owner_data(w, c)
    points = event_set(w, c)
    gaps = cyclic_gaps(points)
    assert g == 15 and C == 7 and u == 10
    assert len(points) == u and set(gaps) == {F(1, u)}
    print("REDUCED-WINDING EVENT MESH")
    print(
        f"  c={c} w={w} g={g} C=c/g={C} u=w/g={u} "
        f"event_count={len(points)} actual_mesh={gaps[0]} raw_1/w={F(1,w)}"
    )


def main() -> None:
    exact_formula_audit()
    defect_identity_audit()
    print_examples()
    print("SUMMARY: ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
