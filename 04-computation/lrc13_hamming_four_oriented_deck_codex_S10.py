#!/usr/bin/env python3
"""Exact replay for THM-810.

The script uses integer arithmetic for every coverage decision.  Fractions are
used only to print sharp margins and the four small base-packet maxima.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, prod


P = 13
LABELS = tuple(range(1, P))
H = (1, 5, 8, 12)


def inv(x: int) -> int:
    return pow(x, -1, P)


def base_count(d: int, replacement: int, owner: int, side: str = "left") -> int:
    """Count the residue in the short d-window, 1 <= d <= 12."""
    target = d * replacement * inv(owner) % P
    if side == "left":
        window = range(-d + 1, d + 1)  # -d < z <= d
    elif side == "right":
        window = range(-d, d)  # -d <= z < d
    else:
        raise ValueError(side)
    return sum((z - target) % P == 0 for z in window)


def count(D: int, replacement: int, owner: int, side: str = "left") -> int:
    """The exact number N_D(replacement, owner)."""
    q, d = divmod(D, P)
    assert 1 <= d <= 12
    return 2 * q + base_count(d, replacement, owner, side)


def direct_counts(D: int, side: str = "left") -> tuple[int, ...]:
    """Independently count every residue in the full oriented D-window."""
    if side == "left":
        window = range(-D + 1, D + 1)  # -D < z <= D
    elif side == "right":
        window = range(-D, D)  # -D <= z < D
    else:
        raise ValueError(side)
    out = [0] * P
    for z in window:
        out[z % P] += 1
    return tuple(out)


def capacity(D: int, replacement: int, owner: int, side: str = "left") -> Fraction:
    return Fraction(count(D, replacement, owner, side), D)


def f(D: int) -> Fraction:
    return max(capacity(D, 1, owner) for owner in LABELS)


def s4(D: int) -> Fraction:
    return sum(sorted((capacity(D, 1, owner) for owner in LABELS), reverse=True)[:4])


def owner_numerators(labels: tuple[int, ...], orders: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    denominator = prod(orders)
    nums = tuple(
        sum(count(D, r, owner) * (denominator // D) for r, D in zip(labels, orders))
        for owner in labels
    )
    return denominator, nums


def owner_sums(labels: tuple[int, ...], orders: tuple[int, ...]) -> tuple[Fraction, ...]:
    den, nums = owner_numerators(labels, orders)
    return tuple(Fraction(n, den) for n in nums)


def update_best(
    best: tuple[int, int, object] | None,
    labels: tuple[int, ...],
    orders: tuple[int, ...],
) -> tuple[int, int, object]:
    den, nums = owner_numerators(labels, orders)
    num = min(nums)
    data = (labels, orders, tuple(Fraction(n, den) for n in nums))
    if best is None or num * best[1] > best[0] * den:
        return num, den, data
    return best


def normalized_one_d2_three_d3() -> tuple[int, Fraction, object]:
    rows = 0
    survivors = []
    best = None
    for tail in combinations(range(2, P), 3):
        labels = (1,) + tail
        orders = (2, 3, 3, 3)
        rows += 1
        den, nums = owner_numerators(labels, orders)
        if min(nums) >= den:
            survivors.append((labels, orders))
        best = update_best(best, labels, orders)
    assert rows == 165 and not survivors and best is not None
    return rows, Fraction(best[0], best[1]), best[2]


def normalized_two_d2(L16: tuple[int, ...]) -> tuple[int, int, Fraction, object]:
    configurations = 0
    rows = 0
    survivors = []
    best = None
    for second in range(2, P):
        remaining = [x for x in range(2, P) if x != second]
        for u, v in combinations(remaining, 2):
            labels = (1, second, u, v)
            # With no D=2 hit, the two remaining colours contribute at most 2/3.
            if any(count(2, 1, owner) + count(2, second, owner) == 0 for owner in labels):
                continue
            configurations += 1
            for D, E in product(L16, repeat=2):
                orders = (2, 2, D, E)
                rows += 1
                den, nums = owner_numerators(labels, orders)
                if min(nums) >= den:
                    survivors.append((labels, orders))
                best = update_best(best, labels, orders)
    assert configurations == 49 and rows == 78_400
    assert not survivors and best is not None
    return configurations, rows, Fraction(best[0], best[1]), best[2]


def normalized_one_d3_three_large(L23: tuple[int, ...]) -> tuple[int, Fraction, object]:
    rows = 0
    survivors = []
    best = None
    # A lone D=3 colour must hit every owner, so normalize it to 1.
    for tail in combinations((3, 5, 8, 10), 3):
        labels = (1,) + tail
        for Ds in product(L23, repeat=3):
            orders = (3,) + Ds
            rows += 1
            den, nums = owner_numerators(labels, orders)
            if min(nums) >= den:
                survivors.append((labels, orders))
            best = update_best(best, labels, orders)
    assert rows == 62_500 and not survivors and best is not None
    return rows, Fraction(best[0], best[1]), best[2]


def normalized_all_d3() -> tuple[int, list[tuple[tuple[int, ...], tuple[Fraction, ...]]]]:
    rows = 0
    survivors = []
    for tail in combinations(range(2, P), 3):
        labels = (1,) + tail
        orders = (3, 3, 3, 3)
        rows += 1
        sums = owner_sums(labels, orders)
        if min(sums) >= 1:
            survivors.append((labels, sums))
    assert rows == 165
    assert survivors == [(H, (Fraction(1),) * 4)]
    return rows, survivors


def no_d4_clique() -> tuple[tuple[int, ...], tuple[tuple[int, tuple[int, ...]], ...]]:
    mutual = tuple(
        y
        for y in range(2, P)
        if count(4, 1, y) == 1 and count(4, y, 1) == 1
    )
    adjacency = tuple(
        (
            x,
            tuple(
                y
                for y in mutual
                if y != x and count(4, x, y) == 1 and count(4, y, x) == 1
            ),
        )
        for x in mutual
    )
    assert mutual == (3, 4, 9, 10)
    assert not any(
        all(count(4, x, y) == count(4, y, x) == 1 for x, y in combinations(T, 2))
        for T in combinations(mutual, 3)
    )
    return mutual, adjacency


def crt_u(label: int, parity: int) -> int:
    return next(u for u in range(1, 3 * P) if u % P == 3 * label % P and u % 3 == parity)


def signed_mod(x: int, modulus: int) -> int:
    r = x % modulus
    return r if r <= modulus // 2 else r - modulus


def eligible_sheets(label: int, parity: int, owner: int) -> tuple[int, ...]:
    u = crt_u(label, parity)
    a = inv(owner)
    return tuple(
        ell
        for ell in range(3)
        if -3 < signed_mod(u * (a + P * ell), 3 * P) <= 3
    )


def exact_overlap_patterns() -> tuple[list[tuple[int, ...]], list[tuple[tuple[int, ...], object]]]:
    feasible = []
    audit = []
    for signs in product((1, 2), repeat=4):
        tables = []
        ok = True
        for owner in H:
            row = tuple(eligible_sheets(r, e, owner) for r, e in zip(H, signs))
            occupied = [ell for cells in row for ell in cells]
            if sorted(occupied) != [0, 1, 2]:
                ok = False
            tables.append((owner, row))
        audit.append((signs, tuple(tables)))
        if ok:
            feasible.append(signs)
    expected = [(1, 1, 1, 1), (1, 2, 2, 1), (2, 1, 1, 2), (2, 2, 2, 2)]
    assert feasible == expected
    assert all(s[0] == s[3] and s[1] == s[2] for s in feasible)
    return feasible, audit


def coset_clock() -> tuple[tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, int], ...]], ...]:
    """Audit the lift-invariant equality clock on the order-three interface."""
    expected = {
        (1, 5, 8, 12): ((2, 10, 11, 16, 23, 28, 29, 37), ((12, 27), (18, 21))),
        (2, 3, 10, 11): ((1, 5, 8, 14, 25, 31, 34, 38), ((3, 36), (15, 24))),
        (4, 6, 7, 9): ((4, 7, 17, 19, 20, 22, 32, 35), ((6, 33), (9, 30))),
    }
    rows = []
    for multiplier in (1, 2, 4):
        labels = tuple(sorted({multiplier * r % P for r in H}))
        core = tuple(3 * r for r in LABELS if r not in labels)
        negative_pairs = tuple((r, (-r) % P) for r in labels if r < (-r) % P)
        pattern_clocks = []
        active_pairs = set()
        for bits in product((1, 2), repeat=2):
            parity = {
                r: bit
                for pair, bit in zip(negative_pairs, bits)
                for r in pair
            }
            lifts = tuple(crt_u(r, parity[r]) for r in labels)
            speeds = core + lifts
            safe = []
            for a in range(1, 3 * P):
                if gcd(a, 3 * P) != 1:
                    continue
                residues = tuple(signed_mod(a * speed, 3 * P) for speed in speeds)
                margin = min(abs(z) for z in residues)
                if margin >= 3:
                    assert margin == 3
                    safe.append(a)
                    active = tuple(
                        speed for speed, z in zip(speeds, residues) if abs(z) == 3
                    )
                    assert len(active) == 2
                    active_residues = {
                        z for speed, z in zip(speeds, residues) if speed in active
                    }
                    assert active_residues == {-3, 3}
                    assert all(speed in core for speed in active)
                    active_pairs.add(tuple(sorted(active)))
            pattern_clocks.append(tuple(safe))
        clock, expected_pairs = expected[labels]
        assert all(value == clock for value in pattern_clocks)
        assert tuple(sorted(active_pairs)) == expected_pairs
        rows.append((labels, clock, expected_pairs))
    return tuple(rows)


def norm(x: Fraction) -> Fraction:
    r = x % 1
    return min(r, 1 - r)


def exact_M(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    candidates = {Fraction(0), Fraction(1)}
    for a in speeds:
        candidates.update(Fraction(k, 2 * a) for k in range(2 * a + 1))
    for a, b in combinations(speeds, 2):
        for denominator in (a + b, abs(a - b)):
            if denominator:
                candidates.update(Fraction(k, denominator) for k in range(denominator + 1))
    scored = [(min(norm(a * t) for a in speeds), t) for t in sorted(candidates)]
    maximum = max(value for value, _ in scored)
    maximizers = tuple(t for value, t in scored if value == maximum)
    return maximum, maximizers


def base_packets(feasible: list[tuple[int, ...]]) -> list[tuple[object, ...]]:
    rows = []
    core = {3 * r for r in range(1, P) if r not in H}
    expected = [Fraction(6, 43), Fraction(2, 17), Fraction(6, 43), Fraction(1, 8)]
    for signs, target in zip(feasible, expected):
        lifts = tuple(crt_u(r, e) for r, e in zip(H, signs))
        speeds = tuple(sorted(core | set(lifts)))
        assert len(speeds) == 12
        common = 0
        for speed in speeds:
            common = gcd(common, speed)
        assert common == 1
        maximum, maximizers = exact_M(speeds)
        assert maximum == target
        rows.append((signs, lifts, speeds, maximum, maximizers[0]))
    return rows


def tournament_fingerprint(side: str) -> tuple[object, ...]:
    # Pair observable: antisymmetric cross-capacity difference.  Every pair ties
    # on H, so the declared numerical tie path supplies the orientation.
    n = len(H)
    edge = [[False] * n for _ in range(n)]
    for i, r in enumerate(H):
        for j in range(i + 1, n):
            s = H[j]
            delta = count(3, r, s, side) - count(3, s, r, side)
            if delta >= 0:  # delta=0 uses H's numerical tie path
                edge[i][j] = True
            else:
                edge[j][i] = True
    scores = tuple(sum(row) for row in edge)
    triangles = sum(
        edge[i][j] and edge[j][k] and edge[k][i]
        or edge[j][i] and edge[k][j] and edge[i][k]
        for i, j, k in combinations(range(n), 3)
    )
    hpaths = sum(all(edge[p[i]][p[i + 1]] for i in range(n - 1)) for p in permutations(range(n)))
    reach = [[i == j or edge[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or reach[i][k] and reach[k][j]
    scc_sizes = []
    unseen = set(range(n))
    while unseen:
        i = min(unseen)
        component = {j for j in unseen if reach[i][j] and reach[j][i]}
        scc_sizes.append(len(component))
        unseen -= component
    return scores, triangles, tuple(sorted(scc_sizes, reverse=True)), hpaths, edge


def incidence_flip_count() -> int:
    return sum(
        (count(3, r, o, "left") > 0) != (count(3, r, o, "right") > 0)
        for r in H
        for o in H
    )


def main() -> None:
    # Infinite-order attenuation identity and its two finite cutoff lists.
    for D in range(1, 1000):
        if D % P == 0:
            continue
        direct = direct_counts(D)
        for r in LABELS:
            for o in LABELS:
                target = D * r * inv(o) % P
                assert count(D, r, o) == direct[target]

    finite_orders = tuple(D for D in range(3, 79) if D % P)
    assert all(f(D) <= Fraction(1, 3) for D in finite_orders)
    assert tuple(D for D in finite_orders if f(D) == Fraction(1, 3)) == (3,)
    assert all(s4(D) <= 1 for D in finite_orders if D >= 4)
    assert tuple(D for D in finite_orders if D >= 4 and s4(D) == 1) == (4,)

    L16 = tuple(D for D in finite_orders if f(D) >= Fraction(1, 6))
    L23 = tuple(D for D in finite_orders if D >= 4 and s4(D) >= Fraction(2, 3))
    assert max(L16) == 72 and max(L23) == 48
    assert all(f(D) < Fraction(1, 6) for D in range(79, 501) if D % P)
    assert all(s4(D) < Fraction(2, 3) for D in range(79, 501) if D % P)

    one2 = normalized_one_d2_three_d3()
    two2 = normalized_two_d2(L16)
    one3 = normalized_one_d3_three_large(L23)
    all3 = normalized_all_d3()
    mutual4, adjacency4 = no_d4_clique()
    feasible, _ = exact_overlap_patterns()
    clock_rows = coset_clock()
    packets = base_packets(feasible)

    left = tournament_fingerprint("left")
    right = tournament_fingerprint("right")
    edge_flips = sum(
        left[4][i][j] != right[4][i][j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    assert left[:4] == right[:4] == ((3, 2, 1, 0), 0, (1, 1, 1, 1), 1)
    assert edge_flips == 0 and incidence_flip_count() == 8

    lines = [
        "THM-810 exact replay",
        "attenuation_identity.direct_D<=999=PASS",
        f"cutoff.f>=1/6.count={len(L16)} max={max(L16)} values={L16}",
        f"cutoff.S4>=2/3.count={len(L23)} max={max(L23)} values={L23}",
        "sharp_bounds.f.D>=3<=1/3 equality_orders=(3,) S4.D>=4<=1 equality_orders=(4,)",
        f"case.one_D2_three_D3.rows={one2[0]} survivors=0 best_min={one2[1]}",
        f"case.two_D2.configurations={two2[0]} rows={two2[1]} survivors=0 best_min={two2[2]}",
        f"case.one_D3_three_large.rows={one3[0]} survivors=0 best_min={one3[1]}",
        f"case.all_D3.rows={all3[0]} survivors={tuple(x[0] for x in all3[1])}",
        f"case.no_D2_D3.mutual_D4={mutual4} adjacency={adjacency4} clique4=0",
        f"overlap.feasible_patterns={tuple(feasible)} count={len(feasible)}",
        f"coset_clock.mod39={clock_rows} margin=1/13",
    ]
    for signs, lifts, speeds, maximum, witness in packets:
        lines.append(
            f"base.signs={signs} lifts={lifts} speeds={speeds} M={maximum} witness={witness}"
        )
    lines.extend(
        [
            f"tournament.left.scores={left[0]} triangles={left[1]} scc={left[2]} H={left[3]}",
            f"tournament.right.scores={right[0]} triangles={right[1]} scc={right[2]} H={right[3]}",
            f"tournament.gauge_edge_flips={edge_flips} incidence_flips={incidence_flip_count()}",
        ]
    )
    digest = sha256(("\n".join(lines) + "\n").encode()).hexdigest()
    lines.append(f"sha256={digest}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
