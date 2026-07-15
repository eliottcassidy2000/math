#!/usr/bin/env python3
"""Exact replay for the Hamming-one full-residue rigidity lemma.

The theorem proved in the companion write-up is the following uniform statement.
For ``1 <= r <= 12`` and ``k >= 1``, put

    W_k(r) = ({1,...,12} - {r}) union {r + 13 k}.

Then ``M(W_k(r)) > 1/13``.  More generally, a shallow full-residue
Hamming-one perturbation of ``c*{1,...,12}`` can be tight only when it restores
the deleted speed.

This file replays every finite atom in that proof using ``Fraction`` arithmetic:

* the twelve core-safe witnesses and their exact interval radii;
* the resulting exact tooth-width thresholds;
* the fact that precisely seven low replacement atoms are not discharged by
  the width inequality;
* explicit strict witnesses, complete least-residue vectors, and binders for
  those seven atoms.

The two infinite branches are algebraic, not bounded searches.  They are
printed with their proof identities:

* for ``r >= 7``, either ``r`` does not divide ``k`` and ``1/r`` is already a
  strict witness, or ``k >= r`` and the replacement exceeds the width bound;
* in the dilated problem, the ``c`` preimages of a missing protected splice
  form a translated grid of order ``D=c/gcd(c,w)``.  Tightness would place that
  whole grid in a phase arc of length ``2/13``.  For ``D>=2`` its circular span
  is ``1-1/D >= 1/2``, so ``D=1`` and ``c|w``.

The bounded deck census at the end is telemetry only.  It independently checks
the exact grid-capacity inequality through ``c<=64`` but is not used as proof
of the uniform deck statement.

Tournament Analysis / assumption challenge
------------------------------------------
The theorem-facing vertices are the lifted missing-splice sheets, followed
after descent by a core-safe interval and the individual replacement teeth.
A runner tournament forgets both sheet multiplicity and metric width.  For
telemetry, sheets are oriented by the exact danger margin
``1/13-||w*t_l||``; the switch reverses that comparison and sheet order is the
tie Hamiltonian gauge.  Both gauges are transitive.  Their bare fingerprints
therefore cannot record the deciding predicate that *all* sheet margins are
nonnegative.  The exact carrier is the sheet-danger incidence deck plus the
metric interval/tooth sidecar.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import gcd
from typing import Sequence


DELTA = F(1, 13)
BASE = tuple(range(1, 13))


def fmt(x: F) -> str:
    """Stable compact formatting for an exact rational."""

    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def least_absolute_residue(n: int, q: int) -> int:
    z = n % q
    return min(z, q - z)


def norm_fraction(x: F) -> F:
    z = x.numerator % x.denominator
    return F(min(z, x.denominator - z), x.denominator)


def clearance_vector(speeds: Sequence[int], a: int, q: int) -> tuple[int, ...]:
    return tuple(least_absolute_residue(a * v, q) for v in speeds)


def clearance(speeds: Sequence[int], a: int, q: int) -> F:
    return F(min(clearance_vector(speeds, a, q)), q)


def binders(speeds: Sequence[int], vector: Sequence[int]) -> tuple[int, ...]:
    minimum = min(vector)
    return tuple(v for v, z in zip(speeds, vector) if z == minimum)


CORE_EXPECTED: dict[int, tuple[int, int, F, tuple[int, ...]]] = {
    1: (1, 14, F(1, 7), (2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2)),
    2: (7, 15, F(2, 15), (7, 6, 2, 5, 3, 4, 4, 3, 5, 2, 6)),
    3: (5, 16, F(1, 8), (5, 6, 4, 7, 2, 3, 8, 3, 2, 7, 4)),
    4: (4, 17, F(2, 17), (4, 8, 5, 3, 7, 6, 2, 2, 6, 7, 3)),
    5: (7, 18, F(1, 9), (7, 4, 3, 8, 6, 5, 2, 9, 2, 5, 6)),
    6: (3, 19, F(2, 19), (3, 6, 9, 7, 4, 2, 5, 8, 8, 5, 2)),
}


THRESHOLD_EXPECTED: dict[int, F] = {
    1: F(14),
    2: F(180, 11),
    3: F(96, 5),
    4: F(68, 3),
    5: F(27),
    6: F(228, 7),
    7: F(14),
    8: F(96, 5),
    9: F(27),
    10: F(40),
    11: F(66),
    12: F(132),
}


LOW_EXPECTED = {
    (1, 1),
    (2, 1),
    (3, 1),
    (4, 1),
    (5, 1),
    (6, 1),
    (6, 2),
}


ATOM_EXPECTED: dict[
    tuple[int, int], tuple[int, int, F, tuple[int, ...], tuple[int, ...]]
] = {
    (1, 1): (
        1,
        16,
        F(1, 8),
        (2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 2),
        (2, 14),
    ),
    (2, 1): (
        9,
        19,
        F(2, 19),
        (9, 8, 2, 7, 3, 6, 4, 5, 5, 4, 6, 2),
        (4, 15),
    ),
    (3, 1): (
        6,
        17,
        F(2, 17),
        (6, 5, 7, 4, 2, 8, 3, 3, 8, 2, 4, 6),
        (6, 11),
    ),
    (4, 1): (
        5,
        19,
        F(2, 19),
        (5, 9, 4, 6, 8, 3, 2, 7, 7, 2, 3, 9),
        (8, 11),
    ),
    (5, 1): (
        4,
        19,
        F(2, 19),
        (4, 8, 7, 3, 5, 9, 6, 2, 2, 6, 9, 4),
        (9, 10),
    ),
    (6, 1): (
        4,
        23,
        F(2, 23),
        (4, 8, 11, 7, 3, 5, 9, 10, 6, 2, 2, 7),
        (11, 12),
    ),
    (6, 2): (
        7,
        44,
        F(1, 11),
        (7, 14, 21, 16, 9, 5, 12, 19, 18, 11, 4, 4),
        (12, 32),
    ),
}


def core_safe_interval_audit() -> dict[int, F]:
    print("CORE_SAFE_INTERVAL_AUDIT")
    thresholds: dict[int, F] = {}
    for r in BASE:
        core = tuple(v for v in BASE if v != r)
        if r <= 6:
            a, q, mu, expected_vector = CORE_EXPECTED[r]
        else:
            a, q, mu = 1, r, F(1, r)
            expected_vector = clearance_vector(core, a, q)

        vector = clearance_vector(core, a, q)
        assert vector == expected_vector
        assert F(min(vector), q) == mu
        assert mu > DELTA

        maximum = max(core)
        rho = (mu - DELTA) / maximum
        threshold = DELTA / rho
        assert threshold == THRESHOLD_EXPECTED[r]
        thresholds[r] = threshold

        print(
            f"r={r:2d} B={maximum:2d} t={a}/{q} mu={fmt(mu):>6s} "
            f"rho={fmt(rho):>8s} T={fmt(threshold):>6s} "
            f"binders={list(binders(core, vector))} vector={list(vector)}"
        )
    print()
    return thresholds


def low_atom_audit(thresholds: dict[int, F]) -> None:
    low: set[tuple[int, int]] = set()
    for r in range(1, 7):
        k = 1
        while r + 13 * k <= thresholds[r]:
            low.add((r, k))
            k += 1
    assert low == LOW_EXPECTED

    print("WIDTH_RESIDUAL")
    print(f"low_atoms={sorted(low)} count={len(low)}")
    print("all other r<=6 rows satisfy w>T and escape by interval/tooth width")
    print()

    print("LOW_ATOM_EXACT_WITNESSES")
    for r, k in sorted(low):
        a, q, expected_mu, expected_vector, expected_binders = ATOM_EXPECTED[(r, k)]
        replacement = r + 13 * k
        speeds = tuple(sorted(v for v in BASE if v != r) + [replacement])
        vector = clearance_vector(speeds, a, q)
        mu = F(min(vector), q)
        actual_binders = binders(speeds, vector)
        assert vector == expected_vector
        assert mu == expected_mu
        assert actual_binders == expected_binders
        assert mu > DELTA
        print(
            f"r={r:2d} k={k:2d} w={replacement:3d} t={a}/{q} "
            f"clearance={fmt(mu):>5s} binders={list(actual_binders)} "
            f"vector={list(vector)}"
        )
    print()


def symbolic_branch_audit(thresholds: dict[int, F]) -> None:
    print("UNIFORM_SYMBOLIC_BRANCHES")
    print(
        "SAFE_INTERVAL: rho=(mu-1/13)/B and w>T=1/(13*rho) imply "
        "2*rho>2/(13*w), so a connected core-safe interval cannot lie in one w-tooth"
    )

    for r in range(7, 13):
        assert gcd(r, 13) == 1
        assert F(1, r) > DELTA
        assert 14 * r > thresholds[r]
        print(
            f"R_GE_7 r={r:2d}: r does not divide k -> t=1/{r} has clearance >=1/{r}; "
            f"r divides k -> w>=14r={14*r}>T={fmt(thresholds[r])}"
        )

    assert F(1, 2) > F(2, 13)
    print(
        "DILATED_DECK: D=c/gcd(c,w)>=2 -> circular_span=1-1/D>=1/2>2/13; "
        "therefore full closed-danger incidence forces D=1 and c|w"
    )
    print(
        "GERM_CLOSURE: the surviving complementary owner is -1/13 at every lifted splice; "
        "a left half-germ is core-safe, and closedness of {||wt||<=1/13} carries danger to the splice"
    )
    print(
        "DESCENT: 13 does not divide c and w congruent cr (mod 13), together with c|w, "
        "give w/c=r+13k; k=0 restores the AP and k>=1 is the loose W_k(r) branch"
    )
    print()


def deck_phases(c: int, r: int, h: int) -> tuple[int, int, tuple[F, ...]]:
    """Return ``(w,a,phase norms)`` for the c lifted splice preimages.

    Here ``w=c*r+13*h`` is the residue-preserving replacement and
    ``a*r=1 mod 13``.  Only ``h mod c`` matters for the deck.
    """

    assert c % 13 != 0
    a = pow(r, -1, 13)
    w = c * r + 13 * h
    phases = tuple(
        norm_fraction(F(w * (a + 13 * ell), 13 * c)) for ell in range(c)
    )
    return w, a, phases


def bounded_deck_telemetry(limit: int = 64) -> tuple[int, int, int, int, tuple[int, int, int, int]]:
    """Finite telemetry only; the printed symbolic span proves the general step."""

    checked = 0
    false_full_covers = 0
    checksum = 0
    order_histogram: Counter[int] = Counter()
    danger_histogram: Counter[int] = Counter()
    best: tuple[F, int, int, int, int] | None = None

    for c in range(2, limit + 1):
        if c % 13 == 0:
            continue
        for r in BASE:
            for h in range(1, c):
                w, a, phases = deck_phases(c, r, h)
                assert w % 13 == (c * r) % 13
                assert w % c != 0
                order = c // gcd(c, w)
                repetition = c // order
                span = F(order - 1, order)
                assert order >= 2
                assert span >= F(1, 2) > F(2, 13)

                dangerous = sum(z <= DELTA for z in phases)
                capacity = repetition * ((2 * order) // 13 + 1)
                assert dangerous <= capacity
                if dangerous == c:
                    false_full_covers += 1

                checked += 1
                order_histogram[order] += 1
                danger_histogram[dangerous] += 1
                checksum += (
                    1000003 * c
                    + 1000033 * r
                    + 1000037 * h
                    + 1000039 * order
                    + 1000081 * dangerous
                    + 1000099 * a
                )

                ratio = F(dangerous, c)
                candidate = (ratio, -c, -r, -h, dangerous)
                if best is None or candidate > best:
                    best = candidate

    assert false_full_covers == 0
    assert best is not None
    _, neg_c, neg_r, neg_h, best_dangerous = best
    best_row = (-neg_c, -neg_r, -neg_h, best_dangerous)

    print("BOUNDED_DECK_TELEMETRY_NOT_A_PROOF")
    print(f"c_range=2..{limit} excluding multiples of 13")
    print(f"nondivisible_rows={checked} false_full_covers={false_full_covers}")
    print(f"deck_order_histogram={dict(sorted(order_histogram.items()))}")
    print(f"dangerous_sheet_count_histogram={dict(sorted(danger_histogram.items()))}")
    c, r, h, dangerous = best_row
    print(
        f"max_danger_fraction={dangerous}/{c} first_row=(c={c},r={r},h={h},w={c*r+13*h})"
    )
    print(f"telemetry_checksum={checksum}")
    print()
    return checked, false_full_covers, checksum, limit, best_row


def orient_by_margin(margins: Sequence[F], larger_wins: bool) -> list[list[bool]]:
    n = len(margins)
    adjacency = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if margins[i] == margins[j]:
                winner, loser = i, j  # sheet order is the tie gauge
            elif (margins[i] > margins[j]) == larger_wins:
                winner, loser = i, j
            else:
                winner, loser = j, i
            adjacency[winner][loser] = True
    return adjacency


def tournament_fingerprint(adjacency: Sequence[Sequence[bool]]) -> tuple[dict[int, int], int, list[int], int]:
    n = len(adjacency)
    scores = [sum(row) for row in adjacency]
    score_histogram = dict(sorted(Counter(scores).items()))

    triangles = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if (
                    adjacency[i][j]
                    and adjacency[j][k]
                    and adjacency[k][i]
                ) or (
                    adjacency[i][k]
                    and adjacency[k][j]
                    and adjacency[j][i]
                ):
                    triangles += 1

    reach = [[i == j or adjacency[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc_sizes: list[int] = []
    while unseen:
        root = min(unseen)
        component = {v for v in unseen if reach[root][v] and reach[v][root]}
        scc_sizes.append(len(component))
        unseen -= component
    scc_sizes.sort(reverse=True)

    # Exact directed Hamiltonian-path count by subset DP; n=12 in the report.
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if count == 0:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adjacency[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    hamiltonian_paths = sum(dp[-1])
    return score_histogram, triangles, scc_sizes, hamiltonian_paths


def tournament_telemetry() -> None:
    c, r, h = 12, 6, 5
    w, a, phases = deck_phases(c, r, h)
    margins = tuple(DELTA - z for z in phases)
    danger_first = orient_by_margin(margins, larger_wins=True)
    escape_first = orient_by_margin(margins, larger_wins=False)
    fingerprint = tournament_fingerprint(danger_first)
    reverse_fingerprint = tournament_fingerprint(escape_first)
    edge_flips = sum(
        danger_first[i][j] != escape_first[i][j]
        for i in range(c)
        for j in range(i + 1, c)
    )
    path = tuple(sorted(range(c), key=lambda ell: (-margins[ell], ell)))
    assert all(danger_first[path[i]][path[i + 1]] for i in range(c - 1))
    assert fingerprint[1] == 0
    assert fingerprint[2] == [1] * c
    assert fingerprint[3] == 1
    assert reverse_fingerprint[1] == 0
    assert reverse_fingerprint[2] == [1] * c
    assert reverse_fingerprint[3] == 1

    print("SHEET_MARGIN_TOURNAMENT_TELEMETRY")
    print(
        f"representative=(c={c},r={r},h={h},w={w},inverse_a={a}) "
        f"deck_order={c//gcd(c,w)}"
    )
    print(f"phase_norms={[fmt(z) for z in phases]}")
    print(f"danger_margins={[fmt(z) for z in margins]}")
    print(f"closed_danger_mask={[int(z >= 0) for z in margins]}")
    print(f"danger_first_score_histogram={fingerprint[0]}")
    print(
        f"danger_first_directed_triangles={fingerprint[1]} "
        f"scc_sizes={fingerprint[2]} hamiltonian_paths={fingerprint[3]}"
    )
    print(f"tie_gauge_hamiltonian_path={list(path)}")
    print(f"edge_flips_under_margin_switch={edge_flips}")
    print(
        "interpretation=both scalar gauges are transitive telemetry; the exact all-sheets-dangerous "
        "bit and the interval/tooth metric are external sidecars"
    )
    print()


def main() -> None:
    print("LRC13 HAMMING-ONE FULL-RESIDUE RIGIDITY — EXACT REPLAY")
    print("arithmetic=fractions.Fraction floating_point=none global_M_enumeration=none")
    print()

    thresholds = core_safe_interval_audit()
    low_atom_audit(thresholds)
    symbolic_branch_audit(thresholds)
    bounded_deck_telemetry()
    tournament_telemetry()

    print("PASS: all finite atoms and exact telemetry reproduced; uniform branches are algebraic as labelled.")


if __name__ == "__main__":
    main()
