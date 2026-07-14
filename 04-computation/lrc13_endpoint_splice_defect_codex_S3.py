#!/usr/bin/env python3
"""Exact endpoint/splice audit for the residue-pinned LRC(13) rigidity leaf.

For a target 1/q and speed w, split the open danger set into teeth

    I(w,m) = {t in R/Z : |t-m/w|_circle < 1/(q w)},  m in Z/wZ.

The exact primary object is the intersection graph of the individual teeth.
Assume their open union is not the whole circle, and let kappa be its number
of connected components.  A protected splice is a point where a tooth ending
from the left meets a tooth starting to the right and no open tooth contains
the point.  Let P be the number of such splices.  Then kappa-P is exactly the
number of open target-safe components.  The non-full-union hypothesis is
automatic in the twelve-speed q=13 bank by the settled LRC lower bound.  The
sweep implementation counts component starts; this agrees with literal
kappa under that hypothesis and is zero in the excluded full-circle case.

For a complete nonzero residue system modulo 13, nominal splices can only
come from the six complementary residue pairs.  Their multiplicity capacity
is

    Pstar = 2 * sum_{r=1}^6 gcd(w_r, w_{13-r}).

Thus the endpoint Euler defect decomposes as

    kappa-P = (kappa-Pstar) + (Pstar-P),

separating overlap-rank shortage from third-runner blocker debt.

All arithmetic below is exact (fractions and integers); no floating point is
used.  The exhaustive bank is the 4095 nonzero height-one lift patterns
w_r = r + 13*k_r, k_r in {0,1}.

Tournament Analysis guardrail
-----------------------------
Vertices in the diagnostic tournament are the six complementary-pair splice
obligations, not runners.  The pairwise observable is directed blocker count:
A -> B when runners owned by pair A block more raw splice events owned by B
than conversely.  The switch/gauge is the sign of that difference.  Ties use
the Hamiltonian path 1 -> 2 -> ... -> 6.  Scores, 3-cycles, SCC sizes, and edge
flips are reported.  This tournament preserves blocker direction but destroys
tooth-component rank and simultaneous coverage.  In particular, positive-
defect examples below can have the same tournament as the AP.  The theorem-
facing carrier is therefore the endpoint graph plus protected-splice sidecar.

This is a reproducible audit, not a proof of uniform strict lift-rigidity.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations, product
from math import gcd


Q13 = 13


def circle_distance(x: Fraction) -> Fraction:
    """Distance of x modulo one to the nearest integer."""
    phase = x % 1
    return min(phase, 1 - phase)


def endpoint_events(speeds: tuple[int, ...], q: int):
    """Return exact tooth endpoint events, with start/end owner lists."""
    events: dict[Fraction, dict[str, list[int]]] = {}
    for w in speeds:
        for m in range(w):
            start = Fraction(q * m - 1, q * w) % 1
            end = Fraction(q * m + 1, q * w) % 1
            events.setdefault(start, {"start": [], "end": []})["start"].append(w)
            events.setdefault(end, {"start": [], "end": []})["end"].append(w)
    return events


def endpoint_stats(speeds: tuple[int, ...], q: int):
    """Exact component-start, protected-splice, and gap counts by endpoint sweep.

    If the open danger union has nonempty complement, ``kappa`` is its literal
    component count.  In the full-circle edge case it is instead zero, the
    circular component-start count; that edge case is absent from the q=13
    twelve-speed bank audited here.
    """
    events = endpoint_events(speeds, q)

    # At t=0 one zero-centred tooth from each runner is active.
    depth = len(speeds)
    kappa = 0
    protected = 0
    raw_splice_times = 0
    gap_starts = 0

    for _t, event in sorted(events.items()):
        starts = len(event["start"])
        ends = len(event["end"])
        before = depth
        continuing = before - ends
        after = continuing + starts

        if starts and ends:
            raw_splice_times += 1

        # No continuing tooth means the point separates open components.
        if continuing == 0 and after > 0:
            kappa += 1
            if before > 0:
                protected += 1

        # Coverage disappears to the right: a genuine open safe gap starts.
        if before > 0 and after == 0:
            gap_starts += 1

        depth = after

    assert depth == len(speeds)
    defect = kappa - protected
    assert defect == gap_starts
    assert defect >= 0
    return {
        "kappa": kappa,
        "protected": protected,
        "defect": defect,
        "raw_splice_times": raw_splice_times,
    }


def brute_pair_overlap_count(u: int, v: int, q: int) -> int:
    """Count intersecting open-tooth pairs directly by the determinant test."""
    modulus = u * v
    count = 0
    for m in range(u):
        for n in range(v):
            determinant = (m * v - n * u) % modulus
            determinant = min(determinant, modulus - determinant)
            if q * determinant < u + v:
                count += 1
    return count


def pair_overlap_formula(u: int, v: int, q: int) -> int:
    """Closed count formula for the q=5,13,14 banks audited here."""
    g = gcd(u, v)
    return g * (1 + 2 * ((u + v - 1) // (q * g)))


def pair_edge_total(speeds: tuple[int, ...], q: int) -> int:
    return sum(pair_overlap_formula(u, v, q) for u, v in combinations(speeds, 2))


def splice_capacity(speeds: tuple[int, ...], q: int) -> int:
    """Raw oriented end/start capacity from q-divisible pair sums."""
    capacity = 0
    for u, v in combinations(speeds, 2):
        if (u + v) % q != 0:
            continue
        h = (u + v) // q
        g = gcd(u, v)
        if h % g == 0:
            capacity += 2 * g
    return capacity


def first_window_route(speeds: tuple[int, ...]) -> str:
    """Report the THM-766 projective branch for the representative tuple."""
    speeds = tuple(sorted(speeds))
    n = len(speeds)
    m, b, w = speeds[0], speeds[-2], speeds[-1]
    if b >= n * m:
        relation = "=" if b == n * m else ">"
        return f"boundary/super: b {relation} n*m ({b} {relation} {n*m})"
    threshold_ok = b * (n + 2) >= n * n * m
    ks = [
        k
        for k in range(1, n)
        if ((n + 1) * k - 1) * m <= w
        and n * w <= ((n + 1) * k + 1) * b
    ]
    return f"first-window: span_ok={threshold_ok}, compatible_k={ks}"


def pair_label(w: int, q: int = Q13) -> int:
    residue = w % q
    assert residue != 0
    return min(residue, q - residue)


def blocker_counts(speeds: tuple[int, ...], q: int = Q13):
    """Directed complementary-pair blocker counts for raw splice obligations."""
    assert q == Q13
    events = endpoint_events(speeds, q)
    counts: dict[tuple[int, int], int] = defaultdict(int)

    for t, event in events.items():
        if not event["start"] or not event["end"]:
            continue
        targets = {pair_label(w, q) for w in event["start"] + event["end"]}
        blocker_sources = {
            pair_label(w, q)
            for w in speeds
            if circle_distance(w * t) < Fraction(1, q)
        }
        for source in blocker_sources:
            for target in targets:
                if source != target:
                    counts[(source, target)] += 1
    return counts


def tournament_fingerprint(speeds: tuple[int, ...]):
    """Gauge blocker counts into a six-vertex tournament and report fingerprints."""
    n = 6
    blocks = blocker_counts(speeds)
    adjacency = [[False] * (n + 1) for _ in range(n + 1)]
    flips = 0

    for i, j in combinations(range(1, n + 1), 2):
        ij = blocks.get((i, j), 0)
        ji = blocks.get((j, i), 0)
        if ij >= ji:  # the equality case is the tie Hamiltonian path i -> j
            winner, loser = i, j
        else:
            winner, loser = j, i
            flips += 1
        adjacency[winner][loser] = True

    scores = [sum(adjacency[i][j] for j in range(1, n + 1)) for i in range(1, n + 1)]
    cycles = 0
    for a, b, c in combinations(range(1, n + 1), 3):
        if (
            adjacency[a][b]
            and adjacency[b][c]
            and adjacency[c][a]
        ) or (
            adjacency[a][c]
            and adjacency[c][b]
            and adjacency[b][a]
        ):
            cycles += 1

    reach = [[adjacency[i][j] or i == j for j in range(n + 1)] for i in range(n + 1)]
    for k in range(1, n + 1):
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])

    unseen = set(range(1, n + 1))
    scc_sizes = []
    while unseen:
        i = min(unseen)
        component = {j for j in unseen if reach[i][j] and reach[j][i]}
        unseen -= component
        scc_sizes.append(len(component))

    blocker_arcs = [
        f"{i}->{j}:{count}"
        for (i, j), count in sorted(blocks.items())
        if count
    ]
    return {
        "scores": scores,
        "cycles3": cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "edge_flips": flips,
        "blocker_arcs": blocker_arcs,
    }


def print_representative_bank() -> None:
    representatives = [
        ("AP12_q13", tuple(range(1, 13)), 13),
        ("2AP_q13", tuple(2 * r for r in range(1, 13)), 13),
        ("q5_perm_lift_beater", (1, 3, 4, 7), 5),
        ("GW13_q14", tuple(range(1, 12)) + (13, 24), 14),
        ("lift_1_to_14", tuple(range(2, 13)) + (14,), 13),
        ("lift_2_to_15", (1,) + tuple(range(3, 13)) + (15,), 13),
        ("lift_12_to_25", tuple(range(1, 12)) + (25,), 13),
        ("lift_7_to_98", tuple(range(1, 7)) + tuple(range(8, 13)) + (98,), 13),
        ("all_odd_lift", (1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25), 13),
    ]

    print("REPRESENTATIVE ENDPOINT EULER DEFECTS")
    for name, speeds, q in representatives:
        speeds = tuple(sorted(speeds))
        stats = endpoint_stats(speeds, q)
        capacity = splice_capacity(speeds, q)
        edges = pair_edge_total(speeds, q)
        print(
            f"{name}: q={q}, V={sum(speeds)}, pair_edges={edges}, "
            f"kappa={stats['kappa']}, Pstar={capacity}, "
            f"P={stats['protected']}, defect={stats['defect']}, "
            f"raw_splice_times={stats['raw_splice_times']}, "
            f"target_tight={stats['defect'] == 0}"
        )
        print(f"  THM766={first_window_route(speeds)}")


def print_pair_formula_audit() -> None:
    banks = [
        ("q13_height1_speed_universe", 13, tuple(range(1, 13)) + tuple(range(14, 26))),
        ("q5_perm_lift_beater", 5, (1, 3, 4, 7)),
        ("q14_GW", 14, tuple(range(1, 12)) + (13, 24)),
    ]
    total = 0
    print("PAIR OVERLAP FORMULA AUDIT")
    for name, q, speeds in banks:
        checked = 0
        mismatches = []
        for u, v in combinations(speeds, 2):
            brute = brute_pair_overlap_count(u, v, q)
            formula = pair_overlap_formula(u, v, q)
            checked += 1
            if brute != formula:
                mismatches.append((u, v, brute, formula))
        total += checked
        print(f"{name}: q={q}, checked_pairs={checked}, mismatches={len(mismatches)}")
        assert not mismatches
    print(f"formula_checks_total={total}, formula_mismatches_total=0")


def print_height_one_audit() -> None:
    categories = Counter()
    defect_histogram = Counter()
    gcd_histogram = Counter()
    zero_rows = []
    defect_two_rows = []
    primitive_defects = []

    for bits in product((0, 1), repeat=12):
        if not any(bits):
            continue
        labelled = tuple(r + 13 * bits[r - 1] for r in range(1, 13))
        speeds = tuple(sorted(labelled))
        stats = endpoint_stats(speeds, 13)
        capacity = splice_capacity(speeds, 13)
        common_gcd = reduce(gcd, speeds)
        defect = stats["defect"]

        if stats["kappa"] > capacity:
            category = "rank_gt_splice_capacity"
        elif defect > 0:
            category = "capacity_survives_but_blocked"
        else:
            category = "zero_defect"

        categories[category] += 1
        defect_histogram[defect] += 1
        gcd_histogram[common_gcd] += 1
        if common_gcd == 1:
            primitive_defects.append(defect)
        if defect == 0:
            zero_rows.append((bits, speeds, common_gcd, stats, capacity))
        if defect == 2:
            defect_two_rows.append((bits, speeds, stats, capacity))

    assert sum(categories.values()) == 4095
    assert categories == Counter(
        {
            "rank_gt_splice_capacity": 4085,
            "capacity_survives_but_blocked": 9,
            "zero_defect": 1,
        }
    )
    assert len(zero_rows) == 1
    assert zero_rows[0][1] == tuple(2 * r for r in range(1, 13))
    assert zero_rows[0][2] == 2
    assert min(primitive_defects) == 2
    assert len(defect_two_rows) == 4

    print("HEIGHT-ONE q=13 LIFT CUBE")
    print("patterns_nonzero=4095")
    print("categories=" + repr(dict(sorted(categories.items()))))
    print("gcd_histogram=" + repr(dict(sorted(gcd_histogram.items()))))
    print("defect_histogram=" + repr(dict(sorted(defect_histogram.items()))))
    print(f"primitive_patterns={len(primitive_defects)}, primitive_min_defect={min(primitive_defects)}")
    print("zero_defect_rows=1")
    bits, speeds, common_gcd, stats, capacity = zero_rows[0]
    print(
        "  bits=" + "".join(map(str, bits))
        + f", speeds={list(speeds)}, gcd={common_gcd}, "
        + f"kappa={stats['kappa']}, Pstar={capacity}, P={stats['protected']}"
    )
    print("defect_two_rows=4")
    for bits, speeds, stats, capacity in defect_two_rows:
        print(
            "  bits=" + "".join(map(str, bits))
            + f", speeds={list(speeds)}, kappa={stats['kappa']}, "
            + f"Pstar={capacity}, P={stats['protected']}"
        )


def print_tournament_audit() -> None:
    rows = [
        ("AP12", tuple(range(1, 13))),
        ("lift_1_to_14", tuple(range(2, 13)) + (14,)),
        ("lift_2_to_15", (1,) + tuple(range(3, 13)) + (15,)),
        ("lift_7_to_98", tuple(range(1, 7)) + tuple(range(8, 13)) + (98,)),
        ("height1_blocker_residual", (1, 2, 4, 6, 7, 10, 11, 12, 16, 18, 21, 22)),
    ]
    print("TOURNAMENT ANALYSIS (DIAGNOSTIC QUOTIENT)")
    print("vertices=complementary-pair splice obligations 1..6")
    print("observable=directed exclusive blocker count")
    print("switch/gauge=sign(count(A blocks B)-count(B blocks A))")
    print("tie_hamiltonian_path=1->2->3->4->5->6")
    for name, speeds in rows:
        stats = endpoint_stats(tuple(sorted(speeds)), 13)
        fingerprint = tournament_fingerprint(tuple(sorted(speeds)))
        print(
            f"{name}: defect={stats['defect']}, scores={fingerprint['scores']}, "
            f"cycles3={fingerprint['cycles3']}, scc_sizes={fingerprint['scc_sizes']}, "
            f"edge_flips={fingerprint['edge_flips']}, "
            f"blocker_arcs={fingerprint['blocker_arcs']}"
        )
    print("guardrail: lift_2_to_15 and lift_7_to_98 have the AP tournament fingerprint")
    print("           but positive endpoint defects; tournament data alone is not theorem-safe.")


def main() -> None:
    print("=" * 88)
    print("LRC(13) ENDPOINT / SPLICE DEFECT AUDIT (codex-2026-07-14-S3)")
    print("=" * 88)
    print("method=exact Fraction endpoint sweep; no floating point")
    print("primary_vertices=individual danger teeth; secondary_vertices=splice obligations")
    print("identity=endpoint_defect = kappa(open tooth graph) - protected_splices")
    print("q13_capacity=Pstar=2*sum_{r=1}^6 gcd(w_r,w_{13-r})")
    print()
    print_pair_formula_audit()
    print()
    print_representative_bank()
    print()
    print_height_one_audit()
    print()
    print_tournament_audit()
    print()
    print("PRESERVE / DESTROY")
    print("unit fractions: preserve forced equality witnesses; destroy lift geometry and gap interiors")
    print("binding pairs: preserve raw splice capacity; destroy third-runner blockers and overlap rank")
    print("splice obligations: preserve blocker incidence; destroy unspliced tooth components")
    print("individual teeth + endpoint events: preserve the exact target-cover predicate")
    print("tournament quotient: preserves blocker direction; destroys simultaneous connectivity")
    print()
    print("SCOPE")
    print("THM-763 makes the primitive tight branch finite, not empty.")
    print("THM-766 supplies local component-tooth cones, not global splice coherence.")
    print("The remaining lemma is: zero endpoint defect at q=13 forces one common AP scale.")


if __name__ == "__main__":
    main()
