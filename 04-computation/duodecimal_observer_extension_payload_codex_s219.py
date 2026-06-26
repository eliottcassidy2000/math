#!/usr/bin/env python3
"""S219: duodecimal observer-extension/cut-payload synthesis.

This is a small exact ledger, not a new brute-force tournament enumerator.  It
reuses the verified S211-S217 counts and asks where the user's recurring 12
comes from around the first shifted A000568/rooted-perspective failure.

Main guardrail:

    P(5) = 48, U(6) = 56, U(5) = 12, U(4) = 4.

So 48 + 12 is 60, not 56.  The useful exact identity is instead

    U(6) = P(5) + U(5) - U(4) = 48 + 12 - 4.

The dozen is real, but it overlaps by the four 4-tournament classes.  This
script records that as an observer-extension/cut-payload inclusion-exclusion
window and then connects it to incident words, ordered-pair edge sectors,
deletion fibers, and S217 rectangle/hourglass residues.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, gcd
from typing import Iterable


U = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
P = {1: 1, 2: 2, 3: 4, 4: 12, 5: 48, 6: 296, 7: 3040}
SELF_CONVERSE = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}

S213_EDGE_SECTORS = [
    ("rooted_5_perspective_plus_incident_word", 1408, "ordered-pair 6-perspectives", "1408/1408 exact"),
    ("forget_endpoint_role", 704, "directed-edge/unordered-pair perspectives", "role quotient"),
    ("sector_size_deck", 35, "individual signatures", "55/56 class decks"),
    ("internal_sector_deck", 59, "individual signatures", "55/56 class decks"),
    ("cross_sector_orientation_deck", 362, "individual signatures", "56/56 class decks"),
    ("full_ordered_pair_sector_deck", 422, "individual signatures", "56/56 class decks"),
]

S216_DELETION_FIBERS = {
    "parent_aut_sizes_on_5_classes": {1: 7, 3: 4, 5: 1},
    "diagonal_word_orbit_sizes": {1: 258, 3: 32, 5: 6},
    "rooted_extension_orbits_per_6_class": {2: 5, 4: 10, 6: 41},
    "raw_labelled_extensions_per_6_class": {2: 5, 4: 10, 6: 15, 8: 12, 10: 11, 12: 2, 14: 1},
    "unique_deletion_parent_classes_per_6_class": {1: 3, 2: 6, 3: 13, 4: 13, 5: 17, 6: 4},
    "max_repeated_deletion_parent_multiplicity_per_6_class": {1: 4, 2: 27, 3: 17, 4: 3, 5: 2, 6: 3},
}


def frac(numer: int, denom: int) -> str:
    f = Fraction(numer, denom)
    return f"{f.numerator}/{f.denominator}" if f.denominator != 1 else str(f.numerator)


def ratio_text(numer: int, denom: int) -> str:
    f = Fraction(numer, denom)
    return f"{numer}/{denom} = {frac(numer, denom)} = {float(f):.6g}"


def weighted_average(hist: dict[int, int]) -> Fraction:
    total = sum(hist.values())
    return Fraction(sum(k * v for k, v in hist.items()), total)


def local_layer(k: int) -> dict[str, int]:
    lines = k * (k + 1)
    rank = 2 * k
    return {
        "k": k,
        "lines": lines,
        "rank": rank,
        "rectangles": lines - rank,
    }


def full_flow(n: int) -> dict[str, int]:
    lines = 2 * comb(n, 3)
    rank = comb(n, 2) - 1 if n >= 2 else 0
    rectangles = 2 * comb(n - 1, 3) if n >= 4 else 0
    hourglass = comb(n - 2, 2) if n >= 4 else 0
    return {
        "n": n,
        "lines": lines,
        "rank": rank,
        "redundancy": lines - rank,
        "rectangles": rectangles,
        "hourglass": hourglass,
    }


def fixed_path_flow(n: int) -> dict[str, int]:
    lines = 2 * comb(n - 1, 3) if n >= 4 else 0
    rank = comb(n - 1, 2) - 1 if n >= 3 else 0
    rectangles = 2 * comb(n - 2, 3) if n >= 5 else 0
    hourglass = comb(n - 3, 2) if n >= 5 else 0
    return {
        "n": n,
        "lines": lines,
        "rank": rank,
        "redundancy": lines - rank,
        "rectangles": rectangles,
        "hourglass": hourglass,
    }


def count_directed_3cycles(vertices: list[str], beats: dict[tuple[str, str], bool]) -> int:
    total = 0
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            for k in range(j + 1, len(vertices)):
                a, b, c = vertices[i], vertices[j], vertices[k]
                ab = beats[(a, b)]
                bc = beats[(b, c)]
                ac = beats[(a, c)]
                # A 3-vertex tournament is cyclic iff it is not transitive.
                if (ab and bc and not ac) or (not ab and not bc and ac):
                    total += 1
    return total


def strongly_connected_components(vertices: list[str], beats: dict[tuple[str, str], bool]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    rev = {v: [] for v in vertices}
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            if beats[(a, b)]:
                adj[a].append(b)
                rev[b].append(a)
            else:
                adj[b].append(a)
                rev[a].append(b)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(comp)
    return comps


def hamiltonian_path_count(vertices: list[str], beats: dict[tuple[str, str], bool]) -> tuple[int, list[str]]:
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    parent: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                a, b = vertices[last], vertices[nxt]
                forward = beats[(a, b)] if (a, b) in beats else not beats[(b, a)]
                if forward:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + count
                    parent.setdefault(key, last)

    full = (1 << n) - 1
    total = sum(dp.get((full, last), 0) for last in range(n))
    if not total:
        return 0, []
    end = next(last for last in range(n) if dp.get((full, last), 0))
    path = [end]
    mask = full
    while mask != (1 << path[-1]):
        key = (mask, path[-1])
        prev = parent[key]
        path.append(prev)
        mask ^= 1 << key[1]
    path.reverse()
    return total, [vertices[i] for i in path]


def carrier_tournament() -> dict[str, object]:
    # Vertices are proof carriers, not tournament runners.  The observable is:
    # how much observer/cut/cycle payload survives, penalized by quotient risk.
    features = {
        "endpoint_owner_packet_sheaf": (5, 5, 4, 5, 0),
        "cross_sector_orientation_word": (5, 4, 3, 4, 1),
        "observer_extension_cut_payload": (5, 5, 2, 4, 1),
        "deletion_parent_fiber_profile": (4, 4, 3, 4, 2),
        "rectangle_hourglass_cycle_residue": (3, 3, 5, 4, 2),
        "incident_word_orbit_under_aut": (4, 4, 1, 3, 2),
        "ordered_pair_edge_sector_deck": (4, 3, 2, 3, 2),
        "rooted_node_perspective_cache": (2, 1, 1, 2, 3),
        "fixed_path_half_tiling_shadow": (2, 1, 2, 2, 4),
        "raw_A000568_class_count": (1, 0, 0, 1, 5),
        "raw_labelled_word_count": (1, 0, 0, 0, 6),
    }
    # feature tuple = observer, cut, cycle residue, LRC applicability, quotient risk
    scores = {
        name: 4 * a + 3 * b + 2 * c + 2 * d - 2 * e
        for name, (a, b, c, d, e) in features.items()
    }
    vertices = list(features)
    beats: dict[tuple[str, str], bool] = {}
    flips: list[tuple[str, str]] = []
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            if scores[a] == scores[b]:
                # Tie Hamiltonian path: prefer the carrier with less quotient risk,
                # then lexicographic order for reproducibility.
                orient = (features[a][4], a) < (features[b][4], b)
            else:
                orient = scores[a] > scores[b]
            beats[(a, b)] = orient
            if not orient:
                flips.append((b, a))
    outdegrees = Counter()
    for v in vertices:
        degree = 0
        for w in vertices:
            if v == w:
                continue
            if (v, w) in beats:
                degree += int(beats[(v, w)])
            else:
                degree += int(not beats[(w, v)])
        outdegrees[degree] += 1
    h_count, h_path = hamiltonian_path_count(vertices, beats)
    return {
        "scores": scores,
        "score_hist": dict(sorted(outdegrees.items())),
        "directed_3cycles": count_directed_3cycles(vertices, beats),
        "scc_sizes": sorted(len(c) for c in strongly_connected_components(vertices, beats)),
        "hamiltonian_path_count": h_count,
        "one_hamiltonian_path": h_path,
        "edge_flips_from_input_order": flips,
    }


def print_shift_ladder() -> None:
    print("1. SHIFTED ROOTED-PERSPECTIVE LADDER")
    print("   m  U(m)  P(m)  U(m+1)  defect  P/U(m+1)  defect/U(m+1)  root_tax=mU-P")
    for m in range(1, 7):
        defect = U[m + 1] - P[m]
        tax = m * U[m] - P[m]
        print(
            f"   {m}  {U[m]:4d}  {P[m]:4d}  {U[m + 1]:6d}  {defect:6d}"
            f"  {frac(P[m], U[m + 1]):>10s}  {frac(defect, U[m + 1]):>14s}  {tax:14d}"
        )
    print()
    print("   First failure window:")
    print(f"     P(5) = {P[5]} = 4 * {U[5]}")
    print(f"     U(6) = {U[6]} = 4 * {U[5]} + {U[6] - 4 * U[5]}")
    print(f"     48 + 12 = {P[5] + U[5]}, so the exact count is not 48+12.")
    print(f"     U(6) = P(5) + U(5) - U(4) = {P[5]} + {U[5]} - {U[4]}")
    print(f"     defect = U(6)-P(5) = {U[6] - P[5]} = U(5)-U(4)")
    print(f"     P(5)/U(6) = {frac(P[5], U[6])}; defect/U(6) = {frac(U[6] - P[5], U[6])}")
    print(f"     U(5)/U(6) = {frac(U[5], U[6])}; U(4)/U(6) = {frac(U[4], U[6])}")
    print("     Denominator-14 form: 12/14 + 3/14 - 1/14 = 1.")
    print()
    print("   Recurrence guardrail:")
    print("     The identity U(6)=P(5)+U(5)-U(4) is exact at the first failure.")
    print(
        f"     It is not a stable recurrence: U(7)-P(6)={U[7]-P[6]}, "
        f"but U(6)-U(5)={U[6]-U[5]}."
    )


def print_duodecimal_bridge() -> None:
    print()
    print("2. THE DOZEN BRIDGE P(4) = U(5) = SC(6)")
    print("   object                                      value  fraction_of_U6")
    rows = [
        ("P(4): rooted 4-perspectives", P[4]),
        ("U(5): 5-tournament isomorphism classes", U[5]),
        ("SC(6): self-converse 6-classes", SELF_CONVERSE[6]),
        ("P(5): rooted 5-perspective cache", P[5]),
        ("U(6)-P(5): first observer-cut defect", U[6] - P[5]),
        ("U(4): overlap subtracted by inclusion-exclusion", U[4]),
        ("U(6): full 6-class quotient", U[6]),
    ]
    for label, value in rows:
        print(f"   {label:48s} {value:5d}  {frac(value, U[6]):>12s}")
    print()
    print("   Exact relations:")
    print(f"     P(4)=U(5)=SC(6)={U[5]}")
    print(f"     P(5)=5*U(5)-SC(6)={5*U[5]}-{SELF_CONVERSE[6]}={P[5]}")
    print(f"     U(6)=5*U(5)-U(4)={5*U[5]}-{U[4]}={U[6]}")
    print(f"     U(6)-P(5)=SC(6)-U(4)={SELF_CONVERSE[6]}-{U[4]}={U[6]-P[5]}")
    print("     Reading: the dozen is a control/fold slice; the additive defect is the")
    print("     dozen after quotienting out the four-class overlap.")


def print_incident_sector_fibers() -> None:
    print()
    print("3. INCIDENT WORDS, EDGE SECTORS, AND DELETION FIBERS")
    print("   Ordered-pair / edge-sector readout from S213:")
    for name, value, kind, note in S213_EDGE_SECTORS:
        print(f"     {name:42s} {value:5d}  {kind:40s} {note}")
    print("     The only 55/56 collision is a converse pair; cross-sector orientation separates it.")
    print()
    print("   Deletion-fiber readout from S216 at 5 -> 6:")
    for name, hist in S216_DELETION_FIBERS.items():
        avg = weighted_average(hist)
        avg_text = frac(avg.numerator, avg.denominator)
        print(f"     {name}: {hist}; weighted_avg={avg_text}")
    print("     Rooted child presentations average 296/56 = 37/7 per unrooted 6-class.")
    print("     Unique deletion-parent classes average 215/56; the sink quotient is not")
    print("     just source class plus a scalar word count.")


def print_rectangle_hourglass() -> None:
    print()
    print("4. RECTANGLE / HOURGLASS CYCLE RESIDUES")
    print("   Local K_{k,k+1} lines:")
    print("   k  lines  rank  rectangle_residue_dim")
    for k in range(1, 8):
        row = local_layer(k)
        print(f"   {k}  {row['lines']:5d}  {row['rank']:4d}  {row['rectangles']:21d}")
    print("   The first dozen line block is K_{3,4}: 12 lines = 6 rank + 6 rectangles.")
    print()
    print("   Global full flow and fixed-path flow:")
    print("   n  full_lines full_rank full_red rect hour  fixed_lines fixed_rank fixed_red rect hour")
    for n in range(4, 8):
        f = full_flow(n)
        pth = fixed_path_flow(n)
        print(
            f"   {n}  {f['lines']:10d} {f['rank']:9d} {f['redundancy']:8d}"
            f" {f['rectangles']:4d} {f['hourglass']:4d}"
            f"  {pth['lines']:11d} {pth['rank']:10d} {pth['redundancy']:9d}"
            f" {pth['rectangles']:4d} {pth['hourglass']:4d}"
        )
    print("   Full redundancy = local rectangles + hourglass cycles.")
    print("   Fixed-path redundancy = 2*C(n-2,3) + C(n-3,2).")
    print("   These residues are the layer-flow version of observer-cut payload:")
    print("   zero residue descends to potentials; nonzero residue names the hidden sidecar.")


def print_abstract_payload() -> None:
    print()
    print("5. OBSERVER-EXTENSION / CUT-PAYLOAD ABSTRACT")
    applications = [
        (
            "A000568 perspectives",
            "parent class + incident word + deletion fiber",
            "rooted cache alone loses the 8 first-failure states",
        ),
        (
            "ordered-pair edge sectors",
            "tip/tail endpoint roles plus cross-sector orientation",
            "sector sizes/internal decks collide at 55/56; orientation gives 56/56",
        ),
        (
            "fixed-path diagonal flow",
            "line potentials plus rectangle/hourglass residues",
            "raw k(k+1) lines duplicate; residues are the certificate",
        ),
        (
            "LRC endpoint-owner packets",
            "observer source, owner strip, route/status payload",
            "apex/residue quotients are unsafe without magnitude and owner data",
        ),
        (
            "Haar/fixed-margin rectangles",
            "2-by-2 zeta residue and owner-strip handoff",
            "same rectangle law as the K_{k,k+1} parity carrier",
        ),
        (
            "pair-good decoy/residual capacitor",
            "blocker teeth, edge-tail/tip sectors, residual current",
            "counts are scouts; the generator key is a sidecar",
        ),
        (
            "matrix sidecar observability",
            "rows are coarse-fiber pairs, columns are hidden coordinates",
            "a quotient is safe only when changing pairs are separated or discharged",
        ),
        (
            "tournament spectrum for LRC",
            "set of winding iso classes across phases plus binding scale",
            "single apex class forgets magnitude; the spectrum keeps the cut",
        ),
    ]
    print("   domain                         retained cut payload                         warning")
    for domain, payload, warning in applications:
        print(f"   {domain:30s} {payload:44s} {warning}")
    print()
    print("   Controlled-forgetting rule:")
    print("     A quotient may forget the observer-extension cut only if the target")
    print("     predicate is fiber-constant, reconstructible from retained sidecars,")
    print("     annihilated by a dual certificate, descended by a family lemma, or")
    print("     routed to a named residual sector.  Otherwise keep the cut payload.")


def print_tournament_analysis() -> None:
    print()
    print("6. TOURNAMENT ANALYSIS OVER CARRIERS")
    ta = carrier_tournament()
    print("   Pairwise observable: retained observer/cut/cycle payload minus quotient risk.")
    print("   Switch/gauge: higher observer+cut+residue+LRC score orients the edge.")
    print("   Vertices challenged: runners, arcs, gaps, layer boundaries, wall events,")
    print("   residues, cover arcs, Fourier modes, matroid circuits, and proof obligations.")
    print("   Chosen vertices are proof carriers, not raw runners.")
    print(f"   score_hist={ta['score_hist']}")
    print(f"   directed_3cycles={ta['directed_3cycles']}")
    print(f"   SCC_sizes={ta['scc_sizes']}")
    print(f"   Hamiltonian_path_count={ta['hamiltonian_path_count']}")
    print("   one Hamiltonian path:")
    print("     " + " > ".join(ta["one_hamiltonian_path"]))  # type: ignore[index]
    print("   edge_flips_from_input_order:")
    for a, b in ta["edge_flips_from_input_order"]:  # type: ignore[index]
        print(f"     {a} > {b}")


def main() -> None:
    print("=" * 80)
    print("S219: duodecimal observer-extension / cut-payload synthesis")
    print("=" * 80)
    print_shift_ladder()
    print_duodecimal_bridge()
    print_incident_sector_fibers()
    print_rectangle_hourglass()
    print_abstract_payload()
    print_tournament_analysis()
    print()
    print("BOTTOM LINE")
    print("  The dozen is not the first-failure additive defect.  It is a reusable")
    print("  control slice: P(4)=U(5)=SC(6)=12, source/sink diagonal words each give")
    print("  12 slices at 5->6, and K_{3,4} is the first 12-line rectangle carrier.")
    print("  The exact first-failure arithmetic is duodecimal inclusion-exclusion:")
    print("  U(6)=48+12-4, or 12/14 + 3/14 - 1/14 = 1.  The missing 8 is the dozen")
    print("  minus the four-class overlap, i.e. observer-extension cut payload after")
    print("  controlled forgetting.")


if __name__ == "__main__":
    main()
