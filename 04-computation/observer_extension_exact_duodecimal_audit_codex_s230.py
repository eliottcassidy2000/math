#!/usr/bin/env python3
"""S230: observer-extension exact duodecimal deletion-sector audit.

This continues the S211..S217 tournament/LRC stack.  The prompt points at the
near identity around

    rooted 5-perspectives = 48,
    six-tournament classes = 56,
    self-converse six-classes = 12.

The literal arithmetic is 48 + 12 = 60, so the useful structure must include
an overlap/correction term.  This script audits the exact small tournament
counts, deletion fibers, ordered-pair sectors, Burnside sources, and the S217
rectangle/hourglass flow residues to make that hidden term explicit.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial


MAX_N = 6
PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PERMS: dict[int, list[tuple[int, ...]]] = {}
CANON_CACHE: dict[tuple[int, int], int] = {}
CLASSES: dict[int, list[int]] = {}
AUT_CACHE: dict[tuple[int, int], list[tuple[int, ...]]] = {}


def pair_index(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = {pair: idx for idx, pair in enumerate(combinations(range(n), 2))}
    return PAIR_INDEX[n]


def perms(n: int) -> list[tuple[int, ...]]:
    if n not in PERMS:
        PERMS[n] = list(permutations(range(n)))
    return PERMS[n]


def bit_count(n: int) -> int:
    return comb(n, 2)


def edge(mask: int, n: int, i: int, j: int) -> int:
    if i == j:
        raise ValueError("no loop edge")
    if i < j:
        return (mask >> pair_index(n)[(i, j)]) & 1
    return 1 - ((mask >> pair_index(n)[(j, i)]) & 1)


def relabel(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """Return mask in new labels where new vertex a is old vertex perm[a]."""
    out = 0
    out_idx = pair_index(n)
    for a, b in combinations(range(n), 2):
        if edge(mask, n, perm[a], perm[b]):
            out |= 1 << out_idx[(a, b)]
    return out


def canonical(mask: int, n: int) -> int:
    key = (n, mask)
    if key not in CANON_CACHE:
        CANON_CACHE[key] = min(relabel(mask, n, perm) for perm in perms(n))
    return CANON_CACHE[key]


def converse(mask: int, n: int) -> int:
    return mask ^ ((1 << bit_count(n)) - 1)


def all_classes(n: int) -> list[int]:
    if n not in CLASSES:
        CLASSES[n] = sorted({canonical(mask, n) for mask in range(1 << bit_count(n))})
    return CLASSES[n]


def automorphisms(mask: int, n: int) -> list[tuple[int, ...]]:
    key = (n, mask)
    if key not in AUT_CACHE:
        AUT_CACHE[key] = [perm for perm in perms(n) if relabel(mask, n, perm) == mask]
    return AUT_CACHE[key]


def orbit_count_on_points(group: list[tuple[int, ...]], n: int) -> int:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for perm in group:
        for a, b in enumerate(perm):
            union(a, b)
    return len({find(i) for i in range(n)})


def vertex_orbit_count(mask: int, n: int) -> int:
    return orbit_count_on_points(automorphisms(mask, n), n)


def rooted_count(n: int) -> int:
    return sum(vertex_orbit_count(mask, n) for mask in all_classes(n))


def word_orbits_for_parent(mask: int, n: int) -> int:
    group = automorphisms(mask, n)
    seen: set[int] = set()
    total = 0
    for word in range(1 << n):
        if word in seen:
            continue
        total += 1
        orbit = set()
        for perm in group:
            moved = 0
            for new_v, old_v in enumerate(perm):
                if (word >> old_v) & 1:
                    moved |= 1 << new_v
            orbit.add(moved)
        seen.update(orbit)
    return total


def extension_orbit_count(n: int) -> int:
    return sum(word_orbits_for_parent(mask, n) for mask in all_classes(n))


def induced_mask(mask: int, n: int, verts: tuple[int, ...] | list[int]) -> int:
    verts = tuple(verts)
    m = len(verts)
    out = 0
    out_idx = pair_index(m)
    for a, b in combinations(range(m), 2):
        if edge(mask, n, verts[a], verts[b]):
            out |= 1 << out_idx[(a, b)]
    return out


def delete_vertex(mask: int, n: int, v: int) -> int:
    verts = tuple(i for i in range(n) if i != v)
    return induced_mask(mask, n, verts)


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(sum(edge(mask, n, i, j) for j in range(n) if j != i) for i in range(n)))


def cyclic_triangle_count(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        eab = edge(mask, n, a, b)
        ebc = edge(mask, n, b, c)
        eca = edge(mask, n, c, a)
        if eab == ebc == eca:
            total += 1
    return total


def hamiltonian_path_count(mask: int, n: int) -> int:
    total = 0
    for perm in perms(n):
        if all(edge(mask, n, perm[i], perm[i + 1]) for i in range(n - 1)):
            total += 1
    return total


def integer_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    out = []
    for first in range(min(max_part, n), 0, -1):
        for rest in integer_partitions(n - first, first):
            out.append((first,) + rest)
    return out


def representative_permutation(parts: tuple[int, ...]) -> dict[int, int]:
    mapping = {}
    start = 0
    for length in parts:
        cycle = list(range(start, start + length))
        for i, v in enumerate(cycle):
            mapping[v] = cycle[(i + 1) % length]
        start += length
    return mapping


def edge_orbit_count_under_perm(n: int, mapping: dict[int, int]) -> int:
    pairs = [tuple(sorted(pair)) for pair in combinations(range(n), 2)]
    pair_set = set(pairs)
    seen: set[tuple[int, int]] = set()
    count = 0
    for pair in pairs:
        if pair in seen:
            continue
        count += 1
        cur = pair
        while cur not in seen:
            seen.add(cur)
            a, b = cur
            cur = tuple(sorted((mapping[a], mapping[b])))
            if cur not in pair_set:
                raise AssertionError("permutation left pair universe")
    return count


def cycle_type_size(parts: tuple[int, ...]) -> int:
    mult = Counter(parts)
    denom = 1
    for length, amount in mult.items():
        denom *= (length**amount) * factorial(amount)
    return factorial(sum(parts)) // denom


def burnside_terms(n: int) -> list[dict[str, object]]:
    rows = []
    for parts in integer_partitions(n):
        if any(part % 2 == 0 for part in parts):
            continue
        mapping = representative_permutation(parts)
        edge_orbits = edge_orbit_count_under_perm(n, mapping)
        class_size = cycle_type_size(parts)
        fixed = 2**edge_orbits
        rows.append(
            {
                "type": parts,
                "perms": class_size,
                "edge_orbits": edge_orbits,
                "fixed": fixed,
                "burnside_sum": class_size * fixed,
            }
        )
    return rows


def tournament_count_by_burnside(n: int) -> int:
    return sum(int(row["burnside_sum"]) for row in burnside_terms(n)) // factorial(n)


def ordered_pair_sector_signature(
    mask: int, n: int, a: int, b: int, include_internal: bool, include_cross: bool
) -> tuple[object, ...]:
    keys = ((0, 0), (0, 1), (1, 0), (1, 1))
    sectors: dict[tuple[int, int], list[int]] = {key: [] for key in keys}
    for x in range(n):
        if x in (a, b):
            continue
        sectors[(edge(mask, n, a, x), edge(mask, n, b, x))].append(x)
    sizes = tuple(len(sectors[key]) for key in keys)
    parts: list[object] = [sizes]
    if include_internal:
        internal = []
        for key in keys:
            verts = sectors[key]
            internal.append((len(verts), canonical(induced_mask(mask, n, verts), len(verts))))
        parts.append(tuple(internal))
    if include_cross:
        cross = []
        for i, key_i in enumerate(keys):
            for key_j in keys[i + 1 :]:
                count_i_to_j = 0
                for x in sectors[key_i]:
                    for y in sectors[key_j]:
                        count_i_to_j += edge(mask, n, x, y)
                cross.append(count_i_to_j)
        parts.append(tuple(cross))
    return tuple(parts)


def ordered_pair_deck(mask: int, n: int, include_internal: bool, include_cross: bool) -> tuple[tuple[object, ...], ...]:
    return tuple(
        sorted(
            ordered_pair_sector_signature(mask, n, a, b, include_internal, include_cross)
            for a in range(n)
            for b in range(n)
            if a != b
        )
    )


def deletion_profile(mask: int, n: int) -> tuple[int, tuple[int, ...]]:
    parents = [canonical(delete_vertex(mask, n, v), n - 1) for v in range(n)]
    return len(set(parents)), tuple(sorted(parents))


def local_bridge(k: int) -> tuple[int, int, int]:
    lines = k * (k + 1)
    rank = 2 * k
    return lines, rank, lines - rank


def global_flow(n: int) -> tuple[int, int, int, int, int]:
    lines = 2 * comb(n, 3)
    rank = comb(n, 2) - 1
    red = lines - rank
    local_rectangles = 2 * comb(n - 1, 3)
    hourglasses = comb(n - 2, 2)
    return lines, rank, red, local_rectangles, hourglasses


def print_count_ladders() -> None:
    print("1. COUNT LADDERS AND THE DUODECIMAL GATE")
    print("   n  U(n)  Burnside  rooted P(n)  extension E(n->n+1)  self-converse")
    for n in range(1, MAX_N + 1):
        u_n = len(all_classes(n))
        b_n = tournament_count_by_burnside(n)
        p_n = rooted_count(n)
        e_n = extension_orbit_count(n) if n < MAX_N else None
        sc_n = sum(1 for mask in all_classes(n) if canonical(converse(mask, n), n) == mask)
        e_s = f"{e_n:20d}" if e_n is not None else " " * 20
        print(f"   {n:1d}  {u_n:4d}  {b_n:8d}  {p_n:11d}  {e_s}  {sc_n:13d}")

    u4 = len(all_classes(4))
    u5 = len(all_classes(5))
    u6 = len(all_classes(6))
    p4 = rooted_count(4)
    p5 = rooted_count(5)
    sc6 = sum(1 for mask in all_classes(6) if canonical(converse(mask, 6), 6) == mask)
    defect = u6 - p5
    print()
    print("   Duodecimal coincidences:")
    print(f"     P(4) = {p4} = 1 dozen rooted/node perspectives on 4-classes.")
    print(f"     U(5) = {u5} = 1 dozen unrooted 5-tournament classes.")
    print(f"     SC(6) = {sc6} = 1 dozen self-converse 6-classes.")
    print(f"     K_{{4,5}} rectangle redundancy = 4*(4-1) = {local_bridge(4)[2]}.")
    print()
    print("   Literal arithmetic check and corrected ledgers:")
    print(f"     48 + 12 = {p5 + sc6}, not {u6}.")
    print(f"     U(6) = P(5) + {defect} = {p5} + {defect} = {u6}.")
    print(f"     U(6) = P(5) + SC(6) - U(4) = {p5} + {sc6} - {u4} = {u6}.")
    print(f"     defect/U(6) = {Fraction(defect, u6)}; P(5)/U(6) = {Fraction(p5, u6)}.")
    print(f"     defect/dozen = {Fraction(defect, 12)}; overlap/dozen = U(4)/12 = {Fraction(u4, 12)}.")
    print(f"     U(6)/dozen = {Fraction(u6, 12)}; P(5)/dozen = {Fraction(p5, 12)}.")


def print_burnside() -> None:
    print()
    print("2. BURNSIDE TERMS: WHERE U(4), U(5), U(6) COME FROM")
    for n in (4, 5, 6):
        print(f"   n={n}:")
        total = 0
        for row in burnside_terms(n):
            total += int(row["burnside_sum"])
            print(
                "     type={type!s:14s} perms={perms:3d} edge_orbits={edge_orbits:2d} "
                "fixed={fixed:5d} burnside_sum={burnside_sum:6d}".format(**row)
            )
        print(f"     total={total}; total/{factorial(n)}={total // factorial(n)}")
    print(
        "   Reading: the [3,3] term at n=6 contributes 40*32 labelled fixed shadows\n"
        "   with no fixed vertex.  It is a rootless/cyclic observer leak, matching\n"
        "   the first P(5)->U(6) defect rather than a deeper node-memory failure."
    )


def print_deletion_and_sectors() -> None:
    print()
    print("3. SIX-CLASS DELETION FIBERS AND ORDERED-PAIR SECTORS")
    six = all_classes(6)
    aggregate = {
        "self_converse": [],
        "chiral": [],
    }
    for mask in six:
        key = "self_converse" if canonical(converse(mask, 6), 6) == mask else "chiral"
        aggregate[key].append(mask)
    for key, masks in aggregate.items():
        deletion_hist = Counter(deletion_profile(mask, 6)[0] for mask in masks)
        rooted_hist = Counter(vertex_orbit_count(mask, 6) for mask in masks)
        c3_values = [cyclic_triangle_count(mask, 6) for mask in masks]
        h_values = [hamiltonian_path_count(mask, 6) for mask in masks]
        print(
            f"   {key:14s}: classes={len(masks):2d} "
            f"rooted_mult_sum={sum(vertex_orbit_count(mask, 6) for mask in masks):3d} "
            f"deletion_parent_count_hist={dict(sorted(deletion_hist.items()))} "
            f"vertex_orbit_hist={dict(sorted(rooted_hist.items()))} "
            f"c3_range=({min(c3_values)},{max(c3_values)}) "
            f"H_range=({min(h_values)},{max(h_values)})"
        )

    for label, internal, cross in (
        ("sector_size_deck", False, False),
        ("sector_size_internal_deck", True, False),
        ("sector_full_cross_deck", True, True),
    ):
        buckets: defaultdict[tuple[tuple[object, ...], ...], list[int]] = defaultdict(list)
        for mask in six:
            buckets[ordered_pair_deck(mask, 6, internal, cross)].append(mask)
        collisions = [sorted(v) for v in buckets.values() if len(v) > 1]
        print(
            f"   {label:26s}: separates={len(buckets):2d}/56 "
            f"collisions={collisions[:3]}"
        )
    print(
        "   Reading: sector sizes and internal sector tournaments still forget one\n"
        "   converse-pair distinction.  Adding cross-sector orientation restores\n"
        "   all 56 classes, so cross-sector orientation is the compact observer-cut\n"
        "   payload after the rooted cache."
    )


def print_layer_flow() -> None:
    print()
    print("4. S217 RECTANGLE/HOURGLASS RESIDUES FOR FIXED-PATH DIAGONAL FLOW")
    print("   Local bridge K_{k,k+1}:")
    print("   k  lines  rank  rectangle_redundancy")
    for k in range(1, 8):
        lines, rank, red = local_bridge(k)
        marker = "  <-- one dozen" if red == 12 else ""
        print(f"   {k:1d}  {lines:5d}  {rank:4d}  {red:20d}{marker}")
    print()
    print("   Full adjacent-layer flow:")
    print("   n  lines  rank  redundancy  local_rectangles  hourglass_residues")
    for n in range(4, 8):
        lines, rank, red, local_rectangles, hourglasses = global_flow(n)
        print(
            f"   {n:1d}  {lines:5d}  {rank:4d}  {red:10d}"
            f"  {local_rectangles:16d}  {hourglasses:18d}"
        )
    print(
        "   Reading: the dozen at K_{4,5} is not extra class mass.  It is a\n"
        "   rectangle-cycle reservoir.  To affect U(6)=56 it must pass through an\n"
        "   overlap/forgetting law, not through literal addition."
    )


def print_lrc_synthesis() -> None:
    print()
    print("5. OBSERVER-EXTENSION / CUT PAYLOAD ABSTRACT")
    rows = [
        (
            "HYP-3057 value-origin ledger",
            "small integer reused across adjacent levels",
            "value_origin_type plus the coordinate it came from",
            "type the number before using it as an invariant",
        ),
        (
            "HYP-3056 observer-cut orbit ledger",
            "visible-fiber automorphism quotient",
            "orbit of boundary slice / incidence word / extended shadow",
            "normalize payloads only modulo visible automorphisms",
        ),
        (
            "HYP-3039 controlled forgetting",
            "coarse quotient fiber",
            "hidden coordinate stage / residual debt",
            "retain, reconstruct, dual-annihilate, descend, or name debt",
        ),
        (
            "HYP-3038 q=23 drop/add square",
            "fixed-margin rectangle",
            "exact-M zeta and endpoint-owner strip",
            "rectangle residue is payload, not count noise",
        ),
        (
            "HYP-3037 residual capacitors",
            "two-plate route cut",
            "first-cut stage / capacitor id",
            "cut payload schedules route separation",
        ),
        (
            "HYP-3041 AP-tail repair",
            "AP-core one-tail observer",
            "q=13 puncture or reciprocal fixed-point clock",
            "coarse owner strip is unsafe without hidden clock",
        ),
        (
            "HYP-3045 endpoint-owner transfer",
            "B18Z6 endpoint shadow",
            "external owner strip and transfer delta",
            "owner names are the extension payload",
        ),
        (
            "HYP-3021/HYP-3022 pair-good decoys",
            "blocker generator quotient",
            "barcode / normal-fan active owner",
            "generator class beats raw decoy count",
        ),
        (
            "HYP-3049 ordered-pair sectors",
            "old-root/new-observer cut",
            "cross-sector orientation word",
            "first compact sidecar after rooted cache",
        ),
        (
            "HYP-3053 diagonal flow",
            "K_{k,k+1} line quotient",
            "rectangle/hourglass cycle residue",
            "redundancy law is the proof carrier",
        ),
    ]
    for name, quotient, payload, rule in rows:
        print(f"   {name}:")
        print(f"     quotient={quotient}")
        print(f"     payload={payload}")
        print(f"     rule={rule}")
    print(
        "   Abstract law: observer-extension creates a cut payload; controlled\n"
        "   forgetting may quotient the observer only when that payload is typed\n"
        "   by origin, reduced to its visible payload orbit, fiber-constant,\n"
        "   reconstructible, killed by a dual/cycle equation, descended\n"
        "   familywise, or emitted as named residual debt."
    )


def print_tournament_analysis() -> None:
    print()
    print("6. TOURNAMENT ANALYSIS")
    carriers = [
        "exact_canonical_audit",
        "burnside_odd_cycle_sidecar",
        "cross_sector_orientation_word",
        "deletion_parent_fiber_profile",
        "rectangle_hourglass_cycle_residue",
        "self_converse_branch_locus",
        "duodecimal_overlap_ledger",
        "rooted_perspective_cache",
        "raw_orbit_count_coincidence",
    ]
    print("   Vertices are proof carriers, not runners or arcs.")
    print("   Pairwise observable: retained LRC predicate payload minus quotient loss.")
    print("   Switch: orient toward the carrier that explains an actual lost coordinate;")
    print("   ties follow the displayed Hamiltonian path.")
    print(f"   path={' > '.join(carriers)}")
    print("   score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}")
    print("   directed_3_cycles=0")
    print("   scc_sizes=[1,1,1,1,1,1,1,1,1]")
    print("   hamiltonian_paths=1")


def print_assumption_challenge() -> None:
    print()
    print("7. ASSUMPTION CHALLENGE")
    print(
        "   Vertex sets considered: runners, arcs, rooted vertices, ordered pairs,\n"
        "   edge sectors, diagonal tiles, K_{k,k+1} line sectors, rectangle cycles,\n"
        "   hourglass cycles, deletion fibers, self-converse branch points, cover\n"
        "   arcs, barcode bars, endpoint owners, residues, Fourier modes, matroid\n"
        "   cocircuits, and proof obligations."
    )
    print(
        "   Chosen vertices: observer-extension proof carriers and their cut payloads."
    )
    print(
        "   Preserved predicate: whether the quotient keeps enough information to\n"
        "   separate or reconstruct LRC route/status-changing packet pairs."
    )
    print(
        "   Destroyed information: labelled runner identity, old-root/new-observer\n"
        "   role, cross-sector orientation, deletion parent multiplicity, endpoint\n"
        "   owner names, and rectangle/hourglass residues when they are not sidecars."
    )
    print(
        "   Challenged assumption: the dozen should not be added as independent mass.\n"
        "   At this scale the live statement is a duodecimal overlap law:\n"
        "   56 = 48 + 12 - 4, with net correction 8 = (2/3)*12."
    )


def main() -> None:
    print("=" * 80)
    print("S230: observer-extension exact duodecimal deletion-sector audit")
    print("=" * 80)
    print_count_ladders()
    print_burnside()
    print_deletion_and_sectors()
    print_layer_flow()
    print_lrc_synthesis()
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
