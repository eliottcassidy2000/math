#!/usr/bin/env python3
"""S221: value-origin ledger for the 12/48/56 tournament carrier knot.

This scout pulls together the count shadows that appeared across S211, S213,
S216, and S217:

* U(n), A000568 unlabelled tournament classes.
* R(n), rooted/node perspectives, i.e. sum of vertex orbits over U(n).
* SC(n), self-converse classes.
* incident-word extension orbits, ordered-pair/edge sectors, deletion fibers.
* fixed-path rectangle/hourglass redundancy dimensions.

The point is not speed.  It is an exact small-n ledger for deciding which
numbers are class counts, which are rooted fibers, and which are presentation
duplication laws.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

import tournament_rigidity_cascade_s589 as tour


KNOWN_U = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}


def converse(mask: int, n: int) -> int:
    return mask ^ ((1 << (n * (n - 1) // 2)) - 1)


def is_self_converse(mask: int, n: int) -> bool:
    return tour.canonical(converse(mask, n), n) == mask


def rooted_count(n: int) -> int:
    return sum(len(tour.vertex_orbits(rep, n)) for rep in tour.classes(n))


def self_converse_count(n: int) -> int:
    return sum(1 for rep in tour.classes(n) if is_self_converse(rep, n))


def lift_mask(mask: int, n: int) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        if tour.edge(mask, n, i, j):
            out = tour.set_edge(out, n + 1, i, j, True)
    return out


def extend_by_word(mask: int, n: int, word: int) -> int:
    """Append vertex n.  Bit i means old vertex i beats the new vertex."""
    out = lift_mask(mask, n)
    for i in range(n):
        out = tour.set_edge(out, n + 1, i, n, bool((word >> i) & 1))
    return out


def act_on_word(word: tuple[int, ...], aut: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(word[aut[i]] for i in range(len(word)))


def word_orbit_reps(parent: int, n: int) -> list[tuple[int, int]]:
    """Return (representative word-as-int, orbit size) under Aut(parent)."""
    auts = tour.automorphisms(parent, n)
    unseen = {tuple((word >> i) & 1 for i in range(n)) for word in range(1 << n)}
    out: list[tuple[int, int]] = []
    while unseen:
        word = min(unseen)
        orbit = {act_on_word(word, aut) for aut in auts}
        rep = min(orbit)
        rep_int = sum(bit << i for i, bit in enumerate(rep))
        out.append((rep_int, len(orbit)))
        unseen -= orbit
    return out


def ordered_pair_canonical(mask: int, n: int, first: int, second: int) -> int:
    others = tuple(v for v in range(n) if v not in (first, second))
    return min(tour.relabel(mask, n, (first, second) + p) for p in permutations(others))


def unordered_pair_canonical(mask: int, n: int, a: int, b: int) -> int:
    return min(ordered_pair_canonical(mask, n, a, b), ordered_pair_canonical(mask, n, b, a))


def oriented_edges(mask: int, n: int) -> list[tuple[int, int]]:
    out = []
    for a, b in combinations(range(n), 2):
        out.append((a, b) if tour.edge(mask, n, a, b) else (b, a))
    return out


def rooted_parent_incident_state(parent: int, root: int, word: int) -> int:
    child = extend_by_word(parent, 5, word)
    old_vertices = tuple(v for v in range(5) if v != root)
    return min(tour.relabel(child, 6, (root, 5) + p) for p in permutations(old_vertices))


SECTORS = ((0, 0), (0, 1), (1, 0), (1, 1))


def induced_canonical(mask: int, n: int, vertices: list[int]) -> int:
    if len(vertices) <= 1:
        return 0
    return tour.canonical(tour.induced(mask, n, tuple(vertices)), len(vertices))


def sector_signature(mask: int, n: int, first: int, second: int, mode: str) -> tuple[object, ...]:
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for x in range(n):
        if x in (first, second):
            continue
        key = (int(tour.edge(mask, n, first, x)), int(tour.edge(mask, n, second, x)))
        groups[key].append(x)

    sizes = tuple((sector, len(groups[sector])) for sector in SECTORS)
    internal = tuple((sector, induced_canonical(mask, n, groups[sector])) for sector in SECTORS)
    cross = []
    for i, a in enumerate(SECTORS):
        for b in SECTORS[i + 1 :]:
            wins = sum(1 for u in groups[a] for v in groups[b] if tour.edge(mask, n, u, v))
            cross.append((a, b, wins, len(groups[a]) * len(groups[b])))
    cross_tuple = tuple(cross)

    if mode == "size":
        return sizes
    if mode == "internal":
        return sizes, internal
    if mode == "cross":
        return sizes, cross_tuple
    if mode == "full":
        return sizes, internal, cross_tuple
    raise ValueError(mode)


def class_sector_deck(mask: int, n: int, mode: str) -> tuple[tuple[object, ...], ...]:
    return tuple(
        sorted(
            sector_signature(mask, n, a, b, mode)
            for a in range(n)
            for b in range(n)
            if a != b
        )
    )


def deletion_profile(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(tour.canonical(tour.delete_vertex(mask, n, v), n - 1) for v in range(n)))


def fixed_path_tournament(n: int, tile_bits: int) -> int:
    mask = 0
    bit = 0
    for i, j in combinations(range(n), 2):
        if j == i + 1:
            value = True
        else:
            value = bool((tile_bits >> bit) & 1)
            bit += 1
        if value:
            mask = tour.set_edge(mask, n, i, j, True)
    return mask


def fixed_path_fibers(n: int) -> Counter[int]:
    fibers: Counter[int] = Counter()
    for bits in range(1 << comb(n - 1, 2)):
        fibers[tour.canonical(fixed_path_tournament(n, bits), n)] += 1
    return fibers


def hamiltonian_paths(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if tour.edge(mask, n, v, u):
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def split_hist_by_self_converse(values: dict[int, int], n: int) -> tuple[dict[int, int], dict[int, int]]:
    sc: Counter[int] = Counter()
    non: Counter[int] = Counter()
    for rep, value in values.items():
        (sc if is_self_converse(rep, n) else non)[value] += 1
    return dict(sorted(sc.items())), dict(sorted(non.items()))


def local_rectangle_dim(k: int) -> int:
    return k * (k - 1)


def fixed_path_flow_row(n: int) -> dict[str, int]:
    tiles = comb(n - 1, 2)
    lines = 2 * comb(n - 1, 3)
    rank = max(0, tiles - 1)
    red = lines - rank
    local = 2 * comb(n - 2, 3) if n >= 5 else 0
    hourglass = comb(n - 3, 2) if n >= 5 else 0
    return {
        "n": n,
        "tiles": tiles,
        "lines": lines,
        "rank": rank,
        "red": red,
        "local_rectangles": local,
        "hourglass": hourglass,
    }


def print_shift_table() -> None:
    print("1. SHIFT / SELF-CONVERSE / FRACTION TABLE")
    print("   m  U(m)  R(m)  U(m+1)  U-R  R/U(m+1)  SC(m)  SC(m+1)")
    for m in range(1, 7):
        u = KNOWN_U[m]
        r = rooted_count(m)
        next_u = KNOWN_U[m + 1]
        gap = next_u - r
        frac = Fraction(r, next_u)
        next_sc = self_converse_count(m + 1) if m + 1 <= 6 else None
        next_sc_text = f"{next_sc:7d}" if next_sc is not None else "      -"
        print(
            f"   {m}  {u:4d}  {r:4d}  {next_u:6d}  {gap:3d}"
            f"  {frac.numerator}/{frac.denominator:<7}  {self_converse_count(m):5d}"
            f"  {next_sc_text}"
        )
    print()
    print("   Arithmetic correction: the first shifted failure is 48 + 8 = 56, not 48 + 12 = 56.")
    print("   The 8 is also SC(5) at this small boundary; that equality does not persist one step later.")
    print("   The 12-thread is diagonal: R(4)=12, U(5)=12, and SC(6)=12.")


def print_twelve_ledger() -> None:
    classes5 = tour.classes(5)
    classes6 = tour.classes(6)
    source_children = {
        tour.canonical(extend_by_word(parent, 5, 0), 6)
        for parent in classes5
    }
    sink_children = {
        tour.canonical(extend_by_word(parent, 5, (1 << 5) - 1), 6)
        for parent in classes5
    }
    print()
    print("2. WHERE THE VALUE 12 APPEARS")
    rows = [
        ("R(4)", rooted_count(4), "unique node perspectives on 4 vertices"),
        ("U(5)", KNOWN_U[5], "unlabelled 5-tournament classes"),
        ("SC(6)", self_converse_count(6), "self-converse 6-tournament classes"),
        ("parents 5->6", len(classes5), "parent classes for the first failing extension"),
        ("source children", len(source_children), "all-zero incident words, exact source deletion slice"),
        ("sink children", len(sink_children), "all-one incident words, exact sink deletion slice"),
        ("triple carriers at n=4", 12, "exact triple perspectives from S214"),
    ]
    for name, value, meaning in rows:
        print(f"   {name:22s} = {value:2d}  {meaning}")
    print(f"   U(6) classes checked = {len(classes6)}")


def print_first_failure_ledger() -> None:
    classes5 = tour.classes(5)
    classes6 = tour.classes(6)
    rooted5 = {
        tour.rooted_canonical(rep, 5, root)
        for rep in classes5
        for root in range(5)
    }

    child_sinks: set[int] = set()
    word_orbit_total = 0
    orbit_size_hist: Counter[int] = Counter()
    for parent in classes5:
        for rep_word, orbit_size in word_orbit_reps(parent, 5):
            word_orbit_total += 1
            orbit_size_hist[orbit_size] += 1
            child_sinks.add(tour.canonical(extend_by_word(parent, 5, rep_word), 6))

    extension_states: set[int] = set()
    states_by_root: dict[int, set[int]] = defaultdict(set)
    targets_by_root: dict[int, set[int]] = defaultdict(set)
    roots_by_target: dict[int, set[int]] = defaultdict(set)
    for parent in classes5:
        for orbit in tour.vertex_orbits(parent, 5):
            root = orbit[0]
            root_key = tour.rooted_canonical(parent, 5, root)
            for word in range(1 << 5):
                state = rooted_parent_incident_state(parent, root, word)
                target = tour.canonical(extend_by_word(parent, 5, word), 6)
                extension_states.add(state)
                states_by_root[root_key].add(state)
                targets_by_root[root_key].add(target)
                roots_by_target[target].add(root_key)

    ordered_pairs = {
        ordered_pair_canonical(rep, 6, a, b)
        for rep in classes6
        for a in range(6)
        for b in range(6)
        if a != b
    }
    directed_edges = {
        ordered_pair_canonical(rep, 6, tail, tip)
        for rep in classes6
        for tail, tip in oriented_edges(rep, 6)
    }
    unordered_pairs = {
        unordered_pair_canonical(rep, 6, a, b)
        for rep in classes6
        for a, b in combinations(range(6), 2)
    }

    print()
    print("3. FIRST FAILURE LEDGER, 5 -> 6")
    print(f"   U(5)={len(classes5)}  R(5)={len(rooted5)}  U(6)={len(classes6)}  defect={len(classes6)-len(rooted5)}")
    print(f"   raw incident words: U(5)*2^5 = {len(classes5) * 32}")
    print(f"   parent-Aut word orbits = rooted R(6) = {word_orbit_total}")
    print(f"   word-orbit size histogram = {dict(sorted(orbit_size_hist.items()))}")
    print(f"   unrooted child sinks reached by word orbits = {len(child_sinks)}")
    print(f"   rooted 5-perspective + incident-word states = {len(extension_states)}")
    print(f"   ordered-pair perspectives on U(6) = {len(ordered_pairs)}")
    print(f"   directed-edge perspectives = {len(directed_edges)}")
    print(f"   unordered-pair perspectives = {len(unordered_pairs)}")
    print(
        "   extension-state histogram per rooted 5-parent = "
        f"{dict(sorted(Counter(len(v) for v in states_by_root.values()).items()))}"
    )
    print(
        "   target classes per rooted 5-parent = "
        f"{dict(sorted(Counter(len(v) for v in targets_by_root.values()).items()))}"
    )
    print(
        "   rooted 5-parents per U(6) target = "
        f"{dict(sorted(Counter(len(v) for v in roots_by_target.values()).items()))}"
    )


def print_edge_sector_ledger() -> None:
    classes6 = tour.classes(6)
    print()
    print("4. ORDERED-PAIR / EDGE-SECTOR SEPARATION")
    for mode in ("size", "internal", "cross", "full"):
        individual = {
            sector_signature(rep, 6, a, b, mode)
            for rep in classes6
            for a in range(6)
            for b in range(6)
            if a != b
        }
        deck_to_classes: dict[tuple[tuple[object, ...], ...], list[int]] = defaultdict(list)
        for rep in classes6:
            deck_to_classes[class_sector_deck(rep, 6, mode)].append(rep)
        unique = len(deck_to_classes)
        collision_sizes = sorted(len(v) for v in deck_to_classes.values() if len(v) > 1)
        print(
            f"   {mode:8s}: individual_sigs={len(individual):3d} "
            f"class_decks={unique:2d}/56 collision_sizes={collision_sizes}"
        )

    internal_decks: dict[tuple[tuple[object, ...], ...], list[int]] = defaultdict(list)
    for rep in classes6:
        internal_decks[class_sector_deck(rep, 6, "internal")].append(rep)
    collision = [v for v in internal_decks.values() if len(v) > 1][0]
    print(f"   only size/internal collision masks = {collision}")
    print(f"   collision is converse pair = {tour.canonical(converse(collision[0], 6), 6) == collision[1]}")
    print("   cross-sector orientation is the first compact repair column.")


def print_deletion_and_path_fibers() -> None:
    classes6 = tour.classes(6)
    vertex_orbit_values = {rep: len(tour.vertex_orbits(rep, 6)) for rep in classes6}
    unique_parent_values = {rep: len(set(deletion_profile(rep, 6))) for rep in classes6}
    max_repeat_values = {
        rep: max(Counter(deletion_profile(rep, 6)).values()) for rep in classes6
    }
    fixed_fibers = fixed_path_fibers(6)

    print()
    print("5. DELETION FIBERS, SELF-CONVERSE BRANCH, FIXED-PATH FIBERS AT n=6")
    for label, values in (
        ("vertex/root orbit count", vertex_orbit_values),
        ("unique deletion parents", unique_parent_values),
        ("max repeated deletion parent", max_repeat_values),
        ("fixed-path fiber H(T)/Aut(T)", fixed_fibers),
    ):
        all_hist = dict(sorted(Counter(values.values()).items()))
        sc_hist, non_hist = split_hist_by_self_converse(values, 6)
        print(f"   {label}:")
        print(f"      all = {all_hist}")
        print(f"      self_converse = {sc_hist}")
        print(f"      non_self_converse = {non_hist}")


def print_rectangle_hourglass_sequences() -> None:
    print()
    print("6. FIXED-PATH RECTANGLE / HOURGLASS RESIDUE SEQUENCE")
    print("   n  free_tiles  lines  rank  red  local_rectangles  hourglass  red=sum")
    for n in range(3, 9):
        row = fixed_path_flow_row(n)
        print(
            f"   {n}  {row['tiles']:10d}  {row['lines']:5d}  {row['rank']:4d}"
            f"  {row['red']:3d}  {row['local_rectangles']:16d}"
            f"  {row['hourglass']:9d}  {row['local_rectangles'] + row['hourglass']:7d}"
        )
    print()
    print("   recurrence:")
    print("      free_tiles(n)=C(n-1,2)")
    print("      fixed_lines(n)=2*C(n-1,3)")
    print("      fixed_rank(n)=C(n-1,2)-1")
    print("      fixed_red(n)=2*C(n-2,3)+C(n-3,2)")
    print("   rectangle residues are local 4-cycles; hourglass residues link adjacent bridges.")


def print_reading() -> None:
    print()
    print("READING")
    print(
        "  The value 12 is a three-level alignment: R(4)=U(5)=SC(6)=12.  It is "
        "not the additive correction at the first shifted failure; that correction "
        "is 8, with R(5)=48 and U(6)=56."
    )
    print(
        "  The clean recurrence object is not U(n) alone.  It is the transport "
        "ladder parent class + incident word orbit -> rooted child -> unrooted "
        "sink, with deletion fibers and edge-sector orientation as sidecars."
    )
    print(
        "  Cross-sector orientation repairs the only size/internal edge-sector "
        "collision at U(6), a converse pair.  S217's rectangle/hourglass residues "
        "are the analogous repair coordinates for fixed-path diagonal flow."
    )
    print(
        "  For LRC, the useful moral is: when a quotient produces a magical small "
        "integer, ask whether it is a class count, a rooted-perspective count, a "
        "self-converse fixed branch, an incident-word orbit count, a deletion-fiber "
        "multiplicity, or a cycle-space residue."
    )


def main() -> None:
    print("=" * 80)
    print("S221: tournament value-origin ledger for 12 / 48 / 56")
    print("=" * 80)
    print_shift_table()
    print_twelve_ledger()
    print_first_failure_ledger()
    print_edge_sector_ledger()
    print_deletion_and_path_fibers()
    print_rectangle_hourglass_sequences()
    print_reading()


if __name__ == "__main__":
    main()
