#!/usr/bin/env python3
"""S629: SC perspective flips and cyclotomic unit-distance carrier checks.

The script has two deliberately separate halves.

1.  Exact tournament enumeration through n=6.  For self-converse tournaments,
    compute automorphisms, anti-automorphisms, rooted vertex perspectives, and
    the induced edge-flip action on those perspectives.
2.  A unit-distance carrier table that keeps the H=7/H=21 guardrail distinct
    from unit-edge counts and from the S628 Moser spine-ladder rows.

The main theorem audited here is not that every anti-automorphism is an
involution on vertices.  That can fail once sigma^2 is a nontrivial
automorphism.  The stable object is the induced involution on Aut(T)-orbits of
vertices, i.e. on rooted perspectives.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations
from math import floor, isqrt, sqrt


PERMS: dict[int, tuple[tuple[int, ...], ...]] = {}
PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}

KNOWN_UNROOTED = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
KNOWN_SELF_CONVERSE = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12}

U_EXACT = {
    3: 3,
    4: 5,
    5: 7,
    6: 9,
    7: 12,
    8: 14,
    9: 18,
    10: 20,
    11: 23,
    12: 27,
    13: 30,
    14: 33,
    15: 37,
    16: 41,
    17: 44,
    18: 48,
    19: 51,
    20: 54,
    21: 57,
}

# S599z literal unit-tournament values in the Harborth/triangular spine gauge.
UNIT_TOURNAMENT_H = {
    3: 1,
    4: 5,
    5: 15,
    6: 43,
    7: 141,
    8: 513,
    9: 1605,
    10: 4915,
}


def pair_index(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = {
            pair: index for index, pair in enumerate(combinations(range(n), 2))
        }
    return PAIR_INDEX[n]


def perms(n: int) -> tuple[tuple[int, ...], ...]:
    if n not in PERMS:
        PERMS[n] = tuple(permutations(range(n)))
    return PERMS[n]


def edge(mask: int, n: int, i: int, j: int) -> bool:
    """Return True iff i beats j."""
    if i == j:
        raise ValueError("no loops")
    if i < j:
        return bool((mask >> pair_index(n)[(i, j)]) & 1)
    return not edge(mask, n, j, i)


def set_edge(mask: int, n: int, i: int, j: int, i_beats_j: bool) -> int:
    if i > j:
        return set_edge(mask, n, j, i, not i_beats_j)
    bit = 1 << pair_index(n)[(i, j)]
    return mask | bit if i_beats_j else mask & ~bit


def relabel(mask: int, n: int, old_for_new: tuple[int, ...]) -> int:
    """Relabel by new vertex i = old vertex old_for_new[i]."""
    out = 0
    for i, j in combinations(range(n), 2):
        out = set_edge(out, n, i, j, edge(mask, n, old_for_new[i], old_for_new[j]))
    return out


def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, p) for p in perms(n))


def converse(mask: int, n: int) -> int:
    """Reverse every arc."""
    m = n * (n - 1) // 2
    return ((1 << m) - 1) ^ mask


def classes(n: int) -> list[int]:
    """Tournament isomorphism representatives.

    By Redei's theorem every tournament has a Hamiltonian path.  We may
    therefore fix the labelled path 0->1->...->n-1 and enumerate only the
    nonconsecutive edges.  This reaches every unlabelled class and is much
    smaller than all labelled tournaments at n=6.
    """
    seen: set[int] = set()
    out: list[int] = []
    base = 0
    for i in range(n - 1):
        base = set_edge(base, n, i, i + 1, True)
    free_pairs = [(i, j) for i, j in combinations(range(n), 2) if j != i + 1]
    for word in range(1 << len(free_pairs)):
        mask = base
        for bit_index, (i, j) in enumerate(free_pairs):
            mask = set_edge(mask, n, i, j, bool((word >> bit_index) & 1))
        c = canonical(mask, n)
        if c not in seen:
            seen.add(c)
            out.append(c)
    return sorted(out)


def automorphisms(mask: int, n: int) -> list[tuple[int, ...]]:
    return [p for p in perms(n) if relabel(mask, n, p) == mask]


def anti_automorphisms(mask: int, n: int) -> list[tuple[int, ...]]:
    """Permutations sigma with sigma(T)=T^op."""
    return [p for p in perms(n) if relabel(mask, n, p) == converse(mask, n)]


def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """Composition in old_for_new convention: relabel(relabel(T,q),p)."""
    return tuple(q[p[i]] for i in range(len(p)))


def cycle_type(p: tuple[int, ...]) -> tuple[int, ...]:
    seen = [False] * len(p)
    sizes = []
    for i in range(len(p)):
        if seen[i]:
            continue
        cur = i
        size = 0
        while not seen[cur]:
            seen[cur] = True
            size += 1
            cur = p[cur]
        sizes.append(size)
    return tuple(sorted(sizes))


def vertex_orbits_from_auts(auts: list[tuple[int, ...]], n: int) -> tuple[tuple[int, ...], ...]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in auts:
        for v in range(n):
            union(v, p[v])

    grouped: dict[int, list[int]] = defaultdict(list)
    for v in range(n):
        grouped[find(v)].append(v)
    return tuple(tuple(vs) for vs in sorted(grouped.values(), key=lambda xs: (len(xs), xs)))


def perspective_flip(
    orbits: tuple[tuple[int, ...], ...], sigma: tuple[int, ...]
) -> tuple[int, ...]:
    orbit_of = {}
    for idx, orbit in enumerate(orbits):
        for v in orbit:
            orbit_of[v] = idx
    return tuple(orbit_of[sigma[min(orbit)]] for orbit in orbits)


def is_involution(mapping: tuple[int, ...]) -> bool:
    return all(mapping[mapping[i]] == i for i in range(len(mapping)))


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum(edge(mask, n, i, j) for j in range(n) if i != j) for i in range(n))
    )


def ham_paths(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            count = dp[seen][v]
            if not count:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if edge(mask, n, v, u):
                    dp[seen | (1 << u)][u] += count
    return sum(dp[(1 << n) - 1])


def directed_three_cycles(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        if edge(mask, n, a, b) and edge(mask, n, b, c) and edge(mask, n, c, a):
            total += 1
        elif edge(mask, n, a, c) and edge(mask, n, c, b) and edge(mask, n, b, a):
            total += 1
    return total


def scc_sizes(mask: int, n: int) -> list[int]:
    adj = [[i != j and edge(mask, n, i, j) for j in range(n)] for i in range(n)]
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    q.append(u)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        v = min(remaining)
        comp = reach(v, adj) & reach(v, rev)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def flip_cycle_summary(mapping: tuple[int, ...]) -> str:
    if not mapping:
        return "()"
    seen = [False] * len(mapping)
    cycles: list[tuple[int, ...]] = []
    for i in range(len(mapping)):
        if seen[i]:
            continue
        cur = i
        cyc = []
        while not seen[cur]:
            seen[cur] = True
            cyc.append(cur)
            cur = mapping[cur]
        cycles.append(tuple(cyc))
    return " ".join(str(c) for c in cycles)


def phi3(x: int) -> int:
    return x * x + x + 1


def harborth_lattice(n: int) -> int:
    return floor(3 * n - sqrt(12 * n - 3))


def moser_plus_edges(m: int) -> int:
    return 8 if m == 0 else 27 * m + 6


def moser_minus_edges(m: int) -> int:
    return 6 if m == 0 else 27 * m + 3


def tournament_route_analysis() -> dict[str, object]:
    """Tournament Analysis over proof routes in this session."""
    routes = [
        {
            "name": "SC anti-coset theorem",
            "traceability": 5,
            "hgap": 4,
            "computability": 5,
            "side_channels": 4,
            "risk": 1,
        },
        {
            "name": "rooted perspective flip atlas",
            "traceability": 4,
            "hgap": 3,
            "computability": 5,
            "side_channels": 5,
            "risk": 1,
        },
        {
            "name": "unit-spine carrier section",
            "traceability": 5,
            "hgap": 3,
            "computability": 4,
            "side_channels": 5,
            "risk": 1,
        },
        {
            "name": "raw H=7/H=21 scalar match",
            "traceability": 1,
            "hgap": 2,
            "computability": 5,
            "side_channels": 0,
            "risk": 5,
        },
        {
            "name": "Eisenstein edge-count echo",
            "traceability": 3,
            "hgap": 3,
            "computability": 4,
            "side_channels": 2,
            "risk": 3,
        },
        {
            "name": "n=21 Moser slab recursion",
            "traceability": 5,
            "hgap": 3,
            "computability": 4,
            "side_channels": 4,
            "risk": 2,
        },
    ]

    def strength(route: dict[str, object]) -> int:
        return (
            3 * int(route["traceability"])
            + 3 * int(route["hgap"])
            + 2 * int(route["side_channels"])
            + int(route["computability"])
            - 3 * int(route["risk"])
        )

    strengths = [strength(route) for route in routes]
    n = len(routes)
    mask = 0
    for i, j in combinations(range(n), 2):
        if strengths[i] > strengths[j]:
            mask = set_edge(mask, n, i, j, True)
        elif strengths[j] > strengths[i]:
            mask = set_edge(mask, n, j, i, True)
        else:
            # Tie Hamiltonian path: keep the richer side-channel object.
            i_wins = int(routes[i]["side_channels"]) >= int(routes[j]["side_channels"])
            mask = set_edge(mask, n, i, j, i_wins)

    return {
        "routes": [r["name"] for r in routes],
        "strengths": dict(zip((r["name"] for r in routes), strengths)),
        "score_histogram": dict(
            sorted(
                Counter(
                    sum(edge(mask, n, i, j) for j in range(n) if i != j)
                    for i in range(n)
                ).items()
            )
        ),
        "directed_3_cycles": directed_three_cycles(mask, n),
        "scc_sizes": scc_sizes(mask, n),
        "hamiltonian_paths": ham_paths(mask, n),
        "ranking": sorted(zip(strengths, [r["name"] for r in routes]), reverse=True),
    }


def print_sc_perspective_section() -> None:
    print("S629 SC perspective flip atlas")
    print("=" * 72)
    print("Exact enumeration through n=6; all checks use reverse-all-arcs converse.")
    print()

    print("Unrooted, rooted, and self-converse census")
    print("n | U(n) | rooted P(n) | rooted orbit distribution | SC classes | SC rooted")
    print("-- | --: | --: | -- | --: | --:")
    examples: dict[int, list[dict[str, object]]] = {}
    sc_flip_hist: dict[int, Counter[tuple[int, int]]] = {}
    anti_cycle_hist: dict[int, Counter[tuple[int, ...]]] = {}
    theorem_failures: list[tuple[int, int, str]] = []

    for n in range(1, 7):
        reps = classes(n)
        if len(reps) != KNOWN_UNROOTED[n]:
            raise AssertionError(("U", n, len(reps), KNOWN_UNROOTED[n]))
        rooted_total = 0
        rooted_dist: Counter[int] = Counter()
        sc_count = 0
        sc_rooted = 0
        sc_flip_hist[n] = Counter()
        anti_cycle_hist[n] = Counter()
        examples[n] = []

        for idx, rep in enumerate(reps):
            auts = automorphisms(rep, n)
            antis = anti_automorphisms(rep, n)
            orbits = vertex_orbits_from_auts(auts, n)
            rooted_total += len(orbits)
            rooted_dist[len(orbits)] += 1
            if antis:
                sc_count += 1
                sc_rooted += len(orbits)
                if len(antis) != len(auts):
                    theorem_failures.append((n, idx, "|Anti| != |Aut|"))
                flips = {perspective_flip(orbits, sigma) for sigma in antis}
                if len(flips) != 1:
                    theorem_failures.append((n, idx, "anti maps disagree on perspectives"))
                flip = next(iter(flips))
                if not is_involution(flip):
                    theorem_failures.append((n, idx, "perspective flip not involutive"))
                fixed = sum(1 for i, j in enumerate(flip) if i == j)
                transpositions = (len(flip) - fixed) // 2
                sc_flip_hist[n][(fixed, transpositions)] += 1
                anti_cycle_hist[n].update(cycle_type(sigma) for sigma in antis)

            if n in (3, 4) or (antis and n <= 5 and len(examples[n]) < 6):
                examples[n].append(
                    {
                        "idx": idx,
                        "H": ham_paths(rep, n),
                        "scores": score_sequence(rep, n),
                        "orbits": orbits,
                        "is_sc": bool(antis),
                        "anti_count": len(antis),
                        "flip": None
                        if not antis
                        else perspective_flip(orbits, antis[0]),
                        "anti_cycle_types": sorted(set(cycle_type(a) for a in antis)),
                    }
                )

        if sc_count != KNOWN_SELF_CONVERSE[n]:
            raise AssertionError(("SC", n, sc_count, KNOWN_SELF_CONVERSE[n]))

        print(
            f"{n} | {len(reps)} | {rooted_total} | {dict(sorted(rooted_dist.items()))} "
            f"| {sc_count} | {sc_rooted}"
        )

    if theorem_failures:
        raise AssertionError(theorem_failures[:5])

    print()
    print("User's small perspective counts, recomputed")
    for n in (3, 4):
        counts = [len(vertex_orbits_from_auts(automorphisms(rep, n), n)) for rep in classes(n)]
        print(f"  n={n}: perspective counts {counts}, total={sum(counts)}, U({n+1})={KNOWN_UNROOTED[n+1]}")

    print()
    print("Class details at n=3 and n=4")
    print("n | idx | H | scores | SC | Aut-orbits | anti cycle types | perspective flip")
    print("-- | --: | --: | -- | -- | -- | -- | --")
    for n in (3, 4):
        for row in examples[n]:
            flip = row["flip"]
            flip_text = "-" if flip is None else flip_cycle_summary(flip)
            print(
                f"{n} | {row['idx']} | {row['H']} | {row['scores']} | "
                f"{'Y' if row['is_sc'] else 'N'} | {row['orbits']} | "
                f"{row['anti_cycle_types']} | {flip_text}"
            )

    print()
    print("SC flip summary through n=6")
    print("n | anti permutation cycle-type histogram | perspective fixed/transposed histogram")
    print("-- | -- | --")
    for n in range(1, 7):
        print(f"{n} | {dict(sorted(anti_cycle_hist[n].items()))} | {dict(sorted(sc_flip_hist[n].items()))}")

    print()
    print("Theorem audit")
    print("  For every self-converse class through n=6:")
    print("  * |Anti(T)|=|Aut(T)|.")
    print("  * every anti-automorphism induces the same map on Aut(T)-vertex orbits.")
    print("  * that rooted-perspective map is an involution.")
    print("  Individual vertex cycles are not canonical; the perspective involution is.")


def print_unit_distance_section() -> None:
    print()
    print("Cyclotomic H-gap and unit-distance carrier table")
    print("=" * 72)
    print(f"Phi_3(2)={phi3(2)}, Phi_3(4)={phi3(4)}, 3*Phi_3(2)={3 * phi3(2)}")
    print("The repo convention treats 7 and 21 as forbidden H values, not forbidden unit-edge counts.")
    print()
    print("n | exact u(n) | Harborth/Eisenstein echo | spine n-1 | bulk exact | literal unit-tournament H | scalar note")
    print("-- | --: | --: | --: | --: | --: | --")
    for n in range(3, 15):
        exact = U_EXACT[n]
        lattice = harborth_lattice(n)
        spine = n - 1
        bulk = exact - spine
        h_val = UNIT_TOURNAMENT_H.get(n)
        notes = []
        if exact in (7, 21):
            notes.append(f"exact edge scalar {exact}")
        if lattice in (7, 21):
            notes.append(f"lattice echo {lattice}")
        if h_val in (7, 21):
            notes.append(f"literal H {h_val}")
        if not notes:
            notes.append("no forbidden-H scalar")
        print(
            f"{n} | {exact} | {lattice} | {spine} | {bulk} | "
            f"{'-' if h_val is None else h_val} | {', '.join(notes)}"
        )

    print()
    print("S628 Moser spine-ladder rows")
    print("family | m | vertices | unit edges | spine | bulk | recursion note")
    print("-- | --: | --: | --: | --: | --: | --")
    rows = [
        ("P_m^+", 0, 8, "cap"),
        ("P_m^-", 0, 6, "cap"),
        ("P_m^+", 1, moser_plus_edges(1), "S626 n=14"),
        ("P_m^-", 1, moser_minus_edges(1), "n=13 companion"),
        ("P_m^-", 2, moser_minus_edges(2), "S626 exact n=21"),
        ("P_m^+", 2, moser_plus_edges(2), "S626 n=22 lane"),
        ("P_m^-", 3, moser_minus_edges(3), "next minus slab"),
        ("P_m^+", 3, moser_plus_edges(3), "next plus slab"),
    ]
    for family, m, edges, note in rows:
        vertices = 8 * m + (6 if family == "P_m^+" else 5)
        spine = vertices - 1
        print(
            f"{family} | {m} | {vertices} | {edges} | {spine} | "
            f"{edges - spine} | {note}"
        )
    print("  Each added full slab contributes +8 vertices, +27 unit edges,")
    print("  i.e. +8 spine steps and +19 bulk edges.")
    print("  Thus n=21 in S628 is P_2^- with 57 unit edges, not H=21.")


def print_route_analysis() -> None:
    print()
    print("Tournament Analysis over S629 proof routes")
    print("=" * 72)
    ta = tournament_route_analysis()
    print(f"route strengths={ta['strengths']}")
    print(f"ranking={ta['ranking']}")
    print(f"score_histogram={ta['score_histogram']}")
    print(f"directed_3_cycles={ta['directed_3_cycles']}")
    print(f"scc_sizes={ta['scc_sizes']}")
    print(f"hamiltonian_paths={ta['hamiltonian_paths']}")
    print()
    print("Assumption challenge")
    print("  Vertices used above were proof routes and carrier obligations, not points.")
    print("  Alternate vertices considered: tournament vertices, rooted perspectives,")
    print("  anti-automorphism coset elements, unit spines, unit-edge packets,")
    print("  Eisenstein shells, Moser slabs, and forbidden-H proof obligations.")
    print("  The chosen quotient preserves the SC flip payload and the H-gap guardrail;")
    print("  it destroys continuous embeddings and raw point coordinates.")


def main() -> None:
    print_sc_perspective_section()
    print_unit_distance_section()
    print_route_analysis()


if __name__ == "__main__":
    main()
