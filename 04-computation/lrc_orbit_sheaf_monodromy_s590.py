#!/usr/bin/env python3
"""S590: orbit-sheaf monodromy probes for LRC/tournament rigidity.

The guiding slogan is:

    rigidity = an orbit quotient with predicate labels that glue;
    flexibility/defect = monodromy in a forgotten fiber.

The script runs three small exact audits:

1. AP unit-clock witnesses under the doubling map.  For odd n, doubling is an
   automorphism of the unit sheaf.  For even n, it becomes a boundary morphism
   into nonunit residues.
2. Dihedral point-sets in Z/N.  Cyclic vertex-transitivity is primitive
   polygon rigidity; full dihedral transitivity also permits two-block
   bracelets.
3. Rooted tournament quotient fibers through n=6.  A quotient lens is useful
   only when labels such as source, fixed-root, and parent class remain pure on
   its fibers.

Tournament Analysis is included over rigidity lenses.  The pairwise observable
is the defect tuple (boundary loss, projection collisions, max fiber, label
mixing); the switch orients toward lower defect, with name order as a tie path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd


def units(n: int) -> list[int]:
    return [u for u in range(1, n) if gcd(u, n) == 1]


def witness_clock_ap(n: int) -> list[int]:
    """AP {1,...,n-1} lonely clock residues, by the HYP-2124 criterion."""
    return [j for j in range(1, n) if gcd(j, n) == 1]


def v2(n: int) -> int:
    out = 0
    while n and n % 2 == 0:
        out += 1
        n //= 2
    return out


def preimages_of_doubling(y: int, n: int) -> list[int]:
    return [x for x in range(n) if (2 * x) % n == y]


def doubling_sheaf_row(n: int) -> dict[str, object]:
    u = units(n)
    witness = witness_clock_ap(n)
    images = [(a, (2 * a) % n) for a in witness]
    kept = [a for a, b in images if gcd(b, n) == 1]
    boundary = [(a, b, gcd(b, n)) for a, b in images if gcd(b, n) != 1]
    image_hist = Counter(g for _, _, g in boundary)
    collapse = Counter(b for _, b in images)
    return {
        "n": n,
        "parity": "odd" if n % 2 else "even",
        "v2": v2(n),
        "units": len(u),
        "witness_is_units": witness == u,
        "kept_in_unit_sheaf": len(kept),
        "boundary_count": len(boundary),
        "image_size": len(collapse),
        "image_gcd_hist": dict(sorted(image_hist.items())),
        "collapse_hist": dict(sorted(Counter(collapse.values()).items())),
        "status": "automorphism" if len(kept) == len(witness) else "boundary",
        "images": images,
        "boundary": boundary,
    }


def circular_gaps(p: tuple[int, ...], n: int) -> tuple[int, ...]:
    return tuple((p[(i + 1) % len(p)] - p[i]) % n for i in range(len(p)))


def min_period(word: tuple[int, ...]) -> int:
    for d in range(1, len(word) + 1):
        if len(word) % d == 0 and all(word[i] == word[i % d] for i in range(len(word))):
            return d
    raise AssertionError("word always has a period")


def apply_dihedral(x: int, n: int, shift: int, reflect: bool) -> int:
    return (shift - x if reflect else shift + x) % n


def dihedral_stabilizer(p: tuple[int, ...], n: int) -> list[tuple[int, bool]]:
    pset = set(p)
    out = []
    for shift in range(n):
        for reflect in (False, True):
            image = {apply_dihedral(x, n, shift, reflect) for x in p}
            if image == pset:
                out.append((shift, reflect))
    return out


def orbit_under_group(
    start: int, pset: set[int], n: int, group: list[tuple[int, bool]]
) -> set[int]:
    seen = {start}
    frontier = [start]
    while frontier:
        x = frontier.pop()
        for shift, reflect in group:
            y = apply_dihedral(x, n, shift, reflect)
            if y in pset and y not in seen:
                seen.add(y)
                frontier.append(y)
    return seen


def dihedral_pointset_audit(max_n: int = 18) -> dict[str, object]:
    total = regular = bracelet = 0
    gap_period_hist: Counter[int] = Counter()
    first_bracelet = None
    first_period_gt2 = None
    for n in range(3, max_n + 1):
        rest = list(range(1, n))
        for mask in range(1 << (n - 1)):
            p = (0,) + tuple(rest[i] for i in range(n - 1) if (mask >> i) & 1)
            if len(p) < 3:
                continue
            group = dihedral_stabilizer(p, n)
            if not group:
                continue
            pset = set(p)
            if orbit_under_group(0, pset, n, group) != pset:
                continue
            gaps = circular_gaps(p, n)
            period = min_period(gaps)
            is_regular = period == 1
            total += 1
            regular += int(is_regular)
            bracelet += int(not is_regular)
            gap_period_hist[period] += 1
            if not is_regular and first_bracelet is None:
                first_bracelet = {
                    "n": n,
                    "p": p,
                    "gaps": gaps,
                    "period": period,
                    "stabilizer_size": len(group),
                }
            if period > 2 and first_period_gt2 is None:
                first_period_gt2 = {"n": n, "p": p, "gaps": gaps, "period": period}
    return {
        "max_n": max_n,
        "dihedral_vt": total,
        "regular": regular,
        "bracelet": bracelet,
        "gap_period_hist": dict(sorted(gap_period_hist.items())),
        "first_bracelet": first_bracelet,
        "first_period_gt2": first_period_gt2,
    }


PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PERMS: dict[int, tuple[tuple[int, ...], ...]] = {}


def pairs(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = {
            (i, j): k for k, (i, j) in enumerate(combinations(range(n), 2))
        }
    return PAIR_INDEX[n]


def perms(n: int) -> tuple[tuple[int, ...], ...]:
    if n not in PERMS:
        PERMS[n] = tuple(permutations(range(n)))
    return PERMS[n]


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("no loops")
    if i < j:
        return bool((mask >> pairs(n)[(i, j)]) & 1)
    return not edge(mask, n, j, i)


def set_edge(mask: int, n: int, i: int, j: int, i_beats_j: bool) -> int:
    if i > j:
        return set_edge(mask, n, j, i, not i_beats_j)
    bit = 1 << pairs(n)[(i, j)]
    return (mask | bit) if i_beats_j else (mask & ~bit)


def relabel(mask: int, n: int, old_for_new: tuple[int, ...]) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        if edge(mask, n, old_for_new[i], old_for_new[j]):
            out = set_edge(out, n, i, j, True)
    return out


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, p) for p in perms(n))


@lru_cache(maxsize=None)
def rooted_canonical(mask: int, n: int, root: int) -> int:
    others = tuple(v for v in range(n) if v != root)
    return min(relabel(mask, n, (root,) + p) for p in permutations(others))


@lru_cache(maxsize=None)
def automorphisms(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(p for p in perms(n) if relabel(mask, n, p) == mask)


def vertex_orbits(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
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

    for p in automorphisms(mask, n):
        for v in range(n):
            union(v, p[v])
    groups: dict[int, list[int]] = defaultdict(list)
    for v in range(n):
        groups[find(v)].append(v)
    return tuple(tuple(vs) for vs in sorted(groups.values(), key=lambda x: (len(x), x)))


def induced(mask: int, n: int, verts: tuple[int, ...]) -> int:
    out = 0
    for i, j in combinations(range(len(verts)), 2):
        if edge(mask, n, verts[i], verts[j]):
            out = set_edge(out, len(verts), i, j, True)
    return out


def delete_vertex(mask: int, n: int, victim: int) -> int:
    return induced(mask, n, tuple(v for v in range(n) if v != victim))


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum(edge(mask, n, i, j) for j in range(n) if j != i) for i in range(n))
    )


@lru_cache(maxsize=None)
def classes(n: int) -> tuple[int, ...]:
    reps = []
    seen = set()
    for mask in range(1 << (n * (n - 1) // 2)):
        c = canonical(mask, n)
        if c not in seen:
            seen.add(c)
            reps.append(c)
    return tuple(sorted(reps))


def rooted_records(n: int) -> list[dict[str, object]]:
    out = []
    for rep in classes(n):
        for orbit in vertex_orbits(rep, n):
            root = min(orbit)
            out_neighbors = tuple(j for j in range(n) if j != root and edge(rep, n, root, j))
            in_neighbors = tuple(j for j in range(n) if j != root and edge(rep, n, j, root))
            score = len(out_neighbors)
            out_class = canonical(induced(rep, n, out_neighbors), score)
            in_class = canonical(induced(rep, n, in_neighbors), n - 1 - score)
            out.append(
                {
                    "rooted": rooted_canonical(rep, n, root),
                    "unrooted": rep,
                    "score": score,
                    "source": score == n - 1,
                    "sink": score == 0,
                    "fixed": len(orbit) == 1,
                    "orbit_size": len(orbit),
                    "parent": canonical(delete_vertex(rep, n, root), n - 1),
                    "split": (score, out_class, in_class),
                    "score_sequence": score_sequence(rep, n),
                }
            )
    return out


def lens_value(record: dict[str, object], lens: str) -> object:
    if lens == "full_rooted":
        return record["rooted"]
    if lens == "unrooted_plus_score":
        return (record["unrooted"], record["score"])
    if lens == "split_no_cross":
        return record["split"]
    if lens == "unrooted":
        return record["unrooted"]
    if lens == "delete_parent":
        return record["parent"]
    if lens == "score_sequence":
        return record["score_sequence"]
    if lens == "root_score":
        return record["score"]
    raise KeyError(lens)


def quotient_stats(records: list[dict[str, object]], lens: str) -> dict[str, object]:
    fibers: dict[object, list[dict[str, object]]] = defaultdict(list)
    for rec in records:
        fibers[lens_value(rec, lens)].append(rec)
    labels = ("source", "sink", "fixed", "orbit_size", "parent")
    mixed = {}
    for label in labels:
        mixed[label] = sum(
            1 for fiber in fibers.values() if len({rec[label] for rec in fiber}) > 1
        )
    fiber_sizes = Counter(len(fiber) for fiber in fibers.values())
    return {
        "lens": lens,
        "unique": len(fibers),
        "collisions": len(records) - len(fibers),
        "max_fiber": max(fiber_sizes) if fiber_sizes else 0,
        "fiber_hist": dict(sorted(fiber_sizes.items())),
        "mixed": mixed,
    }


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            x = stack.pop()
            for y in graph[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen

    rev = [[] for _ in range(n)]
    for i, row in enumerate(adj):
        for j in row:
            rev[j].append(i)
    unused = set(range(n))
    sizes = []
    while unused:
        x = min(unused)
        comp = reach(x, adj) & reach(x, rev)
        sizes.append(len(comp))
        unused -= comp
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(1, n):
        ndp = dict(dp)
        for (mask, last), count in dp.items():
            if mask.bit_count() != size:
                continue
            for nxt in adj[last]:
                if not (mask >> nxt) & 1:
                    ndp[(mask | (1 << nxt), nxt)] = ndp.get((mask | (1 << nxt), nxt), 0) + count
        dp = ndp
    full = (1 << n) - 1
    return sum(count for (mask, _), count in dp.items() if mask == full)


def tournament_analysis(defects: dict[str, tuple[int, int, int, int]]) -> dict[str, object]:
    names = sorted(defects)
    n = len(names)
    adj = [[] for _ in names]
    ties = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i >= j:
                continue
            da, db = defects[a], defects[b]
            if da == db:
                ties += 1
                winner, loser = (a, b) if a < b else (b, a)
            else:
                winner, loser = (a, b) if da < db else (b, a)
            wi, li = names.index(winner), names.index(loser)
            adj[wi].append(li)
    scores = {names[i]: len(adj[i]) for i in range(n)}
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        edges = {(i, j) for i in (a, b, c) for j in adj[i] if j in (a, b, c)}
        if len(edges) == 3:
            out = Counter(i for i, _ in edges)
            if sorted(out.values()) == [1, 1, 1]:
                c3 += 1
    order = sorted(names, key=lambda x: (defects[x], x))
    return {
        "vertices": names,
        "defects": defects,
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3_cycles": c3,
        "scc_sizes": scc_sizes(adj),
        "hamiltonian_paths": count_hamiltonian_paths(adj),
        "tie_edge_count": ties,
        "one_hamiltonian_order": order,
    }


def main() -> None:
    print("S590 ORBIT-SHEAF MONODROMY RIGIDITY AUDIT")
    print("=" * 48)
    print()
    print("1. AP unit-clock sheaf under doubling")
    for n in range(5, 19):
        row = doubling_sheaf_row(n)
        print(
            f"  n={n:2d} {row['parity']:4s} v2={row['v2']} |U|={row['units']:2d} "
            f"kept={row['kept_in_unit_sheaf']:2d} boundary={row['boundary_count']:2d} "
            f"image_size={row['image_size']:2d} gcd_hist={row['image_gcd_hist']} "
            f"collapse={row['collapse_hist']} status={row['status']}"
        )
    n14 = doubling_sheaf_row(14)
    print("  n=14 boundary map:")
    for a, b, g in n14["boundary"]:
        lifts = [(x, gcd(x, 14)) for x in preimages_of_doubling(b, 14)]
        print(f"    unit {a:2d} -> {b:2d} (gcd={g}); all lifts={lifts}")
    print()

    print("2. Dihedral point-set trienerment sheaf")
    pointset = dihedral_pointset_audit(18)
    print(
        f"  searched N<=18; dihedral_VT={pointset['dihedral_vt']} "
        f"regular={pointset['regular']} bracelet={pointset['bracelet']} "
        f"gap_period_hist={pointset['gap_period_hist']}"
    )
    print(f"  first_bracelet={pointset['first_bracelet']}")
    print(f"  first_period_gt2={pointset['first_period_gt2']}")
    print()

    print("3. Rooted tournament quotient sheaf through n=6")
    records = rooted_records(6)
    lenses = [
        "full_rooted",
        "unrooted_plus_score",
        "split_no_cross",
        "unrooted",
        "delete_parent",
        "score_sequence",
        "root_score",
    ]
    stats = [quotient_stats(records, lens) for lens in lenses]
    for st in stats:
        print(
            f"  {st['lens']:<21s} unique={st['unique']:3d} "
            f"collisions={st['collisions']:3d} max_fiber={st['max_fiber']:3d} "
            f"mixed={st['mixed']}"
        )
    print()

    print("4. Tournament Analysis on rigidity lenses")
    stat_by_lens = {st["lens"]: st for st in stats}
    defects = {
        "static_unit_sheaf": (0, 0, 1, 0),
        "source_root_sheaf": (0, 0, 1, 0),
        "cyclic_polygon_spine": (0, 0, 1, 0),
        "dihedral_bracelet_monodromy": (0, pointset["bracelet"], 2, 1),
        "doubling_even_boundary_n14": (
            n14["boundary_count"],
            n14["boundary_count"],
            max(n14["collapse_hist"]) if n14["collapse_hist"] else 0,
            1,
        ),
        "split_profile_no_cross": (
            0,
            stat_by_lens["split_no_cross"]["collisions"],
            stat_by_lens["split_no_cross"]["max_fiber"],
            sum(stat_by_lens["split_no_cross"]["mixed"].values()),
        ),
        "unmarked_shadow": (
            0,
            stat_by_lens["unrooted"]["collisions"],
            stat_by_lens["unrooted"]["max_fiber"],
            sum(stat_by_lens["unrooted"]["mixed"].values()),
        ),
        "root_score_shadow": (
            0,
            stat_by_lens["root_score"]["collisions"],
            stat_by_lens["root_score"]["max_fiber"],
            sum(stat_by_lens["root_score"]["mixed"].values()),
        ),
    }
    ta = tournament_analysis(defects)
    print(f"  pairwise observable: defect tuple = (boundary_loss, collisions, max_fiber, label_mix)")
    print(f"  score_hist={ta['score_hist']}")
    print(f"  directed_3_cycles={ta['directed_3_cycles']}")
    print(f"  scc_sizes={ta['scc_sizes']}")
    print(f"  tie_edge_count={ta['tie_edge_count']}")
    print(f"  hamiltonian_paths={ta['hamiltonian_paths']}")
    print("  one_hamiltonian_order=" + " -> ".join(ta["one_hamiltonian_order"]))
    print()

    print("READING")
    print(
        "  Odd n: doubling is an automorphism of the AP unit-clock sheaf. "
        "Even n: it is a boundary morphism from units to nonunits."
    )
    print(
        "  At n=14 every unit witness doubles into gcd-2 residue; each image has "
        "one unit lift and one nonunit lift, exposing the 2-block seam."
    )
    print(
        "  Dihedral VT point-sets glue local vertex profiles with either trivial "
        "cyclic monodromy or a two-block bracelet monodromy."
    )
    print(
        "  Rooted tournament quotients show the sheaf warning: source is rigid, "
        "but fixed-root, parent, and orbit-size labels become mixed once cross "
        "fibers are forgotten."
    )


if __name__ == "__main__":
    main()
