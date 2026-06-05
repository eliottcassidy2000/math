#!/usr/bin/env python3
"""S672: Paris-Harrington colorings as finite ultrafilter recursion.

This is not an attempt to compute the real Paris-Harrington growth function.
It builds small exact and sampled shadows of the user's picture:

  colorings = side choices on tuple atoms
  bad colorings = an outer-extension tree
  relative largeness = an initial-segment gate that blocks tail escape
  least forced stage = rank at which every coherent bad branch dies

The exact pair-coloring miniature is small enough to exhaust.  The ternary
branch is only a search-pressure scout, used to keep the analogy honest.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from math import comb
from random import Random


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def tuple_atoms(n: int, arity: int) -> tuple[tuple[int, ...], ...]:
    return tuple(combinations(range(n), arity))


def ph_targets(n: int, arity: int, h_size: int) -> tuple[tuple[int, ...], ...]:
    """Return atom-index sets whose monochromaticity would witness PH."""
    atoms = tuple_atoms(n, arity)
    index = {atom: i for i, atom in enumerate(atoms)}
    targets = []
    for h in combinations(range(n), h_size):
        if h_size > min(h):
            targets.append(tuple(index[a] for a in combinations(h, arity)))
    return tuple(targets)


def is_bad_coloring(colors: tuple[int, ...], targets: tuple[tuple[int, ...], ...]) -> bool:
    """A bad coloring has no monochromatic relatively-large target."""
    for target in targets:
        c0 = colors[target[0]]
        if all(colors[i] == c0 for i in target[1:]):
            return False
    return True


def exact_bad_colorings(n: int, arity: int, num_colors: int, h_size: int) -> list[tuple[int, ...]]:
    atoms = tuple_atoms(n, arity)
    targets = ph_targets(n, arity, h_size)
    out = []

    def rec(prefix: list[int]) -> None:
        if len(prefix) == len(atoms):
            colors = tuple(prefix)
            if is_bad_coloring(colors, targets):
                out.append(colors)
            return
        for c in range(num_colors):
            prefix.append(c)
            still_possible = True
            for target in targets:
                if all(i < len(prefix) for i in target):
                    c0 = prefix[target[0]]
                    if all(prefix[i] == c0 for i in target[1:]):
                        still_possible = False
                        break
            if still_possible:
                rec(prefix)
            prefix.pop()

    rec([])
    return out


def extend_pair_coloring(colors: tuple[int, ...], n: int, pattern: int) -> tuple[int, ...]:
    """Extend a pair-coloring on n vertices by colors on edges to vertex n."""
    old_index = {p: i for i, p in enumerate(tuple_atoms(n, 2))}
    new_atoms = tuple_atoms(n + 1, 2)
    out = []
    for a, b in new_atoms:
        if b < n:
            out.append(colors[old_index[(a, b)]])
        else:
            out.append((pattern >> a) & 1)
    return tuple(out)


def edge_count(colors: tuple[int, ...]) -> int:
    return sum(colors)


def degree_sequence(colors: tuple[int, ...], n: int) -> tuple[int, ...]:
    deg = [0] * n
    for color, (a, b) in zip(colors, tuple_atoms(n, 2)):
        if color:
            deg[a] += 1
            deg[b] += 1
    return tuple(deg)


def pair_tree_audit(max_n: int = 6) -> list[dict[str, object]]:
    levels: dict[int, list[tuple[int, ...]]] = {}
    rows = []
    for n in range(1, max_n + 1):
        bad = exact_bad_colorings(n, arity=2, num_colors=2, h_size=3)
        levels[n] = bad
        ext_dist: Counter[int] = Counter()
        by_edge: defaultdict[int, Counter[bool]] = defaultdict(Counter)
        if n < max_n:
            child_targets = ph_targets(n + 1, 2, 3)
            for colors in bad:
                ext_count = 0
                for pattern in range(1 << n):
                    child = extend_pair_coloring(colors, n, pattern)
                    if is_bad_coloring(child, child_targets):
                        ext_count += 1
                ext_dist[ext_count] += 1
                by_edge[edge_count(colors)][ext_count > 0] += 1
        rows.append(
            {
                "N": n,
                "atoms": len(tuple_atoms(n, 2)),
                "targets": len(ph_targets(n, 2, 3)),
                "bad_count": len(bad),
                "extension_count_distribution": dict(sorted(ext_dist.items())),
                "edge_shell_extendability": {
                    e: dict(v) for e, v in sorted(by_edge.items())
                },
            }
        )
    return rows


def exact_pair_channel_audit() -> list[dict[str, object]]:
    """Which coarse channels separate extendable/dead pair bad colorings?"""
    n = 4
    bad = exact_bad_colorings(n, 2, 2, 3)
    child_targets = ph_targets(n + 1, 2, 3)

    def ext_count(colors: tuple[int, ...]) -> int:
        return sum(
            is_bad_coloring(extend_pair_coloring(colors, n, pattern), child_targets)
            for pattern in range(1 << n)
        )

    channels = {
        "edge_count": edge_count,
        "degree_histogram": lambda c: tuple(sorted(degree_sequence(c, n))),
        "rooted_first3_degrees": lambda c: degree_sequence(c, n)[:3],
        "full_coloring": lambda c: c,
    }
    rows = []
    for name, fn in channels.items():
        groups: defaultdict[object, set[bool]] = defaultdict(set)
        sizes: Counter[object] = Counter()
        for colors in bad:
            key = fn(colors)
            groups[key].add(ext_count(colors) > 0)
            sizes[key] += 1
        rows.append(
            {
                "channel": name,
                "groups": len(groups),
                "mixed_extendability_buckets": sum(1 for v in groups.values() if len(v) > 1),
                "max_bucket": max(sizes.values()),
            }
        )
    return rows


def random_bad_pressure(
    n_values: range = range(4, 13),
    trials: int = 2000,
    seed: int = 672,
) -> list[dict[str, object]]:
    rng = Random(seed)
    rows = []
    for n in n_values:
        atoms = tuple_atoms(n, 3)
        targets = ph_targets(n, 3, 4)
        avoid = 0
        min_bad = None
        total_bad = 0
        for _ in range(trials):
            colors = tuple(rng.randrange(3) for _ in atoms)
            bad_targets = 0
            for target in targets:
                c0 = colors[target[0]]
                if all(colors[i] == c0 for i in target[1:]):
                    bad_targets += 1
            if bad_targets == 0:
                avoid += 1
            total_bad += bad_targets
            if min_bad is None or bad_targets < min_bad:
                min_bad = bad_targets
        rows.append(
            {
                "N": n,
                "atoms": len(atoms),
                "targets": len(targets),
                "random_avoid_count": avoid,
                "trials": trials,
                "min_bad_targets_seen": min_bad,
                "avg_bad_targets": round(total_bad / trials, 4),
            }
        )
    return rows


def local_repair_search(n: int, seed: int = 672, restarts: int = 12, steps: int = 1200) -> dict[str, object]:
    """Try to find a ternary bad coloring.  This is heuristic, not a proof."""
    rng = Random(seed + n)
    atoms = tuple_atoms(n, 3)
    targets = ph_targets(n, 3, 4)

    def bad_list(colors: list[int]) -> list[int]:
        out = []
        for i, target in enumerate(targets):
            c0 = colors[target[0]]
            if all(colors[j] == c0 for j in target[1:]):
                out.append(i)
        return out

    best_bad = len(targets)
    for restart in range(restarts):
        colors = [rng.randrange(3) for _ in atoms]
        bad = bad_list(colors)
        for step in range(steps):
            if not bad:
                return {
                    "N": n,
                    "found_bad_coloring": True,
                    "restart": restart,
                    "step": step,
                    "best_bad_targets": 0,
                }
            best_bad = min(best_bad, len(bad))
            target = targets[rng.choice(bad)]
            atom = rng.choice(target)
            old = colors[atom]
            candidates = []
            for c in range(3):
                if c == old:
                    continue
                colors[atom] = c
                nb = bad_list(colors)
                candidates.append((len(nb), c, nb))
            colors[atom] = old
            candidates.sort(key=lambda x: x[0])
            if candidates[0][0] <= len(bad) or rng.random() < 0.01:
                _, c, bad = candidates[0]
                colors[atom] = c
            else:
                colors[atom] = rng.randrange(3)
                bad = bad_list(colors)
    return {
        "N": n,
        "found_bad_coloring": False,
        "restart": None,
        "step": None,
        "best_bad_targets": best_bad,
    }


def ternary_local_search_summary(n_values: range = range(10, 21)) -> list[dict[str, object]]:
    return [local_repair_search(n) for n in n_values]


def orient(a: dict[str, object], b: dict[str, object]) -> bool:
    criteria = [
        "exact_small_evidence",
        "retains_address",
        "outer_extension_rank",
        "quotient_descent_fit",
        "lrc_transfer",
        "overclaim_risk_control",
    ]
    aw = sum(int(a[c]) > int(b[c]) for c in criteria)
    bw = sum(int(b[c]) > int(a[c]) for c in criteria)
    if aw != bw:
        return aw > bw
    return str(a["name"]) < str(b["name"])


def tournament_fingerprint(routes: list[dict[str, object]]) -> dict[str, object]:
    n = len(routes)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if orient(routes[i], routes[j]):
            adj[i][j] = True
        else:
            adj[j][i] = True

    scores = [sum(row) for row in adj]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[b][a] and adj[c][b] and adj[a][c]
        ):
            c3 += 1

    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    sccs = []

    def rdfs(v: int, acc: list[int]) -> None:
        seen[v] = True
        acc.append(v)
        for w in rev[v]:
            if not seen[w]:
                rdfs(w, acc)

    for v in reversed(order):
        if not seen[v]:
            acc: list[int] = []
            rdfs(v, acc)
            sccs.append(tuple(routes[i]["name"] for i in acc))

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used in range(1 << n):
        for v in range(n):
            if not dp[used][v]:
                continue
            for w in range(n):
                if ((used >> w) & 1) == 0 and adj[v][w]:
                    dp[used | (1 << w)][w] += dp[used][v]

    ranking = sorted(
        [(scores[i], routes[i]["name"]) for i in range(n)],
        reverse=True,
    )
    return {
        "score_histogram": dict(sorted(Counter(scores).items())),
        "directed_3cycles": c3,
        "scc_sizes": [len(s) for s in sccs],
        "sccs": sccs,
        "hamiltonian_paths": sum(dp[-1]),
        "ranking": ranking,
    }


def route_tournament() -> tuple[list[dict[str, object]], dict[str, object]]:
    routes = [
        {
            "name": "PH_initial_segment_filter",
            "exact_small_evidence": 2,
            "retains_address": 2,
            "outer_extension_rank": 2,
            "quotient_descent_fit": 2,
            "lrc_transfer": 2,
            "overclaim_risk_control": 2,
        },
        {
            "name": "bad_coloring_outer_extension_tree",
            "exact_small_evidence": 2,
            "retains_address": 2,
            "outer_extension_rank": 2,
            "quotient_descent_fit": 1,
            "lrc_transfer": 2,
            "overclaim_risk_control": 2,
        },
        {
            "name": "endpoint_half_filter_trace",
            "exact_small_evidence": 2,
            "retains_address": 2,
            "outer_extension_rank": 1,
            "quotient_descent_fit": 2,
            "lrc_transfer": 2,
            "overclaim_risk_control": 2,
        },
        {
            "name": "metagraph_ultrafilter_sheaf",
            "exact_small_evidence": 2,
            "retains_address": 2,
            "outer_extension_rank": 1,
            "quotient_descent_fit": 2,
            "lrc_transfer": 2,
            "overclaim_risk_control": 1,
        },
        {
            "name": "LRC_owner_carry_filter_rank",
            "exact_small_evidence": 1,
            "retains_address": 2,
            "outer_extension_rank": 2,
            "quotient_descent_fit": 2,
            "lrc_transfer": 2,
            "overclaim_risk_control": 1,
        },
        {
            "name": "unit_distance_spine_bulk_filter",
            "exact_small_evidence": 1,
            "retains_address": 2,
            "outer_extension_rank": 1,
            "quotient_descent_fit": 1,
            "lrc_transfer": 1,
            "overclaim_risk_control": 1,
        },
        {
            "name": "ordinary_Ramsey_shadow",
            "exact_small_evidence": 2,
            "retains_address": 0,
            "outer_extension_rank": 1,
            "quotient_descent_fit": 0,
            "lrc_transfer": 0,
            "overclaim_risk_control": 2,
        },
        {
            "name": "raw_coloring_count",
            "exact_small_evidence": 2,
            "retains_address": 0,
            "outer_extension_rank": 0,
            "quotient_descent_fit": 0,
            "lrc_transfer": 0,
            "overclaim_risk_control": 2,
        },
        {
            "name": "generic_ultrafilter_CH_analogy",
            "exact_small_evidence": 0,
            "retains_address": 2,
            "outer_extension_rank": 1,
            "quotient_descent_fit": 1,
            "lrc_transfer": 0,
            "overclaim_risk_control": 0,
        },
    ]
    return routes, tournament_fingerprint(routes)


def main() -> None:
    banner("S672 Paris-Harrington ultrafilter coloring recursion")
    print("HYP-2247.")
    print(
        "Interpretation: the metagraph is not the ultrafilter itself; it is a "
        "base over which side-choice filters must descend.  Paris-Harrington "
        "adds a recursion-theoretic gate: coherent bad side choices cannot keep "
        "escaping into the unnamed tail once relative largeness names an initial "
        "segment."
    )

    banner("Exact pair-coloring miniature: x=2 colors pairs into 2 colors")
    pair_rows = pair_tree_audit()
    for row in pair_rows:
        print(row)
    print(
        "Least forced N for this miniature is 6: bad counts are "
        f"{[row['bad_count'] for row in pair_rows]}."
    )
    print(
        "At N=4 only the middle edge-count shell e=3 extends; lower/upper "
        "shells e=2 and e=4 are dead.  At N=5 the surviving middle shell e=5 "
        "has no child."
    )

    banner("Exact N=4 quotient channels for extendability")
    for row in exact_pair_channel_audit():
        print(row)
    print(
        "For this tiny pair case, even edge_count separates extendability.  The "
        "lesson is not that scalars always work; it is that the next useful "
        "coordinate is a derivative one: how many bad children remain."
    )

    banner("Ternary scout: color triples into 3 colors, seek large 4-set")
    print("Random-pressure rows; not a proof of existence or nonexistence.")
    for row in random_bad_pressure():
        print(row)
    print("Local-repair search for explicit bad colorings; heuristic only.")
    for row in ternary_local_search_summary():
        print(row)
    print(
        "The ternary branch keeps escaping much farther than the pair branch.  "
        "That is the right small-size analogy for the true PH growth: the "
        "obstruction is not raw coloring count but tail-escape rank under "
        "outer extension."
    )

    banner("Tournament Analysis over transfer routes")
    routes, fp = route_tournament()
    for route in routes:
        print(route)
    print(fp)

    banner("Repo synthesis")
    print(
        "HYP-2245: divisor-210 and fixed-path tiling cubes have literal "
        "principal ultrafilters, but tournament isomorphism quotient leaks."
    )
    print(
        "HYP-2246: endpoint enumeration repairs quotient leaks by attaching a "
        "three-state L/M/U deletion-owner address."
    )
    print(
        "HYP-2247: Paris-Harrington suggests adding an extension-rank address: "
        "not only upper/lower side, but the number/profile of coherent bad "
        "children still available."
    )
    print(
        "LRC14 proposal: define bad nodes as owner/carry fibers with no witness, "
        "outer extensions as coherent +27/carry-owner lifts, and the relative "
        "largeness gate as an early named residue/owner section.  A proof target "
        "is a PH-style rank that strictly decreases after the owner-private "
        "deletion filter of HYP-2241."
    )
    print(
        "Unit-distance proposal: treat the unit-spine/bulk split as the initial "
        "segment gate.  Bad constructions are those that postpone every required "
        "unit-distance Hamiltonian spine into bulk/tail coordinates; the repair "
        "is to retain spine-owner deletion rank rather than just edge count."
    )


if __name__ == "__main__":
    main()
