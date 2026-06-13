#!/usr/bin/env python3
"""S671: endpoint ultrafilter traces as A000568 enumeration speedups.

HYP-2245 says literal filters live on Boolean address cubes, while the
metagraph quotient leaks unless the missing address descends.  This script
turns that into an enumeration experiment.

Given one representative of each tournament class on n vertices, add a new
vertex by an incident-edge bit vector in {0,1}^n.  Quotient the bit vectors by
Aut(parent).  This is the standard one-vertex extension skeleton, but the new
test is:

  Can paired deletion / upper-lower filter traces classify most child
  candidates before full child canonicalization?

If yes, A000568 enumeration can be staged as:

  parent class -> Aut(parent)-orbits in incident cube -> deletion trace
  -> canonical fallback only on trace collisions.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations
from math import comb
import time


A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
    10: 9733056,
}


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


@lru_cache(maxsize=None)
def pair_list(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {p: i for i, p in enumerate(pair_list(n))}


def bit_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return pair_index(n)[(i, j)]


def has_arc(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        return False
    if i < j:
        return ((mask >> bit_index(n, i, j)) & 1) == 1
    return ((mask >> bit_index(n, j, i)) & 1) == 0


def outdegrees(mask: int, n: int) -> tuple[int, ...]:
    scores = [0] * n
    for i, j in pair_list(n):
        if has_arc(mask, n, i, j):
            scores[i] += 1
        else:
            scores[j] += 1
    return tuple(scores)


@lru_cache(maxsize=None)
def grouped_perms_for_scores(scores: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    groups: dict[int, list[int]] = defaultdict(list)
    for v, s in enumerate(scores):
        groups[s].append(v)
    blocks = [tuple(groups[s]) for s in sorted(groups)]

    def rec(k: int):
        if k == len(blocks):
            yield ()
            return
        for p in permutations(blocks[k]):
            for tail in rec(k + 1):
                yield p + tail

    return tuple(rec(0))


def encode_under_perm(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """New label a is old vertex perm[a]."""
    out = 0
    for a, b in pair_list(n):
        if has_arc(mask, n, perm[a], perm[b]):
            out |= 1 << bit_index(n, a, b)
    return out


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    scores = outdegrees(mask, n)
    best = None
    for p in grouped_perms_for_scores(scores):
        code = encode_under_perm(mask, n, p)
        if best is None or code < best:
            best = code
    assert best is not None
    return best


@lru_cache(maxsize=None)
def automorphisms(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
    scores = outdegrees(mask, n)
    out = []
    for p in grouped_perms_for_scores(scores):
        if encode_under_perm(mask, n, p) == mask:
            out.append(p)
    return tuple(out)


def extend_by_pattern(parent_mask: int, n: int, pattern: int) -> int:
    """Add vertex n. Pattern bit i=1 means old vertex i beats the new vertex."""
    child = 0
    for i, j in pair_list(n):
        if has_arc(parent_mask, n, i, j):
            child |= 1 << bit_index(n + 1, i, j)
    for i in range(n):
        if (pattern >> i) & 1:
            child |= 1 << bit_index(n + 1, i, n)
        # bit 0 means new vertex beats i; pair bit remains 0.
    return child


def permute_pattern(pattern: int, n: int, perm: tuple[int, ...]) -> int:
    """Pattern after relabeling old vertices: new bit a is old bit perm[a]."""
    out = 0
    for a in range(n):
        if (pattern >> perm[a]) & 1:
            out |= 1 << a
    return out


def pattern_orbit_reps(n: int, auts: tuple[tuple[int, ...], ...]) -> list[int]:
    seen = set()
    reps = []
    for pattern in range(1 << n):
        if pattern in seen:
            continue
        orb = {permute_pattern(pattern, n, g) for g in auts}
        seen |= orb
        reps.append(min(orb))
    return reps


def delete_vertex(mask: int, n: int, drop: int) -> int:
    old_vertices = [v for v in range(n) if v != drop]
    out = 0
    for a, b in pair_list(n - 1):
        if has_arc(mask, n, old_vertices[a], old_vertices[b]):
            out |= 1 << bit_index(n - 1, a, b)
    return out


def c3_count(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        ab = has_arc(mask, n, a, b)
        bc = has_arc(mask, n, b, c)
        ca = has_arc(mask, n, c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            total += 1
    return total


def scc_sizes(mask: int, n: int) -> tuple[int, ...]:
    graph = [[] for _ in range(n)]
    rev = [[] for _ in range(n)]
    for i, j in pair_list(n):
        if has_arc(mask, n, i, j):
            graph[i].append(j)
            rev[j].append(i)
        else:
            graph[j].append(i)
            rev[i].append(j)

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
    sizes = []

    def rdfs(v: int) -> int:
        seen[v] = True
        size = 1
        for w in rev[v]:
            if not seen[w]:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used in range(1 << n):
        for v in range(n):
            cur = dp[used][v]
            if not cur:
                continue
            for w in range(n):
                if ((used >> w) & 1) == 0 and has_arc(mask, n, v, w):
                    dp[used | (1 << w)][w] += cur
    return sum(dp[-1])


def score_side(score: int, n_minus_1: int) -> str:
    """Coarse upper/lower side of a deleted vertex score."""
    twice = 2 * score
    if twice < n_minus_1:
        return "L"
    if twice > n_minus_1:
        return "U"
    return "M"


def build_next_level(
    n: int, parent_reps: list[int]
) -> tuple[list[int], list[dict[str, int]]]:
    child_set = set()
    candidates = []
    for parent_id, parent in enumerate(parent_reps):
        auts = automorphisms(parent, n)
        reps = pattern_orbit_reps(n, auts)
        for pattern in reps:
            child = extend_by_pattern(parent, n, pattern)
            can = canonical(child, n + 1)
            child_set.add(can)
            candidates.append(
                {
                    "parent_id": parent_id,
                    "pattern": pattern,
                    "pattern_weight": pattern.bit_count(),
                    "new_score": n - pattern.bit_count(),
                    "child_raw": child,
                    "child_can": can,
                }
            )
    return sorted(child_set), candidates


def build_tower(max_n: int = 8):
    levels: dict[int, list[int]] = {1: [0]}
    transitions = {}
    for n in range(1, max_n):
        t0 = time.time()
        child_reps, candidates = build_next_level(n, levels[n])
        elapsed = time.time() - t0
        levels[n + 1] = child_reps
        transitions[(n, n + 1)] = {
            "candidates": candidates,
            "elapsed": elapsed,
        }
    return levels, transitions


def class_maps(levels: dict[int, list[int]]) -> dict[int, dict[int, int]]:
    return {n: {mask: i for i, mask in enumerate(reps)} for n, reps in levels.items()}


def child_signatures(
    mask: int, n: int, parent_id_map: dict[int, int]
) -> dict[str, tuple]:
    scores = outdegrees(mask, n)
    parent_ids = []
    paired = []
    half = []
    delta_c3 = []
    child_c3 = c3_count(mask, n)
    for v in range(n):
        card = canonical(delete_vertex(mask, n, v), n - 1)
        cid = parent_id_map[card]
        parent_ids.append(cid)
        paired.append((cid, scores[v]))
        half.append((cid, score_side(scores[v], n - 1)))
        delta_c3.append((cid, scores[v], child_c3 - c3_count(card, n - 1)))
    return {
        "card_multiset": tuple(sorted(parent_ids)),
        "half_filter_trace": tuple(sorted(half)),
        "paired_score_trace": tuple(sorted(paired)),
        "paired_c3_loss_trace": tuple(sorted(delta_c3)),
        "score_c3_scc": (
            tuple(sorted(Counter(scores).items())),
            child_c3,
            scc_sizes(mask, n),
        ),
    }


def audit_signature_on_classes(
    n: int, reps: list[int], parent_id_map: dict[int, int]
) -> dict[str, dict[str, object]]:
    buckets: dict[str, defaultdict[tuple, list[int]]] = {
        "score_c3_scc": defaultdict(list),
        "card_multiset": defaultdict(list),
        "half_filter_trace": defaultdict(list),
        "paired_score_trace": defaultdict(list),
        "paired_c3_loss_trace": defaultdict(list),
    }
    for cid, mask in enumerate(reps):
        sigs = child_signatures(mask, n, parent_id_map)
        for name, sig in sigs.items():
            buckets[name][sig].append(cid)

    out = {}
    for name, bucket in buckets.items():
        sizes = [len(v) for v in bucket.values()]
        out[name] = {
            "groups": len(bucket),
            "mixed_buckets": sum(1 for s in sizes if s > 1),
            "max_bucket": max(sizes) if sizes else 0,
            "singleton_classes": sum(s for s in sizes if s == 1),
            "colliding_classes": sum(s for s in sizes if s > 1),
        }
    return out


def audit_signature_on_candidates(
    child_n: int,
    candidates: list[dict[str, int]],
    parent_id_map: dict[int, int],
) -> dict[str, dict[str, object]]:
    buckets: dict[str, defaultdict[tuple, list[int]]] = {
        "score_c3_scc": defaultdict(list),
        "card_multiset": defaultdict(list),
        "half_filter_trace": defaultdict(list),
        "paired_score_trace": defaultdict(list),
        "paired_c3_loss_trace": defaultdict(list),
    }
    for i, row in enumerate(candidates):
        mask = row["child_raw"]
        sigs = child_signatures(mask, child_n, parent_id_map)
        for name, sig in sigs.items():
            buckets[name][sig].append(i)

    out = {}
    for name, bucket in buckets.items():
        mixed_buckets = 0
        mixed_candidates = 0
        pure_duplicate_buckets = 0
        canonical_fallback_candidates = 0
        for idxs in bucket.values():
            true_classes = {candidates[i]["child_can"] for i in idxs}
            if len(true_classes) > 1:
                mixed_buckets += 1
                mixed_candidates += len(idxs)
                canonical_fallback_candidates += len(idxs)
            elif len(idxs) > 1:
                pure_duplicate_buckets += 1
        out[name] = {
            "candidate_groups": len(bucket),
            "mixed_buckets": mixed_buckets,
            "mixed_candidates": mixed_candidates,
            "pure_duplicate_buckets": pure_duplicate_buckets,
            "canonical_fallback_candidates": canonical_fallback_candidates,
            "avoid_full_child_canon_after_trace": len(candidates)
            - canonical_fallback_candidates,
        }
    return out


def card_collision_details(
    n: int, reps: list[int], parent_id_map: dict[int, int]
) -> list[dict[str, object]]:
    buckets: defaultdict[tuple, list[int]] = defaultdict(list)
    signatures = {}
    for cid, mask in enumerate(reps):
        sigs = child_signatures(mask, n, parent_id_map)
        signatures[cid] = sigs
        buckets[sigs["card_multiset"]].append(cid)

    out = []
    for cids in buckets.values():
        if len(cids) <= 1:
            continue
        rows = []
        for cid in cids:
            mask = reps[cid]
            rows.append(
                {
                    "cid": cid,
                    "H": hamiltonian_paths(mask, n),
                    "c3": c3_count(mask, n),
                    "score_hist": tuple(sorted(Counter(outdegrees(mask, n)).items())),
                    "scc": scc_sizes(mask, n),
                }
            )
        first, second = cids[0], cids[1]
        diff = Counter(signatures[first]["half_filter_trace"])
        diff.subtract(Counter(signatures[second]["half_filter_trace"]))
        half_delta = tuple(sorted((k, v) for k, v in diff.items() if v))
        out.append({"classes": rows, "half_delta": half_delta})
    return out


def count_pattern_orbits_only(n: int, parent_reps: list[int]) -> int:
    total = 0
    for parent in parent_reps:
        total += len(pattern_orbit_reps(n, automorphisms(parent, n)))
    return total


def count_directed_3cycles(adj: dict[str, dict[str, bool]]) -> int:
    vertices = list(adj)
    total = 0
    for a, b, c in combinations(vertices, 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes_digraph(adj: dict[str, dict[str, bool]]) -> list[int]:
    vertices = list(adj)

    def reach(start: str, reverse: bool = False) -> set[str]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w in vertices:
                edge = adj[w][v] if reverse else adj[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(vertices)
    sizes = []
    while remaining:
        v = next(iter(remaining))
        comp = reach(v) & reach(v, True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths_digraph(adj: dict[str, dict[str, bool]]) -> int:
    vertices = list(adj)
    index = {v: i for i, v in enumerate(vertices)}
    dp = {}
    for v in vertices:
        dp[(1 << index[v], v)] = 1
    full = (1 << len(vertices)) - 1
    for used in range(1 << len(vertices)):
        for v in vertices:
            cur = dp.get((used, v), 0)
            if not cur:
                continue
            for w in vertices:
                bit = 1 << index[w]
                if used & bit:
                    continue
                if adj[v][w]:
                    dp[(used | bit, w)] = dp.get((used | bit, w), 0) + cur
    return sum(dp.get((full, v), 0) for v in vertices)


def method_tournament():
    # Metrics: exactness, canonical-work avoidance, invariant retention,
    # scaling, implementation readiness, LRC transfer.
    channels = {
        "paired_deletion_filter_trace": (5, 5, 5, 4, 4, 5),
        "endpoint_aut_orbit_generation": (5, 4, 4, 5, 5, 4),
        "modular_substitution_dp": (5, 5, 4, 5, 3, 3),
        "THM410_interval_ledger": (4, 4, 3, 5, 5, 4),
        "half_filter_trace": (4, 4, 4, 4, 5, 5),
        "raw_fixed_path_cube": (5, 1, 2, 2, 4, 2),
        "raw_labeled_enumeration": (5, 0, 1, 0, 2, 1),
    }
    tie_order = list(channels)
    adj = {a: {} for a in channels}
    for a in channels:
        for b in channels:
            if a == b:
                adj[a][b] = False
                continue
            wins = sum(x > y for x, y in zip(channels[a], channels[b]))
            losses = sum(x < y for x, y in zip(channels[a], channels[b]))
            if wins == losses:
                adj[a][b] = tie_order.index(a) < tie_order.index(b)
            else:
                adj[a][b] = wins > losses
    scores = {a: sum(adj[a][b] for b in channels if b != a) for a in channels}
    order = sorted(channels, key=lambda x: (-scores[x], tie_order.index(x)))
    return channels, adj, scores, order


def main() -> None:
    max_exact = 8
    banner("S671 endpoint ultrafilter enumeration audit")
    print("Goal: use HYP-2245 side-choice/descent language as an A000568 speedup.")
    print("Build exact vertex-extension tower through n=8 by Aut(parent) pattern orbits.")

    levels, transitions = build_tower(max_exact)
    ids = class_maps(levels)

    banner("A. Aut(parent)-orbit endpoint generation")
    print(
        "child  A000568  orbit_candidates  duplicate  "
        "raw_patterns/orbit  fixed_path/orbit  labeled/orbit  time_s"
    )
    for n in range(1, max_exact):
        child_n = n + 1
        candidates = transitions[(n, child_n)]["candidates"]
        orbit = len(candidates)
        unique = len(levels[child_n])
        raw_patterns = len(levels[n]) * (1 << n)
        fixed_path = 1 << comb(child_n - 1, 2)
        labeled = 1 << comb(child_n, 2)
        print(
            f"{child_n:5d} {unique:8d} {orbit:17d} "
            f"{orbit / unique:9.3f} {raw_patterns / orbit:17.3f} "
            f"{fixed_path / orbit:16.3f} {labeled / orbit:14.3f} "
            f"{transitions[(n, child_n)]['elapsed']:7.3f}"
        )
        if child_n in A000568 and unique != A000568[child_n]:
            print(f"  WARNING expected {A000568[child_n]}")

    next_orbits = count_pattern_orbits_only(max_exact, levels[max_exact])
    print(
        f"\nOne-step estimate without building n={max_exact+1}: "
        f"from A({max_exact})={len(levels[max_exact])}, "
        f"Aut-orbit endpoint candidates={next_orbits}, "
        f"known A({max_exact+1})={A000568.get(max_exact+1, '?')}"
    )
    if max_exact + 1 in A000568:
        print(
            f"candidate/known ratio={next_orbits / A000568[max_exact + 1]:.3f}; "
            f"raw pattern count={len(levels[max_exact]) * (1 << max_exact)}"
        )

    banner("B. Deletion/filter signatures on finished classes")
    for n in range(3, max_exact + 1):
        summary = audit_signature_on_classes(n, levels[n], ids[n - 1])
        print(f"n={n} classes={len(levels[n])}")
        for name, row in summary.items():
            print(
                f"  {name:22s} groups={row['groups']:5d} "
                f"mixed={row['mixed_buckets']:4d} max={row['max_bucket']:3d} "
                f"singletons={row['singleton_classes']:5d} "
                f"colliding_classes={row['colliding_classes']:5d}"
            )

    banner("B2. Card-deck collision anatomy at n=8")
    details = card_collision_details(max_exact, levels[max_exact], ids[max_exact - 1])
    if not details:
        print("No unpaired card-deck collisions at the top exact level.")
    for i, item in enumerate(details, 1):
        print(f"collision bucket {i}:")
        for row in item["classes"]:
            print(
                f"  cid={row['cid']:4d} H={row['H']:4d} c3={row['c3']:2d} "
                f"scores={row['score_hist']} scc={row['scc']}"
            )
        print("  half-filter symmetric delta:")
        for (card_side, delta) in item["half_delta"]:
            print(f"    {card_side}: {delta:+d}")

    banner("C. Candidate pre-canonicalization fallback audit")
    print(
        "For each endpoint-orbit candidate, compute traces directly from the "
        "candidate.  In production, only mixed trace buckets need full child "
        "canonical fallback."
    )
    for n in range(2, max_exact):
        child_n = n + 1
        candidates = transitions[(n, child_n)]["candidates"]
        summary = audit_signature_on_candidates(child_n, candidates, ids[n])
        print(f"transition {n}->{child_n}: candidates={len(candidates)}")
        for name, row in summary.items():
            avoided = row["avoid_full_child_canon_after_trace"]
            print(
                f"  {name:22s} groups={row['candidate_groups']:6d} "
                f"mixed_buckets={row['mixed_buckets']:4d} "
                f"mixed_candidates={row['mixed_candidates']:6d} "
                f"pure_dup={row['pure_duplicate_buckets']:5d} "
                f"avoid={avoided:6d} ({100*avoided/len(candidates):5.1f}%)"
            )

    banner("D. Other property payloads available while enumerating")
    for n in range(3, max_exact + 1):
        h_values = Counter()
        c3_values = Counter()
        scc_values = Counter()
        for mask in levels[n]:
            h_values[hamiltonian_paths(mask, n)] += 1
            c3_values[c3_count(mask, n)] += 1
            scc_values[scc_sizes(mask, n)] += 1
        print(
            f"n={n}: H distinct={len(h_values)} c3 distinct={len(c3_values)} "
            f"SCC profiles={len(scc_values)} H_range=({min(h_values)},{max(h_values)}) "
            f"c3_range=({min(c3_values)},{max(c3_values)})"
        )

    banner("E. Tournament Analysis over enumeration routes")
    channels, adj, scores, order = method_tournament()
    print("vertices=enumeration/speedup routes")
    print(
        "observable=(exactness, canonical-work avoidance, invariant retention, "
        "scaling, readiness, LRC transfer)"
    )
    print("switch=majority; tie Hamiltonian path=listed priority order")
    print(f"score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"directed_3cycles={count_directed_3cycles(adj)}")
    print(f"scc_sizes={scc_sizes_digraph(adj)}")
    print(f"hamiltonian_paths={hamiltonian_paths_digraph(adj)}")
    print("top_order:")
    for name in order:
        print(f"  {name:34s} score={scores[name]} vector={channels[name]}")

    banner("F. Algorithmic theorem shape")
    print("Endpoint-ultrafilter enumerator:")
    print("  1. Keep class reps and Aut(parent).")
    print("  2. Enumerate Aut(parent)-orbits in the incident Boolean cube.")
    print("  3. Attach deletion/filter trace: (deleted-card id, owner score/side).")
    print("  4. If trace bucket is pure, accept one representative without full child")
    print("     canonicalization; use canonical fallback only on mixed trace buckets.")
    print("  5. Cache H, c3, SCC, and OCF payloads at the same time.")
    print("\nInterpretation:")
    print("  The incident cube has literal ultrafilters.  The paired deletion trace is")
    print("  the address that makes much of that side choice descend through the")
    print("  isomorphism quotient.  This is HYP-2245 as an enumerator design rule.")


if __name__ == "__main__":
    main()
