#!/usr/bin/env python3
"""
S658: small carrier atlases for HYP-2233 frontier lanes.

This is a small-first follow-up, not a proof of any external conjecture.
It builds three finite predicate-preserving atlases:

1. tournament reconstruction decks,
2. union-closed frequency carriers,
3. finite-field Kakeya/Falconer incidence toys.

Tournament Analysis:
  Vertices are the three worked frontier lanes.  Pairwise observable is
  (finite_exactness, ambiguity_reduction, side_channel_visibility,
   near_term_actionability).  The tie Hamiltonian path is the ranked lane list.

Assumption challenge:
  For reconstruction, vertices could be graphs, tournaments, cards, invariants,
  or deck buckets; we choose iso-classes because the preserved predicate is
  reconstruction from projections.  For union-closed sets, vertices could be
  elements, sets, closure generators, or frequency vectors; we compute both
  element and family quotients because the Frankl predicate lives on elements
  but closure lives on sets.  For finite-field incidence, vertices could be
  points, directions, lines, or distance fibers; we choose line-choice carriers
  because they preserve the Kakeya predicate while exposing the distance side
  channel.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations, product
from math import ceil


# ---------------------------------------------------------------------------
# Tournament utilities


def edge_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    idx = 0
    for a in range(n):
        for b in range(a + 1, n):
            if a == i and b == j:
                return idx
            idx += 1
    raise ValueError((n, i, j))


EDGE_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PERMS: dict[int, list[tuple[int, ...]]] = {}


def prepare_tournaments(max_n: int = 6) -> None:
    for n in range(1, max_n + 1):
        EDGE_INDEX[n] = {}
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                EDGE_INDEX[n][(i, j)] = idx
                idx += 1
        PERMS[n] = list(permutations(range(n)))


def has_edge(bits: int, n: int, i: int, j: int) -> bool:
    if i == j:
        return False
    if i < j:
        idx = EDGE_INDEX[n][(i, j)]
        return bool((bits >> idx) & 1)
    idx = EDGE_INDEX[n][(j, i)]
    return not bool((bits >> idx) & 1)


def permute_tournament(bits: int, n: int, perm: tuple[int, ...]) -> int:
    out = 0
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if has_edge(bits, n, perm[i], perm[j]):
                out |= 1 << idx
            idx += 1
    return out


def canon_tournament(bits: int, n: int) -> int:
    return min(permute_tournament(bits, n, perm) for perm in PERMS[n])


def delete_vertex(bits: int, n: int, v: int) -> int:
    keep = [u for u in range(n) if u != v]
    out = 0
    idx = 0
    for i in range(n - 1):
        for j in range(i + 1, n - 1):
            if has_edge(bits, n, keep[i], keep[j]):
                out |= 1 << idx
            idx += 1
    return out


def insert_vertex(bits: int, old_n: int, new_out_mask: int) -> int:
    """Insert vertex old_n; mask bit i says the new vertex beats old vertex i."""
    n = old_n + 1
    out = 0
    for i in range(old_n):
        for j in range(i + 1, old_n):
            if has_edge(bits, old_n, i, j):
                out |= 1 << EDGE_INDEX[n][(i, j)]
    for i in range(old_n):
        if not ((new_out_mask >> i) & 1):
            out |= 1 << EDGE_INDEX[n][(i, old_n)]
    return out


def tournament_iso_classes(max_n: int) -> dict[int, list[int]]:
    classes: dict[int, list[int]] = {1: [0]}
    for n in range(2, max_n + 1):
        candidates: set[int] = set()
        for old in classes[n - 1]:
            for new_out_mask in range(1 << (n - 1)):
                candidates.add(canon_tournament(insert_vertex(old, n - 1, new_out_mask), n))
        classes[n] = sorted(candidates)
    return classes


def scores(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(sum(1 for j in range(n) if has_edge(bits, n, i, j)) for i in range(n)))


def c3_count(bits: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        if has_edge(bits, n, a, b) and has_edge(bits, n, b, c) and has_edge(bits, n, c, a):
            total += 1
        elif has_edge(bits, n, b, a) and has_edge(bits, n, c, b) and has_edge(bits, n, a, c):
            total += 1
    return total


def scc_sizes(bits: int, n: int) -> tuple[int, ...]:
    adj = [[has_edge(bits, n, i, j) for j in range(n)] for i in range(n)]
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, ok in enumerate(graph[u]):
                if ok and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    remaining = set(range(n))
    sizes: list[int] = []
    while remaining:
        start = min(remaining)
        comp = reach(start, adj) & reach(start, radj)
        sizes.append(len(comp))
        remaining -= comp
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(bits: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if has_edge(bits, n, last, nxt):
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def scalar_signature(bits: int, n: int) -> tuple[object, ...]:
    return (hamiltonian_paths(bits, n), scores(bits, n), c3_count(bits, n), scc_sizes(bits, n))


def deck_signatures(bits: int, n: int) -> tuple[tuple[int, ...], tuple[tuple[object, ...], ...], tuple[tuple[int, int, int], ...]]:
    full_cards: list[int] = []
    scalar_cards: list[tuple[object, ...]] = []
    losses: list[tuple[int, int, int]] = []
    h0 = hamiltonian_paths(bits, n)
    c30 = c3_count(bits, n)
    out_scores = [sum(1 for j in range(n) if has_edge(bits, n, i, j)) for i in range(n)]
    for v in range(n):
        card = delete_vertex(bits, n, v)
        ccard = canon_tournament(card, n - 1)
        full_cards.append(ccard)
        scalar_cards.append(scalar_signature(ccard, n - 1))
        losses.append((h0 - hamiltonian_paths(card, n - 1), c30 - c3_count(card, n - 1), out_scores[v]))
    return tuple(sorted(full_cards)), tuple(sorted(scalar_cards, key=repr)), tuple(sorted(losses))


def bucket_stats(items: dict[object, list[int]]) -> tuple[int, int, int]:
    sizes = [len(v) for v in items.values()]
    return len(items), max(sizes), sum(1 for s in sizes if s > 1)


def tournament_deck_atlas(max_n: int = 6) -> list[dict[str, object]]:
    prepare_tournaments(max_n)
    classes_by_n = tournament_iso_classes(max_n)
    rows: list[dict[str, object]] = []
    for n in range(3, max_n + 1):
        classes = classes_by_n[n]
        scalar_buckets: dict[object, list[int]] = defaultdict(list)
        full_deck_buckets: dict[object, list[int]] = defaultdict(list)
        scalar_deck_buckets: dict[object, list[int]] = defaultdict(list)
        loss_buckets: dict[object, list[int]] = defaultdict(list)
        for bits in classes:
            scalar_buckets[scalar_signature(bits, n)].append(bits)
            full_deck, scalar_deck, losses = deck_signatures(bits, n)
            full_deck_buckets[full_deck].append(bits)
            scalar_deck_buckets[scalar_deck].append(bits)
            loss_buckets[(scalar_deck, losses)].append(bits)

        worst_scalar = max(scalar_buckets.values(), key=len)
        worst_full_deck = max(full_deck_buckets.values(), key=len)
        worst_scalar_deck = max(scalar_deck_buckets.values(), key=len)
        rows.append(
            {
                "n": n,
                "classes": len(classes),
                "scalar_bucket_stats": bucket_stats(scalar_buckets),
                "full_deck_bucket_stats": bucket_stats(full_deck_buckets),
                "scalar_deck_bucket_stats": bucket_stats(scalar_deck_buckets),
                "loss_deck_bucket_stats": bucket_stats(loss_buckets),
                "worst_scalar_bucket": [hex(x) for x in worst_scalar[:6]],
                "worst_full_deck_bucket": [hex(x) for x in worst_full_deck[:6]],
                "worst_scalar_deck_bucket": [hex(x) for x in worst_scalar_deck[:6]],
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Union-closed utilities


SUBSET_PERMS: dict[int, list[list[int]]] = {}


def prepare_subset_perms(m: int) -> None:
    maps: list[list[int]] = []
    for perm in permutations(range(m)):
        table = []
        for subset in range(1 << m):
            image = 0
            for i in range(m):
                if subset & (1 << i):
                    image |= 1 << perm[i]
            table.append(image)
        maps.append(table)
    SUBSET_PERMS[m] = maps


def permute_family(mask: int, m: int, table: list[int]) -> int:
    out = 0
    for subset in range(1 << m):
        if mask & (1 << subset):
            out |= 1 << table[subset]
    return out


def canon_family(mask: int, m: int) -> int:
    return min(permute_family(mask, m, table) for table in SUBSET_PERMS[m])


def family_sets(mask: int, m: int) -> list[int]:
    return [s for s in range(1 << m) if mask & (1 << s)]


def is_union_closed(mask: int, m: int) -> bool:
    sets = family_sets(mask, m)
    present = set(sets)
    for a in sets:
        for b in sets:
            if (a | b) not in present:
                return False
    return True


def family_frequencies(mask: int, m: int) -> tuple[int, ...]:
    sets = family_sets(mask, m)
    return tuple(sorted((sum(1 for s in sets if s & (1 << i)) for i in range(m)), reverse=True))


def family_size_distribution(mask: int, m: int) -> tuple[int, ...]:
    ctr = Counter(bin(s).count("1") for s in family_sets(mask, m))
    return tuple(ctr.get(k, 0) for k in range(m + 1))


def union_pressure_signature(mask: int, m: int) -> tuple[int, ...]:
    sets = family_sets(mask, m)
    vals = []
    for a in sets:
        vals.append(sum(bin(a | b).count("1") for b in sets))
    return tuple(sorted(vals, reverse=True))


def union_closed_atlas(max_m: int = 4) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for m in range(1, max_m + 1):
        prepare_subset_perms(m)
        canonicals: set[int] = set()
        for mask in range(1, 1 << (1 << m)):
            if not is_union_closed(mask, m):
                continue
            if all(s == 0 for s in family_sets(mask, m)):
                continue
            canonicals.add(canon_family(mask, m))

        freq_buckets: dict[object, list[int]] = defaultdict(list)
        freq_size_buckets: dict[object, list[int]] = defaultdict(list)
        pressure_buckets: dict[object, list[int]] = defaultdict(list)
        tight = 0
        frankl_failures = 0
        witness_counts: Counter[int] = Counter()
        min_ratio_num = 10**9
        min_ratio_den = 1
        min_examples: list[int] = []

        for fam in sorted(canonicals):
            sets = family_sets(fam, m)
            size = len(sets)
            freqs = family_frequencies(fam, m)
            max_freq = max(freqs)
            if 2 * max_freq < size:
                frankl_failures += 1
            witness_counts[sum(1 for f in freqs if 2 * f >= size)] += 1
            if 2 * max_freq == size:
                tight += 1
            if max_freq * min_ratio_den < min_ratio_num * size:
                min_ratio_num, min_ratio_den = max_freq, size
                min_examples = [fam]
            elif max_freq * min_ratio_den == min_ratio_num * size and len(min_examples) < 6:
                min_examples.append(fam)

            freq_sig = (size, freqs)
            freq_buckets[freq_sig].append(fam)
            freq_size_buckets[(freq_sig, family_size_distribution(fam, m))].append(fam)
            pressure_buckets[(freq_sig, union_pressure_signature(fam, m))].append(fam)

        rows.append(
            {
                "m": m,
                "canonical_families": len(canonicals),
                "frankl_failures": frankl_failures,
                "tight_half_frequency": tight,
                "min_max_frequency_ratio": f"{min_ratio_num}/{min_ratio_den}",
                "witness_count_hist": dict(sorted(witness_counts.items())),
                "freq_bucket_stats": bucket_stats(freq_buckets),
                "freq_plus_size_bucket_stats": bucket_stats(freq_size_buckets),
                "freq_plus_pressure_bucket_stats": bucket_stats(pressure_buckets),
                "min_examples": [hex(x) for x in min_examples],
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Finite-field Kakeya/Falconer utilities


def line_points(p: int, direction: int, b: int) -> list[tuple[int, int]]:
    if direction == p:
        return [(b, y) for y in range(p)]
    return [(x, (direction * x + b) % p) for x in range(p)]


def kakeya_choice_points(p: int, choice: tuple[int, ...]) -> tuple[set[tuple[int, int]], Counter[tuple[int, int]]]:
    counts: Counter[tuple[int, int]] = Counter()
    for direction, b in enumerate(choice):
        counts.update(line_points(p, direction, b))
    return set(counts), counts


def distance_stats(p: int, pts: set[tuple[int, int]]) -> tuple[int, int, int, tuple[int, ...]]:
    all_dist: set[int] = set()
    pinned_counts = []
    points = sorted(pts)
    for i, a in enumerate(points):
        pinned: set[int] = set()
        for j, b in enumerate(points):
            if i == j:
                continue
            d = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) % p
            all_dist.add(d)
            pinned.add(d)
        pinned_counts.append(len(pinned))
    return len(all_dist), min(pinned_counts), max(pinned_counts), tuple(sorted(Counter(pinned_counts).items()))


def multiplicity_profile(counts: Counter[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(counts.values()).items()))


def kakeya_atlas(primes: tuple[int, ...] = (3, 5, 7), sample_limit: int = 256) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for p in primes:
        size_hist: Counter[int] = Counter()
        min_size = 10**9
        min_count = 0
        min_choices: list[tuple[int, ...]] = []
        profile_hist: Counter[tuple[tuple[int, int], ...]] = Counter()

        for choice in product(range(p), repeat=p + 1):
            pts, counts = kakeya_choice_points(p, choice)
            size = len(pts)
            size_hist[size] += 1
            if size < min_size:
                min_size = size
                min_count = 1
                min_choices = [choice]
                profile_hist = Counter({multiplicity_profile(counts): 1})
            elif size == min_size:
                min_count += 1
                profile_hist[multiplicity_profile(counts)] += 1
                if len(min_choices) < sample_limit:
                    min_choices.append(choice)

        dist_hist: Counter[tuple[int, int, int, tuple[int, ...]]] = Counter()
        for choice in min_choices:
            pts, _counts = kakeya_choice_points(p, choice)
            dist_hist[distance_stats(p, pts)] += 1

        rows.append(
            {
                "p": p,
                "choices": p ** (p + 1),
                "size_hist": dict(sorted(size_hist.items())),
                "min_size": min_size,
                "min_count": min_count,
                "min_multiplicity_profiles": {str(k): v for k, v in profile_hist.most_common(4)},
                "sampled_min_distance_stats": {str(k): v for k, v in dist_hist.most_common(6)},
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Cross-lane Tournament Analysis


@dataclass(frozen=True)
class Lane:
    slug: str
    metric: tuple[int, int, int, int]
    note: str


def lane_tournament(lanes: list[Lane]) -> dict[str, object]:
    n = len(lanes)
    wins = [0] * n
    edges: set[tuple[int, int]] = set()
    for i, j in combinations(range(n), 2):
        if lanes[i].metric > lanes[j].metric:
            winner = i
        elif lanes[j].metric > lanes[i].metric:
            winner = j
        else:
            winner = min(i, j)
        loser = j if winner == i else i
        wins[winner] += 1
        edges.add((winner, loser))

    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (i, j) in edges and (j, k) in edges and (k, i) in edges:
            c3 += 1
        if (j, i) in edges and (k, j) in edges and (i, k) in edges:
            c3 += 1

    order = sorted(range(n), key=lambda i: (-wins[i], lanes[i].slug))
    hpaths = 0
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edges for i in range(n - 1)):
            hpaths += 1
    return {
        "top_order": [lanes[i].slug for i in order],
        "score_hist": dict(sorted(Counter(wins).items())),
        "directed_3cycles": c3,
        "hamiltonian_paths": hpaths,
        "edges": [(lanes[i].slug, lanes[j].slug) for i, j in sorted(edges)],
    }


def main() -> None:
    print("S658 frontier small-carrier atlases")
    print("=" * 72)

    print("\nA. Tournament reconstruction deck atlas")
    deck_rows = tournament_deck_atlas(6)
    for row in deck_rows:
        print(
            f"  n={row['n']} classes={row['classes']} "
            f"scalar={row['scalar_bucket_stats']} "
            f"full_deck={row['full_deck_bucket_stats']} "
            f"scalar_deck={row['scalar_deck_bucket_stats']} "
            f"loss_deck={row['loss_deck_bucket_stats']}"
        )
        print(f"    worst_scalar_bucket={row['worst_scalar_bucket']}")
        print(f"    worst_full_deck_bucket={row['worst_full_deck_bucket']}")
        print(f"    worst_scalar_deck_bucket={row['worst_scalar_deck_bucket']}")

    print("\nB. Union-closed frequency atlas")
    uc_rows = union_closed_atlas(4)
    for row in uc_rows:
        print(
            f"  m={row['m']} canonical_families={row['canonical_families']} "
            f"failures={row['frankl_failures']} tight={row['tight_half_frequency']} "
            f"min_ratio={row['min_max_frequency_ratio']}"
        )
        print(f"    witness_count_hist={row['witness_count_hist']}")
        print(f"    freq_buckets={row['freq_bucket_stats']} freq+size={row['freq_plus_size_bucket_stats']} freq+pressure={row['freq_plus_pressure_bucket_stats']}")
        print(f"    min_examples={row['min_examples']}")

    print("\nC. Finite-field Kakeya/Falconer toy atlas")
    kakeya_rows = kakeya_atlas((3, 5, 7))
    for row in kakeya_rows:
        print(
            f"  p={row['p']} choices={row['choices']} min_size={row['min_size']} "
            f"min_count={row['min_count']}"
        )
        print(f"    size_hist={row['size_hist']}")
        print(f"    min_multiplicity_profiles={row['min_multiplicity_profiles']}")
        print(f"    sampled_min_distance_stats={row['sampled_min_distance_stats']}")

    print("\nD. Cross-lane Tournament Analysis")
    loss_deck_exact = all(row["loss_deck_bucket_stats"][1] == 1 for row in deck_rows)
    deck_scalar_lift = max(row["scalar_bucket_stats"][1] - row["scalar_deck_bucket_stats"][1] for row in deck_rows)
    uc_loss = max(row["freq_bucket_stats"][1] for row in uc_rows)
    uc_pressure_gain = max(row["freq_bucket_stats"][1] - row["freq_plus_pressure_bucket_stats"][1] for row in uc_rows)
    kakeya_exact = all(row["min_size"] > row["p"] for row in kakeya_rows)
    kakeya_dist_side = int(all(row["sampled_min_distance_stats"] for row in kakeya_rows))
    lanes = [
        Lane(
            "tournament_decks",
            (5 if loss_deck_exact else 3, min(5, deck_scalar_lift), 5, 5),
            "raw vertex-decks collide through n=6, but the H/c3/deleted-score loss lift resolves all checked buckets",
        ),
        Lane(
            "union_closed_frequency",
            (5, min(5, uc_pressure_gain), 4 if uc_loss > 1 else 5, 5),
            "Frankl holds through m=4 but frequency vectors are visibly lossy",
        ),
        Lane(
            "finite_field_incidence",
            (4 if kakeya_exact else 3, 4, 5 if kakeya_dist_side else 4, 4),
            "Kakeya line-choice carriers expose multiplicity and distance side channels",
        ),
    ]
    for lane in lanes:
        print(f"  lane={lane.slug} metric={lane.metric} note={lane.note}")
    print(f"  tournament={lane_tournament(lanes)}")

    print("\nE. Practical synthesis")
    print("  Reconstruction: raw vertex-decks and H/score/c3/SCC scalars collide")
    print("  through n=6; scalar card-decks help, and the H/c3/deleted-score loss")
    print("  profile resolves every checked collision.  That loss profile is the")
    print("  first side-channel lift to test before broader graph decks.")
    print("  Union-closed: the half-frequency witness is easy at m<=4, but sorted")
    print("  frequency vectors forget many closure systems; set-pressure is a cheap")
    print("  reattachment channel.")
    print("  Finite-field incidence: line-choice Kakeya sets give a small exact")
    print("  direction carrier where multiplicity profiles and distance fibers can be")
    print("  tracked together before touching Euclidean Kakeya/Falconer questions.")


if __name__ == "__main__":
    main()
