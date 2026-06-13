#!/usr/bin/env python3
"""
metagraph_ultrafilter_divisor210_s670.py

S670 audit for the user's "tournament metagraph is an ultrafilter" picture.

The finite Boolean lattice behind the divisor poset of 210 is exact: choosing
one prime gives a principal ultrafilter, and the complementary lower ideal is
the other side of every d <-> 210/d pair.

The tournament tiling cube has the same raw shape: choose one tile coordinate,
take the upper half cube, and complement-tiling XORs every tile bit.  The key
question is whether this side choice descends through the tournament
isomorphism quotient.  The computation below measures the quotient leak.
"""

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, gcd


def banner(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ---------------------------------------------------------------------------
# Divisor lattice of 210 = Boolean lattice B_4.
# ---------------------------------------------------------------------------


PRIMES_210 = (2, 3, 5, 7)
N210 = 210


def divisors_of(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def support_210(d):
    return tuple(p for p in PRIMES_210 if d % p == 0)


def is_filter_210(subset, divisors):
    if not subset or 1 in subset:
        return False
    subset = set(subset)
    for d in subset:
        for e in divisors:
            if e % d == 0 and e not in subset:
                return False
    for a in subset:
        for b in subset:
            if gcd(a, b) not in subset:
                return False
    return True


def is_ultrafilter_210(subset, divisors):
    subset = set(subset)
    if not is_filter_210(subset, divisors):
        return False
    for d in divisors:
        comp = N210 // d
        if (d in subset) == (comp in subset):
            return False
    return True


def divisor_ultrafilter_audit():
    divs = divisors_of(N210)
    ultrafilters = []
    for mask in range(1 << len(divs)):
        subset = {divs[i] for i in range(len(divs)) if (mask >> i) & 1}
        if is_ultrafilter_210(subset, divs):
            ultrafilters.append(tuple(sorted(subset)))

    principal = {}
    for p in PRIMES_210:
        upper = tuple(d for d in divs if d % p == 0)
        lower = tuple(d for d in divs if d % p != 0)
        edges = Counter()
        for d in divs:
            for q in PRIMES_210:
                if d % q != 0:
                    e = d * q
                    if e <= N210 and N210 % e == 0:
                        a = "upper" if d in upper else "lower"
                        b = "upper" if e in upper else "lower"
                        edges[(a, b)] += 1
        principal[p] = {
            "upper": upper,
            "lower": lower,
            "hasse_edges": edges,
        }
    return divs, ultrafilters, principal


# ---------------------------------------------------------------------------
# Fixed-base-path tiling cube for tournaments.
# ---------------------------------------------------------------------------


def tile_pairs(n):
    base_edges = {(i, i + 1) for i in range(n - 1)}
    return [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if (i, j) not in base_edges
    ]


def build_tournament(n, bits, pairs):
    a = [[0] * n for _ in range(n)]
    # Fixed base path: n-1 -> n-2 -> ... -> 1 -> 0.
    for i in range(1, n):
        a[i][i - 1] = 1
    for k, (lo, hi) in enumerate(pairs):
        if (bits >> k) & 1:
            a[lo][hi] = 1
        else:
            a[hi][lo] = 1
    return a


@lru_cache(maxsize=None)
def grouped_permutations(score_tuple):
    groups = defaultdict(list)
    for v, s in enumerate(score_tuple):
        groups[s].append(v)
    blocks = [groups[s] for s in sorted(groups)]

    def rec(rest):
        if not rest:
            yield ()
            return
        head, *tail = rest
        for p in permutations(head):
            for q in rec(tail):
                yield p + q

    return tuple(rec(blocks))


def canonical_tournament(a):
    n = len(a)
    scores = tuple(sum(row) for row in a)
    best = None
    for p in grouped_permutations(scores):
        form = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best


def hamiltonian_path_count(a):
    n = len(a)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for s in range(1, 1 << n):
        for v in range(n):
            cur = dp.get((s, v), 0)
            if cur == 0:
                continue
            for w in range(n):
                if s & (1 << w):
                    continue
                if a[v][w]:
                    dp[(s | (1 << w), w)] = dp.get((s | (1 << w), w), 0) + cur
    return sum(dp.get((full, v), 0) for v in range(n))


def build_tiling_space(n):
    pairs = tile_pairs(n)
    m = len(pairs)
    total = 1 << m
    class_of_bits = {}
    class_id = {}
    members = defaultdict(list)
    reps = {}

    for bits in range(total):
        a = build_tournament(n, bits, pairs)
        c = canonical_tournament(a)
        if c not in class_id:
            cid = len(class_id)
            class_id[c] = cid
            reps[cid] = bits
        cid = class_id[c]
        class_of_bits[bits] = cid
        members[cid].append(bits)

    info = {}
    mask = total - 1
    for cid, bits in reps.items():
        a = build_tournament(n, bits, pairs)
        comp_bits = bits ^ mask
        info[cid] = {
            "H": hamiltonian_path_count(a),
            "size": len(members[cid]),
            "comp_tiling_class": class_of_bits[comp_bits],
        }

    return {
        "n": n,
        "pairs": pairs,
        "m": m,
        "total": total,
        "mask": mask,
        "class_of_bits": class_of_bits,
        "members": members,
        "info": info,
        "num_classes": len(class_id),
    }


def color_summary(ts, k):
    summary = {}
    for cid, bits_list in ts["members"].items():
        colors = {(bits >> k) & 1 for bits in bits_list}
        summary[cid] = tuple(sorted(colors))
    return summary


def audit_coordinate(ts, k):
    summary = color_summary(ts, k)
    pure0 = sum(1 for colors in summary.values() if colors == (0,))
    pure1 = sum(1 for colors in summary.values() if colors == (1,))
    mixed = sum(1 for colors in summary.values() if colors == (0, 1))
    mixed_tilings = sum(
        len(ts["members"][cid]) for cid, colors in summary.items() if colors == (0, 1)
    )

    pair_stats = Counter()
    mask = ts["mask"]
    for bits in range(ts["total"]):
        comp_bits = bits ^ mask
        if bits > comp_bits:
            continue
        ci = ts["class_of_bits"][bits]
        cj = ts["class_of_bits"][comp_bits]
        if ci == cj:
            pair_stats["self_class_pair"] += 1
        else:
            pair_stats["cross_class_pair"] += 1

        colors_i = summary[ci]
        colors_j = summary[cj]
        if len(colors_i) == 1 and len(colors_j) == 1 and colors_i != colors_j:
            pair_stats["pure_blue_black_pair"] += 1
        else:
            pair_stats["quotient_leak_pair"] += 1

    return {
        "k": k,
        "tile": ts["pairs"][k],
        "pure0": pure0,
        "pure1": pure1,
        "mixed": mixed,
        "mixed_tilings": mixed_tilings,
        "pair_stats": pair_stats,
    }


def canonical_graph(n, bits, pairs):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            a[i][j] = 1
            a[j][i] = 1
    deg = tuple(sum(row) for row in a)
    best = None
    for p in grouped_permutations(deg):
        form = tuple(a[p[i]][p[j]] for i in range(n) for j in range(i + 1, n))
        if best is None or form < best:
            best = form
    return best


def degree_even_graph_iso_count(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    classes = set()
    labeled = 0
    for bits in range(1 << len(pairs)):
        deg = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                deg[i] += 1
                deg[j] += 1
        if all(d % 2 == 0 for d in deg):
            labeled += 1
            classes.add(canonical_graph(n, bits, pairs))
    return labeled, len(classes)


# ---------------------------------------------------------------------------
# Tournament Analysis over interpretations.
# ---------------------------------------------------------------------------


def count_directed_3cycles(adj):
    vertices = list(adj)
    total = 0
    for a, b, c in combinations(vertices, 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj):
    vertices = list(adj)

    def reach(start, graph):
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w in vertices:
                if graph[v][w] and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    rev = {v: {w: adj[w][v] for w in vertices} for v in vertices}
    remaining = set(vertices)
    sizes = []
    while remaining:
        v = next(iter(remaining))
        comp = reach(v, adj) & reach(v, rev)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths_in_digraph(adj):
    vertices = list(adj)
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp = {}
    for v in vertices:
        dp[(1 << idx[v], v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in vertices:
            cur = dp.get((mask, v), 0)
            if cur == 0:
                continue
            for w in vertices:
                bit = 1 << idx[w]
                if mask & bit:
                    continue
                if adj[v][w]:
                    dp[(mask | bit, w)] = dp.get((mask | bit, w), 0) + cur
    return sum(dp.get((full, v), 0) for v in vertices)


def interpretation_tournament():
    # Metrics: literal filter law, complement decision, quotient descent,
    # retained address profile, LRC actionability.
    channels = {
        "LRC_owner_filter_profile": (4, 4, 4, 5, 5),
        "tiling_cube_principal_ultrafilter": (5, 5, 3, 4, 4),
        "divisor210_principal_ultrafilter": (5, 5, 5, 2, 2),
        "iso_metagraph_filter_sheaf": (3, 4, 2, 4, 5),
        "degree_even_cycle_space": (3, 2, 2, 3, 3),
        "Royle_even_equinumerosity": (2, 1, 1, 2, 3),
        "raw_cardinality": (1, 0, 0, 0, 1),
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


def main():
    banner("S670 metagraph ultrafilter / divisor-210 audit")
    print("Thesis under test:")
    print("  Raw side-choice is an ultrafilter on a Boolean cube.")
    print("  The tournament metagraph is the quotient shadow; it needs an address")
    print("  coordinate or sheaf if iso classes mix upper/lower membership.")

    banner("A. Divisor lattice of 210")
    divs, ultrafilters, principal = divisor_ultrafilter_audit()
    print(f"divisors={divs}")
    print(f"support map: d -> subset of {PRIMES_210}")
    print(f"number of divisors={len(divs)} = 2^4")
    print(f"brute ultrafilters={len(ultrafilters)}")
    for p in PRIMES_210:
        data = principal[p]
        print(f"  p={p}: upper(size={len(data['upper'])})={data['upper']}")
        print(f"       lower(size={len(data['lower'])})={data['lower']}")
        print(f"       Hasse edge split={dict(data['hasse_edges'])}")
    print("Conclusion: the exact finite ultrafilters are principal prime choices.")
    print("Each decides one side of every complement pair d <-> 210/d.")

    banner("B. Tiling cube, even graph count shadow, and quotient leak")
    spaces = {}
    for n in (4, 5, 6):
        ts = build_tiling_space(n)
        spaces[n] = ts
        labeled_even, degree_even_iso = degree_even_graph_iso_count(n)
        print(f"n={n}")
        print(
            "  fixed-path tilings=2^C(n-1,2)="
            f"{ts['total']} ; labeled degree-even graphs={labeled_even}"
        )
        print(
            f"  tournament iso classes in tiling cube={ts['num_classes']} ; "
            f"degree-even graph iso classes={degree_even_iso}"
        )
        print("  coordinate audits (best first by mixed classes, then leaked pairs):")
        audits = [audit_coordinate(ts, k) for k in range(ts["m"])]
        audits.sort(
            key=lambda x: (
                x["mixed"],
                x["pair_stats"]["quotient_leak_pair"],
                x["mixed_tilings"],
            )
        )
        for a in audits[: min(5, len(audits))]:
            ps = a["pair_stats"]
            print(
                f"    k={a['k']} tile={a['tile']}: "
                f"pure0={a['pure0']} pure1={a['pure1']} mixed={a['mixed']} "
                f"mixed_tilings={a['mixed_tilings']} "
                f"pure_pairs={ps['pure_blue_black_pair']} "
                f"leak_pairs={ps['quotient_leak_pair']} "
                f"self_pairs={ps['self_class_pair']} cross_pairs={ps['cross_class_pair']}"
            )
        best = audits[0]
        if best["mixed"] == 0 and best["pair_stats"]["quotient_leak_pair"] == 0:
            print("  descent verdict: a principal cube ultrafilter descends cleanly.")
        else:
            print(
                "  descent verdict: quotient leaks; the metagraph needs the missing "
                "tile/address coordinate."
            )

    banner("C. Blue/black complement-tiling statement")
    for n, ts in spaces.items():
        pairs = ts["total"] // 2
        comp_cross = 0
        comp_self = 0
        for bits in range(ts["total"]):
            comp_bits = bits ^ ts["mask"]
            if bits > comp_bits:
                continue
            if ts["class_of_bits"][bits] == ts["class_of_bits"][comp_bits]:
                comp_self += 1
            else:
                comp_cross += 1
        print(
            f"n={n}: complement-tiling pairs={pairs}, "
            f"same iso class={comp_self}, different iso class={comp_cross}"
        )
    print("This uses MISTAKE-033's corrected rule: blue/black is tiling complement")
    print("inside Q_m, not tournament complement T^op.")

    banner("D. Tournament Analysis over interpretations")
    channels, adj, scores, order = interpretation_tournament()
    print("vertices=interpretations of the ultrafilter/equinumerosity analogy")
    print(
        "observable=(literal filter law, complement decision, quotient descent, "
        "retained address profile, LRC actionability)"
    )
    print("switch=majority; tie Hamiltonian path=listed priority order")
    print(f"score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"directed_3cycles={count_directed_3cycles(adj)}")
    print(f"scc_sizes={scc_sizes(adj)}")
    print(f"hamiltonian_paths={count_hamiltonian_paths_in_digraph(adj)}")
    print("top_order:")
    for name in order:
        print(f"  {name:36s} score={scores[name]} vector={channels[name]}")

    banner("E. Repo theorem shape")
    print("Literal theorem candidate:")
    print("  The divisor lattice of squarefree N and the tiling cube Q_m are Boolean")
    print("  algebras; principal atom filters decide complement pairs.")
    print("Refined metagraph claim:")
    print("  The tournament metagraph is not itself an ultrafilter in general.")
    print("  It is an observer quotient of many cube ultrafilters.  Mixed iso")
    print("  classes are exactly the proof-obligation defects.")
    print("LRC14 use:")
    print("  Replace raw runner/time vertices by proof obligations, residues,")
    print("  cover arcs, owner bits, and carry cocycles.  A usable ultrafilter")
    print("  is a choice-of-side that stays pure after the quotient; if it leaks,")
    print("  add the missing address coordinate and retest.")


if __name__ == "__main__":
    main()
