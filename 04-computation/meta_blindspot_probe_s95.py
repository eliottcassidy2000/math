#!/usr/bin/env python3
"""
meta_blindspot_probe_s95.py

Small measurements for meta-level hypothesis discovery.

This is intentionally narrow:
1. Scan repo prose for topic pairs that are individually common but rarely
   co-mentioned.
2. Measure the relation between Hamiltonian path count H(T) and minimum
   feedback arc set size on all labeled tournaments at n=6, plus a sample
   at n=7.
3. Classify Royle-even graph classes for n<=7 and test acyclic orientations
   as a possible "dark H" analog.

The goal is not to prove a theorem. It is to find places where the repo's
collection of hypotheses has underexplored joints.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb, sqrt
from pathlib import Path
import random
import statistics
import sys

import networkx as nx


ROOT = Path(__file__).resolve().parents[1]
sys.stdout.reconfigure(line_buffering=True)


TOPICS = {
    "dark_even_graph": ["dark tournament", "even graph", "odd graph", "sign representation", "royle"],
    "feedback_arc": ["feedback arc", "min-fas", "minimum feedback", "transitive distance"],
    "krawtchouk": ["krawtchouk", "band-limited", "bandlimited", "hamming scheme"],
    "cartan_attention": ["cartan", "attention", "dark mode", "gl(4"],
    "path_homology": ["path homology", "glmy", "beta_", "betti"],
    "paley_code": ["paley", "quadratic residue", "golay", "qr code"],
    "vitali_gap": ["vitali", "forbidden h", "h=7", "h=21", "permanent gap"],
    "staircase_tiling": ["staircase", "tiling", "waggly", "blueself"],
    "burnside": ["burnside", "cycle type", "all-odd", "automorphism"],
    "chromatic": ["chromatic", "coloring", "tutte polynomial", "acyclic orientation"],
    "magnitude_reachability": ["magnitude", "reachability", "spectral sequence"],
    "omega_conflict": ["omega", "conflict graph", "independence polynomial"],
    "ising_wick": ["ising", "wick", "yang-lee", "partition function"],
}


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def docs_to_scan():
    bases = ["00-navigation", "01-canon", "03-artifacts", "05-knowledge", "06-writeups", "07-reflections"]
    for base in bases:
        for path in (ROOT / base).rglob("*"):
            if path.suffix in {".md", ".tex", ".txt"}:
                yield path


def topic_hits(text: str) -> set[str]:
    low = text.lower()
    hits = set()
    for topic, needles in TOPICS.items():
        if any(n.lower() in low for n in needles):
            hits.add(topic)
    return hits


def scan_topic_blindspots() -> None:
    banner("1. Prose Co-Occurrence Blindspots")
    doc_hits = {}
    for path in docs_to_scan():
        try:
            hits = topic_hits(path.read_text(errors="ignore"))
        except OSError:
            continue
        if hits:
            doc_hits[path] = hits

    freq = Counter()
    pair = Counter()
    for hits in doc_hits.values():
        freq.update(hits)
        for a, b in combinations(sorted(hits), 2):
            pair[(a, b)] += 1

    print("Topic document counts:")
    for topic, count in freq.most_common():
        print(f"  {topic:22s} {count:4d}")

    candidates = []
    topics = sorted(freq)
    for a, b in combinations(topics, 2):
        if freq[a] < 4 or freq[b] < 4:
            continue
        c = pair[(a, b)]
        expected = freq[a] * freq[b] / max(1, len(doc_hits))
        # High score means both are common but surprisingly separate.
        score = (freq[a] + freq[b]) / (1 + c)
        if c <= max(1, expected * 0.35):
            candidates.append((score, a, b, c, expected))

    print("\nUnder-bridged common pairs:")
    for score, a, b, c, expected in sorted(candidates, reverse=True)[:12]:
        print(f"  {a:22s} + {b:22s} co={c:2d}, expected~{expected:4.1f}, score={score:5.1f}")

    # Greedy triads: pick a third topic with low pairwise overlap to both.
    print("\nSuggested triads to inspect:")
    shown = 0
    for _, a, b, _, _ in sorted(candidates, reverse=True):
        thirds = []
        for t in topics:
            if t in {a, b} or freq[t] < 4:
                continue
            overlap = pair[tuple(sorted((a, t)))] + pair[tuple(sorted((b, t)))]
            mass = freq[a] + freq[b] + freq[t]
            thirds.append((overlap, -mass, t))
        if thirds:
            overlap, neg_mass, t = min(thirds)
            print(f"  ({a}, {b}, {t}) -- pair-overlap-to-third={overlap}, total-docs={-neg_mass}")
            shown += 1
            if shown >= 6:
                break


def adj_from_bits(n: int, bits: int):
    adj = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return adj


def hamiltonian_paths(adj) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v, val in enumerate(row):
            if val == 0:
                continue
            av = adj[v]
            remaining = ((1 << n) - 1) ^ mask
            u = 0
            while remaining:
                if remaining & 1 and av[u]:
                    dp[mask | (1 << u)][u] += val
                remaining >>= 1
                u += 1
    return sum(dp[-1])


def min_feedback_arc_set(adj) -> int:
    """Minimum reversals to make a tournament transitive."""
    n = len(adj)
    m = n * (n - 1) // 2
    dp = [0] * (1 << n)
    for mask in range(1, 1 << n):
        best = -1
        x = mask
        while x:
            lsb = x & -x
            v = lsb.bit_length() - 1
            prev = mask ^ lsb
            gain = 0
            y = prev
            while y:
                lsb2 = y & -y
                u = lsb2.bit_length() - 1
                if adj[u][v]:
                    gain += 1
                y ^= lsb2
            cand = dp[prev] + gain
            if cand > best:
                best = cand
            x ^= lsb
        dp[mask] = best
    return m - dp[-1]


def corr(xs, ys) -> float:
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def tournament_h_fas_probe() -> None:
    banner("2. H(T) versus Minimum Feedback Arc Set")
    rng = random.Random(20260529)
    for n, mode, limit in [(6, "exhaustive", 1 << comb(6, 2)), (7, "sample", 5000)]:
        m = comb(n, 2)
        if mode == "exhaustive":
            bit_values = range(limit)
        else:
            bit_values = (rng.randrange(1 << m) for _ in range(limit))

        by_h = defaultdict(set)
        by_fas = defaultdict(list)
        hs = []
        fs = []
        score_variants = defaultdict(set)
        examples = {}

        for idx, bits in enumerate(bit_values):
            adj = adj_from_bits(n, bits)
            h = hamiltonian_paths(adj)
            f = min_feedback_arc_set(adj)
            hs.append(h)
            fs.append(f)
            by_h[h].add(f)
            by_fas[f].append(h)
            scores = tuple(sorted((sum(adj[v][u] for u in range(n) if u != v) for v in range(n))))
            score_variants[scores].add((h, f))
            if len(by_h[h]) > 1 and h not in examples:
                examples[h] = sorted(by_h[h])

        ambiguous_h = {h: sorted(v) for h, v in by_h.items() if len(v) > 1}
        multi_score = {s: v for s, v in score_variants.items() if len(v) > 1}

        print(f"\nn={n} ({mode}, {len(hs)} labeled tournaments)")
        print(f"  corr(H,FAS) = {corr(hs, fs):.4f}")
        print(f"  H range = {min(hs)}..{max(hs)}, FAS range = {min(fs)}..{max(fs)}")
        print(f"  distinct H = {len(by_h)}, H values with multiple FAS = {len(ambiguous_h)}")
        print(f"  score sequences with multiple (H,FAS) pairs = {len(multi_score)} / {len(score_variants)}")
        print("  H distribution by FAS:")
        for f in sorted(by_fas):
            vals = by_fas[f]
            print(f"    FAS={f}: count={len(vals):6d}, H_min={min(vals):5d}, H_med={statistics.median(vals):7.1f}, H_max={max(vals):7d}")
        if ambiguous_h:
            print("  first ambiguous H examples:")
            for h in sorted(ambiguous_h)[:6]:
                print(f"    H={h}: FAS values {ambiguous_h[h]}")


def graph_bits_from_graph(g: nx.Graph, nodes) -> int:
    bits = 0
    k = 0
    node_list = list(nodes)
    for i in range(len(node_list)):
        for j in range(i + 1, len(node_list)):
            if g.has_edge(node_list[i], node_list[j]):
                bits |= 1 << k
            k += 1
    return bits


def edge_sign_under_perm(edges, perm):
    pos = {e: i for i, e in enumerate(edges)}
    image = []
    for u, v in edges:
        a, b = perm[u], perm[v]
        if a > b:
            a, b = b, a
        image.append(pos[(a, b)])
    inv = 0
    for i in range(len(image)):
        for j in range(i + 1, len(image)):
            if image[i] > image[j]:
                inv += 1
    return -1 if inv % 2 else 1


def is_royle_even(g: nx.Graph) -> bool:
    nodes = list(g.nodes())
    edges = sorted((min(u, v), max(u, v)) for u, v in g.edges())
    if not edges:
        return True
    gm = nx.algorithms.isomorphism.GraphMatcher(g, g)
    for mapping in gm.isomorphisms_iter():
        perm = {u: mapping[u] for u in nodes}
        if edge_sign_under_perm(edges, perm) < 0:
            return False
    return True


def acyclic_orientations_count(g: nx.Graph) -> int:
    x = nx.chromatic_polynomial(g)
    val = int(x.subs({x.free_symbols.pop(): -1})) if x.free_symbols else int(x)
    return abs(((-1) ** g.number_of_nodes()) * val)


def royle_dark_probe() -> None:
    banner("3. Royle-Even Graphs and a Candidate Dark-H Analog")
    atlas = nx.graph_atlas_g()
    by_n = defaultdict(list)
    for g in atlas:
        n = g.number_of_nodes()
        if 1 <= n <= 6:
            by_n[n].append(nx.convert_node_labels_to_integers(g))

    for n in range(1, 7):
        even = []
        dark = []
        for g in by_n[n]:
            record = (
                g.number_of_edges(),
                tuple(sorted(dict(g.degree()).values())),
                acyclic_orientations_count(g),
            )
            if is_royle_even(g):
                even.append(record)
            else:
                dark.append(record)

        print(f"\nn={n}: graph classes={len(by_n[n])}, Royle-even={len(even)}, dark={len(dark)}")
        if even:
            vals = [r[2] for r in even]
            print(f"  even acyclic orientations: min={min(vals)}, median={statistics.median(vals):.1f}, max={max(vals)}")
        if dark:
            vals = [r[2] for r in dark]
            print(f"  dark acyclic orientations: min={min(vals)}, median={statistics.median(vals):.1f}, max={max(vals)}")
            print("  smallest dark classes by (edges, degrees, acyclic_orientations):")
            for rec in sorted(dark)[:5]:
                print(f"    {rec}")

        if even and dark:
            even_vals = {r[2] for r in even}
            dark_vals = {r[2] for r in dark}
            overlap = sorted(even_vals & dark_vals)
            print(f"  acyclic-orientation value overlap count = {len(overlap)}")
            if overlap[:8]:
                print(f"  overlap examples = {overlap[:8]}")


def main() -> None:
    scan_topic_blindspots()
    tournament_h_fas_probe()
    royle_dark_probe()

    banner("4. Candidate Meta-Hypotheses")
    print("HYP-A: H and min-FAS are not functionally equivalent, but they share a monotone")
    print("       low-frequency component. FAS is the order parameter; H is the entropy")
    print("       inside each order stratum.")
    print("HYP-B: The graph-side 'dark H' should not be acyclic orientations alone:")
    print("       that statistic overlaps heavily between Royle-even and dark classes.")
    print("       A better analog should include the sign representation, e.g. a signed")
    print("       acyclic-orientation count weighted by automorphism edge sign.")
    print("HYP-C: The blind triads involving dark/even graphs, Krawtchouk/band limits,")
    print("       and feedback arc distance suggest a missing theory of the cycle-space")
    print("       quotient where distance-to-transitive and sign representation are dual")
    print("       projections of the same GF(2) object.")


if __name__ == "__main__":
    main()
