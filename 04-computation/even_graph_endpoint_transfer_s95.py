#!/usr/bin/env python3
"""
even_graph_endpoint_transfer_s95.py

Endpoint-transfer analogue for the even-graph quotient of the same fixed-path
tiling cube.

The tiling cube Q_m maps bijectively to labeled even subgraphs of K_n by XORing
fundamental cycles of the base path.  This script builds the quotient by
undirected graph isomorphism, then applies the same endpoint extension

    Q_{m(n+1)} = Q_{m(n)} x Q_{n-1}.

Unlike tournament fibers, even-graph class fibers are labeled orbit sizes
n!/|Aut(G)|, so their parity is not uniformly odd.  The mod-2 transfer boundary
is the set of child even-graph classes with odd orbit size.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import comb


def tile_pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]


def all_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def extend_bits(parent_n, parent_bits, signature):
    old = {}
    for idx, pair in enumerate(tile_pairs(parent_n)):
        old[pair] = (parent_bits >> idx) & 1

    child_bits = 0
    for idx, (i, j) in enumerate(tile_pairs(parent_n + 1)):
        if j < parent_n:
            bit = old[(i, j)]
        else:
            bit = (signature >> i) & 1
        if bit:
            child_bits |= 1 << idx
    return child_bits


def fundamental_cycles(n):
    edges = all_edges(n)
    edge_idx = {e: i for i, e in enumerate(edges)}
    cycles = []
    for i, j in tile_pairs(n):
        bits = 0
        for k in range(i, j):
            bits |= 1 << edge_idx[(k, k + 1)]
        bits |= 1 << edge_idx[(i, j)]
        cycles.append(bits)
    return edges, edge_idx, cycles


def tiling_to_even(bits, cycles):
    graph = 0
    for idx, cyc in enumerate(cycles):
        if (bits >> idx) & 1:
            graph ^= cyc
    return graph


def graph_degrees(n, graph_bits, edges):
    deg = [0] * n
    for idx, (i, j) in enumerate(edges):
        if (graph_bits >> idx) & 1:
            deg[i] += 1
            deg[j] += 1
    return deg


def graph_canonical(n, graph_bits, edges, edge_idx):
    deg = graph_degrees(n, graph_bits, edges)
    blocks = defaultdict(list)
    for v, d in enumerate(deg):
        blocks[d].append(v)
    groups = [blocks[d] for d in sorted(blocks)]
    best = None

    def rec(k):
        if k == len(groups):
            yield []
            return
        for p in permutations(groups[k]):
            for rest in rec(k + 1):
                yield list(p) + rest

    for order in rec(0):
        nb = 0
        for idx, (i, j) in enumerate(edges):
            a, b = sorted((order[i], order[j]))
            if (graph_bits >> edge_idx[(a, b)]) & 1:
                nb |= 1 << idx
        if best is None or nb < best:
            best = nb
    return best


def build_even_level(n):
    m = comb(n - 1, 2)
    edges, edge_idx, cycles = fundamental_cycles(n)
    class_map = {}
    class_info = []
    tiling_cid = []

    for bits in range(1 << m):
        graph = tiling_to_even(bits, cycles)
        can = graph_canonical(n, graph, edges, edge_idx)
        if can not in class_map:
            class_map[can] = len(class_info)
            class_info.append({
                "canon": can,
                "fiber": 0,
                "edges": can.bit_count(),
                "degree_seq": tuple(sorted(graph_degrees(n, can, edges))),
            })
        cid = class_map[can]
        class_info[cid]["fiber"] += 1
        tiling_cid.append(cid)

    return {
        "n": n,
        "m": m,
        "classes": class_info,
        "tiling_cid": tiling_cid,
    }


def gf2_rank(bit_rows):
    basis = {}
    for row in bit_rows:
        x = row
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def support_stats(rows, cols):
    row_support = [len(row) for row in rows]
    col_support = [0] * cols
    for row in rows:
        for j in row:
            col_support[j] += 1
    return {
        "support": sum(row_support),
        "row_min": min(row_support) if row_support else 0,
        "row_max": max(row_support) if row_support else 0,
        "row_avg": sum(row_support) / len(row_support) if row_support else 0,
        "col_min": min(col_support) if col_support else 0,
        "col_max": max(col_support) if col_support else 0,
        "col_avg": sum(col_support) / cols if cols else 0,
    }


def parity_summary(rows, cols, expected_col_parity):
    bit_rows = []
    xor_rows = 0
    col_parity = [0] * cols
    odd_entries = 0
    row_failures = 0

    for row in rows:
        bits = 0
        row_parity = 0
        for j, weight in row.items():
            if weight % 2:
                bits |= 1 << j
                col_parity[j] ^= 1
                row_parity ^= 1
                odd_entries += 1
        if row_parity:
            row_failures += 1
        bit_rows.append(bits)
        xor_rows ^= bits

    expected_bits = 0
    col_failures = 0
    for j, parity in enumerate(expected_col_parity):
        if parity:
            expected_bits |= 1 << j
        if col_parity[j] != parity:
            col_failures += 1

    return {
        "odd_entries": odd_entries,
        "rank": gf2_rank(bit_rows),
        "row_failures": row_failures,
        "col_failures": col_failures,
        "boundary_weight": sum(expected_col_parity),
        "xor_boundary_ok": xor_rows == expected_bits,
    }


def analyze_transfer(parent, child):
    n = parent["n"]
    rows = [Counter() for _ in parent["classes"]]

    for parent_bits, parent_cid in enumerate(parent["tiling_cid"]):
        for sig in range(1 << (n - 1)):
            child_bits = extend_bits(n, parent_bits, sig)
            child_cid = child["tiling_cid"][child_bits]
            rows[parent_cid][child_cid] += 1

    col_sums = [0] * len(child["classes"])
    for row in rows:
        for j, weight in row.items():
            col_sums[j] += weight

    row_errors = []
    for i, row in enumerate(rows):
        expected = (1 << (n - 1)) * parent["classes"][i]["fiber"]
        actual = sum(row.values())
        if actual != expected:
            row_errors.append((i, actual, expected))

    col_errors = []
    for j, actual in enumerate(col_sums):
        expected = child["classes"][j]["fiber"]
        if actual != expected:
            col_errors.append((j, actual, expected))

    expected_col_parity = [c["fiber"] % 2 for c in child["classes"]]
    return {
        "rows": rows,
        "row_errors": row_errors,
        "col_errors": col_errors,
        "support": support_stats(rows, len(child["classes"])),
        "parity": parity_summary(rows, len(child["classes"]), expected_col_parity),
        "weight_spectrum": Counter(w for row in rows for w in row.values()),
        "child_odd_fiber_spectrum": Counter(
            (c["edges"], c["degree_seq"]) for c in child["classes"] if c["fiber"] % 2
        ),
    }


def compact_counter(counter, limit=12):
    rows = sorted(counter.items())
    if len(rows) <= limit:
        return dict(rows)
    return {
        "head": rows[: limit // 2],
        "tail": rows[-limit // 2 :],
        "distinct": len(rows),
    }


def main():
    print("=" * 78)
    print("EVEN GRAPH ENDPOINT TRANSFER S95")
    print("=" * 78)
    print("Same endpoint cube recursion, but quotienting by even-graph isomorphism.")
    print("Column parity now means odd labeled orbit size, not Redei oddness.")

    levels = {n: build_even_level(n) for n in range(2, 8)}
    summaries = []
    for n in range(2, 7):
        parent = levels[n]
        child = levels[n + 1]
        d = analyze_transfer(parent, child)
        summaries.append((n, d))
        ss = d["support"]
        ps = d["parity"]

        print("\n" + "-" * 78)
        print(
            f"n={n}->{n+1}: parent even classes={len(parent['classes'])}, "
            f"child even classes={len(child['classes'])}, multiplier=2^{n-1}"
        )
        print(
            "  exact sums: "
            f"row/col errors={len(d['row_errors'])}/{len(d['col_errors'])}"
        )
        print(
            "  support: "
            f"nonzero={ss['support']}, row={ss['row_min']}..{ss['row_max']} "
            f"avg={ss['row_avg']:.2f}, col={ss['col_min']}..{ss['col_max']} "
            f"avg={ss['col_avg']:.2f}"
        )
        print(
            "  GF(2) transfer: "
            f"odd_entries={ps['odd_entries']}, rank={ps['rank']}, "
            f"row_fail={ps['row_failures']}, col_fail={ps['col_failures']}, "
            f"odd-orbit boundary={ps['boundary_weight']}, "
            f"xor_boundary_ok={ps['xor_boundary_ok']}"
        )
        print(f"  transfer weight spectrum={compact_counter(d['weight_spectrum'])}")
        print(f"  child odd-fiber class shapes={compact_counter(d['child_odd_fiber_spectrum'])}")

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    fields = [
        ("even class counts", lambda n, d: len(levels[n]["classes"])),
        ("support", lambda n, d: d["support"]["support"]),
        ("GF2 rank", lambda n, d: d["parity"]["rank"]),
        ("parent class count", lambda n, d: len(levels[n]["classes"])),
        ("odd entries", lambda n, d: d["parity"]["odd_entries"]),
        ("odd-orbit child boundary", lambda n, d: d["parity"]["boundary_weight"]),
    ]
    for name, fn in fields:
        print(f"  {name}: {[fn(n, d) for n, d in summaries]}")
    print(f"  child even class counts: {[len(levels[n + 1]['classes']) for n, _ in summaries]}")

    print("\n" + "=" * 78)
    print("CONSTRAINT STATEMENTS")
    print("=" * 78)
    print("""
1. The endpoint transfer row/column theorem is quotient-agnostic.
   Any quotient of the fixed-path tiling cube inherits a transfer matrix whose
   rows sum to 2^(n-1) times parent fiber and whose columns sum to child fiber.

2. The tournament quotient is special because every child fiber is odd.
   For even graphs, child fibers are orbit sizes n!/|Aut(G)| and can be even.
   The mod-2 boundary is therefore the set of even-graph classes whose
   automorphism group absorbs the full 2-adic valuation of n!.

3. Unlike the tournament quotient, the even-graph endpoint transfer does NOT
   have full row rank over GF(2) once n >= 3.  Recursive parity-injectivity is
   therefore not just a property of arbitrary tiling-cube quotients; it is a
   tournament-side feature that the coarser cycle-space quotient forgets.
""")


if __name__ == "__main__":
    main()
