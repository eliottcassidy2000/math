#!/usr/bin/env python3
"""Exact deficit-two test from the complete SC frame catalog, not the tiling cube.

Inputs are the independently audited pairing_n{3,5,7,9}.json files.
Hamiltonian paths are counted by subset DP; possible automorphism orders
are filtered BEFORE graph-isomorphism backtracking is used.
"""
import argparse
from collections import Counter
from itertools import permutations
import hashlib
import json
import math
from pathlib import Path
import time


def adjacency(n, mask):
    k = n // 2
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    a = [0] * n
    def edge(i, j, bit):
        if bit:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    for i in range(k):
        edge(i, k, 1)
        edge(k, k + 1 + i, 1)
        edge(i, k + 1 + i, (mask >> i) & 1)
    for z, (i, j) in enumerate(pairs):
        bit = (mask >> (k + z)) & 1
        edge(i, j, bit)
        edge(k + 1 + j, k + 1 + i, bit)
        bit = (mask >> (k + len(pairs) + z)) & 1
        edge(i, k + 1 + j, bit)
        edge(j, k + 1 + i, bit)
    return a


def paths(a):
    n = len(a)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            out = a[v] & (full ^ mask)
            while out:
                bit = out & -out
                out -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[full])


def strong(a):
    full = (1 << len(a)) - 1
    for start in range(len(a)):
        seen = front = 1 << start
        while front:
            bit = front & -front
            front -= bit
            nxt = a[bit.bit_length() - 1] & ~seen
            seen |= nxt
            front |= nxt
        if seen != full:
            return False
    return True


def refinement(a):
    """Invariant equitable partition, so Aut divides the product of cell factorials."""
    colors = [x.bit_count() for x in a]
    while True:
        labels = sorted(set(colors))
        keys = [(colors[v], tuple(sum((a[v] >> w) & 1
                    for w in range(len(a)) if colors[w] == c)
                    for c in labels)) for v in range(len(a))]
        lookup = {x: i for i, x in enumerate(sorted(set(keys)))}
        new = [lookup[x] for x in keys]
        if len(set(new)) == len(set(colors)):
            return new
        colors = new


def automorphisms(a, colors):
    n = len(a)
    order = sorted(range(n), key=lambda v: (colors.count(colors[v]), v))
    candidates = {v: [w for w in range(n) if colors[w] == colors[v]] for v in order}
    mapping = [-1] * n
    count = 0
    def visit(depth, used):
        nonlocal count
        if depth == n:
            count += 1
            return
        v = order[depth]
        for w in candidates[v]:
            if (used >> w) & 1:
                continue
            if any(((a[v] >> u) & 1) != ((a[w] >> mapping[u]) & 1)
                   for u in order[:depth]):
                continue
            mapping[v] = w
            visit(depth + 1, used | (1 << w))
            mapping[v] = -1
    visit(0, 0)
    return count


def run(source, destination):
    results = []
    for n in (3, 5, 7, 9):
        start = time.monotonic()
        path = source / f"pairing_n{n}.json"
        data = json.loads(path.read_text())
        filters = Counter()
        tested = []
        accepted = []
        rows = []
        for row in data["rows"]:
            a = adjacency(n, row["mask"])
            H = paths(a)
            B = row["blue"]
            assert math.comb(n, 3) - sum(math.comb(x.bit_count(), 2) for x in a) == row["c3"]
            if n <= 5:
                assert H == sum(all((a[u] >> v) & 1 for u,v in zip(p,p[1:]))
                                for p in permutations(range(n)))
            is_strong = strong(a)
            candidate = H // (B + 2) if H % (B + 2) == 0 else None
            record = dict(id=row["id"], mask=row["mask"], H=H, B=B,
                          strong=is_strong, required_aut=candidate)
            rows.append(record)
            filters["all"] += 1
            if is_strong:
                filters["strong"] += 1
            if candidate is None:
                continue
            filters["integral_required_aut"] += 1
            if candidate % 2 == 0 or math.factorial(n) % candidate:
                continue
            filters["odd_divisor_n_factorial"] += 1
            if candidate % row["centralizer"]:
                continue
            filters["centralizer_divides"] += 1
            colors = refinement(a)
            bound = math.prod(math.factorial(s) for s in Counter(colors).values())
            if bound % candidate:
                continue
            filters["refined_cell_group_divides"] += 1
            aut = automorphisms(a, colors)
            record["aut"] = aut
            record["mass"] = H // aut
            record["deficit"] = H // aut - B
            assert H % aut == 0
            tested.append(record)
            if aut == candidate:
                accepted.append(record)
        result = dict(n=n, source=str(path), source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      classes=len(data["rows"]), filters=dict(filters),
                      automorphism_checks=tested, deficit_two=accepted,
                      elapsed_seconds=round(time.monotonic()-start, 3), rows=rows)
        results.append(result)
        print(json.dumps({k: v for k,v in result.items() if k != "rows"}), flush=True)
    # Independent direct permutation controls, never applied to all order-nine classes.
    controls = []
    s4 = [2, 4, 9, 3]  # 0->1->2->3->0, 2->0, 3->1.
    c5 = [sum(1 << ((v+d) % 5) for d in (1,2)) for v in range(5)]
    for name, a, B in (("S4", s4, 3), ("regular_C5", c5, 3)):
        n=len(a)
        H=sum(all((a[u]>>v)&1 for u,v in zip(p,p[1:])) for p in permutations(range(n)))
        aut=sum(all(((a[u]>>v)&1)==((a[p[u]]>>p[v])&1)
                    for u in range(n) for v in range(n) if u!=v)
                for p in permutations(range(n)))
        assert H == paths(a)
        assert aut == automorphisms(a, refinement(a))
        controls.append(dict(name=name, adjacency=a, H=H, aut=aut, mass=H//aut,
                             B=B, deficit=H//aut-B, strong=strong(a)))
    assert controls[0]["deficit"] == 2
    assert controls[1]["deficit"] == 0
    destination.write_text(json.dumps(dict(results=results, controls=controls), indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    here=Path(__file__).resolve().parent
    p.add_argument("--source", type=Path, default=here.parent / "odd_blue_spectrum_hollow_cubes_20261020")
    p.add_argument("--output", type=Path, default=here / "verification.json")
    args=p.parse_args()
    run(args.source, args.output)
