#!/usr/bin/env python3
"""Standalone exact replay. Python 3 standard library only."""
from itertools import combinations, permutations
from math import comb
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path


def ensure(condition, message):
    if not condition:
        raise RuntimeError(message)


def arcs(n, mask):
    out = [0] * n
    for b, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (j, i) if (mask >> b) & 1 else (i, j)
        out[u] |= 1 << v
    return out


def hamilton(out):
    n = len(out)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v]:
                nxt = out[v] & ~mask
                while nxt:
                    bit = nxt & -nxt
                    nxt -= bit
                    dp[mask | bit][bit.bit_length() - 1] += dp[mask][v]
    return sum(dp[-1])


def potential_cycles(n):
    edge = {pair: bit for bit, pair in enumerate(combinations(range(n), 2))}
    answer = []
    for length in range(3, n + 1, 2):
        for vertices in combinations(range(n), length):
            for tail in permutations(vertices[1:]):
                word = (vertices[0],) + tail
                positive = negative = 0
                for u, v in zip(word, word[1:] + word[:1]):
                    bit = 1 << edge[tuple(sorted((u, v))) ]
                    if u > v:
                        positive |= bit
                    else:
                        negative |= bit
                support = sum(1 << v for v in vertices)
                answer.append((positive, negative, support, word))
    return answer


def cycles_from_mask(n, mask, potential):
    return [(support, word) for positive, negative, support, word in potential
            if (mask & positive) == positive and not (mask & negative)]


def small_packing(cycles):
    # Used only at n <= 6, when no three nontrivial cycles can be disjoint.
    supports = [s for s, _ in cycles]
    compatible = [(i, j) for i, j in combinations(range(len(cycles)), 2)
                  if not supports[i] & supports[j]]
    c = [1, len(cycles), len(compatible)]
    while len(c) > 1 and not c[-1]:
        c.pop()
    return c, compatible


def graph_betti(vertices, edges):
    if not vertices:
        return [0, 0]
    parent = list(range(vertices))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in edges:
        parent[root(i)] = root(j)
    components = len({root(i) for i in range(vertices)})
    # Each compatible pair becomes two distinct parallel edges.
    return [components, 2 * len(edges) - vertices + components]


def join(*blocks):
    result = []
    offset = 0
    n = sum(map(len, blocks))
    for block in blocks:
        after = ((1 << n) - 1) ^ ((1 << (offset + len(block))) - 1)
        result.extend((out << offset) | after for out in block)
        offset += len(block)
    return result


def reverse_outer(n):
    out = arcs(n, 0)
    out[0] ^= 1 << (n - 1)
    out[-1] |= 1
    return out


def encode(out):
    return sum(1 << b for b, (i, j) in enumerate(combinations(range(len(out)), 2))
               if (out[j] >> i) & 1)


def gf2_rank(columns):
    basis = {}
    for value in columns:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def cells_from_supports(supports):
    # Every face is a compatible cycle subset, signed modulo all-sign reversal.
    faces = [[]]
    for r in range(1, len(supports) + 1):
        level = []
        for face in combinations(range(len(supports)), r):
            occupied = 0
            for i in face:
                if occupied & supports[i]:
                    break
                occupied |= supports[i]
            else:
                for signs in range(1 << (r - 1)):
                    level.append((face, signs << 1))
        if not level:
            break
        faces.append(level)
    return faces


def boundary_gf2(faces):
    ranks = [0] * len(faces)
    for r in range(2, len(faces)):
        previous = {cell: i for i, cell in enumerate(faces[r-1])}
        columns = []
        for face, signs in faces[r]:
            column = 0
            for j in range(r):
                reduced = face[:j] + face[j+1:]
                newsigns = (signs & ((1 << j) - 1)) | ((signs >> (j+1)) << j)
                if newsigns & 1:
                    newsigns ^= (1 << (r - 1)) - 1
                column ^= 1 << previous[(reduced, newsigns)]
            columns.append(column)
        ranks[r] = gf2_rank(columns)
    betti = []
    for r in range(1, len(faces)):
        incoming = ranks[r+1] if r+1 < len(ranks) else 0
        betti.append(len(faces[r]) - ranks[r] - incoming)
    return betti


def main():
    report = {"scope": "all labelled tournaments through n=6; separate exact SCC controls", "enumeration": []}
    digest = sha256()
    for n in range(7):
        potential = potential_cycles(n)
        spectrum = set()
        gradings = defaultdict(set)
        total = 1 << comb(n, 2)
        for mask in range(total):
            h = hamilton(arcs(n, mask))
            cycles = cycles_from_mask(n, mask, potential)
            c, pairs = small_packing(cycles)
            ensure(sum(value * (1 << r) for r, value in enumerate(c)) == h,
                   f"OCF failed at n={n}, mask={mask}")
            k = sum(value * (1 << (r-1)) for r, value in enumerate(c) if r)
            ensure(2 * k + 1 == h, "cell total mismatch")
            spectrum.add(h)
            gradings[h].add(tuple(c))
            digest.update(json.dumps([n, mask, h, c], separators=(",", ":")).encode() + b"\n")
        report["enumeration"].append({"n": n, "labelled_tournaments": total,
            "distinct_H": len(spectrum), "distinct_packing_polynomials": len({c for cs in gradings.values() for c in cs}),
            "H_with_multiple_gradings": sorted(h for h, cs in gradings.items() if len(cs) > 1)})
    report["exhaustive_rows_sha256"] = digest.hexdigest()
    c3 = reverse_outer(3)
    s4 = reverse_outer(4)
    a = join(reverse_outer(5), [0])
    b = join(c3, c3)
    controls = []
    for name, out in [("strong5_with_sink", a), ("two_triangle_SCCs", b)]:
        mask = encode(out)
        cyc = cycles_from_mask(6, mask, potential_cycles(6))
        c, pairs = small_packing(cyc)
        f = [2**(r-1)*c[r] for r in range(1, len(c))]
        beta = graph_betti(len(cyc), pairs)
        controls.append({"name": name, "n": 6, "mask": mask,
                         "H": hamilton(out), "packing_coefficients": c,
                         "cell_counts": f, "betti_Q_and_F2": beta,
                         "cycles": [list(word) for _, word in cyc]})
    ensure([x["H"] for x in controls] == [9, 9], "hostile counts")
    ensure([x["packing_coefficients"] for x in controls] == [[1, 4], [1, 2, 1]], "hostile grading")
    ensure([x["betti_Q_and_F2"] for x in controls] == [[4, 0], [1, 1]], "hostile topology")
    report["same_order_same_H_different_topology"] = controls
    grammar = []
    for r in range(1, 7):
        supports = [1 << i for i in range(r)]
        cells = cells_from_supports(supports)
        expected = [comb(r, j) * 2**(j-1) for j in range(1, r+1)]
        ensure([len(x) for x in cells[1:]] == expected, "RP face vector")
        beta = boundary_gf2(cells)
        ensure(beta == [1]*r, "RP mod-two homology")
        if r <= 3:
            ensure(hamilton(join(*([c3]*r))) == 3**r, "triangle join H")
        grammar.append({"triangle_SCCs": r, "S4_centers": 0,
                        "cell_counts": expected, "F2_betti": beta,
                        "cell_total": (3**r-1)//2})
    for r in range(0, 6):
        # The two S4 cycles have the same abstract resource; every C3 has its own.
        supports = [1, 1] + [1 << (i+1) for i in range(r)]
        cells = cells_from_supports(supports)
        coeff = [comb(r, j) + 2*comb(r, j-1) if 1 <= j <= r else 2 if j == r+1 else 0
                 for j in range(1, r+2)]
        expected = [2**(j-1)*coeff[j-1] for j in range(1, r+2)]
        ensure([len(x) for x in cells[1:]] == expected, "mixed face vector")
        beta = boundary_gf2(cells)
        ensure(beta == [1]*r + [2], "mixed mod-two homology")
        if r <= 3:
            ensure(hamilton(join(*([c3]*r), s4)) == 5*3**r, "mixed join H")
        grammar.append({"triangle_SCCs": r, "S4_centers": 1,
                        "cell_counts": expected, "F2_betti": beta,
                        "cell_total": (5*3**r-1)//2})
    report["grammar_controls"] = grammar
    # Hostile inverse controls: an unrestricted discrete K with m vertices gives m cells.
    report["unrestricted_complex_controls"] = []
    for m in (3, 10):
        cells = cells_from_supports([1]*m)
        ensure(sum(map(len, cells)) == m, "unrestricted inverse hostile")
        report["unrestricted_complex_controls"].append({"vertices": m, "cell_total": m})
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    Path(__file__).with_name("RESULT.json").write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
