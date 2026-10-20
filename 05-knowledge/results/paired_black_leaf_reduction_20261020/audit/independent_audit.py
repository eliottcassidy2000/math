#!/usr/bin/env python3
"""Independent SC9 leaf audit. No producer imports; exact Python integers."""
import argparse
import hashlib
import itertools as it
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / "odd_blue_spectrum_hollow_cubes_20261020"
PRODUCER = HERE.parent


def frame(n, code):
    k = n // 2
    A = [[0] * n for _ in range(n)]
    def set_arc(i, j, direction):
        A[i][j] = direction
        A[j][i] = 1 - direction
    positions = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for i in range(k):
        set_arc(i, k, 1)
        set_arc(k, i + k + 1, 1)
        set_arc(i, i + k + 1, (code >> i) & 1)
    for t, (i, j) in enumerate(positions):
        x = (code >> (k + t)) & 1
        set_arc(i, j, x)
        set_arc(i + k + 1, j + k + 1, 1 - x)
        y = (code >> (k + len(positions) + t)) & 1
        set_arc(i, j + k + 1, y)
        set_arc(j, i + k + 1, y)
    rho = list(range(k + 1, n)) + [k] + list(range(k))
    assert all(A[i][j] == A[rho[j]][rho[i]]
               for i in range(n) for j in range(n))
    return A, rho


def encode_reindexed_frame(A, p):
    k = len(p)
    values = [A[p[i]][p[i] + k + 1] for i in range(k)]
    values += [A[p[i]][p[j]] for i in range(k) for j in range(i + 1, k)]
    values += [A[p[i]][p[j] + k + 1]
               for i in range(k) for j in range(i + 1, k)]
    return sum(x << j for j, x in enumerate(values))


def count_paths(A):
    """Top-down backwards recurrence, unlike the producer's forward DP."""
    n = len(A)
    @lru_cache(None)
    def ending_at(vertices, last):
        prev = vertices ^ (1 << last)
        if not prev:
            return 1
        return sum(ending_at(prev, j) for j in range(n)
                   if (prev >> j) & 1 and A[j][last])
    return sum(ending_at((1 << n) - 1, j) for j in range(n))


def symmetric_words(A, rho):
    pairs = [(i, rho[i]) for i in range(len(A)) if i < rho[i]]
    fixed = [i for i in range(len(A)) if i == rho[i]]
    count = 0
    for p in it.permutations(range(len(pairs))):
        for choices in it.product((0, 1), repeat=len(pairs)):
            left = [pairs[j][c] for j, c in zip(p, choices)]
            word = left + fixed + [rho[j] for j in left[::-1]]
            count += all(A[x][y] for x, y in zip(word, word[1:]))
    return count


def is_strong(A):
    # Floyd-Warshall transitive closure, independently of traversal code.
    reach = [[bool(x) for x in row] for row in A]
    for k in range(len(A)):
        for i in range(len(A)):
            for j in range(len(A)):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    return all(reach[i][j] for i in range(len(A)) for j in range(len(A))
               if i != j)


def invariant_cells(A):
    labels = [sum(row) for row in A]
    while True:
        signatures = [(labels[i], tuple(sorted(
            (labels[j], A[i][j]) for j in range(len(A)) if j != i)))
            for i in range(len(A))]
        names = {s: k for k, s in enumerate(sorted(set(signatures)))}
        new = [names[s] for s in signatures]
        if len(set(new)) == len(set(labels)):
            return [[i for i, c in enumerate(new) if c == d]
                    for d in sorted(set(new))]
        labels = new


def restricted_automorphisms(A, cells=None):
    if cells is None:
        cells = [[i for i in range(len(A)) if sum(A[i]) == d]
                 for d in sorted({sum(row) for row in A})]
    result = []
    edges = [(i, j) for i in range(len(A)) for j in range(i + 1, len(A))]
    for images in it.product(*(it.permutations(c) for c in cells)):
        p = [None] * len(A)
        for c, image in zip(cells, images):
            for x, y in zip(c, image):
                p[x] = y
        if all(A[i][j] == A[p[i]][p[j]] for i, j in edges):
            result.append(p)
    return result


def exhaustive_control(name, A):
    n = len(A)
    paths = [p for p in it.permutations(range(n))
             if all(A[i][j] for i, j in zip(p, p[1:]))]
    aut = len(restricted_automorphisms(A))
    blue_paths = 0
    closing = 0
    cycle_types = Counter()
    for p in paths:
        r = [None] * n
        for x, y in zip(p, p[::-1]):
            r[x] = y
        is_blue = all(A[i][j] == A[r[j]][r[i]]
                      for i in range(n) for j in range(n))
        blue_paths += is_blue
        if A[p[-1]][p[0]]:
            closing += 1
            if p[0] == 0:
                rotations = []
                blue_cuts = 0
                for shift in range(n):
                    c = dict(zip(p, p[shift:] + p[:shift]))
                    if all(A[i][j] == A[c[i]][c[j]]
                           for i in range(n) for j in range(n)):
                        rotations.append(shift)
                    cut = p[shift:] + p[:shift]
                    rev = dict(zip(cut, cut[::-1]))
                    blue_cuts += all(A[i][j] == A[rev[j]][rev[i]]
                                     for i in range(n) for j in range(n))
                d = len(rotations)
                q = n // d
                assert blue_cuts % d == 0
                blue_orbits = blue_cuts // d
                assert blue_orbits in (0, math.gcd(2, q))
                cycle_types[str((q, blue_orbits))] += 1
    assert len(paths) == count_paths(A)
    assert len(paths) % aut == blue_paths % aut == 0
    return dict(name=name, H=len(paths), aut=aut,
                blue_paths=blue_paths, B=blue_paths // aut,
                mass=len(paths) // aut,
                deficit=(len(paths) - blue_paths) // aut,
                strong=is_strong(A), HC=closing // n,
                cycle_types=dict(cycle_types))


def audit_catalog(n):
    source = SOURCE / f"pairing_n{n}.json"
    data = json.loads(source.read_text())
    k = n // 2
    pair_permutations = list(it.permutations(range(k)))
    covered = set()
    records = []
    counts = Counter()
    for row in data["rows"]:
        A, rho = frame(n, row["mask"])
        codes = [encode_reindexed_frame(A, p) for p in pair_permutations]
        orbit = set(codes)
        assert not (covered & orbit), ("duplicate frame class", n, row["id"])
        covered |= orbit
        centralizer = codes.count(row["mask"])
        assert centralizer * len(orbit) == math.factorial(k)
        h_rho = symmetric_words(A, rho)
        assert h_rho % centralizer == 0
        B = h_rho // centralizer
        assert (B, centralizer, h_rho) == (
            row["blue"], row["centralizer"], row["symmetric_HP"])
        c3 = sum(A[i][j] * A[j][t] * A[t][i]
                 + A[j][i] * A[t][j] * A[i][t]
                 for i in range(n) for j in range(i+1,n) for t in range(j+1,n))
        assert c3 == row["c3"]
        H = count_paths(A)
        record = dict(id=row["id"], mask=row["mask"], H=H, B=B,
                      strong=is_strong(A), required_aut=None)
        records.append(record)
        counts["all"] += 1
        counts["strong"] += record["strong"]
        if H % (B + 2):
            continue
        counts["integral_required_aut"] += 1
        candidate = H // (B + 2)
        record["required_aut"] = candidate
        if candidate % 2 == 0 or math.factorial(n) % candidate:
            continue
        counts["odd_divisor_n_factorial"] += 1
        if candidate % centralizer:
            continue
        counts["centralizer_divides"] += 1
        cells = invariant_cells(A)
        upper_group = math.prod(math.factorial(len(c)) for c in cells)
        if upper_group % candidate:
            continue
        counts["refined_cell_group_divides"] += 1
        # Intentionally use only outdegree classes, a coarser independent search.
        aut = len(restricted_automorphisms(A))
        record.update(aut=aut, mass=H//aut, deficit=H//aut-B)
        assert H % aut == 0
        assert aut != candidate
    assert covered == set(range(1 << (k*k)))
    # Compare producer output only after independent calculations have finished.
    produced = json.loads((PRODUCER / "verification.json").read_text())
    prod = next(x for x in produced["results"] if x["n"] == n)
    assert dict(counts) == prod["filters"]
    for a, b in zip(records, prod["rows"]):
        assert a == b
    result = dict(n=n, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  covered_frames=len(covered), filters=dict(counts),
                  survivors=[r for r in records if "aut" in r],
                  all_rows_agree_with_producer=True)
    print(json.dumps(result), flush=True)
    return result


def main():
    catalogs = [audit_catalog(n) for n in (3, 5, 7, 9)]
    s4_bits = [2,4,9,3]
    s4 = [[(x >> j) & 1 for j in range(4)] for x in s4_bits]
    c3 = [[int((j-i)%3 == 1) for j in range(3)] for i in range(3)]
    c5 = [[int((j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]
    paley7 = [[int((j-i)%7 in (1,2,4)) for j in range(7)] for i in range(7)]
    controls = [exhaustive_control(name, A) for name,A in (
        ("C3", c3), ("S4",s4), ("regular_C5",c5), ("Paley7",paley7))]
    assert [(r["mass"],r["B"]) for r in controls] == [(1,1),(5,3),(3,3),(9,3)]
    # C5 followed by its converse: direct DP plus degree-preserving permutations.
    stack = [[0]*10 for _ in range(10)]
    for i in range(10):
        for j in range(10):
            if i < 5 and j < 5: stack[i][j] = c5[i][j]
            elif i >= 5 and j >= 5: stack[i][j] = c5[j-5][i-5]
            else: stack[i][j] = int(i < j)
    auts = restricted_automorphisms(stack)
    rho = list(range(5,10)) + list(range(5))
    centralizer = sum(all(p[rho[i]] == rho[p[i]] for i in range(10)) for p in auts)
    H = count_paths(stack)
    B = symmetric_words(stack,rho) // centralizer
    assert (H, len(auts), B) == (225,25,3)
    controls.append(dict(name="C5_to_C5op", H=H, aut=len(auts), B=B,
                         mass=9, deficit=6, paired_black_degree=3))
    for row in controls: print(json.dumps(row),flush=True)
    (HERE / "verification.json").write_text(
        json.dumps(dict(catalogs=catalogs, controls=controls),indent=2)+"\n")
    print("ALL INDEPENDENT CHECKS PASSED",flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    args = parser.parse_args()
    SOURCE = args.source
    main()
