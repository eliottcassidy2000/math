#!/usr/bin/env python3
"""Independent direct-adjacency audit of normalized path relocations.
No imported census, Burnside routine or production classifier is used.
The n<=6 oracle enumerates every vertex permutation and checks the full path.
"""
import argparse
import itertools
import json
import time
from collections import Counter
from pathlib import Path

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

class DSU:

    def __init__(self, n):
        self.p = list(range(n))
        self.sz = [1] * n

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        a, b = (self.find(a), self.find(b))
        if a == b:
            return
        if self.sz[a] < self.sz[b]:
            a, b = (b, a)
        self.p[b] = a
        self.sz[a] += self.sz[b]

def audit(n, oracle=True, compare_dir=None):
    t = time.perf_counter()
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    m = len(pairs)
    N = 1 << m
    ix = {p: k for k, p in enumerate(pairs)}

    def decode(x):
        a = [[False] * n for _ in range(n)]
        for i in range(n - 1):
            a[i][i + 1] = True
        for k, (i, j) in enumerate(pairs):
            a[i][j] = bool(x >> k & 1)
        for i in range(n):
            for j in range(i + 1, n):
                a[j][i] = not a[i][j]
        return a

    def encode(a, p):
        return sum((int(a[p[i]][p[j]]) << k for k, (i, j) in enumerate(pairs)))

    def ispath(a, p):
        return all((a[p[i]][p[i + 1]] for i in range(n - 1)))

    def reflect(x):
        return sum(((x >> ix[n - 1 - j, n - 1 - i] & 1) << k for k, (i, j) in enumerate(pairs)))
    ordinary = DSU(N)
    rightinc = 0
    collapsed = 0
    allmovechecked = 0
    perms = list(itertools.permutations(range(n))) if oracle else None
    oracle_keys = []
    reflection_keys = []
    hp_total = 0
    moves = []
    fullindex = list(range(n))
    for x in range(N):
        a = decode(x)
        actual_right = set()
        for i in range(n):
            for j in range(i + 2, n):
                p = fullindex[:i] + fullindex[i + 1:j + 1] + [i] + fullindex[j + 1:]
                direct = ispath(a, p)
                predicate = a[j][i] and (i == 0 or a[i - 1][i + 1]) and (j == n - 1 or a[i][j + 1])
                require(predicate == direct, (n, x, i, j, 'right legality'))
                allmovechecked += 1
                if direct:
                    y = encode(a, p)
                    ordinary.union(x, y)
                    rightinc += 1
                    collapsed += x == y
                    actual_right.add(y)
                    moves.append((x, y))
        if oracle:
            paths = [p for p in perms if ispath(a, p)]
            keys = [encode(a, p) for p in paths]
            require(x in keys, 'verification failed')
            oracle_keys.append(min(keys))
            hp_total += len(paths)
            reflection_keys.append(min(min(keys), min((reflect(k) for k in keys))))
            for i in range(n):
                rest = fullindex[:i] + fullindex[i + 1:]
                for gap in range(n):
                    p = rest[:gap] + [i] + rest[gap:]
                    if p == fullindex:
                        continue
                    direct = ispath(a, p)
                    if abs(gap - i) == 1:
                        require(not direct, 'verification failed')
                    if direct:
                        y = encode(a, p)
                        if gap > i:
                            require(y in actual_right, 'verification failed')
                        else:
                            ay = decode(y)
                            inv = list(range(n))
                            u = inv.pop(gap)
                            inv.insert(i, u)
                            require(ispath(ay, inv) and encode(ay, inv) == x, 'verification failed')
    roots = [ordinary.find(x) for x in range(N)]
    if oracle:
        d = {}
        for x in range(N):
            r = roots[x]
            d.setdefault(r, oracle_keys[x])
            require(d[r] == oracle_keys[x], 'verification failed')
        require(len(d) == len(set(oracle_keys)), 'verification failed')
        for key, group in itertools.groupby(sorted(range(N), key=lambda x: oracle_keys[x]), key=lambda x: oracle_keys[x]):
            require(len({roots[x] for x in group}) == 1, 'verification failed')
    ordinary_count = len(set(roots))
    self_roots = {roots[x] for x in range(N) if roots[x] == roots[reflect(x)]}
    merged = DSU(N)
    for x, y in moves:
        merged.union(x, y)
    for x in range(N):
        merged.union(x, reflect(x))
    merged_roots = [merged.find(x) for x in range(N)]
    if oracle:
        d = {}
        for x in range(N):
            r = merged_roots[x]
            d.setdefault(r, reflection_keys[x])
            require(d[r] == reflection_keys[x], 'verification failed')
        require(len(d) == len(set(reflection_keys)), 'verification failed')
    merged_count = len(set(merged_roots))
    require(2 * merged_count == ordinary_count + len(self_roots), 'verification failed')
    fixed = [x for x in range(N) if reflect(x) == x]
    sigma_cycles = 0
    seen = set()
    for pair in pairs:
        if pair not in seen:
            sigma_cycles += 1
            seen.add(pair)
            seen.add((n - 1 - pair[1], n - 1 - pair[0]))
    require(len(fixed) == 1 << sigma_cycles, 'verification failed')
    hist = Counter()
    line_count = blue = loops = 0
    if n >= 3:
        for x in range(N):
            y = N - 1 ^ x
            require(y != x and reflect(y) == N - 1 ^ reflect(x), 'verification failed')
            require((reflect(x) == x) == (reflect(y) == y), 'verification failed')
            if x < y:
                line_count += 1
                colour = reflect(x) == x
                blue += colour
                u, v = sorted((merged_roots[x], merged_roots[y]))
                hist[u, v, colour] += 1
                loops += u == v
        require(line_count == N // 2 and blue == len(fixed) // 2, 'verification failed')
    if compare_dir is not None:
        import csv
        ordinary_min = {}
        merged_min = {}
        for x in range(N):
            ordinary_min.setdefault(roots[x], x)
            merged_min.setdefault(merged_roots[x], x)
        with (Path(compare_dir) / f'nodes_{n}.tsv').open() as f:
            nodes = list(csv.DictReader(f, delimiter='\t'))
        require(len(nodes) == ordinary_count, 'verification failed')
        by_representative = {int(r['ordinary_representative']): r for r in nodes}
        require(set(by_representative) == set(ordinary_min.values()), 'verification failed')
        counts = Counter(roots)
        for x, r in by_representative.items():
            require(int(r['tilings_tc']) == counts[roots[x]], 'verification failed')
            require(merged_roots[int(r['merged_representative'])] == merged_roots[x], 'verification failed')
            require(bool(int(r['self_converse'])) == (roots[x] == roots[reflect(x)]), 'verification failed')
        if n >= 3:
            actual = Counter()
            with (Path(compare_dir) / f'lines_{n}.tsv').open() as f:
                for r in csv.DictReader(f, delimiter='\t'):
                    u = merged_min[merged_roots[int(r['u'])]]
                    v = merged_min[merged_roots[int(r['v'])]]
                    u, v = sorted((u, v))
                    actual[u, v, bool(int(r['blue']))] += int(r['multiplicity'])
            expected = Counter()
            for (r, s, b), c in hist.items():
                u, v = sorted((merged_min[r], merged_min[s]))
                expected[u, v, b] += c
            require(actual == expected, (n, 'complete line tables mismatch'))
    out = {'n': n, 'free_bits': m, 'tilings': N, 'ordinary_classes': ordinary_count, 'merged_classes': merged_count, 'self_converse_classes': len(self_roots), 'right_legal_incidences_including_collapses': rightinc, 'collapsed_right_incidences': collapsed, 'right_legality_predicates_checked': allmovechecked, 'all_permutation_partition_oracle': oracle, 'production_node_and_full_line_tables_match': compare_dir is not None, 'permutation_hamilton_words_examined': hp_total if oracle else None, 'reflection_fixed_tilings': len(fixed), 'complement_lines': line_count if n >= 3 else None, 'blue_lines': blue if n >= 3 else None, 'merged_loops_including_multiplicity': loops if n >= 3 else None, 'ordinary_fibre_sizes': dict(sorted(Counter(Counter(roots).values()).items())), 'seconds': round(time.perf_counter() - t, 4)}
    print(json.dumps(out, sort_keys=True), flush=True)
    return out
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-n', type=int, default=7)
    ap.add_argument('--oracle-max-n', type=int, default=6)
    ap.add_argument('--output', default='independent_results.json')
    ap.add_argument('--compare-dir', default=None)
    args = ap.parse_args()
    rows = [audit(n, n <= args.oracle_max_n, args.compare_dir) for n in range(1, args.max_n + 1)]
    Path(args.output).write_text(json.dumps({'status': 'PASS', 'rows': rows}, indent=2) + '\n')
