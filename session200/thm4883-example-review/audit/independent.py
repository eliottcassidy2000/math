"""Independent enumeration: no subset DP and no imports from verify.py."""
from itertools import permutations
from math import gcd
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_path(rows, p):
    return all(rows[a] & (1 << b) for a, b in zip(p, p[1:]))


def profile_from_permutations(rows, root):
    n = len(rows)
    paths = [p for p in permutations(range(n)) if is_path(rows, p)]
    covers = 0
    # Each ordered cover is represented uniquely by its vertex order and split.
    # Its first path begins at root; the complementary second path may be empty.
    for p in permutations(range(n)):
        if p[0] != root:
            continue
        for split in range(1, n + 1):
            covers += is_path(rows, p[:split]) and is_path(rows, p[split:])
    deletion = sum(is_path(rows, p) for p in permutations(i for i in range(n) if i != root))
    rooted = [p for p in paths if p[0] == root]
    cycles = sum(bool(rows[p[-1]] & (1 << root)) for p in rooted)
    return [len(paths), sum(p[-1] == root for p in paths), covers, deletion], len(rooted), cycles


def all_path_count(rows):
    # Direct recursive path enumeration, with no memoization or subset tables.
    n = len(rows)
    full = (1 << n) - 1
    def walk(start, last, used):
        if used == full:
            return 1, int(bool(rows[last] & (1 << start)))
        paths = cycles = 0
        for vertex in range(n):
            if not (used >> vertex) & 1 and (rows[last] >> vertex) & 1:
                p, c = walk(start, vertex, used | (1 << vertex))
                paths += p
                cycles += c
        return paths, cycles
    totals = [walk(start, start, 1 << start) for start in range(n)]
    # A directed Hamilton cycle is counted once from each starting vertex.
    cycle_sum = sum(c for _, c in totals)
    require(cycle_sum % n == 0, 'cycle normalization')
    return sum(p for p, _ in totals), cycle_sum // n


def check_tournament(rows):
    n = len(rows)
    for i in range(n):
        require(rows[i] >= 0 and rows[i] < (1 << n), 'out-of-range adjacency')
        require(not (rows[i] >> i) & 1, 'loop')
        for j in range(i):
            require(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1, 'non-tournament')


receipt = json.loads((Path(__file__).resolve().parents[1] / 'VALIDATION_NORMAL.json').read_text())
result = {'status': 'PASS', 'method': 'Permutation source profiles and unmemoized direct path enumeration'}
for name, expected in [('ZOJ', [9, 3, 9, 3]), ('ZOJJ', [15, 6, 18, 9])]:
    item = receipt[name]
    check_tournament(item['rows'])
    profile, starts, cycles = profile_from_permutations(item['rows'], item['root'])
    require(profile == expected and starts == cycles == 1, ('source', name))
    result[name] = {'profile': profile, 'root_start_paths': starts, 'hamilton_cycles': cycles}
for name, expected in [('N', 417), ('M', 867)]:
    rows = receipt[name]['rows']
    check_tournament(rows)
    paths, cycles = all_path_count(rows)
    require(paths == expected and cycles == 2, ('graft', name))
    result[name] = {'hamilton_paths': paths, 'hamilton_cycles': cycles, 'order': len(rows)}
require(gcd(69, 3) == 3 and 867 - 2 * 417 == 33, 'offset arithmetic')
result['universal_identity'] = '69(2k+1)-66U+24 = 6(23k-11U+15)+3, for all integers k,U'
result['scope'] = 'Literal counts and universal mod-6 obstruction only; no full THM-4883 audit'
print(json.dumps(result, sort_keys=True, indent=2))
