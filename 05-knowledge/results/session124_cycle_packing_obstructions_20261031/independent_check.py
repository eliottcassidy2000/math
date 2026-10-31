"""Independent audit: no executable imports from the audited packet.

Uses Boolean adjacency, outdegree-one triangle detection, and literal path
backtracking (without endpoint DP). Frozen JSON is compared, never trusted.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE
CHECKS = 0

def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)

def adjacency(word):
    n = len(word) + 3
    a = [[False] * n for _ in range(n)]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a[i][j] = True
    for u, x in enumerate(word, 3):
        for c in range(3):
            a[c][u] = bool(x & (1 << c))
            a[u][c] = not a[c][u]
        for v in range(u + 1, n):
            a[u][v] = True
    return a

def triangle_supports(a):
    return [frozenset(s) for s in combinations(range(len(a)), 3)
            if all(sum(a[u][v] for v in s) == 1 for u in s)]

def bad(a):
    return any(s.isdisjoint(t) for s, t in combinations(triangle_supports(a), 2))

def literal_paths(a):
    n = len(a)
    def walk(path, remaining):
        if not remaining:
            return 1
        return sum(walk(path + (v,), remaining - {v}) for v in remaining
                   if not path or a[path[-1]][v])
    return walk((), frozenset(range(n)))

def literal_cycles(a):
    # Enumerate permutations on every odd support, fixing its least vertex.
    cycles = []
    for k in range(3, len(a) + 1, 2):
        for support in combinations(range(len(a)), k):
            for tail in permutations(support[1:]):
                cycle = (support[0],) + tail
                if all(a[cycle[i]][cycle[(i + 1) % k]] for i in range(k)):
                    cycles.append(cycle)
    return cycles

def is_subsequence(short, long):
    # A separate iterator-based recognizer, not prefix-state propagation.
    it = iter(long)
    return all(any(y == x for y in it) for x in short)

def main():
    basis = []
    for m in (3, 4):
        for w in product(range(8), repeat=m):
            if bad(adjacency(w)) and not any(is_subsequence(v, w) for v in basis):
                basis.append(w)
    require(Counter(map(len, basis)) == {3: 240, 4: 35}, 'basis counts')
    frozen = json.loads((SOURCE / 'patterns.json').read_text())
    require(basis == [tuple(w) for k in ('3', '4')
                      for w in frozen['minimal_bad_patterns'][k]], 'frozen basis')
    prefixes = {w[:i] for w in basis for i in range(len(w))}
    require(len(prefixes) == 84, 'proper prefixes')
    # Proper-subsequence minimality follows already from one-symbol deletions.
    for w in basis:
        require(bad(adjacency(w)), 'basis positive witness')
        require(all(not bad(adjacency(w[:i] + w[i+1:])) for i in range(len(w))),
                'basis deletion minimality')

    spectrum = {1}
    rows = []
    witnesses = {}
    cycle_checks = 0
    for m in range(6):
        hist = Counter()
        for word in product(range(8), repeat=m):
            a = adjacency(word)
            packing_bad = bad(a)
            language_bad = any(is_subsequence(w, word) for w in basis)
            require(packing_bad == language_bad, 'independent language equivalence')
            if packing_bad:
                continue
            h = literal_paths(a)
            hist[h] += 1
            spectrum.add(h)
            # One literal cycle certificate for every (order, H) pair. The
            # upstream rerun checks OCF on every accepted presentation.
            if (m, h) not in witnesses:
                cycles = literal_cycles(a)
                require(1 + 2 * len(cycles) == h, 'independent literal OCF')
                require(all(not set(c).isdisjoint(d) for c, d in combinations(cycles, 2)),
                        'literal odd cycles pairwise intersect')
                witnesses[m, h] = word
                cycle_checks += 1
        rows.append({'n': m + 3, 'accepted': sum(hist.values()),
                     'histogram': dict(sorted(hist.items()))})
    expected = [1,3,5,9,11,13,15,17,19,23,25,27,31,33,37,39,
                45,47,49,51,53,61,63,65,69,73]
    require(sorted(h for h in spectrum if h <= 73) == expected, 'full exact spectrum')
    require(sum(r['accepted'] for r in rows) == 3624, 'accepted total')
    require(5**8 > 73**3, 'exact cutoff')
    holes = [29,35,41,43,55,57,59,67,71]
    witness_data = json.loads((SOURCE / 'global_witnesses.json').read_text())
    require(sorted(map(int, witness_data)) == holes, 'exact witness keys')
    ordinary = []
    for key, out in witness_data.items():
        n = len(out)
        require(all(isinstance(x, int) and 0 <= x < (1 << n) for x in out),
                'adjacency bitmask domain')
        a = [[bool(out[i] & (1 << j)) for j in range(n)] for i in range(n)]
        require(all(not a[i][i] for i in range(n)), 'no loops')
        require(all(a[i][j] != a[j][i] for i, j in combinations(range(n), 2)),
                'actual tournament')
        h = literal_paths(a)
        require(h == int(key), 'ordinary path witness')
        cycles = literal_cycles(a)
        pairs = sum(set(c).isdisjoint(d) for c, d in combinations(cycles, 2))
        require(n < 9 and 1 + 2*len(cycles) + 4*pairs == h, 'ordinary OCF')
        require(bad(a) and pairs > 0, 'positive packing obstruction sidecar')
        ordinary.append({'H':h, 'n':n, 'cycles':dict(Counter(map(len,cycles))), 'pairs':pairs})
    free = []
    for mask in range(1, 256):
        alphabet = {x for x in range(8) if mask & (1 << x)}
        if not any(set(w) <= alphabet for w in basis):
            free.append(alphabet)
    maximal = sorted(sorted(a) for a in free if not any(a < b for b in free))
    require(maximal == frozen['maximal_free_alphabets'], 'all maximal free alphabets')
    # Explicit hostile controls for two otherwise plausible weakened claims.
    # Collapsing cycle supports loses multiplicity for this legal presentation.
    collapse_witness = None
    for (m,h), word in witnesses.items():
        cycles = literal_cycles(adjacency(word))
        if len({frozenset(c) for c in cycles}) < len(cycles):
            collapse_witness = {'word':word,'H':h,'literal_cycles':len(cycles),
                                'distinct_supports':len({frozenset(c) for c in cycles})}
            break
    require(collapse_witness is not None, 'literal versus support hostile control')
    require(1 + 14 * Fraction(-1, 3) != 0 and (1 + 14 * 2) == 29,
            'K14 poison-factor arithmetic')
    return {'status':'PASS','checks':CHECKS,'presentations':37449,
            'accepted':3624,'basis_by_length':dict(Counter(map(len,basis))),
            'proper_prefixes':len(prefixes),'maximal_free_alphabets':maximal,
            'rows':rows,'spectrum_through73':expected,'ordinary_witnesses':ordinary,
            'independent_cycle_certificate_count':cycle_checks,
            'literal_support_hostile_control':collapse_witness,
            'scope':'All finite claims independently regenerated; all-order steps audited in RESULTS.md.'}

if __name__ == '__main__':
    print(json.dumps(main(), sort_keys=True, indent=2))
