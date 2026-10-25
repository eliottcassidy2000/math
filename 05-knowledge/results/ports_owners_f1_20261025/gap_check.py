"""Exact next-band weighted-core records; explicit quotient and clone checks.

No inherited producer imports. The all-order theorem is in gap.md; these
bounded checks do not prove its quantifiers. Output only, no receipt writes.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import comb
import json


CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for b in xs[1:]:
        for rest in matchings(tuple(x for x in xs[1:] if x != b)):
            yield ((a, b),) + rest


def df(k):
    answer = 1
    for a in range(1, 2*k, 2):
        answer *= a
    return answer


def choose(n, r):
    return comb(n, r) if n >= r >= 0 else 0


def indecomp(edges, n):
    return all(any(a < cut <= b for a, b in edges) for cut in range(1, n))


def overlap_graphs(k):
    xs = tuple(range(2*k-1))
    for c in xs:
        for a, b in combinations((v for v in xs if v != c), 2):
            rem = tuple(v for v in xs if v not in (a, b, c))
            for m in matchings(rem):
                yield tuple(sorted(m + (tuple(sorted((a, c))), tuple(sorted((b, c))))))


def compositions(total, parts, cap):
    if parts == 1:
        if 0 <= total <= cap:
            yield (total,)
        return
    for a in range(max(0, total-cap*(parts-1)), min(cap, total)+1):
        for rest in compositions(total-a, parts-1, cap):
            yield (a,) + rest


def make_record(k, edges, deficits, strong, overlap=False):
    t = 2*k-1 if overlap else 2*k
    gaps = tuple(k-d for d in deficits)
    if strong:
        gaps = (0,) + gaps + (0,)
    need(len(gaps) == t+1, 'full gap arity')
    slots = []
    weights = []
    for i, g in enumerate(gaps):
        need(0 <= g <= k, 'gap weight in cap')
        if g:
            weights.append(g)
        if i < t:
            slots.append(len(weights))
            weights.append(1)
    n = len(weights)
    rows = [(1 << n) - (1 << (i+1)) for i in range(n)]
    for a, b in edges:
        a, b = slots[a], slots[b]
        rows[a] ^= 1 << b
        rows[b] ^= 1 << a
    return tuple(rows), tuple(weights), gaps


def clones(rows):
    mask = (1 << len(rows))-1
    return tuple((a, b) for a, b in combinations(range(len(rows)), 2)
                 if not ((rows[a] ^ rows[b]) & (mask ^ (1 << a) ^ (1 << b))))


def strong_matrix(rows):
    n = len(rows)
    reach = [[bool(rows[i] >> j & 1) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    return all(all(row) for row in reach)


def relabel(rows, weights, order):
    return (tuple(sum(((rows[a] >> b) & 1) << j for j, b in enumerate(order))
                  for a in order), tuple(weights[a] for a in order))


def serial(rows, weights, order):
    rr, ww = relabel(rows, weights, order)
    size = (len(rows)+7)//8
    return bytes(ww) + b''.join(r.to_bytes(size, 'little') for r in rr)


def canonical(rows, weights, independent=False):
    cells = defaultdict(list)
    if independent:
        for i, w in enumerate(weights):
            cells[(w, rows[i].bit_count())].append(i)
        prefix = ()
    else:
        large = tuple(i for i, w in enumerate(weights) if w > 1)
        prefix = tuple(sorted(large, key=lambda i: -sum(rows[i] >> j & 1 for j in large)))
        need(all(rows[a] >> b & 1 for j, a in enumerate(prefix) for b in prefix[j+1:]),
             'large induced tournament transitive')
        for i, w in enumerate(weights):
            if w == 1:
                cells[tuple(rows[i] >> j & 1 for j in prefix)].append(i)
    best = None
    for ps in product(*(permutations(cells[c]) for c in sorted(cells))):
        order = prefix + tuple(x for p in ps for x in p)
        key = serial(rows, weights, order)
        if best is None or key < best:
            best = key
    return best, max((len(v) for v in cells.values()), default=0)


def optimum(rows, weights):
    n = len(rows)
    dp = [0] + [10**9]*((1 << n)-1)
    for mask in range(1, 1 << n):
        for a in range(n):
            if mask >> a & 1:
                rest = mask ^ (1 << a)
                cost = weights[a]*sum(weights[b] for b in range(n)
                                      if rest >> b & 1 and rows[a] >> b & 1)
                dp[mask] = min(dp[mask], dp[rest]+cost)
    return dp[-1]


def sequences(limit):
    I = [0]*(limit+1)
    J = [0]*(limit+1)
    for k in range(1, limit+1):
        I[k] = df(k)-sum(I[j]*df(k-j) for j in range(1, k))
        J[k] = (k-1)*I[k]+sum(I[j]*I[k-j] for j in range(1, k))
    return I, J


def formula(k, delta, strong):
    I, J = sequences(k)
    if strong:
        return (I[k]*(choose(2*k+delta-2, delta)
                      -(2*k-1)*choose(k+delta-3, delta-k-1))
                -(2*k-3)*I[k-1]*choose(k+delta-3, delta-k)
                +J[k]*choose(k+delta-4, delta-k-1))
    return df(k)*(choose(2*k+delta, delta)
                  -(2*k+1)*choose(k+delta-1, delta-k-1)
                  -choose(k+delta-1, delta-k)
                  +(k-1)*choose(k+delta-2, delta-k-1))


def verify_layer(k, delta, strong):
    seen = set()
    independent = set()
    counts = Counter()
    max_cell = 0
    dp_checks = covariance = 0
    sample_hostile = None
    for overlap in (False, True):
        remaining = delta-(k+1 if overlap else 0)
        if remaining < 0:
            continue
        t = 2*k-1 if overlap else 2*k
        family = overlap_graphs(k) if overlap else matchings(tuple(range(t)))
        g = t-1 if strong else t+1
        for edges in family:
            if strong and not indecomp(edges, t):
                continue
            for deficits in compositions(remaining, g, k-1 if overlap else k):
                rows, weights, gaps = make_record(k, edges, deficits, strong, overlap)
                kind = 'overlap' if overlap else ('missing' if 0 in (gaps[1:-1] if strong else gaps) else 'full')
                counts['raw_'+kind] += 1
                cp = clones(rows)
                bad = any(gaps[i] == 0 and (i-1, i) in edges for i in range(1, t)) if not overlap else False
                need(bool(cp) == bad, 'complete clone iff')
                if cp:
                    counts['rejected'] += 1
                    if sample_hostile is None:
                        sample_hostile = {'edges': edges, 'gaps': gaps, 'clone_pairs': cp}
                    continue
                need(sum(weights) == 2*k*k+(k if strong else 3*k)-delta, 'mass')
                # Every source is retained in the all case; strongness is an independently
                # checked constructor predicate in the strong case, not a hidden filter.
                if strong or len(seen) < 12:
                    need(strong_matrix(rows) == (strong or (gaps[0] == gaps[-1] == 0 and indecomp(edges, t))), 'literal strongness')
                key, cell = canonical(rows, weights)
                max_cell = max(max_cell, cell)
                need(key not in seen, 'distinct records have distinct weighted-isomorphism owners')
                seen.add(key)
                counts[kind] += 1
                if k == 2:
                    other, _ = canonical(rows, weights, independent=True)
                    need(other not in independent, 'independent quotient has no collision')
                    independent.add(other)
                    need(optimum(rows, weights) == k, 'independent exact weighted feedback')
                    dp_checks += 1
                if covariance < 8:
                    rr, ww = relabel(rows, weights, tuple(reversed(range(len(rows)))))
                    need(canonical(rr, ww)[0] == key, 'reversal label covariance')
                    order = tuple(range(1, len(rows)))+(0,)
                    rr, ww = relabel(rows, weights, order)
                    need(canonical(rr, ww)[0] == key, 'cyclic label covariance')
                    covariance += 2
    need(len(seen) == formula(k, delta, strong), 'exact band formula')
    if k == 2:
        need(len(independent) == len(seen), 'independent class total')
    return {'k': k, 'delta': delta, 'strong': strong, 'classes': len(seen),
            'sectors': dict(sorted(counts.items())), 'max_unit_cell': max_cell,
            'feedback_checks': dp_checks, 'relabel_checks': covariance,
            'clone_hostile': sample_hostile}


def main():
    # All existing k2 layers, its exceptional delta3, and both genuinely new
    # band levels at k3/4. The formula's k2 delta3 is separately finite-exact.
    pairs = [(2, d) for d in range(4)] + [(3, 4), (4, 5), (4, 6)]
    layers = [verify_layer(k, d, s) for k, d in pairs for s in (False, True)]
    I, J = sequences(8)
    graph_counts = []
    for k in range(2, 6):
        graphs = tuple(overlap_graphs(k))
        need(len(set(graphs)) == len(graphs) == (k-1)*df(k), 'all overlap defect graphs')
        count = sum(indecomp(e, 2*k-1) for e in graphs)
        need(count == J[k], 'independent overlap indecomposability recurrence')
        graph_counts.append({'k': k, 'all': len(graphs), 'indecomposable': count})
    # A tempting single-P3 multiplier misses the two-sided placement term.
    need(J[2] == 3 and (2-1)*I[2] == 2, 'J is not merely (k-1)I')
    # Independent older complete compiler coefficients, cited as receipt controls.
    small = {(x['delta'], x['strong']): x['classes'] for x in layers if x['k'] == 2}
    need([small[d, False] for d in range(4)] == [3, 15, 42, 81], 'inherited full k2 mass coefficients')
    need([small[d, True] for d in range(4)] == [2, 6, 11, 15], 'inherited strong k2 mass coefficients')
    print(json.dumps({'status': 'FINITE-EXACT controls; proofs in gap.md',
                      'layers': layers, 'overlap_graph_counts': graph_counts,
                      'I_1_to_8': I[1:], 'J_1_to_8': J[1:], 'checks': CHECKS}, indent=2))


if __name__ == '__main__':
    main()
