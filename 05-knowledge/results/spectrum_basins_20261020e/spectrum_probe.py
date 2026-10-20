"""Raw Hamiltonian-count controls: complete small spectrum and apex recurrence.

Universe one: every presentation with standard directed Hamiltonian path,
orders one through six. This covers every tournament isomorphism class,
without claiming uniform multiplicities. Canonicalization uses ALL actual
Hamiltonian words, so its canonical form stays inside this same universe.
Universe two: all binary one-apex cuts over transitive bases through order8.
Every raw count comes from an independent subset DP, never an orbit count.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(t, u, v):
    return bool(t[u] >> v & 1)


def count_paths(t):
    n = len(t)
    if not n:
        return 1
    full = (1 << n)-1
    table = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask in range(1, full+1):
        for u, value in enumerate(table[mask]):
            if not value:
                continue
            unused_out = t[u] & (full ^ mask)
            while unused_out:
                bit = unused_out & -unused_out
                v = bit.bit_length()-1
                table[mask | bit][v] += value
                unused_out ^= bit
    return sum(table[full])


def paths(t):
    n = len(t)
    def visit(prefix, unused):
        if not unused:
            yield prefix
        for v in sorted(unused):
            if not prefix or arc(t, prefix[-1], v):
                yield from visit(prefix+(v,), unused-{v})
    yield from visit((), set(range(n)))


def canonical(t):
    return min(tuple(sum(1 << j for j, v in enumerate(p) if arc(t, u, v))
                     for u in p) for p in paths(t))


def presentations(n):
    free = [(u, v) for u, v in combinations(range(n), 2) if v != u+1]
    for mask in range(1 << len(free)):
        t = [0]*n
        for u in range(n-1):
            t[u] |= 1 << (u+1)
        for j, (u, v) in enumerate(free):
            if mask >> j & 1:
                t[u] |= 1 << v
            else:
                t[v] |= 1 << u
        yield tuple(t)


def components(t):
    n = len(t)
    reach = [t[u] | (1 << u) for u in range(n)]
    for v in range(n):
        for u in range(n):
            if reach[u] >> v & 1:
                reach[u] |= reach[v]
    unused = set(range(n))
    result = []
    while unused:
        u = min(unused)
        c = {v for v in unused if reach[u] >> v & 1 and reach[v] >> u & 1}
        unused -= c
        result.append(sorted(c))
    return result


def induced(t, vertices):
    return tuple(sum(1 << j for j, v in enumerate(vertices) if arc(t, u, v))
                 for u in vertices)


def minimum_strong_paths(n):
    need(n >= 3, 'strong order is below three')
    if n % 3 == 0:
        return 3*5**((n-3)//3)
    if n % 3 == 1:
        return 5**((n-1)//3)
    return 9*5**((n-5)//3)


def one_apex(bits):
    m = len(bits)
    t = [sum(1 << j for j in range(i+1, m)) for i in range(m)]+[0]
    for i, bit in enumerate(bits):
        if bit:
            t[m] |= 1 << i
        else:
            t[i] |= 1 << m
    return tuple(t)


def apex_recurrence(bits):
    h, a, last = 1, 0, 0
    for bit in bits:
        if not bit:
            h += a+last
        a, last = 2*a+bit, bit
    return h


def apex_pair_formula(bits):
    return 1+sum((2 if j == i+1 else 1 << (j-i-1))
                 for i, j in combinations(range(len(bits)), 2)
                 if bits[i] and not bits[j])


def small_spectrum():
    rows, strong_polynomials, witnesses = [], {}, {}
    expected_classes = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56}
    for n in range(1, 7):
        raw = Counter()
        strong_values = set()
        reps = {}
        count = 0
        for t in presentations(n):
            h = count_paths(t)
            need(h % 2 == 1, 'Redei parity failed')
            cs = components(t)
            prod = 1
            core_size = 0
            for c in cs:
                hc = count_paths(induced(t, c))
                prod *= hc
                if len(c) > 1:
                    core_size += len(c)
                    need(hc >= minimum_strong_paths(len(c)), 'Busch floor control failed')
            need(prod == h, 'SCC multiplicativity failed')
            need(3**core_size <= h**3, 'logarithmic core bound failed')
            if 3**core_size == h**3:
                need(all(len(c) in (1,3) for c in cs), 'core equality has non-C3 component')
            rep = canonical(t)
            need(rep not in reps or reps[rep] == h, 'canonical count mismatch')
            reps[rep] = h
            raw[h] += 1
            if len(cs) == 1 and n >= 3:
                strong_values.add(h)
                witnesses.setdefault(str(h), dict(order=n, adjacency=t))
            count += 1
        need(len(reps) == expected_classes[n], 'ordinary class coverage mismatch')
        for rep, h in reps.items():
            if n >= 3 and len(components(rep)) == 1:
                key = str(h)
                strong_polynomials.setdefault(key, {})[str(n)] = (
                    strong_polynomials.setdefault(key, {}).get(str(n), 0)+1)
        rows.append(dict(order=n, standard_path_presentations=count,
                         ordinary_isomorphism_classes=len(reps),
                         ordinary_raw_H_counts=dict(sorted(Counter(reps.values()).items())),
                         raw_H_values=sorted(raw), strong_raw_H_values=sorted(strong_values)))
    need(all(7 not in r['raw_H_values'] and 21 not in r['raw_H_values'] for r in rows),
         'excluded value found')
    need('13' in witnesses, 'missing H13 hostile')
    hostile = tuple(witnesses['13']['adjacency'])
    one_deletions = [count_paths(induced(hostile, [u for u in range(5) if u != v]))
                     for v in range(5)]
    two_deletions = [pair for pair in combinations(range(5), 2)
                     if count_paths(induced(hostile, [u for u in range(5) if u not in pair])) == 1]
    need(min(one_deletions) > 1 and bool(two_deletions), 'H13 feedback-vertex hostile differs')
    witnesses['13']['single_deletion_H'] = one_deletions
    witnesses['13']['transitive_two_vertex_deletion'] = two_deletions[0]
    polynomials = {
        k:v for k,v in sorted(strong_polynomials.items(), key=lambda kv:int(kv[0]))
        if int(k) <= 23}
    @lru_cache(None)
    def coefficient(target, order):
        if order < 0:
            return 0
        if target == 1:
            return 1
        total = 0
        for ds, polynomial in polynomials.items():
            d = int(ds)
            if target % d:
                continue
            for degree, multiplicity in polynomial.items():
                for prefix in range(order-int(degree)+1):
                    total += multiplicity*coefficient(target//d, order-int(degree)-prefix)
        return total
    for row in rows:
        for target in range(1, 24, 2):
            need(coefficient(target, row['order']) == row['ordinary_raw_H_counts'].get(target, 0),
                 'fixed-raw-count rational recurrence disagrees with complete class census')
    return dict(census=rows, complete_strong_count_polynomials_through_H23=polynomials,
                fixed_raw_count_series_coefficients_through_order12={
                    str(target):[coefficient(target,n) for n in range(13)]
                    for target in range(1,24,2)}, first_strong_witnesses=witnesses)


def apex_controls():
    rows = []
    for m in range(9):
        vals = set()
        for mask in range(1 << m):
            bits = tuple((mask >> i) & 1 for i in range(m))
            t = one_apex(bits)
            h = count_paths(t)
            need(h == apex_recurrence(bits) == apex_pair_formula(bits), 'apex formulas disagree')
            vals.add(h)
            active_interval = [i for i in range(m) if bits[i]]
            zero_positions = [i for i in range(m) if not bits[i]]
            cs = components(t)
            if active_interval and zero_positions and min(active_interval) < max(zero_positions):
                left, right = min(active_interval), max(zero_positions)
                expected = set(range(left, right+1)) | {m}
                actual = [set(c) for c in cs if len(c)>1]
                need(actual == [expected], 'apex strong-core interval fails')
            else:
                need(h == 1, 'inversion-free apex cut is not transitive')
        rows.append(dict(transitive_base_order=m, cuts=1 << m,
                         raw_H_values=sorted(vals)))
    for a in range(1, 9):
        for b in range(1, 9):
            need(apex_recurrence((1,)*a+(0,)*b) == (2**a-1)*(2**b-1)+2,
                 'two-run family formula fails')
    need(all(13 not in row['raw_H_values'] for row in rows), 'one-apex H13 hostile failed')
    return dict(census=rows, total_literal_cuts=sum(r['cuts'] for r in rows),
                two_run_arithmetic_checks=64,
                minimal_lost_state=dict(old_words=[[0],[1]], old_H=[1,1],
                                       append_zero_H=[apex_recurrence((0,0)),apex_recurrence((1,0))]))


def main():
    result = dict(status='FINITE-EXACT small controls; the unbounded reductions are proved in spectrum.md',
                  count_type='raw directed Hamiltonian paths, not automorphism orbits or blue paths',
                  small_spectrum=small_spectrum(), one_apex=apex_controls(),
                  source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest())
    output = json.dumps(result, indent=2)+'\n'
    Path(__file__).with_suffix('.json').write_text(output, encoding='utf-8')
    print(output, end='')


if __name__ == '__main__':
    main()
