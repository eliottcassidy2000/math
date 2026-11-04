#!/usr/bin/env python3
"""Named old-port obstructions and a fresh one-ready extension.

Independent subset tables and literal fixed-start/permutation controls.
The all-extension assertions are proved in NOTE.md, not inferred from the
one-birth finite check. No removable assertions carry the verification.
"""
from pathlib import Path
from itertools import combinations, permutations
from math import factorial
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def graph(word):
    rows = [2, 0] + [0] * len(word)
    for i, c in enumerate(map(int, word), 2):
        for p in (0, 1):
            if c >> p & 1:
                rows[p] |= 1 << i
            else:
                rows[i] |= 1 << p
        for j in range(2, i):
            rows[j] |= 1 << i
    return tuple(rows)


def extend(rows, out_new):
    n = len(rows)
    return tuple(x | (0 if out_new >> i & 1 else 1 << n)
                 for i, x in enumerate(rows)) + (out_new,)


def birth(rows, root):
    return extend(rows, ((1 << len(rows)) - 1) ^ (1 << root))


def tables(rows):
    n = len(rows)
    end = [[0] * n for _ in range(1 << n)]
    start = [[0] * n for _ in range(1 << n)]
    h = [1] + [0] * ((1 << n) - 1)
    for v in range(n):
        end[1 << v][v] = start[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                if rest:
                    end[mask][v] = sum(end[rest][u] for u in range(n)
                                       if rows[u] >> v & 1)
                    start[mask][v] = sum(start[rest][u] for u in range(n)
                                         if rows[v] >> u & 1)
        h[mask] = sum(end[mask])
    full = (1 << n) - 1
    f = [sum(end[mask][v] * h[full ^ mask]
             for mask in range(1, full + 1)) for v in range(n)]
    return h, end, start, f


def source(m):
    b = graph('2' + '0' * (m - 1) + '32')
    first = birth(b, 1)
    c = birth(first, len(b))
    names = ['c_' + str(m), 'x', 'y']
    names += ['c_' + str(i) for i in range(1, m)]
    names += ['birth_a', 'birth_b', 'reset_root', 'reset_child']
    return b, first, c, names


def literal_start(rows, v):
    return sum(all(rows[u] >> w & 1 for u, w in zip(p, p[1:]))
               for tail in permutations(i for i in range(len(rows)) if i != v)
               for p in [(v,) + tail])


def literal_f(rows, v):
    total = 0
    for p in permutations(range(len(rows))):
        cut = p.index(v)
        total += all(rows[p[i]] >> p[i+1] & 1
                     for i in range(len(rows)-1) if i != cut)
    return total


def surj(k, r):
    # Inclusion-exclusion independently evaluates r! S(k,r).
    from math import comb
    return sum((-1)**(r-j) * comb(r, j) * j**k for j in range(r+1))


def main():
    controls = {'endpoint_subset_checks': 0, 'labelled_tournaments': 0,
                'one_birth_old_port_checks': 0, 'literal_fixed_start_checks': 0}
    # Every labelled tournament through order five, every induced subset,
    # every retained endpoint. The conditional lower bound has a zero case.
    for n in range(1, 6):
        pairs = list(combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            rows = [0] * n
            for k, (u, v) in enumerate(pairs):
                a, b = (u, v) if bits >> k & 1 else (v, u)
                rows[a] |= 1 << b
            _, _, starts, _ = tables(rows)
            for mask in range(1, 1 << n):
                for v in range(n):
                    if mask >> v & 1:
                        need(starts[-1][v] == 0 or starts[-1][v] >= starts[mask][v],
                             'conditional endpoint insertion monotonicity')
                        controls['endpoint_subset_checks'] += 1
            if n <= 4:
                for v in range(n):
                    need(literal_start(rows, v) == starts[-1][v], 'literal endpoint DP')
                    controls['literal_fixed_start_checks'] += 1
            controls['labelled_tournaments'] += 1
    # A new universal source kills the old endpoint count: the zero
    # alternative cannot be omitted from the theorem.
    zero_control = extend((2, 0), 3)
    need(tables(zero_control)[2][-1][0] == 0, 'zero endpoint hostile')

    families = []
    skeleton_values = []
    for m in range(2, 7):
        b, first, c, names = source(m)
        hb, _, sb, _ = tables(b)
        hf, _, sf, _ = tables(first)
        h, _, starts, f = tables(c)
        p, q = len(b), len(b)+1
        value = 15*2**m - 1
        need(h[-1] == value, 'balanced H')
        need(sb[-1][1] == 3, 'old x has three starts')
        need(hb[((1 << len(b))-1) ^ 2] == 3, 'B minus x is TT_m then C3')
        need(sf[-1][p] == 3*2**m + 2, 'first newborn start count')
        indegrees = [sum(row >> v & 1 for row in c) for v in range(len(c))]
        need(indegrees[p] == indegrees[q] == 1, 'the two low-indegree roots')
        need(all(indegrees[v] >= 3 for v in range(len(c)) if v not in (p, q)),
             'all other old roots permanently excluded by indegree')
        triangle = [0, len(b)-2, len(b)-1]
        need(all(c[u] >> v & 1 for u, v in zip(triangle, triangle[1:]+triangle[:1])),
             'fixed old triangle avoids x and both reset names')
        need(f[p] == value and f[q]-value == 3*2**m+2, 'reset gaps')
        need(starts[-1][0] == 1, 'retained c_m has a unique start')
        need(f[0] == 26*3**m + 8*2**m - 25, 'mixed exponential F formula')
        skeleton_values.append(f[0])
        if m <= 4:
            need(literal_start(c, 0) == 1, 'literal retained unique start')
            controls['literal_fixed_start_checks'] += 1
        if m == 2:
            need(literal_f(c, 0) == f[0], 'literal two-bin control')
        d = birth(c, 0)
        hd, _, sd, fd = tables(d)
        predicted = 26*3**m + 23*2**m - 27
        need(hd[-1] == predicted and fd[-1]-hd[-1] == 1, 'fresh ready root')
        need(predicted > 3*value+4, 'fresh ready branch exceeds both target budgets')
        families.append({'m': m, 'H': value, 'rows': c, 'names': names,
                         'profiles': [{'name': names[v], 'indegree': indegrees[v],
                                       'start': starts[-1][v], 'gap': f[v]-value}
                                      for v in range(len(c))],
                         'old_triangle': [names[v] for v in triangle],
                         'fresh_ready': {'H': predicted, 'rows': d,
                                         'names': names + ['fresh_ready_root'],
                                         'root': len(c), 'in_port': 0,
                                         'gap': fd[-1]-hd[-1],
                                         'attachment_out_mask': d[-1]}})
        if m <= 3:
            for out_new in range(1 << len(c)):
                child = extend(c, out_new)
                hc, _, _, fc = tables(child)
                need(all(fc[v]-hc[-1] not in (1, 2) for v in range(len(c))),
                     'complete one-birth old-root obstruction')
                controls['one_birth_old_port_checks'] += len(c)
    coefficients = []
    for k, value in enumerate(skeleton_values, 1):
        residual = value - sum(coefficients[r-1]*surj(k, r) for r in range(1, k))
        need(residual % factorial(k) == 0 and residual >= 0, 'positive skeleton coefficient')
        coefficients.append(residual // factorial(k))
    need(coefficients == [241, 250, 78, 0, 0], 'five-run finite skeleton certificate')
    sources = [
        '05-knowledge/results/orbit_spectrum_session_20261103/coverage/NOTE.md',
        '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ADAPTIVE.md',
        '05-knowledge/results/two_entry_funnels_20261103/READY.md',
        '05-knowledge/results/session124_tournament_cut_response_20261031/README.md']
    print(json.dumps({'status': 'PASS', 'scope': 'All labelled tournaments through5 for endpoint lemma; every one-birth cut for C_2,C_3; five fixed exterior skeleton evaluations; all-size conclusions require NOTE proofs.',
                      'controls': controls, 'families': families,
                      'F_c_m_skeleton_coefficients': coefficients,
                      'zero_endpoint_hostile': zero_control,
                      'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                      'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
