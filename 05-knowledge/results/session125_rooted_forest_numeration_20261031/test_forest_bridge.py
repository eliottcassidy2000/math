"""Exhaustive independent controls. Explicit failures survive Python -O."""
import hashlib
import itertools
import json
import math
import pathlib
import re
import sys
from collections import Counter
from fractions import Fraction
import forest_bridge as b

checks = 0

def check(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)


def tilings(n):
    if n == 0:
        yield ()
    if n >= 1:
        for t in tilings(n - 1):
            yield (1,) + t
    if n >= 2:
        for t in tilings(n - 2):
            yield (2,) + t


def forests(n):
    # Independent cut-mask enumeration, with every possible root in each block.
    for cuts in range(1 << (n - 1)):
        ends = [0] + [i + 1 for i in range(n - 1) if cuts >> i & 1] + [n]
        lengths = [y - x for x, y in zip(ends, ends[1:])]
        for roots in itertools.product(*(range(1, k + 1) for k in lengths)):
            yield tuple(zip(lengths, roots))


def forest_count_dsu(n, edges):
    total = 0
    for mask in range(1 << len(edges)):
        parent = list(range(n))
        sizes = [1] * n
        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for i, (u, v) in enumerate(edges):
            if mask >> i & 1:
                ru, rv = find(u), find(v)
                if ru == rv:
                    break
                parent[ru] = rv
                sizes[rv] += sizes[ru]
        else:
            total += math.prod(sizes[i] for i in range(n) if parent[i] == i)
    return total


def det_leibniz(a):
    total = 0
    for perm in itertools.permutations(range(len(a))):
        inversions = sum(perm[i] > perm[j] for i in range(len(a)) for j in range(i+1, len(a)))
        total += (-1)**inversions * math.prod(a[i][perm[i]] for i in range(len(a)))
    return total


def det_rational(a):
    a = [list(map(Fraction, row)) for row in a]
    result = Fraction(1)
    for k in range(len(a)):
        pivot = next(i for i in range(k, len(a)) if a[i][k])
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            result = -result
        p = a[k][k]
        result *= p
        for i in range(k+1, len(a)):
            q = a[i][k] / p
            for j in range(k+1, len(a)):
                a[i][j] -= q*a[k][j]
    return result


def fib(n):
    x, y = 0, 1
    for _ in range(n):
        x, y = y, x+y
    return x


path_rows = []
word_total = tiling_total = forest_total = 0
for length in range(12):
    n = length + 1
    legal_words = set()
    for word in itertools.product(range(3), repeat=length):
        word_total += 1
        expected = re.search('21*2', ''.join(map(str, word))) is None
        check(b.legal(word) == expected, ('legal', word))
        if expected:
            legal_words.add(word)
            t = b.word_to_tiling(word)
            f = b.word_to_forest(word)
            check(b.tiling_to_word(t) == word, ('word-tiling inverse', word))
            check(b.forest_to_word(f) == word, ('word-forest inverse', word))
            check(len(t) + t.count(2) == 2*n-1, ('tiling length', word))
            check(sum(k for k, r in f) == n, ('forest order', word))
            check(sum(d != 0 for d in word) == n-len(f) == t.count(2), ('edge statistic', word))
            check(b.word_to_forest((0,) + word) == ((1,1),) + f, ('leading zero singleton', word))
            check(b.word_to_tiling(word + (0,)) == t + (1,1), ('right zero monomers', word))
            check(b.word_to_forest(word + (0,)) == f + ((1,1),), ('right zero singleton', word))
            if word:
                check(bool(word[0]) == (f[0][0] != 1), ('canonical first component', word))
            if word and word[0] and sum(d != 0 for d in word) == 1:
                expected_f = ((2, word[0]),) + ((1,1),)*(len(word)-1)
                check(f == expected_f, ('one-edge canonical subclass', word))
            value = sum(d*fib(2*i+2) for i,d in enumerate(reversed(word)))
            shifted = sum(d*fib(2*i+4) for i,d in enumerate(reversed(word)))
            beatty = (3*value + math.isqrt(5*value*value))//2 + 1 if value else 0
            check(shifted == beatty, ('right shift Beatty value', word))
    all_t = set(tilings(2*n-1))
    all_f = set(forests(n))
    tiling_total += len(all_t)
    forest_total += len(all_f)
    mapped_t, mapped_f = set(), set()
    for t in all_t:
        word = b.tiling_to_word(t)
        mapped_t.add(word)
        check(b.word_to_tiling(word) == t, ('tiling inverse', t))
        check(b.forest_to_tiling(b.tiling_to_forest(t)) == t, ('tiling forest inverse', t))
    for f in all_f:
        mapped_f.add(b.forest_to_word(f))
        check(b.tiling_to_forest(b.forest_to_tiling(f)) == f, ('forest tiling inverse', f))
    check(mapped_t == legal_words == mapped_f, ('onto and injection', n))
    check(len(all_t) == len(all_f) == len(legal_words) == fib(2*n), ('counts', n))
    values = [sum(d*fib(2*i+2) for i,d in enumerate(reversed(w))) for w in sorted(legal_words)]
    check(values == list(range(fib(2*n))), ('numeric rank', n))
    hist = Counter(sum(d != 0 for d in w) for w in legal_words)
    check(hist == {k: math.comb(2*n-1-k, k) for k in range(n)}, ('weighted formula', n))
    for t in range(4):
        a = b.identity_laplacian(n, [(i, i+1) for i in range(n-1)])
        a = [[t*(a[i][j]-int(i == j))+int(i == j) for j in range(n)] for i in range(n)]
        check(b.det_int(a) == sum(count*t**k for k, count in hist.items()), ('weighted determinant', n, t))
    path_rows.append({'vertices': n, 'word_length': length, 'count': len(all_f),
                      'edge_coefficients': [hist[k] for k in range(n)]})

# Entire unfiltered universe of labelled simple graphs of order <= 5.
graph_rows = []
all_counts = {}
graph_total = 0
for n in range(6):
    pairs = tuple(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        edges = tuple(e for i, e in enumerate(pairs) if mask >> i & 1)
        dsu = forest_count_dsu(n, edges)
        det = b.det_int(b.identity_laplacian(n, edges))
        check(det == dsu, ('matrix forest', n, mask))
        check(det == det_leibniz(b.identity_laplacian(n, edges)), ('determinant algorithm', n, mask))
        all_counts[(n, edges)] = dsu
        graph_total += 1
    graph_rows.append({'vertices': n, 'all_labelled_graphs': 1 << len(pairs)})

arrival_total = 0
for n in range(5):
    pairs = tuple(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        edges = tuple(e for i, e in enumerate(pairs) if mask >> i & 1)
        old = all_counts[(n, edges)]
        for port_mask in range(1 << n):
            ports = tuple(i for i in range(n) if port_mask >> i & 1)
            daughter = tuple(sorted(edges + tuple((v, n) for v in ports)))
            expected = all_counts[(n+1, daughter)]
            check(b.arrival_count(n, edges, ports) == expected, ('all arrivals', n, mask, ports))
            arrival_total += 1
            a = b.identity_laplacian(n, edges)
            columns = [b.solve(a, [int(i == v) for i in range(n)]) for v in ports]
            h = [[columns[j][v] + int(i == j) for j in range(len(ports))] for i,v in enumerate(ports)]
            h_ones = b.solve(h, [1]*len(ports))
            green_prediction = old * det_rational(h) * (1 + sum(h_ones))
            check(green_prediction == expected, ('boundary Green response', n, mask, ports))
            if len(ports) == 1:
                a = b.identity_laplacian(n, edges)
                v = ports[0]
                minor = [[a[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
                check(expected == 2*old + b.det_int(minor), ('leaf response', n, mask, v))

p3 = ((0, 1), (1, 2))
a = b.identity_laplacian(3, p3)
response = [b.solve(a, [int(i == v) for i in range(3)])[v] for v in range(3)]
check(response == [Fraction(5,8), Fraction(1,2), Fraction(5,8)], 'P3 Green diagonal')
check(b.arrival_count(3, p3, [0]) == 21, 'endpoint leaf')
check(b.arrival_count(3, p3, [1]) == 20, 'middle leaf')
check(b.det_int([]) == 1 and fib(0) == 0, 'empty versus formal initial')

def terminal_lock(word):
    locked = False
    for digit in word:
        if digit == 0:
            locked = False
        elif digit == 2:
            if locked:
                raise ValueError('illegal CG transition')
            locked = True
    return locked

check(not terminal_lock((2, 0, 1)) and terminal_lock((1, 0, 2)),
      '201 MSB/LSD terminal state distinction')
check(b.word_to_tiling((2, 0, 1)) == (2, 1, 1, 1, 2), '201 oriented tiling')
check(b.word_to_tiling((1, 0, 2)) == (1, 2, 1, 2, 1), '102 oriented tiling')


square_counts = {}
for n in (12, 30, 31, 32, 33):
    edges = b.square_sum_edges(n)
    square_counts[n] = b.det_int(b.identity_laplacian(n, edges))
check(square_counts[12] == 9240, 'G12 factorization')
arrivals = []
for n in (31, 32, 33):
    old_edges = b.square_sum_edges(n-1)
    new_edges = b.square_sum_edges(n)
    ports = [u for u, v in new_edges if v == n-1]
    count = b.arrival_count(n-1, old_edges, ports)
    check(count == square_counts[n], ('square arrival', n))
    arrivals.append({'new_label': n, 'ports_one_based': [v+1 for v in ports], 'rooted_forests': count})

for word in ((2,2), (2,1,2), (3,), (-1,), (1,2,1,2,0)):
    try:
        b.word_to_tiling(word)
    except ValueError:
        checks += 1
    else:
        raise RuntimeError(('invalid word accepted', word))
for t in ((), (2,), (1,1), (3,), (1,0)):
    for converter in (b.tiling_to_word, b.tiling_to_forest):
        try:
            converter(t)
        except ValueError:
            checks += 1
        else:
            raise RuntimeError(('invalid tiling accepted', t))

samples = []
for word in itertools.product(range(3), repeat=2):
    if b.legal(word):
        samples.append({'word': ''.join(map(str, word)), 'value': 3*word[0]+word[1],
                        'tiling': ''.join('M' if t == 1 else 'D' for t in b.word_to_tiling(word)),
                        'forest_blocks': b.word_to_forest(word)})

result = {'status': 'PASS', 'checks': checks, 'assertions_removed_under_O': False,
          'all_ternary_words_lengths_0_through_11': word_total,
          'independently_generated_odd_tilings': tiling_total,
          'independently_generated_rooted_path_forests': forest_total,
          'path_rows': path_rows, 'graph_universe': graph_rows,
          'all_labelled_graphs_orders_0_through_5': graph_total,
          'all_arrival_queries_old_orders_0_through_4': arrival_total,
          'square_sum_rooted_forest_counts': square_counts, 'square_sum_arrivals': arrivals,
          'p3_port_hostile': {'before': 8, 'endpoint_response': '5/8', 'middle_response': '1/2',
                             'endpoint_after': 21, 'middle_after': 20},
          'P3_word_table': samples,
          'source_sha256': {name: hashlib.sha256(pathlib.Path(__file__).with_name(name).read_bytes()).hexdigest()
                            for name in ('forest_bridge.py', 'test_forest_bridge.py')}}
output = pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).with_name('verification.normal.json')
output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
