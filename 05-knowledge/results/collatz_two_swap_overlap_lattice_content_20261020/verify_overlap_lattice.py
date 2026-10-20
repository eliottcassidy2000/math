#!/usr/bin/env python3
"""Exact all-r overlap controls, short lattice vectors, and Sylvester resultants.

No floating-point arithmetic, symbolic package, or imported Collatz implementation.
The finite run checks the theorem, not a proposed all-length integer exclusion.
"""
from hashlib import sha256
from math import comb, gcd
from pathlib import Path
import json

checks = 0


def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ArithmeticError(message)


def independent(v, w):
    return any(v[i]*w[j] != v[j]*w[i]
               for i, j in [(0, 1), (0, 2), (1, 2)])


def short_pair(p, r, u):
    first = None
    for height in range(1, p+1):
        for b in range(-height, height+1):
            for c in range(-height, height+1):
                for m in range(-height, height+1):
                    v = (b, c, m)
                    if max(map(abs, v)) != height:
                        continue
                    if (m-b*r-c*u) % p:
                        continue
                    # Exactly one of v,-v is retained.
                    if next(x for x in v if x) < 0:
                        continue
                    if first is None:
                        first = v
                    elif independent(first, v):
                        h1 = max(map(abs, first))
                        need((h1*height)**3 <= p*p,
                             'successive-minimum product bound')
                        return first, v
    raise ArithmeticError('missing independent lattice pair')


def polynomial(p, r, u, family, vector):
    b, c, m = vector
    if c < 0:
        b, c, m = -b, -c, -m
    a, remainder = divmod(m-b*r-c*u, p)
    need(remainder == 0, 'not a lattice relation')
    height = max(map(abs, vector))
    need(abs(a) <= height, 'legal range did not bound a')
    kappa = 1 if family == 'plus' else 2
    e = a-b+kappa*c
    left = 2**max(-e, 0)*3**max(-b, 0)*5**max(c, 0)
    sign = -1 if family == 'plus' and c % 2 else 1
    right = sign*2**max(e, 0)*3**max(b, 0)*5**max(-c, 0)
    # F=left*z^mplus*(z-1)^cplus-right*z^mminus*(z-1)^cminus.
    degree = max(max(m, 0)+max(c, 0), max(-m, 0)+max(-c, 0))
    f = [0]*(degree+1)
    for coefficient, shift, power in [
            (left, max(m, 0), max(c, 0)),
            (-right, max(-m, 0), max(-c, 0))]:
        for j in range(power+1):
            f[shift+j] += coefficient*comb(power, j)*(-1)**(power-j)
    while len(f) > 1 and not f[-1]:
        f.pop()
    need(any(f), 'identically zero cleared polynomial')
    need(f[0] != 0 and sum(f) != 0, 'extraneous root zero or one')
    need((len(f) == 1) == (m == 0 and c == 0),
         'constant-polynomial classification')
    need(len(f)-1 <= 2*height, 'degree bound')
    norm = sum(map(abs, f))
    base = 38 if family == 'plus' else 46
    need(norm <= base**height, 'coefficient norm bound')
    return f, dict(a=a, b=b, c=c, m=m, height=height,
                   degree=len(f)-1, coefficient_l1_norm=norm)


def determinant(matrix):
    n = len(matrix)
    if not n:
        return 1
    a = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for k in range(n-1):
        if not a[k][k]:
            pivot = next((j for j in range(k+1, n) if a[j][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                need(numerator % previous == 0, 'Bareiss nonexact division')
                a[i][j] = numerator//previous
        for i in range(k+1, n):
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def resultant(f, h):
    m, n = len(f)-1, len(h)-1
    need(m+n > 0, 'two constants cannot be an independent pair')
    rows = []
    for i in range(n):
        rows.append([0]*i+f[::-1]+[0]*(n-1-i))
    for i in range(m):
        rows.append([0]*i+h[::-1]+[0]*(m-1-i))
    return determinant(rows)


def carry(word):
    value = 0
    for i, bit in enumerate(word):
        value = (3 if bit else 1)*value+bit*2**i
    positional = sum(2**i*3**sum(word[i+1:])
                     for i, bit in enumerate(word) if bit)
    need(value == positional, 'carry methods disagree')
    return value


def move(source, edges):
    word = source[:]
    for i in edges:
        j = (i+1) % len(word)
        need(word[i] != word[j], 'nonactual adjacent swap')
        word[i], word[j] = word[j], word[i]
    return word


def actual_overlaps(source):
    q = len(source)
    paths = set()
    for i in range(q):
        if source[i] == source[(i+1) % q]:
            continue
        intermediate = source[:]
        intermediate[i], intermediate[(i+1) % q] = (
            intermediate[(i+1) % q], intermediate[i])
        for j in {(i-1) % q, (i+1) % q}:
            if intermediate[j] != intermediate[(j+1) % q]:
                paths.add((i, j))
    return paths


def check_case(p, r, u, family, source, labels, forced_pair=None):
    q = p+r
    edges = ([labels[u], labels[u+r]] if family == 'plus'
             else [(labels[u+r]-1) % q, (labels[u]-1) % q])
    word = move(source, edges)
    D = 2**q-3**p
    B = carry(word)
    g = gcd(abs(D), B)
    need(gcd(g, 30) == 1, 'content has forbidden denominator prime')
    vectors = short_pair(p, r, u) if forced_pair is None else forced_pair
    need(independent(*vectors), 'dependent input vectors')
    f, data1 = polynomial(p, r, u, family, vectors[0])
    h, data2 = polynomial(p, r, u, family, vectors[1])
    R = resultant(f, h)
    need(R != 0, 'zero resultant counterexample')
    need(R % g == 0, 'full content does not divide resultant')
    need(abs(R) <= data1['coefficient_l1_norm']**data2['degree']
         * data2['coefficient_l1_norm']**data1['degree'],
         'Sylvester-Hadamard bound')
    return {'family': family, 'p': p, 'r': r, 'u': u, 'q': q,
            'D': D, 'B': B, 'content': g, 'edges': edges,
            'word': ''.join(map(str, word)), 'vectors': [data1, data2],
            'polynomials': [f, h], 'resultant': R,
            'clock_resultant_gcd': gcd(abs(D), abs(R))}


trace = sha256()
shapes = 0
paths = 0
counts = {'plus': 0, 'minus': 0}
clock_signs = {'positive': 0, 'negative': 0}
max_content = None
integer_hits = []
extra_clock_factors = []
dependent_projection_count = 0
small_r_paths = 0
cut_crossings = 0
for p in range(2, 41):
    for r in range(1, p):
        if gcd(p, r) != 1:
            continue
        q = p+r
        source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
        labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q
                  for j in range(p)}
        predicted = set()
        shapes += 1
        for family, residues in [
                ('plus', range(min(r, p-r))),
                ('minus', range(max(0, p-2*r), p-r))]:
            for u in residues:
                row = check_case(p, r, u, family, source, labels)
                predicted.add(tuple(row['edges']))
                trace.update(json.dumps(row, sort_keys=True,
                                        separators=(',', ':')).encode()+b'\n')
                paths += 1
                counts[family] += 1
                clock_signs['positive' if row['D'] > 0 else 'negative'] += 1
                small_r_paths += 2*r < p
                cut_crossings += q-1 in row['edges']
                x, y = row['vectors']
                dependent_projection_count += x['m']*y['c'] == y['m']*x['c']
                if max_content is None or row['content'] > max_content['content']:
                    max_content = row
                if row['content'] == abs(row['D']):
                    integer_hits.append(row)
                if row['clock_resultant_gcd'] > row['content'] and len(extra_clock_factors) < 5:
                    extra_clock_factors.append(row)
        need(predicted == actual_overlaps(source),
             'all-r legality omits or adds an actual overlapping path')

# Hostile dependent-projection and constant-polynomial cases. They must not
# be silently discarded by an assertion that projected vectors are independent.
controls = []
p, r, u = 5, 3, 1
q = p+r
source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q for j in range(p)}
for pair in [((0, 1, 1), (5, 1, 1)), ((5, 0, 0), (0, 1, 1))]:
    controls.append(check_case(p, r, u, 'plus', source, labels, pair))

primepowers = []
for family, p, r, u, expected in [
        ('plus', 418, 237, 180, 11849), ('minus', 203, 165, 0, 169)]:
    q = p+r
    source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
    labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q for j in range(p)}
    row = check_case(p, r, u, family, source, labels)
    need(row['content'] == expected, 'prime-power depth control')
    primepowers.append(row)

need(len(integer_hits) == 2 and all(x['p'] == 2 for x in integer_hits),
     'finite integer hits changed; this is not an all-length theorem')
print(json.dumps({
    'status': 'FINITE-EXACT support for the all-r overlap content theorem',
    'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks, 'universe': '2<=p<=40, 1<=r<p, gcd(p,r)=1; '
    'every actual nonundo overlapping adjacent-swap path of the source',
    'shapes': shapes, 'paths': paths, 'family_path_counts': counts,
    'path_clock_signs': clock_signs, 'paths_with_r_below_p_over_2': small_r_paths,
    'original_cut_crossing_paths': cut_crossings,
    'dependent_projection_short_pairs': dependent_projection_count,
    'maximum_content': max_content, 'integer_hits': integer_hits,
    'extra_clock_resultant_factor_controls': extra_clock_factors,
    'forced_degeneracy_controls': controls, 'primepower_controls': primepowers,
    'all_case_trace_sha256': trace.hexdigest(),
    'scope_warning': 'The theorem bounds actual content subexponentially. '
                     'This finite census is not an all-length overlap exclusion.'
}, indent=2))
