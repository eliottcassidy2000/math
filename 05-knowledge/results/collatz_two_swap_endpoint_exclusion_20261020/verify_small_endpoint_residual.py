#!/usr/bin/env python3
"""Complete all-clock endpoint residual p<100; no floating-point decisions."""
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

checks = 0


def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ArithmeticError(message)


def affine_carry(word):
    b = 0
    for i, bit in enumerate(word):
        b = (3 if bit else 1) * b + bit * 2**i
    return b


def positional_carry(word):
    positions = [i for i, bit in enumerate(word) if bit]
    return sum(2**i * 3**(len(positions)-j-1)
               for j, i in enumerate(positions))


def apply_path(source, edges):
    word = source[:]
    stages = [''.join(map(str, word))]
    for edge in edges:
        other = (edge + 1) % len(word)
        need(word[edge] != word[other], 'nonactual swap')
        word[edge], word[other] = word[other], word[edge]
        stages.append(''.join(map(str, word)))
    return word, stages


def verify_integer_orbit(word, fixed):
    value = fixed
    orbit = []
    for bit in word:
        need(value % 2 == bit, 'integer branch parity mismatch')
        orbit.append(value)
        value = (3*value+1)//2 if bit else value//2
    need(value == fixed, 'integer orbit fails to close')
    return orbit


trace = sha256()
shapes = 0
paths = 0
signs = {'positive': 0, 'negative': 0}
cut_crossings = 0
nonunit = {'plus': 0, 'minus': 0}
integer_hits = []
shape_counts_by_weight = {}

for p in range(2, 100):
    residues = [1] if p == 2 else range(p//2+1, p)
    for r in residues:
        if gcd(p, r) != 1:
            continue
        q = p+r
        D = 2**q-3**p
        need(D != 0, 'impossible zero clock')
        source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
        need(len(source) == q and sum(source) == p, 'source shape')
        # Label 0 is the final one at q-1, not the first one in this cut.
        labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q
                  for j in range(p)}
        need(len(labels) == p and all(source[i] for i in labels.values()),
             'source labelling')
        plus = [labels[p-1-r], labels[p-1]]
        minus = [(labels[r]-1) % q, (labels[0]-1) % q]
        shapes += 1
        shape_counts_by_weight[p] = shape_counts_by_weight.get(p, 0)+1
        signs['positive' if D > 0 else 'negative'] += 1
        for family, edges in [('plus', plus), ('minus', minus)]:
            word, stages = apply_path(source, edges)
            need(sum(word) == p, 'weight changed')
            B = affine_carry(word)
            need(B == positional_carry(word), 'two carry methods disagree')
            g = gcd(abs(D), B)
            formula = (gcd(2*13**p-10**p, 3*13**r-2*10**r)
                       if family == 'plus'
                       else gcd(9**p-2*5**p, 2*9**r-3*5**r))
            need(g == formula, 'content formula disagreement')
            paths += 1
            cut_crossings += q-1 in edges
            nonunit[family] += g > 1
            row = {'family': family, 'p': p, 'r': r, 'q': q,
                   'D': D, 'B': B, 'content': g, 'edges': edges,
                   'word': stages[-1]}
            trace.update(json.dumps(row, sort_keys=True,
                                    separators=(',', ':')).encode()+b'\n')
            if B % abs(D) == 0:
                fixed = B//D
                orbit = verify_integer_orbit(word, fixed)
                integer_hits.append(dict(row, path=stages,
                                         fixed_point=fixed, orbit=orbit))
            if p >= 3:
                need(g < abs(D), 'unexpected integer endpoint at p>=3')

need(shapes == 1502 and paths == 3004, 'complete small universe count')
need([(x['family'], x['p'], x['r'], x['word'], x['fixed_point'])
      for x in integer_hits]
     == [('plus', 2, 1, '101', -7), ('minus', 2, 1, '110', -5)],
     'small integer exception classification')

print(json.dumps({
    'status': 'FINITE-EXACT complete all-clock endpoint residual',
    'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks,
    'universe': '2<=p<100, q=p+r, gcd(p,r)=1, p/2<r<p for p>=3; '
                'the legal boundary (p,r)=(2,1) is explicitly included',
    'shapes': shapes, 'paths': paths, 'shape_clock_signs': signs,
    'original_cut_crossing_paths': cut_crossings,
    'nonunit_path_counts': nonunit,
    'shape_counts_by_weight': shape_counts_by_weight,
    'integer_hits': integer_hits,
    'all_case_trace_sha256': trace.hexdigest(),
    'scope_warning': 'Only the two specified overlapping endpoint chains; '
                     'not all primitive two-edit words and not only critical clocks.'
}, indent=2))
