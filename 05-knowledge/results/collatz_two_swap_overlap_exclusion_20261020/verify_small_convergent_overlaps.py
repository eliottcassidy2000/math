#!/usr/bin/env python3
"""Exact log intervals and all actual overlap candidates below weight 4096."""
from fractions import Fraction as F
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


def logarithm(t, terms=80):
    lower = 2*sum((t**(2*j+1)/F(2*j+1)
                   for j in range(terms)), F())
    upper = lower+2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
    return lower, upper


L2, U2 = logarithm(F(1, 3))
L3, U3 = logarithm(F(1, 2))
need(L2 > F(2, 3), 'log2 strict lower bound')
N = 2**24
need(N == 4096**2 and 4095**2 < N, 'Legendre magnitude/weight relation')
need(F(3, 2*(3*N+2)) < F(1, 2*N), 'negative integer-minimum bound')
lo, hi = L3/U2, U3/L2
hprev, h, kprev, k = 0, 1, 1, 0
convergents = []
prefix = []
for step in range(50):
    a = lo.numerator//lo.denominator
    need(a == hi.numerator//hi.denominator, 'ambiguous partial quotient')
    prefix.append(a)
    hprev, h = h, a*h+hprev
    kprev, k = k, a*k+kprev
    q, p = h, k
    if p >= 4096:
        terminal = {'q': q, 'p': p}
        break
    if p >= 2 and p < q < 2*p:
        convergents.append((p, q))
    lo, hi = 1/(hi-a), 1/(lo-a)
else:
    raise ArithmeticError('no terminal denominator')
need(convergents == [(2, 3), (5, 8), (12, 19), (41, 65),
                     (53, 84), (306, 485), (665, 1054)],
     'small convergent universe')
need(terminal == {'q': 24727, 'p': 15601}, 'terminal convergent')


def direct_carry(word):
    b = 0
    for i, bit in enumerate(word):
        b = (3 if bit else 1)*b+bit*2**i
    positions = [i for i, bit in enumerate(word) if bit]
    positional = sum(2**i*3**(len(positions)-j-1)
                     for j, i in enumerate(positions))
    need(b == positional, 'carry algorithms disagree')
    return b


def path(source, edges):
    word = source[:]
    for i in edges:
        j = (i+1) % len(word)
        need(word[i] != word[j], 'nonactual edge')
        word[i], word[j] = word[j], word[i]
    return word


def all_actual_overlaps(source):
    result = set()
    q = len(source)
    for i in range(q):
        if source[i] == source[(i+1) % q]:
            continue
        intermediate = source[:]
        intermediate[i], intermediate[(i+1) % q] = (
            intermediate[(i+1) % q], intermediate[i])
        for j in {(i-1) % q, (i+1) % q}:
            if intermediate[j] != intermediate[(j+1) % q]:
                result.add((i, j))
    return result


rows = []
hits = []
trace = sha256()
# The known -17 cycle is not subject to the new-cycle minimum bound.
shapes = sorted(set(convergents+[(7, 11)]))
for p, q in shapes:
    r = q-p
    need(gcd(p, r) == 1 and 1 <= r < p, 'candidate shape hypothesis')
    D = 2**q-3**p
    source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
    labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q for j in range(p)}
    predicted = set()
    count = 0
    maximum_content = 0
    for family, residues in [('plus', range(min(r, p-r))),
                             ('minus', range(max(0, p-2*r), p-r))]:
        for u in residues:
            edges = ([labels[u], labels[u+r]] if family == 'plus'
                     else [(labels[u+r]-1) % q, (labels[u]-1) % q])
            predicted.add(tuple(edges))
            word = path(source, edges)
            B = direct_carry(word)
            content = gcd(abs(D), B)
            maximum_content = max(maximum_content, content)
            record = {'p': p, 'q': q, 'r': r, 'family': family,
                      'u': u, 'edges': edges, 'B': B, 'D': D,
                      'content': content, 'word': ''.join(map(str, word))}
            trace.update(json.dumps(record, sort_keys=True,
                                    separators=(',', ':')).encode()+b'\n')
            if content == abs(D):
                fixed = B//D
                value = fixed
                orbit = []
                for bit in word:
                    need(value % 2 == bit, 'integral branch parity')
                    orbit.append(value)
                    value = (3*value+1)//2 if bit else value//2
                need(value == fixed, 'integral orbit closure')
                hits.append(dict(record, fixed_point=fixed, orbit=orbit))
            if p >= 3:
                need(content < abs(D), 'unexcluded integral overlap')
            count += 1
    need(predicted == all_actual_overlaps(source),
         'candidate legality not exhaustive')
    rows.append({'p': p, 'q': q, 'r': r, 'clock': D, 'paths': count,
                 'maximum_content': maximum_content,
                 'reason_included': 'known negative 17-cycle' if (p, q) == (7, 11)
                                    else 'small ordinary convergent'})
need(sum(row['paths'] for row in rows) == 906, 'candidate path count')
need([(x['p'], x['q'], x['family'], x['fixed_point']) for x in hits]
     == [(2, 3, 'plus', -7), (2, 3, 'minus', -5)],
     'integer exception classification')
print(json.dumps({
    'status': 'FINITE-EXACT complete candidate check after signed descent reduction',
    'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks, 'log_series_terms': 80,
    'continued_fraction_prefix': prefix, 'terminal_convergent': terminal,
    'small_convergent_shapes': [{'p': p, 'q': q} for p, q in convergents],
    'candidate_shapes': rows, 'total_paths': 906, 'integer_hits': hits,
    'all_case_trace_sha256': trace.hexdigest(),
    'dependency': 'The separate complete signed-descent certificate through 2^24 '
                  'is required to make this reduced candidate universe complete.',
    'scope': 'Primitive k=1 overlapping two-one chains, not arbitrary two edits.'
}, indent=2))
