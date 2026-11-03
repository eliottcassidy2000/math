#!/usr/bin/env python3
"""Independent arithmetic/tape execution and literal adjacency-matrix replay.
Does not import check.py. All checks survive python -O.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
C = json.loads((ROOT / 'certificate.json').read_text())


def check(value, message):
    if not value:
        raise ValueError(message)


def tournament(size, mask):
    matrix = [[False for _ in range(size)] for _ in range(size)]
    offset = 0
    for u in range(size):
        for v in range(u + 1, size):
            matrix[u][v] = bool(mask & (2**offset))
            matrix[v][u] = not matrix[u][v]
            offset += 1
    return matrix


def join(blocks):
    size = sum(len(block) for block in blocks)
    matrix = [[u < v for v in range(size)] for u in range(size)]
    offset = 0
    for block in blocks:
        for u in range(len(block)):
            for v in range(len(block)):
                matrix[offset + u][offset + v] = block[u][v]
        offset += len(block)
    return matrix


def graph(bits):
    return join([tournament(5, 8 if bit == 0 else 10) for bit in bits])


def digest(matrix):
    row_integers = [sum((2**v) for v, edge in enumerate(row) if edge) for row in matrix]
    return hashlib.sha256(json.dumps(row_integers, separators=(',', ':')).encode()).hexdigest()


def replay(matrix, change):
    u, v = change['arc']
    check(u != v and matrix[u][v] != matrix[v][u], 'invalid arc')
    matrix[u][v], matrix[v][u] = matrix[v][u], matrix[u][v]
    check(digest(matrix) == change['graph_sha256'], 'literal arc replay differs from certificate')


def direct_profile(matrix):
    n = len(matrix)
    counts = [0] * (n + 1)
    for word in itertools.permutations(range(n)):
        mandatory = sum(not matrix[word[i]][word[i + 1]] for i in range(n - 1))
        free = n - 1 - mandatory
        for number_paths in range(mandatory + 1, n + 1):
            counts[number_paths] += math.comb(free, number_paths - 1 - mandatory)
    for k in range(1, n + 1):
        check(counts[k] % math.factorial(k) == 0, 'nonintegral profile')
        counts[k] //= math.factorial(k)
    return counts


for seed in C['seed']:
    check(direct_profile(tournament(5, seed['mask'])) == seed['profile'], 'seed profile disagreement')
triangle = [[False, True, False], [False, False, True], [True, False, False]]
check(direct_profile(triangle) == [0, 3, 3, 1], 'triangle separator profile')

# Arithmetic reference implementation never uses the exponent-step code.
fractions = [(3, 2), (5, 3)]
value = 12
observed = [value]
fr = C['fractran_demo']
primes, capacities = fr['primes'], fr['capacities']
exponents = fr['trace'][0]['source']
bits = [bit for e, cap in zip(exponents, capacities) for bit in ([1]*e + [0]*(cap-e))]
matrix = graph(bits)
for step in fr['trace']:
    choices = [i for i, (a, b) in enumerate(fractions) if value % b == 0]
    if not choices:
        check(step['status'] == 'halt', 'reference arithmetic halt')
        break
    i = choices[0]
    check(step['status'] == 'step' and step['rule'] == i, 'reference arithmetic priority')
    a, b = fractions[i]
    value = value // b * a
    observed.append(value)
    for change in step['flips']:
        replay(matrix, change)
    target_bits = []
    for p, cap in zip(primes, capacities):
        remaining, e = value, 0
        while remaining % p == 0:
            e += 1
            remaining //= p
        target_bits += [1]*e + [0]*(cap-e)
    check(matrix == graph(target_bits), 'arithmetic matrix endpoint')
check(observed == fr['integers'], 'arithmetic macrotrace')

# Independent dictionary tape with signed coordinates; no use of encoded cell fields.
tm = C['tm_demo']
tape, cursor, state = {}, 0, 'A'
coordinates = tm['tape_coordinates']
transitions = {('A', 0): ('B', 1, 1), ('A', 1): ('B', 1, -1),
               ('B', 0): ('A', 1, -1), ('B', 1): ('H', 1, 1)}
states = ['A', 'B', 'H']
initial_bits = []
for position in coordinates:
    initial_bits += [0] + [int(position == cursor and q == state) for q in states]
matrix = graph(initial_bits)
check(digest(matrix) == tm['trace'][0]['graph_sha256'], 'TM initial matrix')
tape_flips = 0
for record in tm['trace'][1:]:
    new_state, symbol, move = transitions[(state, tape.get(cursor, 0))]
    tape[cursor] = symbol
    cursor += move
    state = new_state
    check(record['tape'] == [tape.get(i, 0) for i in coordinates], 'reference TM tape')
    check(coordinates[record['head']] == cursor and record['state'] == state, 'reference TM control')
    for change in record['flips']:
        replay(matrix, change)
        tape_flips += 1
    bits = []
    for position in coordinates:
        bits += [tape.get(position, 0)] + [int(position == cursor and q == state) for q in states]
    check(matrix == graph(bits), 'TM matrix endpoint')
check(state == 'H' and sum(tape.values()) == 4 and tape_flips == 16, 'TM result')

# Explicit growth: insert a blank B_0 before the final separator of one full register.
before = join([triangle, tournament(5, 10), triangle])
padded = join([triangle, tournament(5, 10), tournament(5, 8), triangle])
old_to_new = [i if i < 8 else i + 5 for i in range(11)]
check(all(before[u][v] == padded[old_to_new[u]][old_to_new[v]]
          for u in range(11) for v in range(11)), 'growth changed an old arc')
after_increment = join([triangle, tournament(5, 10), tournament(5, 10), triangle])
changes = [(u, v) for u in range(16) for v in range(u+1, 16)
           if padded[u][v] != after_increment[u][v]]
check(changes == [(8, 10)], 'growth increment not one literal flip')

# Check all four n=10 profiles against the exact integer falling-basis product,
# independently reconstructed from the five-vertex polynomial.
for item in C['two_slot_products']:
    for r in range(11):
        falling = 1
        value = item['profile'][0]
        for j in range(1, 11):
            falling *= r-j+1
            value += item['profile'][j] * falling
        check(value == (r**5 + 6*r**3 + 2*r)**2, 'two-slot polynomial')

result = {'seed_profiles_recomputed': 2, 'separator_profile_recomputed': True,
          'fractran_arcs_replayed': 10, 'tm_arcs_replayed': tape_flips,
          'growth_old_order': 11, 'growth_new_order': 16,
          'growth_old_arcs_unchanged': True, 'growth_increment_arc': [8, 10],
          'reference_arithmetic_trace': observed, 'reference_tm_steps': 6}
(ROOT / 'independent_certificate.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
