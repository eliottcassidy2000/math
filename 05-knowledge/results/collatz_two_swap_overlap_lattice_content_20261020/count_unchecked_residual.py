#!/usr/bin/env python3
"""Count, but do not test carries in, the finite overlap residual."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, gcd
from pathlib import Path
import json

checks = 0


def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ArithmeticError(message)


limit = 4096
need(sum((F(23, 6)**j/factorial(j) for j in range(20)), F()) > 46,
     'positive Taylor certificate for log46')
need(F(13, 12)-F(23, 24) == F(1, 8) > F(1, 9),
     'exponential margin')
need(16**3 == limit, 'cube threshold')

# Directly count the exact two legality intervals.
direct_shapes = direct_paths = 0
for p in range(2, limit):
    for r in range(1, p):
        if gcd(p, r) == 1:
            direct_shapes += 1
            direct_paths += 2*min(r, p-r)

# Independently use Euler phi and Mobius inversion of weighted coprime
# residue sums. Neither count generates a source word or branch carry.
phi = list(range(limit))
mu = [1]*limit
for prime in range(2, limit):
    if phi[prime] != prime:
        continue
    for n in range(prime, limit, prime):
        phi[n] -= phi[n]//prime
        mu[n] *= -1
    for n in range(prime*prime, limit, prime*prime):
        mu[n] = 0
weighted_half_sums = [0]*limit
for d in range(1, limit):
    if not mu[d]:
        continue
    for n in range(d, limit, d):
        k = n//(2*d)
        weighted_half_sums[n] += mu[d]*d*k*(k+1)//2
independent_shapes = sum(phi[2:])
independent_paths = 2+4*sum(weighted_half_sums[3:])
need(direct_shapes == independent_shapes == 5097971, 'shape count')
need(direct_paths == independent_paths == 6959303782, 'ordered path count')

print(json.dumps({
    'status': 'EXACT COUNTS ONLY, NOT AN INTEGRALITY CENSUS',
    'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks,
    'universe': '2<=p<4096; every coprime 1<=r<p; both exact legal '
                'overlapping two-one parameter ranges',
    'source_shapes': direct_shapes,
    'ordered_paths': direct_paths,
    'source_shapes_after_omitting_p2_exception': direct_shapes-1,
    'ordered_paths_after_omitting_p2_exception': direct_paths-2,
    'independent_count': 'Euler phi plus Mobius inversion of coprime half sums',
    'not_done': 'No branch carry or integer-cycle check in this script. '
                'The large finite residual is not claimed to have been excluded.'
}, indent=2))
