#!/usr/bin/env python3
"""Exact 64-state norm-observer control; no producer imports or dependencies."""
import hashlib
import json
from pathlib import Path

checks = 0

def require(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)


def mul(a, b):
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        a <<= 1
        if a & 64:
            a ^= 73  # X^6+X^3+1
        b >>= 1
    return answer


def power(a, n):
    result = 1
    while n:
        if n & 1:
            result = mul(result, a)
        a = mul(a, a)
        n >>= 1
    return result


def labels(keys):
    seen = {}
    return [seen.setdefault(k, len(seen)) for k in keys]


def refine(partition, operations):
    return labels([tuple([partition[x]] + [partition[f[x]] for f in operations])
                   for x in range(64)])


def trajectory(partition, operations):
    result = [partition]
    while True:
        nxt = refine(result[-1], operations)
        if nxt == result[-1]:
            return result
        result.append(nxt)


def classes(partition):
    return [[x for x in range(64) if partition[x] == i]
            for i in range(len(set(partition))) ]


b = 3  # 1+t in the declared polynomial basis
phi = power(b, 21)
F4 = [x for x in range(64) if power(x, 4) == x]
F8 = [x for x in range(64) if power(x, 8) == x]
require(len(F4) == 4 and len(F8) == 8, 'subfield sizes')
require(len({power(b, k) for k in range(63)}) == 63, 'primitive 63-clock')
require(power(2, 9) == 1 and power(2, 3) != 1, 'hostile t has order 9')
require(F4 == sorted([0, 1, phi, phi ^ 1]), 'F4 explicit')
require(b not in F4, 'multiplier is not scalar')
for x in range(64):
    require(power(x, 64) == x, 'F64 Frobenius')
    if x:
        require(mul(x, power(x, 62)) == 1, 'every nonzero element invertible')
    for y in range(64):
        require(mul(x, y) == mul(y, x), 'commutativity')
        require(mul(x, y ^ 1) == (mul(x, y) ^ x), 'distributivity control')

norms = [(power(x, 21), power(x, 9)) for x in range(64)]
p0 = labels(norms)
require(len(set(norms)) == 22, '21 norm values plus zero')
require(sorted(map(len, classes(p0))) == [1] + [3]*21, 'norm fibres')
for x in range(64):
    for y in range(64):
        scalar_equiv = (x == y == 0) or (x != 0 and y in [mul(c, x) for c in F4 if c])
        require((norms[x] == norms[y]) == scalar_equiv, 'norm fibres equal scalar orbits')

M = [mul(b, x) for x in range(64)]
A = [x ^ 1 for x in range(64)]
m_path = trajectory(p0, [M])
a_path = trajectory(p0, [A])
joint_path = trajectory(p0, [M, A])
require([len(set(p)) for p in m_path] == [22], 'M quotient closed')
require([len(set(p)) for p in a_path] == [22, 63], 'A repair has one strict round')
require([len(set(p)) for p in joint_path] == [22, 63, 64], 'joint repair has two strict rounds')
remaining = [c for c in classes(joint_path[1]) if len(c) > 1]
require(remaining == [sorted([phi, mul(phi, phi)])], 'unique unresolved pair')
require([c for c in classes(a_path[-1]) if len(c) > 1] == remaining, 'A stable residual pair')
require(norms[M[phi]] == norms[M[mul(phi, phi)]], 'one multiplication alone does not distinguish')
require(norms[A[M[phi]]] != norms[A[M[mul(phi, phi)]]], 'future word A after M distinguishes')
require(len({norms[power(b, k)] for k in range(63)}) == 21, 'visible multiplicative period 21')
require(norms[1] == norms[phi] and 1 != phi, 'hostile return observer loses phase')
# Direct bounded-word observation is a second construction of each partition.
word_functions = [list(range(64))]
for depth, p in enumerate(joint_path):
    direct = labels([tuple(norms[f[x]] for f in word_functions) for x in range(64)])
    require(direct == p, 'refinement equals all observed futures through depth')
    word_functions += [[op[f[x]] for x in range(64)] for f in list(word_functions) for op in (M, A)]

out = {
    'status': 'FINITE-EXACT', 'universe': 'all 64 residues modulo X^6+X^3+1 over F2',
    'checks': checks, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'encoding': 'six low-to-high polynomial coefficient bits; t=2, b=1+t=3',
    'phi': phi, 'F4': F4, 'F8': F8, 'norm_fibre_sizes': [1] + [3]*21,
    'multiplication_counts': [len(set(p)) for p in m_path],
    'addition_counts': [len(set(p)) for p in a_path],
    'joint_counts': [len(set(p)) for p in joint_path],
    'unique_pair_after_first_round': remaining[0],
    'distinguishing_word': 'x -> b*x -> b*x+1; then observe both norms',
    'partitions': joint_path,
    'limitations': ['not a native tournament action', 'finite result does not bound general refinement depth']
}
print(json.dumps(out, sort_keys=True, indent=2))
