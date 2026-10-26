#!/usr/bin/env python3
"""Independent exact finite-group audit by all projective matrices and all permutations."""
from itertools import combinations, permutations, product
from collections import Counter
import json

checks = 0


def require(v, message):
    global checks
    checks += 1
    if not v:
        raise RuntimeError(message)


def normalize(v):
    v = tuple(x % 7 for x in v)
    z = next(x for x in v if x)
    return tuple(x * pow(z, -1, 7) % 7 for x in v)


def mobius(v):
    a, b, c, d = v
    out = []
    for t in range(8):
        numerator, denominator = ((a, c) if t == 7 else ((a * t + b) % 7, (c * t + d) % 7))
        out.append(numerator * pow(denominator, -1, 7) % 7 if denominator else 7)
    return tuple(out)


def compose(g, h):
    return tuple(g[h[x]] for x in range(8))


def closure(gens):
    identity = tuple(range(8))
    seen = {identity}
    todo = [identity]
    while todo:
        h = todo.pop()
        for g in gens:
            gh = compose(g, h)
            if gh not in seen:
                seen.add(gh)
                todo.append(gh)
    return seen


def image(block, g):
    return tuple(sorted(g[x] for x in block))


matrices = {normalize(v) for v in product(range(7), repeat=4)
            if (v[0] * v[3] - v[1] * v[2]) % 7}
group = {}
for v in matrices:
    det = (v[0] * v[3] - v[1] * v[2]) % 7
    g = mobius(v)
    require(len(set(g)) == 8, 'projective map permutes eight points')
    require(g not in group, 'normalized matrix action faithful')
    group[g] = 1 if det in (1, 2, 4) else -1
square = {g for g, sign in group.items() if sign == 1}
require(len(group) == 336 and len(square) == 168, 'projective group sizes')

U, A, D = [mobius(v) for v in [(2, -1, 1, 0), (2, 1, 1, 0), (1, 2, 0, 1)]]
require(closure([U, D]) == square, 'U,D generate exactly PSL')
require(closure([U, A, D]) == set(group), 'all branches generate exactly PGL')
require(compose(U, U)[2] == A[2] == 6 and group[compose(U, U)] != group[A],
        'same endpoint different determinant bit')

blocks = set(combinations(range(8), 4))
orbits = []
while blocks:
    seed = min(blocks)
    orbit = {image(seed, g) for g in square}
    orbits.append(orbit)
    blocks -= orbit
require(sorted(map(len, orbits)) == [14, 14, 42], 'all four-subset orbits')
designs = sorted([o for o in orbits if len(o) == 14], key=lambda x: sorted(x))
plus, minus = designs
require(not plus & minus, 'design block sets disjoint')
for design in designs:
    counts = Counter(triple for block in design for triple in combinations(block, 3))
    require(set(counts) == set(combinations(range(8), 3)) and set(counts.values()) == {1},
            'every triple belongs to exactly one block')
for g, sign in group.items():
    require({image(b, g) for b in plus} == (plus if sign == 1 else minus),
            'determinant exactly controls design transport')
    require({image(b, g) for b in minus} == (minus if sign == 1 else plus),
            'transport both designs')

# Full permutation oracle is independent of matrix generation.
single = ordered = unordered = 0
pair = frozenset([frozenset(plus), frozenset(minus)])
for g in permutations(range(8)):
    gp = frozenset(image(b, g) for b in plus)
    gm = frozenset(image(b, g) for b in minus)
    single += gp == plus
    ordered += gp == plus and gm == minus
    unordered += frozenset([gp, gm]) == pair
require((single, ordered, unordered) == (1344, 168, 336), 'all-permutation automorphism oracle')

# Construct an affine F2^3 marking of one design, without importing a producer table.
affine = {b for b in combinations(range(8), 4) if b[0] ^ b[1] ^ b[2] ^ b[3] == 0}
marking = next((0,) + tail for tail in permutations(range(1, 8))
               if {image(b, (0,) + tail) for b in plus} == affine)
inverse = tuple(marking.index(x) for x in range(8))
for g in square:
    f = compose(marking, compose(g, inverse))
    require(all(f[x ^ y] == (f[x] ^ f[y] ^ f[0]) for x in range(8) for y in range(8)),
            'marked PSL action is affine over F2')

# Enumerate the entire ternary homogeneous quadric, not just its parameter image.
conic = {normalize(v) for v in product(range(7), repeat=3)
         if any(v) and (v[0] * v[0] + v[1] * v[1] - v[2] * v[2]) % 7 == 0}
parametric = {normalize((t * t - 1, 2 * t, t * t + 1)) for t in range(7)}
parametric.add((1, 0, 1))
require(conic == parametric and len(conic) == 8, 'full projective conic equals parameter chart')


def field_mul(x, y, degree, modulus):
    # Literal coefficient convolution, then polynomial long division.
    raw = 0
    for i in range(degree):
        for j in range(degree):
            if (x >> i) & 1 and (y >> j) & 1:
                raw ^= 1 << (i + j)
    while raw.bit_length() > degree:
        raw ^= modulus << (raw.bit_length() - degree - 1)
    return raw


marking_census = {}
for p, d, modulus in [(3, 2, 7), (7, 3, 11), (31, 5, 37)]:
    powers = [1]
    for i in range(p - 1):
        powers.append(field_mul(powers[-1], 2, d, modulus))
    require(set(powers) == set(range(1, p + 1))
            and field_mul(powers[-1], 2, d, modulus) == 1,
            'direct field powers exhaust all nonzero elements')
    inverse_map = {t: next(u for u in range(1, p) if t * u % p == p - 1)
                   for t in range(1, p)}
    affine_count = 0
    for k in range(1, p):
        E = [powers[k * t % p] for t in range(p)] + [0]
        transport = [None] * (p + 1)
        for t in range(p + 1):
            s = 0 if t == p else p if t == 0 else inverse_map[t]
            transport[E[t]] = E[s]
        affine = all(transport[x ^ y] == (transport[x] ^ transport[y] ^ transport[0])
                     for x in range(p + 1) for y in range(p + 1))
        require(affine == (p in (3, 7)), 'all primitive markings at the three specified fields')
        affine_count += affine
        if p == 31 and k == 1:
            require([transport[x] for x in range(4)] == [1, 0, 18, 14],
                    'independent literal field32 plane image')
            require(transport[0] ^ transport[1] ^ transport[2] ^ transport[3] == 29,
                    'independent field32 XOR defect')
    marking_census[str(p)] = {'markings': p - 1, 'affine_markings': affine_count}

# Identify the author's named designs, not only the unordered orbit pair.
E = [1, 2, 4, 3, 6, 7, 5, 0]
cube = [field_mul(field_mul(x, x, 3, 11), x, 3, 11) for x in range(8)]
named_plus = {b for b in combinations(range(8), 4)
              if E[b[0]] ^ E[b[1]] ^ E[b[2]] ^ E[b[3]] == 0}
named_minus = {b for b in combinations(range(8), 4)
               if cube[E[b[0]]] ^ cube[E[b[1]]] ^ cube[E[b[2]]] ^ cube[E[b[3]]] == 0}
require({frozenset(named_plus), frozenset(named_minus)} == pair,
        'exact named exponential designs equal independent PSL orbit pair')
for g, sign in group.items():
    require({image(b, g) for b in named_plus} == (named_plus if sign == 1 else named_minus),
            'named signs and determinant transport agree')

print(json.dumps({'status': 'PASS', 'checks': checks,
                  'groups': {'PGL': len(group), 'PSL': len(square)},
                  'four_subset_orbits': sorted(map(len, orbits)),
                  'automorphisms': {'single_design': single, 'ordered_pair': ordered,
                                    'unordered_pair': unordered},
                  'design_plus': sorted(plus), 'design_minus': sorted(minus),
                  'one_affine_marking': marking,
                  'endpoint_hostile': {'root': 2, 'UU': compose(U, U)[2], 'A': A[2]},
                  'conic_points': sorted(conic),
                  'primitive_marking_census': marking_census,
                  'named_design_positive': sorted(named_plus),
                  'named_design_negative': sorted(named_minus),
                  'scope': 'Exact F7 finite universe, all 8! permutations; not an all-prime theorem.'},
                 indent=2, sort_keys=True))
