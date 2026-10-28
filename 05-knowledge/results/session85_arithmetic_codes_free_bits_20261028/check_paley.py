#!/usr/bin/env python3
"""Exhaustive labelled Fano/Paley7 sidecar checks; Python standard library only.

No assert statements: all checks remain active under python -O.
Universes: all 35 triples, all STS(7) by exact cover and relabelling,
all 128 cyclic orientations of a fixed STS, all 5040 vertex permutations.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
import json
from functools import reduce
from math import gcd, isqrt, prod


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


V = tuple(range(1, 8))
PAIRS = tuple(combinations(V, 2))
TRIPLES = tuple(combinations(V, 3))
PERMS = tuple(permutations(V))
F = frozenset(t for t in TRIPLES if t[0] ^ t[1] ^ t[2] == 0)


def image_plane(G, p):
    return frozenset(tuple(sorted(p[x - 1] for x in t)) for t in G)


def image_arcs(A, p):
    return frozenset((p[x - 1], p[y - 1]) for x, y in A)


def cyclic(A):
    return frozenset((a, b, c) for a, b, c in TRIPLES
                     if ((a, b) in A and (b, c) in A and (c, a) in A)
                     or ((b, a) in A and (c, b) in A and (a, c) in A))


def reverse(A):
    return frozenset((b, a) for a, b in A)


def is_sts(G):
    return len(G) == 7 and Counter(p for t in G for p in combinations(t, 2)) == Counter(PAIRS)


def sts_exact_cover():
    pairsets = {t: frozenset(combinations(t, 2)) for t in TRIPLES}
    by_pair = {p: [t for t in TRIPLES if p in pairsets[t]] for p in PAIRS}
    answer = set()

    def visit(remaining, selected):
        if not remaining:
            answer.add(frozenset(selected))
            return
        p = min(remaining)
        for t in by_pair[p]:
            if pairsets[t] <= remaining:
                visit(remaining - pairsets[t], selected + [t])

    visit(frozenset(PAIRS), [])
    return answer


def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 0b1011  # x^3+x+1
    return result


def field_power(a, n):
    out = 1
    for _ in range(n):
        out = mul(out, a)
    return out


def h_dp(A):
    dp = {(1 << (x - 1), x): 1 for x in V}
    for mask in range(1, 128):
        for y in V:
            if mask & (1 << (y - 1)):
                smaller = mask ^ (1 << (y - 1))
                if smaller:
                    dp[mask, y] = sum(dp.get((smaller, x), 0) for x in V if (x, y) in A)
    return sum(dp[127, y] for y in V)


def h_brute(A):
    return sum(all((p[i], p[i + 1]) in A for i in range(6)) for p in PERMS)


def dot(a, x):
    return (a & x).bit_count() % 2


def fano_orientations(G):
    lines = sorted(G)
    for mask in range(128):
        arcs = set()
        for j, (a, b, c) in enumerate(lines):
            tri = ((a, b), (b, c), (c, a))
            arcs.update((y, x) if mask & (1 << j) else (x, y) for x, y in tri)
        yield frozenset(arcs)


powers = tuple(field_power(2, i) for i in range(7))
need(powers == (1, 2, 4, 3, 6, 7, 5), 'field power chart')
log = {x: i for i, x in enumerate(powers)}
T = frozenset((x, y) for x in V for y in V if (log[y] - log[x]) % 7 in (1, 2, 4))
CYC = cyclic(T)
inv = tuple(field_power(x, 6) for x in V)
frob = tuple(mul(x, x) for x in V)
G0 = image_plane(F, inv)
need(is_sts(F) and is_sts(G0) and not F & G0 and CYC == F | G0, 'reciprocal Fano split')
need(image_arcs(T, inv) == reverse(T), 'inversion reverses tournament')
need(image_plane(F, frob) == F and image_arcs(T, frob) == T, 'Frobenius incidence/orientation')
need(all(sum((x, y) in T for y in V) == 3 for x in V), 'regularity')
need(Counter(p for t in CYC for p in combinations(t, 2)) == Counter({p: 2 for p in PAIRS}), '2-(7,3,2)')

planes = sts_exact_cover()
need(planes == {image_plane(F, p) for p in PERMS}, 'exact-cover/relabel STS agreement')
need(len(planes) == 30 and all(is_sts(G) for G in planes), '30 STS')
intersection_hist = Counter(len(F & G) for G in planes)
need(intersection_hist == {0: 8, 1: 14, 3: 7, 7: 1}, 'STS intersection distribution')
paley = {image_arcs(T, p) for p in PERMS}
need(len(paley) == 240, '240 labelled Paley tournaments')
cycles = {A: cyclic(A) for A in paley}
need(all(sum(G <= C for G in planes) == 2 for C in cycles.values()), 'exactly two Fano planes in every Paley')
partners = {G for G in planes if not F & G}
compatible = {A for A, C in cycles.items() if F <= C}
fibres = defaultdict(set)
for A in compatible:
    fibres[cyclic(A) - F].add(A)
need(len(compatible) == 16 and set(fibres) == partners, '16 orientations/eight partners')
need(all(len(As) == 2 and reverse(next(iter(As))) in As for As in fibres.values()), 'partner fibres are reversals')

# Complete independent-line universe, not only permutations of a chosen Paley.
orientations = set(fano_orientations(F))
need(len(orientations) == 128, '128 independent Fano-line orientations')
for G in partners:
    need({A for A in orientations if cyclic(A) == F | G} == fibres[G] and len(fibres[G]) == 2,
         'entire cyclic-triple fibre has exactly two tournaments in all128 universe')
need({A for A in orientations if cyclic(A) - F in partners} == compatible, 'no non-Paley completion of disjoint union')
hist = Counter()
for A in orientations:
    hd = h_dp(A)
    need(hd == h_brute(A), 'Hamilton DP/permutation disagreement')
    hist[hd] += 1
    need((hd == 189) == (A in compatible), 'H=189 iff Paley in bounded family')
need(hist == {189: 16, 171: 112}, '128-orientation Hamilton census')

# All 120 unordered disjoint pairs, and every labelled cyclic orientation.
union_fibres = defaultdict(set)
for A, C in cycles.items():
    union_fibres[C].add(A)
disjoint_pairs = [(G, H) for G, H in combinations(planes, 2) if not G & H]
need(len(disjoint_pairs) == 120 and len(union_fibres) == 120, '120 pairs/unions')
for G, H in disjoint_pairs:
    As = {A for A in fano_orientations(G) if H <= cyclic(A)}
    need(As == union_fibres[G | H] and len(As) == 2, 'unique reversed completion for every disjoint pair')

# Point-code and line-code have different coordinate types; a dot polarity is explicit.
normal_lines = {a: tuple(x for x in V if dot(a, x) == 0) for a in V}
need(set(normal_lines.values()) == F, 'line-normal identification')
dual_frob = tuple(next(b for b in V if set(normal_lines[b]) ==
                      {frob[x - 1] for x in normal_lines[a]}) for a in V)
need(frob == (1, 4, 5, 6, 7, 2, 3) and dual_frob == (1, 6, 7, 2, 3, 4, 5)
     and frob != dual_frob, 'Frobenius is contragredient on line normals, not same-label equivariance')

line_code = {m for m in range(128) if reduce(
    int.__xor__, (a for a in V if m & (1 << (a - 1))), 0) == 0}
need(Counter(m.bit_count() for m in line_code) == {0: 1, 3: 7, 4: 7, 7: 1}, 'Hamming weights')
encoded = {}
for v in range(8):
    for s in range(2):
        m = sum((s ^ dot(a, v)) << (a - 1) for a in V)
        A = set(T)
        for a, line in normal_lines.items():
            if m & (1 << (a - 1)):
                for x, y in combinations(line, 2):
                    if (x, y) in A:
                        A.remove((x, y)); A.add((y, x))
                    else:
                        A.remove((y, x)); A.add((x, y))
        encoded[s, v] = (m, frozenset(A))
need({m for m, A in encoded.values()} == line_code, 'affine evaluation parametrization')
need({A for m, A in encoded.values()} == compatible, 'exact Hamming orientation locus')
need(all(encoded[0, v][0] ^ encoded[1, v][0] == 127 and
         reverse(encoded[0, v][1]) == encoded[1, v][1] and
         cyclic(encoded[0, v][1]) == cyclic(encoded[1, v][1]) for v in range(8)), 'one reversal bit')
need(len({cyclic(encoded[0, v][1]) - F for v in range(8)}) == 8, 'three partner bits')

GL = [p for p in PERMS if image_plane(F, p) == F]
autT = [p for p in PERMS if image_arcs(T, p) == T]
need(len(GL) == 168 and len(autT) == 21, 'automorphism sizes')
need(all(p in GL for p in autT), 'Paley automorphisms preserve chosen Fano')
orbit = {image_arcs(T, p) for p in GL}
need(len(orbit) == 8 and reverse(T) not in orbit and compatible == orbit | {reverse(A) for A in orbit}, 'two GL orbits')
need(orbit == {A for (s, v), (m, A) in encoded.items() if s == 0}, 'even-codeword GL orbit/chirality')
need(len({image_plane(G0, p) for p in GL}) == 8, 'transitive partner action/no invariant choice')
need(len([p for p in GL if image_plane(G0, p) == G0]) == 21, 'partner stabilizer')
need(h_dp(T) == 189 and 189 // 21 == 9 and 189 % 21 == 0, 'H/Aut quotient')

# Optional arithmetic transport is a permutation of a declared prime bank only.
prime = {1: 7, 2: 2, 3: 5, 4: 19, 5: 89, 6: 59, 7: 29}
source = {1, 2, 3}
target = {frob[x - 1] for x in source}
n, image_n = prod(prime[x] for x in source), prod(prime[x] for x in target)
sigma_n, sigma_image = prod(prime[x] + 1 for x in source), prod(prime[x] + 1 for x in target)
need((n, image_n, sigma_n, sigma_image) == (70, 11837, 144, 14400), 'inherited Frobenius arithmetic witness')


# Exact arithmetic-to-tournament map. Point columns a become line normals a
# only after the declared dot-product polarity; T is the affine origin.
all_transport = {}
for m in range(128):
    selected = [a for a in V if m & (1 << (a - 1))]
    N = prod(prime[a] for a in selected)
    divisors = [prod(prime[selected[j]] for j in range(len(selected)) if k & (1 << j))
                for k in range(1 << len(selected))]
    sigma_direct = sum(divisors)
    need(len(set(divisors)) == 1 << len(selected) and all(N % d == 0 for d in divisors), 'squarefree divisor list')
    need(sigma_direct == prod(prime[a] + 1 for a in selected), 'direct divisor-sum/multiplicative agreement')
    flips = frozenset(p for a in selected for p in combinations(normal_lines[a], 2))
    A = frozenset((y, x) if tuple(sorted((x, y))) in flips else (x, y) for x, y in T)
    square = isqrt(sigma_direct) ** 2 == sigma_direct
    need(square == (m in line_code) == (A in compatible), 'all128 square-sigma iff Paley bridge')
    all_transport[m] = (N, sigma_direct, A, flips)
need({row[2] for row in all_transport.values()} == orientations, 'arithmetic/line orientation bijection')
for m in range(128):
    for k in range(128):
        N, _, A, fm = all_transport[m]
        M, _, B, fk = all_transport[k]
        NM, _, C, fc = all_transport[m ^ k]
        need(N * M // gcd(N, M) ** 2 == NM and fm ^ fk == fc, 'star/line-flip composition')
minimal = [m for m in line_code if m and not any(k and k != m and k & m == k for k in line_code)]
need(Counter(m.bit_count() for m in minimal) == {3: 7, 4: 7}, '14 minimal zero patterns')
need(127 in line_code and 127 not in minimal and all_transport[127][2] == reverse(T), 'all-seven nonminimal reversal')
# Chirality recovery from observed prime/line-flip coordinates.
recovery_counts = Counter()
for observed in range(128):
    observed_points = {a for a in V if observed & (1 << (a - 1))}
    classes = defaultdict(set)
    for m in line_code:
        classes[m & observed].add(m.bit_count() % 2)
    recoverable = all(len(values) == 1 for values in classes.values())
    contains_line = any(set(line) <= observed_points for line in F)
    need(recoverable == contains_line, 'chirality recoverability iff contains normal-space Fano line')
    if recoverable:
        recovery_counts[observed.bit_count()] += 1
need(recovery_counts == {3: 7, 4: 28, 5: 21, 6: 7, 7: 1}, 'chirality recovery census')

need(all(all(p % d for d in range(2, isqrt(p) + 1)) for p in prime.values()), 'bank primality')

out = {
    'status': 'all explicit checks passed',
    'universes': {'points': 7, 'triples': 35, 'vertex_permutations': 5040,
                  'STS7_exact_cover_and_relabelling': 30, 'independent_Fano_line_orientations': 128,
                  'labelled_Paley_tournaments': 240, 'unordered_disjoint_STS_pairs': 120},
    'field': {'polynomial': 'x^3+x+1', 'alpha_powers': powers, 'Fano_lines': sorted(F),
              'reciprocal_lines': sorted(G0), 'Frobenius_on_labels_1_to_7': frob},
    'fixed_fano': {'intersection_histogram': dict(sorted(intersection_hist.items())),
                   'disjoint_partners': len(partners), 'compatible_Paley_orientations': len(compatible),
                   'orientation_fibre_per_partner': 2, 'GL_orbit_sizes': [8, 8],
                   'hamming_line_code_weights': dict(sorted(Counter(m.bit_count() for m in line_code).items()))},
    'Hamilton_paths': {'Paley': 189, 'aut_group': 21, 'ratio': 9,
                      'independent_line_orientation_histogram': dict(sorted(hist.items()))},
    'sidecar_chart': [{'partner_coordinate_v': v, 'relative_mask_s0': encoded[0, v][0],
                      'relative_mask_s1': encoded[1, v][0],
                      'partner_lines': sorted(cyclic(encoded[0, v][1]) - F)} for v in range(8)],
    'Frobenius_prime_bank_witness': {'n': n, 'image_n': image_n, 'sigma': sigma_n, 'image_sigma': sigma_image},
    'arithmetic_orientation_bridge': {'all_subsets_checked': 128, 'star_pairs_checked': 16384,
        'square_sigma_and_Paley': 16, 'nonsquare_sigma_and_nonPaley': 112,
        'minimal_square_sigma_flip_patterns_by_weight': dict(sorted(Counter(m.bit_count() for m in minimal).items())),
        'line_normal_to_prime': prime, 'GL_chirality_equals_selection_parity_on_Paley_locus': True,
        'chirality_recovering_observation_sets_by_size': dict(sorted(recovery_counts.items())),
        'all_seven_flips': 'global reversal, nonminimal codeword'},
    'scope': 'Labelled order-seven finite structures. Uniform-count information needs chosen charts; no natural-number dynamics or idoneality transfer.'
}
print(json.dumps(out, indent=2, sort_keys=True))
