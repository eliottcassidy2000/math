#!/usr/bin/env python3
"""Bounded exact native-word, CRT fusion, and ternary-rank checks.
Standard library only. The 51-term EM prefix is imported OEIS data, not refactored.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod
import hashlib
import json

CHECKS = 0

def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(label)

def primes_to(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** .5) + 1))]

def valuation(n):
    return (n & -n).bit_length() - 1

def step(n, sign):
    x = 3 * n + sign
    a = valuation(x)
    return x >> a, a

def compile_word(word):
    """Seed P=2n. Return exact allowed n residues at every relevant odd prime."""
    wheel = primes_to(max(word))
    prefix = [1]
    for p in word:
        prefix.append(prefix[-1] * p)
    obstruction = []
    for j, r in enumerate(word):
        for i in range(j):
            if word[i] > r and (prefix[j] // prefix[i]) % r == 1:
                obstruction.append([i + 1, j + 1, r])
    allowed = {}
    for r in wheel:
        if r == 2:
            continue
        if r in word:
            j = word.index(r)
            a = -pow(2 * prefix[j], -1, r) % r
            allowed[r] = {a} if not any(z[2] == r for z in obstruction) else set()
        else:
            forbidden = {-pow(2 * prefix[i], -1, r) % r
                         for i, p in enumerate(word) if p > r}
            allowed[r] = set(range(r)) - forbidden
    R = prod(allowed)
    density = Fraction(prod(map(len, allowed.values())), R)
    return allowed, R, density, obstruction

def direct_word(n, word):
    P = 2 * n
    for p in word:
        least = next((r for r in primes_to(p) if (P + 1) % r == 0), None)
        if least != p:
            return False
        P *= p
    return True

def collatz_cylinder(signs, exponents):
    A, B = 0, 0
    for e, a in zip(signs, exponents):
        B = 3 * B + e * (1 << A)
        A += a
    mod = 1 << (A + 1)
    root = ((1 << A) - B) * pow(3 ** len(signs), -1, mod) % mod
    return root, mod, A, B

def crt_pair(a, m, b, n):
    return (a + m * ((b - a) * pow(m, -1, n) % n)) % (m * n)

def exact_rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    rank = 0
    for col in range(len(a[0])):
        piv = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        v = a[rank][col]
        a[rank] = [x / v for x in a[rank]]
        for i in range(rank + 1, len(a)):
            if a[i][col]:
                v = a[i][col]
                a[i] = [x - v * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank

def conditional_kernel(j, h, source_set):
    """Unnormalised exact weighted numerator via source CRT, preserving branch overlaps."""
    q, t = 3 ** j, 3 ** (j - 1)
    d = 2 * t
    source_mod = 3 ** max(h, j - 1)
    units = [u for u in range(q) if u % 3]
    ix = {u: i for i, u in enumerate(units)}
    M = [[0] * d for _ in range(d)]
    for k in range(2, d + 2):
        weight = 1 << (d + 1 - k)
        for sign in (-1, 1):
            dyad = 1 << (k + 1)
            root = ((1 << k) - sign) * pow(3, -1, dyad) % dyad
            for z in range(source_mod):
                n = root + z * dyad
                if n % (3 ** h) not in source_set:
                    continue
                u, a = step(n, 1)
                v, b = step(n, -1)
                require(min(a, b) == 1 and max(a, b) == k, 'native paired valuations')
                M[ix[u % q]][ix[v % q]] += weight
    return M

# Complete finite universe of all distinct words over four odd primes.
word_rows = []
fusion_examples = []
for size in range(1, 5):
    for word in permutations((3, 5, 7, 11), size):
        allowed, R, density, obstruction = compile_word(word)
        direct = []
        for n in range(1, 2 * R, 2):
            compiler = all(n % r in values for r, values in allowed.items())
            native = direct_word(n, word)
            require(compiler == native, 'word compiler agrees with direct integer least-factor gates')
            if native:
                direct.append(n % R)
        require(len(direct) == density * R, 'exact complete-wheel density')
        require(bool(direct) == (not obstruction), 'internal priority obstruction iff unrealizable')
        word_rows.append({'word': word, 'admissible': bool(direct),
                          'density_among_odd_sources': str(density), 'wheel': R,
                          'priority_obstructions': obstruction})
        if not direct:
            continue
        # Every signed length-two valuation word, not a random sampling.
        for signs in product((-1, 1), repeat=2):
            for exponents in product(range(1, 4), repeat=2):
                c, mod, A, B = collatz_cylinder(signs, exponents)
                n = crt_pair(c, mod, direct[0], R)
                if n == 0:
                    n += mod * R
                require(direct_word(n, word), 'CRT fusion retains native Mullin gate')
                x = n
                got = []
                for e in signs:
                    x, a = step(x, e)
                    got.append(a)
                require(tuple(got) == exponents, 'CRT fusion retains exact native Collatz word')
                require((1 << A) * x == 9 * n + B, 'native affine carry endpoint')
                if word in ((7, 5, 3), (3, 5, 7)) and signs == (1, -1) and exponents == (2, 3):
                    fusion_examples.append({'word': word, 'signs': signs, 'exponents': exponents,
                                            'odd_source': n, 'Mullin_seed': 2 * n,
                                            'Collatz_endpoint': x, 'combined_density': str(density / (1 << A))})
require(compile_word((7, 3))[2] == 0, 'hostile word 7 then 3 impossible')
require(compile_word((7, 5, 3))[2] > 0, 'positive control 7 then 5 then 3')

rank_rows = []
for j in range(1, 4):
    base = conditional_kernel(j, 0, {0})
    require(exact_rank(base) == 2 * 3 ** (j - 1), 'unconditional kernel full rank')
    for size in range(1, 4):
        for values in combinations(range(3), size):
            M = conditional_kernel(j, 1, set(values))
            expected = 2 if j == 1 else 2 * size * 3 ** (j - 2)
            rank = exact_rank(M)
            require(rank == expected, 'ternary source subset exact rank')
            rank_rows.append({'level': j, 'source_modulus': 3, 'source_set': values,
                              'rank': rank, 'expected_rank': expected})
    if j >= 2:
        for values in ({3}, {6}, {3, 6}):
            M = conditional_kernel(j, 2, values)
            expected = 2 if j == 2 else 2 * len(values) * 3 ** (j - 3)
            rank = exact_rank(M)
            require(rank == expected, 'squarefree ternary-shell / birth-cylinder exact rank')
            rank_rows.append({'level': j, 'source_modulus': 9, 'source_set': sorted(values),
                              'rank': rank, 'expected_rank': expected})

# Independent composition: impose a real Mullin word on the source, then
# compute paired Collatz outputs at the multiplicatively transported endpoint.
transport_rows = []
for word in ((3, 5, 7), (7, 5, 3), (7, 13), (5, 7)):
    allowed, R, density, obstruction = compile_word(word)
    require(not obstruction, 'transport test word admissible')
    A = prod(word)
    h = int(A % 3 == 0)
    S = allowed.get(3, {0, 1, 2})
    target_h = h + 1
    target_set = {A * s % (3 ** target_h) for s in S}
    for j in (2, 3):
        q, t = 3 ** j, 3 ** (j - 1)
        d = 2 * t
        units = [u for u in range(q) if u % 3]
        ix = {u: i for i, u in enumerate(units)}
        M = [[0] * d for _ in range(d)]
        source_mod = R * 3 ** max(0, j - 2 - h)
        for k in range(2, d + 2):
            weight = 1 << (d + 1 - k)
            dyad = 1 << (k + 1)
            for sign in (-1, 1):
                root_endpoint = ((1 << k) - sign) * pow(3, -1, dyad) % dyad
                root = root_endpoint * pow(A, -1, dyad) % dyad
                for z in range(source_mod):
                    n = root + z * dyad
                    if not all(n % r in values for r, values in allowed.items()):
                        continue
                    require(direct_word(n, word), 'transport retains direct native Mullin gate')
                    u, a = step(A * n, 1)
                    v, b = step(A * n, -1)
                    require(min(a, b) == 1 and max(a, b) == k,
                            'transport retains endpoint native valuations')
                    M[ix[u % q]][ix[v % q]] += weight
        expected = conditional_kernel(j, target_h, target_set)
        total, expected_total = sum(map(sum, M)), sum(map(sum, expected))
        require(total > 0, 'transport matrix nonempty')
        require(all(M[i][z] * expected_total == expected[i][z] * total
                    for i in range(d) for z in range(d)),
                'native Mullin endpoint law equals predicted ternary-cylinder kernel')
        transport_rows.append({'word': word, 'level': j,
                               'endpoint_source_modulus': 3 ** target_h,
                               'endpoint_source_set': sorted(target_set),
                               'endpoint_kernel_rank': exact_rank(M)})

# Imported primary data; checks below are conditional on the published least-factor labels.
EM = [2,3,7,43,13,53,5,6221671,38709183810571,139,2801,11,17,5471,
52662739,23003,30693651606209,37,1741,1313797957,887,71,7127,109,23,97,
159227,643679794963466223081509857,103,1079990819,9539,3143065813,29,
3847,89,19,577,223,139703,457,9649,61,4357,
87991098722552272708281251793312351581099392851768893748012603709343,
107,127,3313,
227432689108589532754984915075774848386671439568260420754414940780761245893,
59,31,211]
P, seen, records, shell = 1, set(), [], Counter()
frontier_q, frontier_start = 2, 0
for index, p in enumerate(EM, 1):
    require((P + 1) % p == 0, 'imported EM term divides current product plus one')
    require(p not in seen, 'imported EM terms distinct')
    seen.add(p)
    P *= p
    if p == frontier_q:
        wait = index - frontier_start
        records.append({'frontier_prime': p, 'activation_index': frontier_start,
                        'selected_index': index, 'waiting_steps': wait,
                        'conjectured_quadratic_bound': p * p})
        require(wait <= p * p, 'finite canonical frontier quadratic-deadline probe')
        frontier_q = next(q for q in primes_to(1000) if q not in seen)
        frontier_start = index
    if index >= 2:
        n = P // 2
        require(n % 2 == 1 and n % 3 == 0 and n % 9 != 0, 'native squarefree 3-shell')
        shell[n % 72] += 1
require(frontier_q == 41 and frontier_start == 50, 'known prefix frontier')
# Canonical-orbit rank ceiling is structural, without any mixing assumption.
empirical_rank_rows = []
for j in range(1, 5):
    q = 3 ** j
    units = [u for u in range(q) if u % 3]
    ix = {u: i for i, u in enumerate(units)}
    M = [[0] * len(units) for _ in units]
    source_mod = 3 ** (j - 1)
    depth = min(2, j - 1)
    seed_classes = {3 % (3 ** depth), 6 % (3 ** depth)}
    source_classes = [s for s in range(source_mod) if s % (3 ** depth) in seed_classes]
    row_strip = {(3 * s + 1) * pow(2, -1, q) % q for s in source_classes}
    column_strip = {(3 * s - 1) * pow(2, -1, q) % q for s in source_classes}
    P = 1
    for index, p in enumerate(EM, 1):
        P *= p
        if index < 2:
            continue
        u, a = step(P // 2, 1)
        v, b = step(P // 2, -1)
        require(u % q in row_strip or v % q in column_strip,
                'actual canonical sample obeys unconditional strip support')
        M[ix[u % q]][ix[v % q]] += 1
    ceiling = 2 if j <= 2 else 4 * 3 ** (j - 3)
    rank = exact_rank(M)
    require(len(row_strip) + len(column_strip) == ceiling,
            'canonical shell strip size gives rank ceiling')
    require(rank <= ceiling, 'canonical empirical matrix structural rank ceiling')
    require(rank == ceiling, 'published 51-term prefix saturates observed rank ceiling')
    empirical_rank_rows.append({'level': j, 'modulus': q, 'sample_indices': [2, 51],
                                'empirical_rank': rank, 'unconditional_rank_ceiling': ceiling,
                                'allowed_low_valuation_rows': sorted(row_strip),
                                'allowed_low_valuation_columns': sorted(column_strip)})

# Pin every veto to a distinct smaller selected prime.
vetoes = []
for q in primes_to(100):
    P = 1
    blockers = []
    for index, p in enumerate(EM, 1):
        if P % q == 0:
            break
        if (P + 1) % q == 0 and p != q:
            require(p < q and p not in blockers, 'finite-veto charge')
            blockers.append(p)
        P *= p
    if blockers:
        vetoes.append({'target_prime': q, 'distinct_blockers': blockers})
result = {
 'status': 'FINITE-EXACT compiler, CRT-fusion, and rank checks; sourced-prefix probes are not refactorizations',
 'checks': CHECKS,
 'universe': {'word_alphabet': [3,5,7,11], 'word_lengths': [1,2,3,4],
              'signed_Collatz_prefix_length': 2, 'valuation_alphabet': [1,2,3],
              'ternary_levels': [1,2,3]},
 'word_count': len(word_rows), 'admissible_word_count': sum(x['admissible'] for x in word_rows),
 'words': word_rows, 'fusion_examples': fusion_examples, 'rank_checks': rank_rows,
 'native_endpoint_transport_checks': transport_rows,
 'known_EM_prefix_source': 'https://oeis.org/A000945/b000945.txt',
 'known_EM_prefix_length': len(EM), 'canonical_frontier_records': records,
 'open_frontier': {'prime': frontier_q, 'activation_index': frontier_start,
                   'observed_unselected_steps': len(EM) - frontier_start,
                   'conjectured_latest_selection_index': frontier_start + frontier_q ** 2},
 'finite_veto_examples': vetoes,
 'canonical_empirical_rank_checks': empirical_rank_rows,
 'shell_mixing_probe': {'modulus': 72, 'sample_indices': [2,51],
                       'counts': dict(sorted(shell.items())),
                       'uniform_comparison_count_per_shell_residue': '25/4',
                       'warning': '50 terms, descriptive only; no convergence evidence or test significance claimed'},
 'scope': 'Finite cylinders over varying seeds, not canonical-orbit coverage or fixed-sign Collatz convergence.'}
result['payload_sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
print(json.dumps(result, indent=2, sort_keys=True))
