#!/usr/bin/env python3
"""Exact finite controls; standard library, no inherited producer imported."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import comb, prod
from pathlib import Path

checks = 0


def require(ok, name):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(name)


def factor(n):
    out = []
    p = 2
    while p*p <= n:
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        if a:
            out.append((p, a))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def aliquot(n):
    return sum(d for d in range(1, n) if n % d == 0)


def divisor_multicomplex(n):
    F = factor(n)
    top = tuple(a for _, a in F)
    D = set(product(*(range(a+1) for _, a in F))) - {top}
    return F, D


def polynomial(adj):
    """Independence polynomial via deletion/link, not subset enumeration."""
    @lru_cache(None)
    def rec(S):
        if not S:
            return (1,)
        b = S & -S
        v = b.bit_length()-1
        a = rec(S ^ b)
        c = rec((S ^ b) & ~adj[v])
        out = [0]*max(len(a), len(c)+1)
        for j, x in enumerate(a):
            out[j] += x
        for j, x in enumerate(c):
            out[j+1] += x
        return tuple(out)
    return rec((1 << len(adj))-1)


def evalpoly(a, x=2):
    return sum(c*x**j for j, c in enumerate(a))


def tournament(n, mask):
    rows = [0]*n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if mask >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def path_count(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v, c in enumerate(dp[S]):
            if c:
                A = rows[v] & ~S
                while A:
                    b = A & -A
                    dp[S | b][b.bit_length()-1] += c
                    A ^= b
    return sum(dp[-1])


def is_auto(rows, p):
    return all(bool(rows[i] >> j & 1) == bool(rows[p[i]] >> p[j] & 1)
               for i, j in combinations(range(len(rows)), 2))


def direct_counts(rows, feedback=False):
    n = len(rows)
    H = A = 0
    fa = n*(n-1)//2
    for p in permutations(range(n)):
        H += all(rows[p[i]] >> p[i+1] & 1 for i in range(n-1))
        A += is_auto(rows, p)
        if feedback:
            fa = min(fa, sum(rows[p[j]] >> p[i] & 1
                             for i in range(n) for j in range(i+1, n)))
    return H, A, fa


def odd_cycles(rows):
    """Literal directed cycles, modulo rotation only; repeated supports retained."""
    out = []
    for k in range(3, len(rows)+1, 2):
        for S in combinations(range(len(rows)), k):
            for tail in permutations(S[1:]):
                c = (S[0],) + tail
                if all(rows[c[j]] >> c[(j+1) % k] & 1 for j in range(k)):
                    out.append(c)
    return out


def ocf(rows):
    C = odd_cycles(rows)
    supports = [sum(1 << v for v in c) for c in C]
    adj = [0]*len(C)
    for i, j in combinations(range(len(C)), 2):
        if supports[i] & supports[j]:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return polynomial(adj), dict(sorted(Counter(map(len, C)).items())), supports


def order_join(parts):
    n = sum(map(len, parts))
    rows = []
    start = 0
    for B in parts:
        end = start + len(B)
        later = sum(1 << j for j in range(end, n))
        rows.extend((r << start) | later for r in B)
        start = end
    return tuple(rows)


# Reconstruct the candidate graph from pairwise faces on actual degree-one
# monomials, then compare its entire independence complex with the divisor set.
flag_counts = Counter()
for n in range(2, 1001):
    F, D = divisor_multicomplex(n)
    r = len(F)
    units = [tuple(int(j == i) for j in range(r)) for i in range(r)]
    active = [i for i in range(r) if units[i] in D]
    edges = {(i, j) for i, j in combinations(active, 2)
             if tuple(int(t in (i, j)) for t in range(r)) not in D}
    recovered = set()
    for bits in product((0, 1), repeat=len(active)):
        S = {active[j] for j, bit in enumerate(bits) if bit}
        if not any(i in S and j in S for i, j in edges):
            recovered.add(tuple(int(i in S) for i in range(r)))
    flag = recovered == D
    require(flag == (sum(a for _, a in F) <= 2), 'flag iff total multiplicity <=2')
    weighted = sum(prod(p**e for (p, _), e in zip(F, E)) for E in D)
    require(weighted == aliquot(n), 'weighted complex equals direct divisor sum')
    require(len(D) == prod(a+1 for _, a in F)-1, 'F1 shadow tau minus1')
    require(weighted == prod(sum(p**e for e in range(a+1)) for p, a in F)-n,
            'product of chains minus top')
    flag_counts['flag' if flag else 'not_flag'] += 1
require(aliquot(1) == 0, 'void n1 boundary')
for n, expected in [(8, 7), (18, 21), (64, 63)]:
    require(aliquot(n) == expected, 'specific touched values')
require(aliquot(64) == 3*(1+4+16), 'even and odd exponent grouping')
for p in (2, 3, 5, 7):
    for r in range(1, 6):
        require(sum(p**j for j in range(3*r)) == sum(p**j for j in range(r))*(1+p**r+p**(2*r)),
                'field-size grouping retains inner exponent sidecar')
F8, D8 = divisor_multicomplex(8)
require(sorted(sum(E) for E in D8) == [0, 1, 2], '1 X X2 graded divisor profile')
F18, D18 = divisor_multicomplex(18)
require(sum(2**sum(E) for E in D18) == 13, 'uniform2 shadow differs from prime weights')
require(aliquot(8) != evalpoly((1, 2, 1)), 'two independent cycles give9 not7')
require(comb(1, 2) < 1, 'one degree1 face cannot support degree2 face')
for value, expected in [(7, [(3, 0)]), (21, [(4, 3), (6, 2), (8, 1), (10, 0)])]:
    pairs = [(a, b) for a in range(value) for b in range(value)
             if 1+2*a+4*b == value and b <= comb(a, 2)]
    require(pairs == expected, 'small independence coefficient constraints')
require(3**3 > 21, 'three independent vertices already contribute27')
K3 = (6, 5, 3)
P4 = (2, 5, 10, 4)
require(polynomial(K3) == (1, 3) and evalpoly(polynomial(K3)) == 7, 'graph7 positive')
require(polynomial(P4) == (1, 4, 3) and evalpoly(polynomial(P4)) == 21, 'graph21 positive')

# Full small labelled tournament universe, with independent vertex-word oracle.
small = {}
for n in range(1, 6):
    Hs, tcs, profiles = set(), set(), set()
    for mask in range(1 << (n*(n-1)//2)):
        rows = tournament(n, mask)
        H, A, _ = direct_counts(rows)
        alpha, _, _ = ocf(rows)
        require(path_count(rows) == H, 'DP versus literal Hamiltonian words')
        require(evalpoly(alpha) == H, 'literal directed-cycle OCF')
        require(H % 2 == 1 and H % A == 0, 'oddness and free-action divisibility')
        require(H not in (7, 21) and H//A not in (7, 21), 'finite small gap control')
        Hs.add(H); tcs.add(H//A); profiles.add(alpha)
    small[n] = {'labelled_tournaments': 1 << (n*(n-1)//2),
                'H_spectrum': sorted(Hs), 'tc_spectrum': sorted(tcs),
                'odd_cycle_profiles': [list(a) for a in sorted(profiles)]}

C3 = (2, 4, 1)
stack_receipts = []
for k in (1, 2, 3):
    rows = order_join([C3]*k)
    alpha, _, _ = ocf(rows)
    require(alpha == tuple(comb(k, j) for j in range(k+1)), 'disjoint cycle binomial profile')
    require(path_count(rows) == evalpoly(alpha) == 3**k, 'SCC path multiplication')
    stack_receipts.append({'triangles': k, 'alpha': alpha, 'H': 3**k})

witness_rows = {
    'H31': (42, 60, 25, 48, 33, 4),
    'H63': (102, 252, 248, 241, 225, 192, 128, 1),
    'H127': (106, 60, 89, 48, 97, 68, 10),
    'Paley7': tuple(sum(1 << j for j in range(7) if (j-i) % 7 in (1, 2, 4)) for i in range(7)),
}
witnesses = {}
for name, rows in witness_rows.items():
    H, A, fa = direct_counts(rows, feedback=(name == 'H127'))
    alpha, lengths, supports = ocf(rows)
    require(path_count(rows) == H == evalpoly(alpha), 'witness DP words OCF agreement')
    require(H == {'H31': 31, 'H63': 63, 'H127': 127, 'Paley7': 189}[name], 'named count')
    require(H % A == 0, 'witness orbit quotient integral')
    if name == 'H63':
        require(alpha == (1, 31) and A == 1, 'H63 complete conflict rigid')
        require(lengths == {3: 8, 5: 17, 7: 6}, 'H63 directed length multiplicities')
    if name == 'H127':
        require(fa == 5, 'fixed H127 class has feedback5 not3')
    if name == 'Paley7':
        require(A == 21 and H//A == 9, 'Paley H Aut tc are different coordinates')
        require(len(set(supports)) < len(supports), 'support collapse loses directed cycles')
    witnesses[name] = {'rows': rows, 'H': H, 'Aut': A, 'tc': H//A,
                       'alpha': alpha, 'cycle_lengths': lengths,
                       'min_feedback': fa if name == 'H127' else None}

# Only the three-point ordinary Collatz cycle descends to multiplier4 here.
def collatz(n):
    return 3*n+1 if n % 2 else n//2

for x in (1, 2, 4):
    require(collatz(x) == (4*x) % 7, 'cycle Collatz equals multiplier4')
paley = witness_rows['Paley7']
for i, j in combinations(range(7), 2):
    require(bool(paley[i] >> j & 1) == bool(paley[(4*i) % 7] >> ((4*j) % 7) & 1),
            'multiplier4 preserves Paley arcs')
require(all(paley[a] >> b & 1 for a, b in [(1, 2), (2, 4), (4, 1)]), 'Paley induced triangle')
require(2 % 7 == 9 % 7 and collatz(2) % 7 != collatz(9) % 7, 'parity prevents global mod7 quotient')

# Mark F8 multiplicative exponents by x^3+x+1; explicitly test additive incidence.
def mul8(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return z

labels = [1]
for _ in range(6):
    labels.append(mul8(labels[-1], 2))
require(labels == [1, 2, 4, 3, 6, 7, 5] and mul8(labels[-1], 2) == 1, 'explicit F8 marking')
for r in range(7):
    require(mul8(labels[r], labels[r]) == labels[2*r % 7], 'doubling exponent is Frobenius')
    require(mul8(mul8(labels[r], labels[r]), mul8(labels[r], labels[r])) == labels[4*r % 7],
            'multiplier4 is Frobenius squared')
require(labels[1] ^ labels[2] ^ labels[4] == 0, 'marked124 is a Fano line')
changed = [labels[3*r % 7] for r in (1, 2, 4)]
require(changed == [3, 5, 7] and prod(changed) != 0 and (changed[0] ^ changed[1] ^ changed[2]) != 0,
        'different primitive generator does not preserve this Fano line')

def rank2(vectors):
    basis = {}
    for x in vectors:
        while x:
            j = x.bit_length()-1
            if j in basis:
                x ^= basis[j]
            else:
                basis[j] = x
                break
    return len(basis)

require(rank2([labels[r] for r in (1, 2, 4)]) == 2, 'Fano line orbit is not a normal basis')
require(rank2(changed) == 3, 'other Frobenius orbit is a normal basis')

out = {'status': 'PASS', 'checks': checks,
       'divisor_universe': {'n_min': 2, 'n_max': 1000, 'flag_counts': dict(flag_counts),
                            'criterion': 'total prime-factor multiplicity <=2; ghost coordinates removed'},
       'small_tournaments': small, 'triangle_stacks': stack_receipts,
       'literal_witnesses': witnesses,
       'graph_controls': {'K3': [1, 3], 'P4': [1, 4, 3]},
       'collatz_cycle': {'cycle': [1, 4, 2], 'modulus': 7, 'multiplier': 4,
                         'no_quotient_hostile': [2, 9]},
       'F8_marking': {'power_table': labels, 'line_image': [labels[r] for r in (1, 2, 4)],
                      'changed_generator_image': changed, 'line_rank': 2, 'normal_orbit_rank': 3},
       'scope': 'Finite independent controls; general flag criterion proved in note; H exclusions inherited; tc all-order exclusion OPEN'}
parser = argparse.ArgumentParser()
parser.add_argument('--output')
args = parser.parse_args()
encoded = json.dumps(out, sort_keys=True, indent=2) + '\n'
if args.output:
    Path(args.output).write_text(encoded)
else:
    print(encoded, end='')
