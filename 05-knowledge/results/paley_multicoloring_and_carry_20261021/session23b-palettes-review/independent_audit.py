#!/usr/bin/env python3
"""Independent small exact controls. No large lexicographic power enumeration."""
from collections import deque, Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
import json
import hashlib

ROOT = Path(__file__).resolve().parent
CHECKS = 0

def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(label)

PROOF = ROOT.parent / 'session23b-palettes' / 'PALEY_MULTICOLORING_AND_CARRY.md'
PROOF_SHA256 = hashlib.sha256(PROOF.read_bytes()).hexdigest()
EXPECTED_PROOF_SHA256 = 'fe009bd8cdb5a7ab03fb8758a7c969594295ff403e7c6da910378847499afa2d'
need(PROOF_SHA256 == EXPECTED_PROOF_SHA256, 'final audited proof hash')


def tournament(n, signature):
    out = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if signature >> k & 1 else (j, i)
        out[a] |= 1 << b
    return tuple(out)


def paley7():
    return tuple(sum(1 << ((i+d) % 7) for d in (1, 2, 4)) for i in range(7))


def acyclic_masks(out):
    # Independent source-removal test, not the block support criterion.
    n = len(out)
    good = [False] * (1 << n)
    good[0] = True
    for mask in range(1, 1 << n):
        for v in range(n):
            bit = 1 << v
            if mask & bit and out[v] & (mask ^ bit) == mask ^ bit:
                good[mask] = good[mask ^ bit]
                break
    return good


def exact_dichromatic(out):
    good = acyclic_masks(out)
    dp = [0] + [len(out) + 1] * ((1 << len(out)) - 1)
    for mask in range(1, len(dp)):
        pivot = mask & -mask
        sub = mask
        while sub:
            if sub & pivot and good[sub]:
                dp[mask] = min(dp[mask], 1 + dp[mask ^ sub])
            sub = (sub - 1) & mask
    return dp[-1]


def exact_kfold(out, k):
    # Exhaustive BFS of exact vertex-incidence vectors, independently finding
    # the minimum number of legal nonempty colors, with repeats allowed.
    if k == 0:
        return 0
    n = len(out)
    good = acyclic_masks(out)
    supports = [tuple(i for i in range(n) if mask >> i & 1)
                for mask in range(1, 1 << n) if good[mask]]
    start, goal = (0,) * n, (k,) * n
    distance = {start: 0}
    queue = deque([start])
    while queue:
        state = queue.popleft()
        nxt_distance = distance[state] + 1
        for support in supports:
            if any(state[i] == k for i in support):
                continue
            nxt = list(state)
            for i in support:
                nxt[i] += 1
            nxt = tuple(nxt)
            if nxt in distance:
                continue
            if nxt == goal:
                return nxt_distance
            distance[nxt] = nxt_distance
            queue.append(nxt)
    raise AssertionError('singletons make every target reachable')


def lex_product(t, u):
    p, q = len(t), len(u)
    return tuple(sum(1 << (j*q+b)
                     for j in range(p) for b in range(q)
                     if (i != j and t[i] >> j & 1)
                     or (i == j and u[a] >> b & 1))
                 for i in range(p) for a in range(q))


def check_supports(t, u):
    p, q = len(t), len(u)
    actual = acyclic_masks(lex_product(t, u))
    t_good, u_good = acyclic_masks(t), acyclic_masks(u)
    for mask, expected in enumerate(actual):
        parts = [(mask >> (q*i)) & ((1 << q) - 1) for i in range(p)]
        support = sum(1 << i for i, part in enumerate(parts) if part)
        via_blocks = t_good[support] and all(u_good[part] for part in parts)
        need((expected == via_blocks), (t, u, mask))
    return len(actual)


def explicit_palette(k):
    palettes = {
        0: [],
        1: [(0,1,2), (3,4), (5,6)],
        2: [(0,1,2), (0,1,2), (3,4,5), (4,5,6), (3,6)],
        3: [tuple((i+d) % 7 for d in (0,1,2)) for i in range(7)]
    }
    m, r = divmod(k, 3)
    return palettes[3] * m + palettes[r]


p = paley7()
good = acyclic_masks(p)
legal_hist = Counter(mask.bit_count() for mask, value in enumerate(good) if value)
need((legal_hist == {0: 1, 1: 7, 2: 21, 3: 21}), 'legal_hist == {0: 1, 1: 7, 2: 21, 3: 21}')
need((exact_dichromatic(p) == 3), 'exact_dichromatic(p) == 3')
for k in range(101):
    palette = explicit_palette(k)
    incidences = Counter(i for color in palette for i in color)
    need((all(incidences[i] == k for i in range(7))), 'all(incidences[i] == k for i in range(7))')
    need((len(palette) == (7*k + 2)//3), 'len(palette) == (7*k + 2)//3')
    need((all(len(set(color)) == len(color) and good[sum(1 << i for i in color)]
               for color in palette)), 'all(len(set(color)) == len(color) and good[sum(1 << i for i in color)]\n               for color in palette)')
exact_budgets = {k: exact_kfold(p, k) for k in range(5)}
need((exact_budgets == {0: 0, 1: 3, 2: 5, 3: 7, 4: 10}), 'exact_budgets == {0: 0, 1: 3, 2: 5, 3: 7, 4: 10}')

# All labeled tournaments of orders 1, 2, 3, all 121 ordered pairs.
small = [tournament(n, signature) for n in range(1, 4)
         for signature in range(1 << (n*(n-1)//2))]
support_states = 0
for t in small:
    for u in small:
        support_states += check_supports(t, u)
        need((exact_dichromatic(lex_product(t, u)) == exact_kfold(t, exact_dichromatic(u))), 'exact_dichromatic(lex_product(t, u)) == exact_kfold(t, exact_dichromatic(u))')
# Full Paley base, but small blocks, for all 16,384 induced subsets.
support_states += check_supports(p, tournament(2, 0))

# The Paley cyclic out-neighborhood triples are a Fano plane.
lines = [tuple((i+d) % 7 for d in (1,2,4)) for i in range(7)]
need((Counter(tuple(sorted(pair)) for line in lines for pair in combinations(line, 2)) == Counter({pair: 1 for pair in combinations(range(7),2)})), 'Counter(tuple(sorted(pair)) for line in lines for pair in combinations(line, 2)) == Counter({pair: 1 for pair in combinations(range(7),2)})')
need((all(not good[sum(1 << i for i in line)] for line in lines)), 'all(not good[sum(1 << i for i in line)] for line in lines)')
monochrome_hist = Counter()
for bits in range(128):
    count = sum(len({bits >> i & 1 for i in line}) == 1 for line in lines)
    need((count >= 1), 'count >= 1')
    monochrome_hist[count] += 1

# Enumerate all 7! relabelings rather than assume the automorphism group.
automorphisms = []
for perm in permutations(range(7)):
    if all(bool(p[i] >> j & 1) == bool(p[perm[i]] >> perm[j] & 1)
           for i in range(7) for j in range(7)):
        automorphisms.append(perm)
need((len(automorphisms) == 21), 'len(automorphisms) == 21')
need((len([perm for perm in automorphisms if {perm[0],perm[1]} == {0,1}]) == 1), 'len([perm for perm in automorphisms if {perm[0],perm[1]} == {0,1}]) == 1')

# Formal recurrence/denominator checks in exact rational polynomial arithmetic.
import sympy as sp
z = sp.symbols('z')
num, den = 3*z + 2*z**2 + 2*z**3, (1-z)*(1-z**3)
need((sp.gcd(num, den) == 1), 'sp.gcd(num, den) == 1')
need((sp.degree(den) == 4), 'sp.degree(den) == 4')
budget = [sp.Rational((7*k+2)//3) for k in range(80)]
need((sp.expand(sum(budget[k]*z**k for k in range(80))*den).as_poly(z).all_coeffs() is not None), 'sp.expand(sum(budget[k]*z**k for k in range(80))*den).as_poly(z).all_coeffs() is not None')
need((all(budget[k] == budget[k-1] + budget[k-3] - budget[k-4] for k in range(4,80))), 'all(budget[k] == budget[k-1] + budget[k-3] - budget[k-4] for k in range(4,80))')
need((sp.series(num/den, z, 0, 20).removeO() == sum(budget[k]*z**k for k in range(20))), 'sp.series(num/den, z, 0, 20).removeO() == sum(budget[k]*z**k for k in range(20))')

# Rounded dynamics exact controls, not a finite experiment claiming aperiodicity.
rounded_cases = 0
for a in range(3, 21):
    for b in range(2, a):
        if gcd(a,b) != 1:
            continue
        for x0 in range(1,5):
            xs, es = [x0], []
            for n in range(80):
                xs.append((a*xs[-1]+b-1)//b)
                es.append(b*xs[-1]-a*xs[-2])
                need((0 <= es[-1] < b), '0 <= es[-1] < b')
                need((xs[-1] > xs[-2]), 'xs[-1] > xs[-2]')
            q = Fraction(a,b)
            norms = [Fraction(x,q**n) for n,x in enumerate(xs)]
            need((all(norms[n+1] >= norms[n] for n in range(80))), 'all(norms[n+1] >= norms[n] for n in range(80))')
            need((norms[-1] <= x0 + Fraction(b-1,a-b)), 'norms[-1] <= x0 + Fraction(b-1,a-b)')
            need((norms[-1] == x0 + sum(Fraction(es[n],b*q**(n+1)) for n in range(80))), 'norms[-1] == x0 + sum(Fraction(es[n],b*q**(n+1)) for n in range(80))')
            # Test every 1..8-step affine identity on real carry blocks.
            for L in range(1,9):
                for n in range(0,20):
                    c = sum(Fraction(es[n+j], b)*q**(L-1-j) for j in range(L))
                    need((xs[n+L] == q**L*xs[n] + c), 'xs[n+L] == q**L*xs[n] + c')
                    alpha = c/(1-q**L)
                    need((alpha <= 0 < xs[n]), 'alpha <= 0 < xs[n]')
                    need((xs[n+L]-alpha == q**L*(xs[n]-alpha)), 'xs[n+L]-alpha == q**L*(xs[n]-alpha)')
            rounded_cases += 1
xs=[1]
es=[]
for n in range(200):
    xs.append((7*xs[-1]+2)//3)
    es.append(3*xs[-1]-7*xs[-2])
q=Fraction(7,3)
lo = Fraction(xs[-1],q**200)
hi = lo + Fraction(1,2*q**200)
need((1 <= lo <= hi <= Fraction(3,2)), '1 <= lo <= hi <= Fraction(3,2)')

result = {
    'status': 'PASS',
    'checks': CHECKS,
    'proof_sha256': PROOF_SHA256,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'paley_legal_subset_histogram': dict(sorted(legal_hist.items())),
    'exact_bfs_palette_budgets': exact_budgets,
    'explicit_palette_certificates_checked': list(range(101)),
    'small_tournament_substitution_pairs': len(small)**2,
    'induced_subset_support_controls': support_states,
    'fano_monochromatic_line_histogram_over_128_colorings': dict(sorted(monochrome_hist.items())),
    'paley_automorphisms': len(automorphisms),
    'two_color_distinguishing_black_set': [0,1],
    'palette_gf_reduced_denominator_degree': 4,
    'rounded_parameter_initial_value_cases': rounded_cases,
    'rounded_steps_per_case': 80,
    'paley_depth_initial_values': xs[:25],
    'paley_depth_carry_initial_values': es[:50],
    'normalized_limit_decimal_approximation_not_exact_claim': float(lo),
    'normalized_limit_certified_interval': {'lower': str(lo), 'upper': str(hi)},
    'large_lexicographic_power_enumeration': False
}
(ROOT/'audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['explicit_palette_certificates_checked','normalized_limit_certified_interval']},indent=2))
