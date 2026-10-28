#!/usr/bin/env python3
"""Exact, bounded audit of the divisor-response future-state theorem.
No external packages; runs identically under python3 and python3 -O.
Universe: all partitions into odd parts >=3 of total mass <=40, empty included.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd
import json

LIMIT = 40

def require(test, message):
    if not test:
        raise RuntimeError(message)

def partitions(n, lo=3):
    if n == 0:
        yield ()
    for p in range(lo, n + 1, 2):
        for tail in partitions(n - p, p):
            yield (p,) + tail

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))

def mobius(n):
    parity = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            parity = -parity
            if n % p == 0:
                return 0
        p += 1
    return -parity if n > 1 else parity

def gcd_exponent(mu):
    return sum((p - 1) // 2 for p in mu) + sum(gcd(p, q) for i, p in enumerate(mu) for q in mu[:i])

def edge_orbit_exponent(mu):
    # Independent definition: construct the permutation, count its unordered-edge orbits.
    n = sum(mu)
    sigma = list(range(n))
    start = 0
    for p in mu:
        for i in range(p):
            sigma[start + i] = start + (i + 1) % p
        start += p
    unseen = {(i, j) for i in range(n) for j in range(i + 1, n)}
    count = 0
    while unseen:
        edge = next(iter(unseen))
        count += 1
        while edge in unseen:
            unseen.remove(edge)
            edge = tuple(sorted((sigma[edge[0]], sigma[edge[1]])))
    return count

def profile(mu):
    return {d: sum(p % d == 0 for p in mu) for d in range(1, LIMIT + 1, 2)}

def response(mu, p):
    return gcd_exponent(mu + (p,)) - gcd_exponent(mu)

rows = []
collision_cells = defaultdict(list)
checked_additions = 0
for mass in range(LIMIT + 1):
    for mu in partitions(mass):
        t = len(mu)
        e = gcd_exponent(mu)
        require(edge_orbit_exponent(mu) == e, f'edge orbit mismatch {mu}')
        M = profile(mu)
        require(2 * e == t * (t - 1) + sum(phi(d) * M[d] ** 2 for d in M if d >= 3), f'quadratic form {mu}')
        A = {1: 0}
        for p in range(3, LIMIT + 1, 2):
            inc = response(mu, p)
            require(inc == (p - 1) // 2 + sum(phi(d) * M[d] for d in divisors(p)), f'add-part {mu,p}')
            A[p] = inc - (p - 1) // 2 - t
            checked_additions += 1
        recovered_M = {1: t}
        for d in range(3, LIMIT + 1, 2):
            numerator = sum(mobius(d // k) * A[k] for k in divisors(d))
            require(numerator % phi(d) == 0, f'nonintegral reconstruction {mu,d}')
            recovered_M[d] = numerator // phi(d)
        require(recovered_M == M, f'divisor inversion {mu}')
        recovered_m = {j: sum(mobius(k) * recovered_M[j * k] for k in range(1, LIMIT // j + 1, 2)) for j in range(3, LIMIT + 1, 2)}
        require(recovered_m == {j: mu.count(j) for j in recovered_m}, f'multiple inversion {mu}')
        collision_cells[(mass,t,e)].append(mu)
        rows.append((mass,t,e,mu))

collisions = [(key, values) for key, values in sorted(collision_cells.items()) if len(values) > 1]
first_key, first_fibre = collisions[0]
require(first_key == (14,2,7), f'first collision changed {first_key}')
require(first_fibre == [(3,11),(5,9)], f'first fibre changed {first_fibre}')
left, right = first_fibre
separating_actions = [p for p in range(3,LIMIT+1,2) if response(left,p) != response(right,p)]
require(separating_actions[0] == 5, 'smallest separator changed')
require(response(left,5) == 4 and response(right,5) == 8, 'hostile responses changed')
require(response(left,3) == response(right,3) == 5, 'positive control changed')
values = [response(left,5),response(right,5)]
mean = sum(Fraction(v,2) for v in values)
variance = sum((Fraction(v)-mean)**2 / 2 for v in values)
require(variance == 4, 'conditional-variance mismatch')
# Every collision pair has a separating allowed one-step continuation.
pairs = 0
for key, fibre in collisions:
    for i, mu in enumerate(fibre):
        for nu in fibre[:i]:
            require(any(response(mu,p) != response(nu,p) for p in range(3,LIMIT+1,2)), f'undetected pair {mu,nu}')
            pairs += 1

# Unit-containing boundary: a different universe, explicitly audited separately.
def odd_partitions_with_units(n, lo=1):
    if n == 0:
        yield ()
    for part in range(lo, n+1, 2):
        for tail in odd_partitions_with_units(n-part, part):
            yield (part,) + tail
unit_cells = defaultdict(list)
unit_count = 0
for mass in range(21):
    for mu in odd_partitions_with_units(mass):
        unit_cells[(mass,len(mu),gcd_exponent(mu))].append(mu)
        unit_count += 1
unit_collisions = [(key,fibre) for key,fibre in sorted(unit_cells.items()) if len(fibre)>1]
unit_key, unit_fibre = unit_collisions[0]
require(unit_key == (8,2,4), 'unit minimality boundary changed')
require(unit_fibre == [(1,7),(3,5)], 'unit hostile changed')
require([response(mu,3) for mu in unit_fibre] == [3,5], 'unit response changed')
require([response(mu,1) for mu in unit_fibre] == [2,2], 'unit query positive control changed')
single_query_pair = [(1,1,1),(3,)]
require([sum(mu) for mu in single_query_pair] == [3,3], 'single-query masses')
require([response(mu,3)-1 for mu in single_query_pair] == [3,3], 'single-query corrected responses')
response_only_pair = [(1,)*7,(3,3,5)]
require([[response(mu,p)-(p-1)//2 for p in (3,5)] for mu in response_only_pair] == [[7,7],[7,7]], 'finite-no-length hostile changed')

result = {
    'status':'FINITE-EXACT PASS',
    'universe':f'all odd-part >=3 cores of mass <={LIMIT}, including empty',
    'core_count':len(rows),
    'independent_edge_orbit_checks':len(rows),
    'add_part_checks':checked_additions,
    'divisor_and_multiplicity_inversions':len(rows),
    'coarse_collision_cells':len(collisions),
    'coarse_collision_pairs_separated':pairs,
    'minimal_same_mass_length_exponent_collision':{
        'mass':14,'length':2,'exponent':7,'cores':[left,right],
        'smallest_allowed_separating_part':5,
        'add_5_increments':values,'resulting_exponents':[11,15],
        'add_3_positive_control_increments':[5,5],
        'equal_fibre_conditional_variance':str(variance),
        'fixed_tournament_counts_after_add_5':[2**11,2**15],
    },
    'unit_containing_boundary':{
        'universe':'all cores with odd parts >=1 of mass <=20',
        'core_count':unit_count,
        'minimal_same_mass_length_exponent_collision':{'mass':8,'length':2,'exponent':4,'cores':unit_fibre,'add_3_increments':[3,5],'add_1_increments':[2,2]},
        'single_query_same_mass_response_only_collision':{'allowed_queries':[3],'cores':single_query_pair,'mass':3,'corrected_responses':[3,3],'lengths':[3,1],'present_exponents':[3,1]},
        'finite_query_response_only_collision':{'allowed_queries':[3,5],'cores':response_only_pair,'corrected_responses':[[7,7],[7,7]],'lengths':[7,3],'present_exponents':[gcd_exponent(mu) for mu in response_only_pair]},
        'scope':'The no-length example concerns responses alone; its current energies differ. With length retained, or an extra prime query beyond support, unit multiplicity is recoverable.'
    },
    'scope':'The all-cutoff theorem is proved in the companion note; this computation is bounded verification. No impossibility of aggregate/analytic recurrence is claimed.'
}
print(json.dumps(result,indent=2))
