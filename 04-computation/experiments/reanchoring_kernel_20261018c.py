#!/usr/bin/env python3
"""Exact binomial reanchoring, first-block carry transport, and S4 hostile.

Reproduce with Python 3 (standard library only):
  python -B -X utf8 04-computation/experiments/reanchoring_kernel_20261018c.py --out 05-knowledge/results/reanchoring_kernel_20261018c.json
Checks use explicit exceptions and remain active under python -O.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def integer(value, name):
    if type(value) is not int:
        raise TypeError(name + ' must be an integer, not a coercible scalar')


def sequence(a, b, c, last):
    for value, name in zip((a, b, c, last), ('a', 'b', 'c', 'last')):
        integer(value, name)
    if last < 0:
        raise ValueError('last must be nonnegative')
    u = [Fraction(1)]
    for n in range(last):
        u.append(((a*n*n+a*n+b)*u[n] - (c*n*n*u[n-1] if n else 0))/(n+1)**2)
    return u


def integral_sequence(a, b, c, last):
    values = sequence(a, b, c, last)
    require(all(x.denominator == 1 for x in values), ('nonintegral', a, b, c))
    return [x.numerator for x in values]


def shifted(u, s, epsilon=1):
    integer(s, 'shift')
    integer(epsilon, 'epsilon')
    if epsilon not in (-1, 1):
        raise ValueError('the integer carry API supports only orientation +1/-1')
    return [epsilon**n * sum(comb(n, k)*s**(n-k)*u[k] for k in range(n+1))
            for n in range(len(u))]


def reanchor(a, b, c, d, epsilon):
    for value, name in zip((a, b, c, d, epsilon), ('a', 'b', 'c', 'd', 'epsilon')):
        integer(value, name)
    if epsilon not in (-1, 1):
        raise ValueError('integer API supports only orientation +1/-1')
    if d*(d*d-a*d+c):
        raise ValueError('root guard P(d)=0 is necessary for the operator family')
    return ((a-3*d)*epsilon, (b-d)*epsilon, 3*d*d-2*a*d+c)


def primes_to(bound):
    return [n for n in range(2, bound+1) if all(n % d for d in range(2, int(n**0.5)+1))]


def carry_correction(u, b, s, p, r):
    """Exact integer C_(p,r), without reducing its numerator before /p."""
    low = sum(comb(p+r, k)*s**(p+r-k)*u[k] for k in range(p))
    high_base = b*sum(comb(p+r, p+j)*s**(r-j)*u[j] for j in range(r+1))
    new_base = (b+s)*sum(comb(r, j)*s**(r-j)*u[j] for j in range(r+1))
    numerator = low+high_base-new_base
    require(numerator % p == 0, ('nonintegral correction', p, r, s))
    return numerator//p


def carry(u, b, p):
    require(len(u) >= 2*p, ('missing full first-two blocks', p))
    result = []
    for r in range(p):
        numerator = u[p+r]-b*u[r]
        require(numerator % p == 0, ('first-block Lucas failure', b, p, r))
        result.append(numerator//p)
    return result


def polynomial_lucas(u, v, s, epsilon, p):
    # Coefficients of (1-s*x/eps)^(p-1) F_u((x/eps)/(1-s*x/eps)).
    for n in range(p):
        rhs = epsilon**n * sum(u[k]*comb(p-1-k, n-k)*(-s)**(n-k)
                               for k in range(n+1))
        require((v[n]-rhs) % p == 0, ('polynomial identity', p, n))


def s4_hostile():
    edges = list(itertools.combinations(range(4), 2))
    matchings = [frozenset((frozenset(edges[i]), frozenset(edges[5-i]))) for i in range(3)]
    require(len(set(matchings)) == 3, 'three perfect matchings')
    matchings = list(dict.fromkeys(matchings))
    def action(perm):
        e = tuple(edges.index(tuple(sorted((perm[a], perm[b])))) for a, b in edges)
        m = tuple(matchings.index(frozenset(frozenset(perm[a] for a in edge) for edge in mat))
                  for mat in matchings)
        return e, m
    permutations = list(itertools.permutations(range(4)))
    actions = [action(p) for p in permutations]
    edge_kernel = [p for p, (e, m) in zip(permutations, actions) if e == tuple(range(6))]
    matching_kernel = [p for p, (e, m) in zip(permutations, actions) if m == tuple(range(3))]
    require(len(edge_kernel) == 1 and len(matching_kernel) == 4, 'different action kernels')
    require(len(set(e for e, m in actions)) == 24 and len(set(m for e, m in actions)) == 6,
            'action image sizes')
    witness = (1, 0, 3, 2)
    e, m = action(witness)
    require(m == tuple(range(3)) and e != tuple(range(6)), 'double-transposition hostile')
    return {'ordered_frames': 24, 'edge_action_kernel_size': 1,
            'crossratio_action_kernel_size': 4, 'witness_permutation': witness,
            'witness_edge_action': e, 'witness_matching_action': m,
            'consequence': 'No S4-equivariant bijection between the six K4 edges and a generic six-value crossratio orbit.'}


def variable_pair_labels():
    edges = list(itertools.combinations(range(4), 2))
    basis = ((0, 1, 2), (0, 1, 3), (0, 2, 3))
    labels = (1, 0, 0, 0)
    require(all(len({labels[v] for v in cycle}) > 1 for cycle in basis), 'all basis cycles mixed')
    require(len({labels[v] for v in (1, 2, 3)}) == 1, 'dependent triangle monochromatic')
    basis_xor = set()
    for cycle in basis:
        for i in range(3):
            edge = tuple(sorted((cycle[i], cycle[(i+1) % 3])))
            basis_xor.symmetric_difference_update((edge,))
    require(basis_xor == {(1, 2), (1, 3), (2, 3)}, 'cycle-space relation')
    require(all((3*s+4) % 6 != 0 for s in range(2)), 'isolated pin needs its prescribed ternary label')
    feasible = 0
    for label in itertools.product(range(3), repeat=4):
        # In K4 the fibre subgraphs are cliques, hence bipartite iff size <=2.
        fibres_bipartite = all(label.count(t) <= 2 for t in range(3))
        actual = any(all((3*shade[a]+4*label[a]) % 6 != (3*shade[b]+4*label[b]) % 6
                         for a, b in edges)
                     for shade in itertools.product(range(2), repeat=4))
        require(actual == fibres_bipartite, ('variable pair label equivalence', label))
        feasible += actual
    return {'graph': 'abstract K4; no plane embedding claimed', 'labels': labels,
            'basis': basis, 'missed_odd_cycle': (1, 2, 3),
            'labelings_tested': 81, 'labelings_with_six_colour_lift': feasible,
            'pin_guard_hostile': {'vertices': 1, 'edges': [], 'pin_colour': 0, 'proposed_label': 1},
            'consequence': 'A cycle basis does not suffice to constrain all monochromatic odd cycles.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    primes = [p for p in primes_to(97) if p >= 7]
    last = 2*max(primes)-1
    A = integral_sequence(7, 2, -8, last)
    # Independent definition of Franel numbers, not a second recurrence evaluator.
    require(A == [sum(comb(n, k)**3 for k in range(n+1)) for n in range(last+1)],
            'Franel binomial-cube oracle')
    transforms = [(-1, 1, (10, 3, 9), 'C'), (8, -1, (17, 6, 72), 'F')]
    records = []
    all_first = []
    for d, epsilon, target, name in transforms:
        require(reanchor(7, 2, -8, d, epsilon) == target, ('root parameters', name))
        s = -d
        v = shifted(A, s, epsilon)
        require(v == integral_sequence(*target, last), ('exact transformed recurrence', name))
        # t=d+eps*z => inverse t'=-d/eps + t/eps; orientation remains eps.
        require(shifted(v, d*epsilon, epsilon) == A, ('inverse map', name))
        failures, p_records = [], []
        for p in primes:
            old = carry(A, 2, p)
            new = carry(v, target[1], p)
            polynomial_lucas(A, v, s, epsilon, p)
            correction_residues = []
            for r in range(p):
                correction = carry_correction(A, 2, s, p, r)
                naive = sum(comb(r, j)*s**(r-j)*old[j] for j in range(r+1))
                predicted = epsilon**(r+1)*(naive+correction)
                require((new[r]-predicted) % p == 0, ('corrected carry', name, p, r))
                correction_residues.append(correction % p)
                if (new[r]-epsilon**(r+1)*naive) % p:
                    failures.append({'p': p, 'r': r, 'actual': new[r] % p,
                                     'naive': epsilon**(r+1)*naive % p,
                                     'correction_before_orientation': correction % p})
            p_records.append({'p': p, 'nonzero_corrections': sum(x != 0 for x in correction_residues),
                              'correction_sha256': hashlib.sha256(bytes(correction_residues)).hexdigest()})
        require(failures, ('a hostile to correction deletion is required', name))
        records.append({'target': name, 'parameters': target, 'root': d,
                        'orientation': epsilon, 'initial_terms': v[:10],
                        'naive_transport_failures': len(failures),
                        'first_naive_failure': failures[0], 'primes': p_records})
        all_first.append(failures[0])
    presentations = []
    for d in (0, -1, 8):
        for epsilon in (-1, 1):
            target = reanchor(7, 2, -8, d, epsilon)
            require(shifted(A[:41], -d, epsilon) == integral_sequence(*target, 40),
                    ('six oriented presentations', d, epsilon))
            presentations.append({'d': d, 'epsilon': epsilon, 'parameters': target})
    try:
        reanchor(7, 2, -8, 1, 1)
    except ValueError:
        pass
    else:
        raise AssertionError('nonroot must be rejected')
    require(1*(1-7-8) == -14, 'nonroot produces nonzero constant in cubic')
    for bad in (True, 1.0):
        try:
            shifted([1, 2], bad)
        except TypeError:
            pass
        else:
            raise AssertionError('scalar coercion erased declared type')
    require(sequence(17, 5, 72, 2)[2] == Fraction(123, 4), 'literal F variant')
    require(sequence(11, 5, 125, 2)[2] == Fraction(5, 2), 'literal D signature variant')
    group = s4_hostile()
    pair_labels = variable_pair_labels()
    payload = {'status': 'PROVED scoped identities / FINITE-EXACT bounded checks',
               'universe': {'sequence_indices': [0, last], 'primes': primes,
                            'carry_indices': 'every 0 <= r < p', 'source': [7, 2, -8]},
               'presentations': presentations, 'transforms': records, 'S4_hostile': group,
               'variable_pair_labels': pair_labels,
               'checks': CHECKS,
               'limitations': ['No Hadwiger-Nelson or Collatz reduction.',
                               'Classical binomial transform is attributed to Zagier and Liu.',
                               'Carry theorem assumes first-block Lucas; bounded checks are not that theorem.',
                               'Integer rescaling API is restricted to epsilon=+1 or -1.']}
    args.out.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'checks': CHECKS, 'out': str(args.out), 'first_hostiles': all_first}, indent=2))


if __name__ == '__main__':
    main()
