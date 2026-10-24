#!/usr/bin/env python3
"""Independent bounded controls for the synthesis proof review.

No producer checker is imported. These finite models exercise proof mechanisms;
they are not a native resource authentication or an all-size proof by testing.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from pathlib import Path

CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def multiplicity_controls():
    cases = 0
    wrong_false_positives = 0
    for a in range(1, 6):
        weight = 2 ** a
        # L=weight*[[1,-1],[-1,1]], K_(2)=C_weight; root is vertex 0.
        for b in range(a + 2):
            multiplicity = 2 ** b
            for coefficient in range(-4, 5):
                actual = multiplicity * coefficient
                ordinary = actual % weight == 0
                paired = coefficient * multiplicity % weight == 0
                require(ordinary == paired, 'ordinary and divided conventions disagree')
                # Directly enumerate every harmonic character at each modulus,
                # rooting h at vertex 0; adding a constant has zero action.
                preserves_all = True
                for k in range(1, a + 3):
                    modulus = 2 ** k
                    harmonics = [h for h in range(modulus) if weight*h % modulus == 0]
                    preserved = all(actual*h % modulus == 0 for h in harmonics)
                    preserves_all &= preserved
                require(preserves_all == ordinary, 'full row-image criterion mismatch')
                wrong = actual * multiplicity % weight == 0
                wrong_false_positives += wrong and not ordinary
                cases += 1
    # Smallest illuminating double-scaling hostile: L=4B, Delta=2B, h=(0,1).
    require(4*1 % 4 == 0, 'hostile character not harmonic')
    require(2*1 % 4 == 2, 'hostile actual perturbation did not fail')
    require(2*2 % 4 == 0, 'hostile wrong double-scaled test did not pass')
    require(1*2 % 4 == 2, 'hostile normalized paired test did not fail')
    require(wrong_false_positives > 0, 'double-scaling hostile absent')
    return {'weighted_two_vertex_cases': cases,
            'wrong_double_scaling_false_positives': wrong_false_positives,
            'hostile': {'L_edge_weight': 4, 'actual_Delta_edge_weight': 2,
                        'b': 1, 'h_mod4': [0, 1], 'actual_Delta_h_mod4': [2, 2]}}


def conformal_vectors(d):
    ranges = [range(0, x+1) if x >= 0 else range(x, 1) for x in d]
    return [e for e in product(*ranges) if any(e)]


def conformal(e, d):
    return all((0 <= x <= y) if y >= 0 else (y <= x <= 0) for x, y in zip(e, d))


def test_exact_fiber(states, homogeneous, name):
    atoms = set()
    routes = 0
    for z, w in product(states, repeat=2):
        d = tuple(b-a for a, b in zip(z, w))
        require(homogeneous(d), name + ': difference not in lattice')
        candidates = [e for e in conformal_vectors(d) if homogeneous(e)]
        local_atoms = [e for e in candidates
                       if not any(f != e and conformal(f, e) for f in candidates)]
        remaining = d
        current = z
        while any(remaining):
            choices = [e for e in local_atoms if conformal(e, remaining)]
            require(bool(choices), name + ': no conformal minimal element')
            e = min(choices)
            current = tuple(a+b for a, b in zip(current, e))
            remaining = tuple(a-b for a, b in zip(remaining, e))
            require(current in states, name + ': atom leaves exact fiber')
            require(all(min(a,b) <= c <= max(a,b) for a,b,c in zip(z,w,current)),
                    name + ': atom leaves endpoint box')
            atoms.add(e)
        require(current == w, name + ': route misses endpoint')
        routes += 1
    return {'name': name, 'states': len(states), 'ordered_endpoint_pairs': routes,
            'atoms_used': len(atoms)}


def graver_controls():
    reports = []
    for n in (2, 3):
        for modulus in (2, 3, 4):
            for labels in product(range(modulus), repeat=n):
                groups = defaultdict(list)
                for x in product(range(3), repeat=n):
                    groups[(sum(x), sum(a*b for a,b in zip(x,labels)) % modulus)].append(x)
                def homogeneous(d):
                    return sum(d) == 0 and sum(a*b for a,b in zip(d,labels)) % modulus == 0
                for (margin, charge), states in groups.items():
                    reports.append(test_exact_fiber(states, homogeneous,
                        f'modular_box_n{n}_m{modulus}_{labels}_{margin}_{charge}'))
    # One-hot cycle with capacity slacks, all physical alternative resources retained.
    n = 4
    def encode(bits):
        records = tuple(v for bit in bits for v in (1-bit, bit))
        occupancy = [sum(records[2*i+b] for i in range(n) for b in (0,1)
                         if (i+b) % n == r) for r in range(n)]
        return records + tuple(1-x for x in occupancy)
    states = [encode(bits) for bits in product((0,1), repeat=n)
              if min(encode(bits)[2*n:]) >= 0]
    require(len(states) == 2, 'cycle capacity encoding not exact')
    def capacity_homogeneous(d):
        return (all(d[2*i]+d[2*i+1] == 0 for i in range(n)) and
                all(sum(d[2*i+b] for i in range(n) for b in (0,1)
                        if (i+b) % n == r) + d[2*n+r] == 0 for r in range(n)))
    reports.append(test_exact_fiber(states, capacity_homogeneous, 'one_hot_capacity_cycle'))
    # A complete ternary parity factor with integral one-hot occurrence consistency.
    tuples = [x for x in product((0,1), repeat=3) if sum(x) % 2 == 0]
    states = []
    integral_pairs = 0
    for bits in product((0,1), repeat=3):
        for k, value in enumerate(tuples):
            integral_pairs += 1
            indicators = tuple(int(j == k) for j in range(len(tuples)))
            consistent = all(bits[i] == sum(indicators[j]*t[i] for j,t in enumerate(tuples))
                             for i in range(3))
            require(consistent == (bits == value), 'factor table occurrence mismatch')
            if consistent:
                states.append(tuple(v for bit in bits for v in (1-bit,bit)) + indicators)
    def factor_homogeneous(d):
        return (all(d[2*i]+d[2*i+1] == 0 for i in range(3)) and sum(d[6:]) == 0 and
                all(d[2*i+1] == sum(d[6+j]*t[i] for j,t in enumerate(tuples))
                    for i in range(3)))
    reports.append(test_exact_fiber(states, factor_homogeneous, 'one_hot_full_parity_factor'))
    # Cyclic local consistency does not authenticate a fractional relaxation.
    neq = [(0,1), (1,0)]
    edges = [(0,1),(1,2),(2,0)]
    integral_solutions = 0
    for bits in product((0,1), repeat=3):
        for choices in product(range(2), repeat=3):
            integral_solutions += all((bits[i],bits[j]) == neq[c]
                                      for (i,j),c in zip(edges,choices))
    require(integral_solutions == 0, 'odd-cycle disequality gained integral solution')
    half = Fraction(1,2)
    require(half+half == 1, 'fractional one-hot normalization')
    for edge in edges:
        for occurrence in range(2):
            for value in (0,1):
                marginal = sum(half for row in neq if row[occurrence] == value)
                require(marginal == half, 'fractional local consistency failed')
    return {'fibers': len(reports),
            'ordered_endpoint_pairs': sum(r['ordered_endpoint_pairs'] for r in reports),
            'complete_one_hot_controls': reports[-2:],
            'factor_assignment_pairs_checked': integral_pairs,
            'fractional_odd_cycle_feasible_integral_empty': True}


def gray_resource_controls():
    reports=[]
    for r in range(2,5):
        m=2**r
        gray=[i^(i>>1) for i in range(m)]
        require(len(set(gray)) == m, 'Gray source choices not distinct')
        require(all((gray[i]^gray[(i+1)%m]).bit_count()==1 for i in range(m)),
                'not a cyclic Gray list')
        accepted=[]
        for bits in product((0,1), repeat=m):
            assignment=[gray[(i+bits[i])%m] for i in range(m)]
            if len(set(assignment))==m:
                accepted.append(bits)
        require(accepted == [tuple([0]*m),tuple([1]*m)], 'Gray resource fiber not two-state')
        # Independent replay of the one-extra-block vacancy schedule, at resource level.
        state=[(g,0) for g in gray]
        route=[tuple(state)]
        state[0]=(gray[0],1);route.append(tuple(state))
        for i in range(m-1,0,-1):
            state[i]=(gray[(i+1)%m],0);route.append(tuple(state))
        state[0]=(gray[1],1);route.append(tuple(state))
        state[0]=(gray[1],0);route.append(tuple(state))
        require(len(route)-1 == m+2, 'Gray schedule event count')
        for before,after in zip(route,route[1:]):
            changed=[i for i in range(m) if before[i]!=after[i]]
            require(len(changed)==1,'Gray event changes multiple slots')
            require(len(set(after))==m,'Gray enlarged resource collision')
            i=changed[0];p,l=before[i];q,k=after[i]
            require((p^q).bit_count()+(l^k).bit_count()==1,'Gray nonlocal cube event')
            require(q in (gray[i],gray[(i+1)%m]),'Gray event left eligible menu')
        require(route[-1]==tuple((gray[(i+1)%m],0) for i in range(m)), 'Gray endpoint mismatch')
        reports.append({'r':r,'slots':m,'complete_binary_assignments':2**m,
                        'accepted_states':2,'enlarged_events':len(route)-1})
    return reports


def matroid_controls():
    nontrivial_fibers=0
    tests=0
    for modulus in range(2,9):
        labels=list(range(modulus))
        for k in range(1,modulus):
            fibers=defaultdict(list)
            for A in combinations(labels,k):
                fibers[sum(A)%modulus].append(frozenset(A))
            for states in fibers.values():
                if len(states)<2:
                    continue
                nontrivial_fibers+=1
                for A,B in permutations(states,2):
                    for e in A-B:
                        require(not any((A-{e})|{f} in states for f in B-A),
                                'distinct-label charge fiber admitted basis exchange')
                        tests+=1
    return {'distinct_label_nontrivial_cyclic_fibers':nontrivial_fibers,
            'failed_basis_exchange_instances':tests,
            'scope':'Fixed-cardinality fixed-charge families on these same labelled slots only'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    args=parser.parse_args()
    report={'multiplicity':multiplicity_controls(), 'graver':graver_controls(),
            'gray_resource_skeleton':gray_resource_controls(), 'matroid':matroid_controls()}
    report.update(status='PASS',checks=CHECKS,
                  scope='Independent bounded abstract proof controls; no native packet checker imported or executed')
    Path(args.out).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','checks':CHECKS,'output':args.out}))

if __name__=='__main__':
    main()
