"""Independent exact audit of Session64's ordinary reverse valuation pass.

The reference is imported only as the implementation under test.  This file
reconstructs unmerged prefix coordinates from literal odd cycle multisets,
checks complete paths against Burnside summands, and uses tournament permutation
Burnside (rather than the implementation's orbit marking) as a class oracle.
All gates are explicit and survive python -O.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
import json
from math import factorial, gcd, prod
from pathlib import Path
import sys


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def v2(x):
    check(x > 0, 'positive valuation input')
    return (x & -x).bit_length() - 1


def odd(x):
    return x >> v2(x)


def partitions(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n + 1, 2):
        for tail in partitions(n - d, d):
            yield (d,) + tail


def raw_weight(part):
    multiplicities = Counter(part)
    exponent = sum((d - 1) // 2 for d in part)
    exponent += sum(gcd(a, b) for a, b in combinations(part, 2))
    denominator = prod(d ** k * factorial(k) for d, k in multiplicities.items())
    return Fraction(1 << exponent, denominator), exponent - sum(v2(factorial(k)) for k in multiplicities.values())


def prefix_key(n, stage, part):
    remaining = n - sum(part)
    defects = tuple(sum(d - gcd(c, d) for c in part)
                    for d in range(1, min(stage, remaining) + 1, 2))
    return stage, remaining, defects


def prefix_value(n, stage, part):
    key = prefix_key(n, stage, part)
    _, remaining, defects = key
    bound = max(((value + d - 2) // (d - 1)
                 for d, value in zip(range(1, min(stage, remaining) + 1, 2), defects)
                 if d > 1), default=0)
    weight, _ = raw_weight(part)
    result = odd(factorial(sum(part))) * weight * (1 << (remaining * (len(part) - bound)))
    check(result.denominator == 1, ('independent prefix integral', n, key, part))
    return key, int(result)


def independent_prefixes(impl, dag, Y, Z):
    if dag.n == 0:
        check(Y == Z == [1], 'empty tournament boundary')
        return 1, 1, None
    n = dag.n
    top = n if n % 2 else n - 1
    histories = [()]
    lookup = {key: i for i, key in enumerate(dag.vertices)}
    total_histories = total_states = 0
    merge_witness = None
    for stage in range(top, -2, -2):
        expected = defaultdict(int)
        history_by_key = defaultdict(list)
        for part in histories:
            key, value = prefix_value(n, stage, part)
            expected[key] += value
            history_by_key[key].append(part)
        actual_keys = {key for key in dag.vertices if key[0] == stage}
        check(set(expected) == actual_keys, ('independent reachable stage graph', n, stage))
        for key, value in expected.items():
            index = lookup[key]
            check(Y[index] == value, ('unmerged prefix Y', n, key))
            z = Fraction(value) * Fraction(2) ** (dag.potential[index] - dag.baseline)
            check(z.denominator == 1 and z == Z[index], ('independent prefix Z', n, key))
            parts = history_by_key[key]
            if merge_witness is None and stage > 0 and key[1] > 0 and len(parts) > 1:
                merge_witness = {'n':n, 'vertex':list(key), 'prefixes':[list(parts[0]), list(parts[1])]}
        total_states += len(expected)
        total_histories += len(histories)
        if stage == -1:
            break
        following = []
        for part in histories:
            s = n - sum(part)
            key = prefix_key(n, stage, part)
            choices = [s] if stage == 1 else range(s // stage + 1)
            expected_edges = []
            for k in choices:
                extended = part + (stage,) * k
                target_key = prefix_key(n, stage - 2, extended)
                expected_edges.append((k, lookup[target_key]))
                following.append(extended)
            check(sorted(expected_edges) == sorted((e.k, e.target) for e in dag.edges[lookup[key]]),
                  ('unmerged exact outgoing choices', n, key, part))
        histories = following
    return total_states, total_histories, merge_witness


def complete_path_checks(impl, dag, Y, Z):
    n = dag.n
    total = Fraction(0)
    minval = None
    partitions_checked = 0
    K_values = (1, 2, 3, 8)
    support = {K:set() for K in K_values}
    edge_by_count = [{e.k:(j,e) for j,e in enumerate(outgoing)} for outgoing in dag.edges]
    for part in partitions(n):
        weight, valuation = raw_weight(part)
        total += weight
        minval = valuation if minval is None else min(minval, valuation)
        source = 0
        path_edges = []
        sumJ = 0
        productT = 1
        counts = Counter(part)
        while source != dag.terminal:
            stage = dag.vertices[source][0]
            edge_number, edge = edge_by_count[source][counts[stage]]
            path_edges.append((source, edge_number))
            sumJ += edge.shift
            t = odd(factorial(edge.u + edge.ell * edge.k)) // (odd(factorial(edge.u)) * edge.ell ** edge.k * odd(factorial(edge.k)))
            productT *= t
            source = edge.target
        directT = Fraction(odd(factorial(n)), prod(d ** k * odd(factorial(k)) for d,k in counts.items()))
        check(directT.denominator == 1 and productT == directT and productT % 2 == 1,
              ('full odd-factor telescope', n, part))
        check(sumJ == valuation, ('path cost equals raw summand valuation', n, part))
        for K in K_values:
            if valuation < dag.baseline + K:
                support[K].update(path_edges)
        partitions_checked += 1
    check(total.denominator == 1, ('partition integrality', n))
    check(Y[-1] == total * odd(factorial(n)), ('independent terminal Y normalization', n))
    check(Z[-1] * (1 << dag.baseline) == Y[-1], ('independent terminal Z normalization', n))
    check(minval == dag.baseline, ('independent minimum valuation', n))
    support_checks = 0
    for K in K_values:
        filtered = {(v,j) for v,outgoing in enumerate(dag.edges) for j,e in enumerate(outgoing)
                    if dag.prefix_cost[v] + dag.reduced_cost(v,e) < K}
        check(filtered == support[K], ('exact complete-path edge support', n, K))
        support_checks += 1
    return int(total), partitions_checked, support_checks


def permutation_burnside(n):
    """Independent unrestricted census by permutation action on unordered pairs."""
    pairs = list(combinations(range(n), 2))
    pair_index = {pair:i for i,pair in enumerate(pairs)}
    fixed_sum = 0
    for permutation in permutations(range(n)):
        target = []
        flipped = []
        for i,j in pairs:
            a,b = permutation[i], permutation[j]
            target.append(pair_index[tuple(sorted((a,b)))])
            flipped.append(a > b)
        visited = set()
        orbits = 0
        impossible = False
        for start in range(len(pairs)):
            if start in visited:
                continue
            orbits += 1
            current = start
            parity = False
            while current not in visited:
                visited.add(current)
                parity ^= flipped[current]
                current = target[current]
            if parity:
                impossible = True
                break
        fixed_sum += 0 if impossible else 1 << orbits
    count, remainder = divmod(fixed_sum, factorial(n))
    check(remainder == 0, ('permutation Burnside integral', n))
    return count


def run(source):
    spec = spec_from_file_location('implementation_under_audit', source / 'reverse_ordinary.py')
    impl = module_from_spec(spec)
    sys.modules[spec.name] = impl
    spec.loader.exec_module(impl)
    checked_edges = checked_coordinates = modular_passes = 0
    prefix_states = prefix_histories = partition_cases = support_checks = 0
    full_census = []
    merge_witness = None
    absolute_zero = None
    missing_stage = None
    costs = []
    orders = list(range(51)) + [60,80,100]
    for n in orders:
        dag = impl.build_dag(n)
        count, normalized, Y, Z, detail = impl.exact_audit(dag)
        checked_coordinates += len(Y)
        check(len(set(dag.vertices)) == len(dag.vertices), ('unique complete vertex keys', n))
        for v, outgoing in enumerate(dag.edges):
            ell,s,D = dag.vertices[v]
            b = max(((D[i]+2*i-1)//(2*i) for i in range(1,len(D))), default=0)
            for e in outgoing:
                ell2,s2,D2 = dag.vertices[e.target]
                check(e.target > v and ell2 == ell-2 and s2 == s-ell*e.k and e.u == n-s,
                      ('edge structure', n, v, e))
                check(D2 == tuple(D[i]+e.k*(d-gcd(ell,d)) for i,d in enumerate(range(1,min(ell-2,s2)+1,2))),
                      ('edge exact defect merge', n, v, e))
                bp = max(((D2[i]+2*i-1)//(2*i) for i in range(1,len(D2))), default=0)
                sell = D[(ell-1)//2] if ell<=s else 0
                expected_shift = e.k*(ell-1)//2+ell*e.k*(e.k-1)//2-v2(factorial(e.k))+e.k*(ell*b-sell)+s2*(e.k+b-bp)
                check(expected_shift == e.shift and e.shift >= 0, ('integer J', n, v))
                numerator = odd(factorial(n-s2))
                denominator = odd(factorial(n-s))*ell**e.k*odd(factorial(e.k))
                t,r = divmod(numerator, denominator)
                check(r==0 and t>0 and t%2==1, ('exact integer odd T', n, v))
                check(impl.edge_factor(e, impl.odd_factorials(n)) == t, ('T implementation', n, v))
                for K in (1,2,8,17):
                    check(impl.edge_factor(e, impl.odd_factorials(n,1<<K),1<<K) == t % (1<<K), ('modular odd factor', n, v,K))
                checked_edges += 1
            if outgoing:
                check(dag.potential[v] == min(e.shift+dag.potential[e.target] for e in outgoing), ('suffix Bellman equality', n,v))
        check(dag.potential[-1]==0 and dag.prefix_cost[0]==0, ('potential boundaries',n))
        predecessors = [[] for _ in dag.vertices]
        for v,outgoing in enumerate(dag.edges):
            for e in outgoing:
                predecessors[e.target].append(dag.prefix_cost[v]+dag.reduced_cost(v,e))
        for v in range(1,len(dag.vertices)):
            check(dag.prefix_cost[v] == min(predecessors[v]), ('prefix Bellman equality',n,v))
        for K in (1,2,3,8,17):
            residue, work, values = impl.normalized_modular(dag,K)
            full, fullwork, fullvalues = impl.normalized_modular(dag,K,False)
            check(values == fullvalues == [z % (1<<K) for z in Z], ('all modular coordinates', n,K))
            check(residue == full == normalized % (1<<K), ('modular terminal',n,K))
            check(all(0<=x<(1<<K) for x in values), ('K-bit stored residues', n,K))
            edge_count = sum(len(x) for x in dag.edges)
            check(work['edge_updates']+work['valuation_filtered_edges']==edge_count==fullwork['edge_updates'], ('full accounting',n,K))
            independent_updates = sum(dag.prefix_cost[v]+dag.reduced_cost(v,e)<K for v,E in enumerate(dag.edges) for e in E)
            check(work['edge_updates']==independent_updates, ('update count',n,K))
            modular_passes += 1
            if absolute_zero is None:
                for v,(y,z) in enumerate(zip(Y,Z)):
                    if y and y%(1<<K)==0 and z%(1<<K)!=0:
                        absolute_zero = {'n':n,'K':K,'vertex':v,'key':dag.vertices[v], 'Y':y,'Z':z,'baseline':dag.baseline,'suffix_potential':dag.potential[v]}
                        break
            if n in (30,50,60,80,100) and K in (2,8):
                costs.append({'n':n,'K':K,'vertices_allocated':len(dag.vertices),'all_edges_preprocessed':edge_count,'modular_edge_updates':independent_updates})
        if n <= 32:
            states,histories,merge = independent_prefixes(impl,dag,Y,Z)
            prefix_states += states
            prefix_histories += histories
            if merge_witness is None and merge is not None:
                merge_witness = merge
            answer,cases,support = complete_path_checks(impl,dag,Y,Z)
            check(answer == count, ('direct partition final',n))
            partition_cases += cases
            support_checks += support
        if n <= 7:
            census = permutation_burnside(n)
            check(census == count, ('unrestricted pair-action Burnside',n))
            full_census.append({'n':n,'classes':census,'permutations':factorial(n),'labelled_tournaments':1<<(n*(n-1)//2)})
        if missing_stage is None:
            for v,E in enumerate(dag.edges):
                for e in E:
                    if dag.vertices[v][1:] == dag.vertices[e.target][1:] and dag.vertices[v][1]>0:
                        missing_stage = {'n':n,'source':dag.vertices[v], 'target':dag.vertices[e.target], 'k':e.k,'source_index':v,'target_index':e.target,'failure':'Erasing stage identifies this edge endpoint with its source, making a self-loop in the purported decreasing DAG.'}
        if n == 50:
            check(dag.baseline==25 and normalized%4==2 and count%4==0,'order50 precision hostile')
            check(any(e.shift-dag.baseline < 0 for e in dag.edges[0]), 'initial-edge-only subtraction hostile')
    check(merge_witness is not None and absolute_zero is not None and missing_stage is not None,'actual hostiles present')
    for n in (-1,True,False,1.5,None,'3'):
        try:
            impl.build_dag(n)
        except ValueError:
            pass
        else:
            raise RuntimeError(('bad order accepted',n))
    for K in (-1,0,True,False,1.5,None,'3'):
        try:
            impl.normalized_modular(dag,K)
        except ValueError:
            pass
        else:
            raise RuntimeError(('bad precision accepted',K))
    normal = (source/'reverse-normal.json').read_bytes()
    optimized = (source/'reverse-optimized.json').read_bytes()
    check(normal==optimized,'producer supplied normal/-O receipt equality')
    return {'status':'PASS', 'implementation_sha256':hashlib.sha256((source/'reverse_ordinary.py').read_bytes()).hexdigest(),
            'orders':{'consecutive':[0,50],'additional':[60,80,100]},
            'all_edge_T_J_checks':checked_edges,'all_exact_prefix_coordinates':checked_coordinates,
            'modular_comparisons':modular_passes,'precision_bits':[1,2,3,8,17],
            'independent_prefix_orders':[0,32],'independent_unmerged_prefix_histories':prefix_histories,
            'independent_merged_prefix_coordinates':prefix_states,'independent_complete_partitions':partition_cases,
            'exact_budget_edge_support_comparisons':support_checks,
            'unrestricted_pair_action_burnside':full_census,'state_merge_witness':merge_witness,
            'missing_stage_hostile':missing_stage,'absolute_zero_hostile':absolute_zero,
            'costs':costs,'producer_receipts_identical':True,
            'scope':'Ordinary e=0, all unrestricted tournament classes; no core/excess, port, connectivity, Newton, runtime, peak-memory, or smaller structural recurrence claim.',
            'precision_qualification':'K-bit stored residues and modular arithmetic; Python temporary products/shifts may exceed K bits before reduction.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'session64-counting')
    args = parser.parse_args()
    print(json.dumps(run(args.source),indent=2,sort_keys=True))
