"""Production reverse-potential controls and operation-count receipts."""
from collections import Counter
from fractions import Fraction
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations
from math import factorial, gcd, prod
from pathlib import Path
import argparse
import hashlib
import json
import sys

from reverse import OrdinaryDAG, baseline, factorial_twos, need, odd_factorials


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
PRIOR = RESULTS / 'port_residual_integral_20261025'


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def partitions(n, least=1):
    if n == 0:
        yield ()
    for first in range(least, n + 1, 2):
        for tail in partitions(n-first, first):
            yield (first,) + tail


def audit_paths(graph):
    n = graph.n
    odd = odd_factorials(n)
    checks = 0
    costs = []
    total = Fraction(0)
    for part in partitions(n):
        mult = Counter(part)
        exponent = sum((d-1)//2 for d in part) + sum(gcd(a,b) for a,b in combinations(part,2))
        t = sum(factorial_twos(k) for k in mult.values())
        endpoint, cost = graph.prefix(mult)
        need(endpoint == graph.terminal and cost == exponent-t, ('partition path valuation',n,part))
        denominator = prod(d**k * factorial(k) for d,k in mult.items())
        total += Fraction(1 << exponent, denominator)
        quotient, remainder = divmod(odd[n], prod(d**k * odd[k] for d,k in mult.items()))
        need(remainder == 0 and quotient % 2 == 1, ('odd path factor',n,part))
        costs.append(cost)
        checks += 1
    need(min(costs) == graph.B, ('minimum path is termwise baseline',n))
    need(total.denominator == 1, ('unmerged integer sum',n))
    return checks, total.numerator


def hostiles():
    merged = OrdinaryDAG(26)
    left, _ = merged.prefix({23:1}, stop_stage=3)
    right, _ = merged.prefix({9:2,5:1}, stop_stage=3)
    need(left == right, 'different past cycle counts truly merge')
    state = merged.states[left]
    need((state.stage,state.remaining,state.defect) == (3,3,(0,2)), 'exact merge sidecars')
    need(merged.structure()['stage_distinctions_lost_if_omitted'] > 0, 'stage is not a redundant key')

    graph = OrdinaryDAG(7)
    skip = next(edge for edge in graph.edges[0] if edge.multiplicity == 0)
    need(skip.cost == 0 and skip.cost-graph.B < 0, 'naive first-edge baseline subtraction fails')
    need(graph.reduced_cost(0,skip) >= 0, 'reverse potential repairs the skip')
    absolute, _ = graph.residue(1, normalized=False)
    normalized, _ = graph.residue(1)
    need(absolute == 0 and normalized == 1, 'absolute low bits do not preserve normalized low bits')

    graph50 = OrdinaryDAG(50)
    low, _ = graph50.residue(2)
    full, _ = graph50.exact(audit_original=True)
    actual_valuation = (full & -full).bit_length()-1
    need(graph50.B == 25 and actual_valuation == 26 and low == 2, 'B is not valuation of the total')
    need(any(graph50.potential[i] < graph50.B for i in range(len(graph50.states))),
         'reachable gauge genuinely includes negative displayed powers')
    return {'exact_merge': {'n':26,'stage':3,'remaining':3,'defect':[0,2],
                             'cycle_histories':[[23],[9,9,5]]},
            'zero_skip': {'n':7,'original_cost':skip.cost,'baseline':graph.B,
                          'reduced_cost':graph.reduced_cost(0,skip)},
            'absolute_truncation': {'n':7,'absolute_mod2':absolute,'normalized_mod2':normalized},
            'summand_vs_sum': {'n':50,'minimum_summand_valuation':graph50.B,
                               'sum_valuation':actual_valuation,'normalized_mod4':low}}


def main(receipt=None):
    old_integer_path = PRIOR / 'counting/check.py'
    precision_path = PRIOR / 'precision/check.py'
    integer = load('old_integer_response_controls',old_integer_path)
    precision = load('bounded_cycle_precision_controls',precision_path)
    path_checks = coordinate_checks = residue_checks = pruning_checks = 0
    samples = []
    for n in range(51):
        graph = OrdinaryDAG(n)
        exact, stats = graph.exact(audit_original=True)
        independent, _ = integer.integral_target(n,core=False)
        need(exact == independent, ('exact integer response agreement',n))
        coordinate_checks += stats['coordinate_checks']
        if n <= 34:
            count, unmerged = audit_paths(graph)
            path_checks += count
            need(exact == unmerged, ('unmerged partition count',n))
        normalized_integer = exact >> graph.B
        for K in range(1,9):
            result, work = graph.residue(K)
            filtered, _ = precision.normalized_residue(n,K)
            need(result == normalized_integer % (1<<K) == filtered, ('three normalized paths',n,K))
            need(work['edge_scans'] == graph.structure()['edges'], 'full structural graph visited before safe pruning')
            residue_checks += 1
            if K in (1,3,8):
                unpruned, uw = graph.residue(K,prune_edges=False,skip_zero_sources=False)
                only_pruned, pw = graph.residue(K,prune_edges=True,skip_zero_sources=False)
                need(result == unpruned == only_pruned, ('edge pruning preserves residue',n,K))
                need(uw['modular_updates'] == graph.structure()['edges'], 'unpruned reference evaluates every edge')
                need(pw['modular_updates']+pw['cost_pruned_edges'] == uw['modular_updates'], 'isolated pruning count')
                pruning_checks += 1

    for n in (50,80,100):
        graph = OrdinaryDAG(n)
        exact, exact_work = graph.exact(audit_original=True)
        independent, _ = integer.integral_target(n,core=False)
        need(exact == independent, ('sample exact engine agreement',n))
        modes = []
        for K in (1,2,6,16):
            normalized, nw = graph.residue(K)
            absolute, aw = graph.residue(graph.B+K,normalized=False)
            need(absolute % (1<<graph.B) == 0 and (absolute>>graph.B) == normalized,
                 ('same target through absolute and normalized gauges',n,K))
            need(normalized == (exact>>graph.B)%(1<<K), ('sample exact residue',n,K))
            if K <= 6:
                filtered, _ = precision.normalized_residue(n,K)
                need(filtered == normalized, ('sample bounded-cycle oracle',n,K))
            modes.append({'normalized_bits':K,'normalized_residue':normalized,
                          'normalized_work':nw,'absolute_work':aw})
        samples.append({'n':n,'structure':graph.structure(),
                        'exact_bit_observations':{k:v for k,v in exact_work.items() if k!='normalized_integer'},
                        'precision_runs':modes})

    for n in (-1,True,1.5,'3'):
        try:
            OrdinaryDAG(n)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid order accepted')
    graph = OrdinaryDAG(4)
    for K in (0,-1,True,1.5):
        try:
            graph.residue(K)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid precision accepted')

    report = {'status':'PASS','exact_and_intermediate_orders':[0,50],
              'normalized_precision_bits':[1,8],
              'normalized_comparisons':residue_checks,
              'pruned_unpruned_comparisons':pruning_checks,
              'intermediate_coordinate_checks':coordinate_checks,
              'unmerged_partition_orders':[0,34],'partition_path_checks':path_checks,
              'hostiles':hostiles(),'samples':samples,
              'source_sha256':{path.relative_to(RESULTS).as_posix():hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for path in (HERE/'reverse.py',HERE/'check.py',old_integer_path,precision_path)},
              'source_hash_semantics':'SHA256 with CRLF normalized to LF; POSIX relative paths',
              'scope':'ordinary structural DAG integration; operation and operand-bit counters, not timings or bit-complexity bounds'}
    text = json.dumps(report,indent=2,sort_keys=True)
    if receipt is not None:
        Path(receipt).write_text(text+'\n',encoding='utf-8')
    print(text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt')
    main(parser.parse_args().receipt)
