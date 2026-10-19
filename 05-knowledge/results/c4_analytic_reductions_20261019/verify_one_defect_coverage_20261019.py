#!/usr/bin/env python3
"""Exact stdlib audit and whole-core coverage for the proved defect/slack test.
Depends only on the bundled verify_saturated_near_twin_lemma_20261019.py.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import gzip
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parent
BASE = ROOT/'verify_saturated_near_twin_lemma_20261019.py'
spec = importlib.util.spec_from_file_location('base', BASE)
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)


def directed_defect(pred, succ, u, v):
    return ((pred[v] & ~pred[u]).bit_count() +
            (succ[u] & ~succ[v]).bit_count())


def audit_injections():
    stats = Counter()
    for n in range(1,7):
        for pred in B.natural_posets(n):
            stats['posets'] += 1
            succ = B.successors(pred)
            pairs = [(u,v) for u in range(n) for v in range(n) if u != v
                     and not ((pred[u] >> v & 1) or (pred[v] >> u & 1))
                     and directed_defect(pred,succ,u,v) <= 1]
            if not pairs:
                continue
            words = B.extensions(pred)
            ranks = {w:{v:i for i,v in enumerate(w)} for w in words}
            for u,v in pairs:
                left = {w for w in words if ranks[w][u] < ranks[w][v]}
                right = set(words) - left
                lower = pred[v] & ~pred[u]
                upper = succ[u] & ~succ[v]
                delta = lower.bit_count() + upper.bit_count()
                if delta == 0:
                    images = {B.relabel(w,{u:v,v:u}) for w in left}
                    B.check(len(images) == len(left) and images <= right, 'Zero-defect injection failed')
                    B.check(len(left) <= len(right), 'Zero-defect bound failed')
                    kind = 'zero_defect'
                else:
                    t = (lower | upper).bit_length() - 1
                    if upper:
                        first = {w for w in left if ranks[w][v] < ranks[w][t]}
                        second = left - first
                        mapping = {u:v,t:u,v:t}
                        kind = 'one_extra_successor'
                    else:
                        first = {w for w in left if ranks[w][t] < ranks[w][u]}
                        second = left - first
                        mapping = {u:t,t:v,v:u}
                        kind = 'one_extra_predecessor'
                    swap = {B.relabel(w,{u:v,v:u}) for w in first}
                    rotate = {B.relabel(w,mapping) for w in second}
                    B.check(len(swap) == len(first) and swap <= right, 'One-defect swap failed')
                    B.check(len(rotate) == len(second) and rotate <= right, 'One-defect rotation failed')
                    B.check(len(left) <= 2*len(right), 'One-defect bound failed')
                    stats['rotation_images'] += len(second)
                stats[kind + '_ordered_pairs'] += 1
                stats['extension_pair_controls'] += len(words)
    return dict(stats)


def build_coverage():
    census = Counter()
    records = []
    catalogue = hashlib.sha256()
    for pred in B.c4_cores():
        succ = B.successors(pred)
        degrees = [8-(pred[v] | succ[v]).bit_count() for v in range(9)]
        if max(degrees) > 7:
            continue
        index = census['degree_valid_cores']
        census['degree_valid_cores'] += 1
        catalogue.update((json.dumps(pred,separators=(',',':'))+'\n').encode())
        epsilon = [0 if v in (0,1) else 7-degrees[v] for v in range(9)]
        candidates = []
        root_candidates = []
        for u,v in combinations(range(9),2):
            if (pred[u] >> v & 1) or (pred[v] >> u & 1):
                continue
            d1 = directed_defect(pred,succ,u,v)
            d2 = directed_defect(pred,succ,v,u)
            if max(d1,d2)+epsilon[u]+epsilon[v] <= 1:
                witness = {'pair':[u,v], 'core_directed_defects':[d1,d2],
                           'exterior_slack_bounds':[epsilon[u],epsilon[v]]}
                candidates.append(witness)
                if u in (0,1) or v in (0,1):
                    root_candidates.append(witness)
        census['root_criterion_cores'] += bool(root_candidates)
        census['any_pair_criterion_cores'] += bool(candidates)
        census['residual_cores'] += not bool(candidates)
        if candidates:
            records.append({'core_index':index, 'predecessor_masks':pred,
                            'core_degrees':degrees, 'witness':candidates[0]})
    B.check(dict(census) == {'degree_valid_cores':79121, 'root_criterion_cores':49317,
                            'any_pair_criterion_cores':49592, 'residual_cores':29529},
            'Defect coverage census mismatch')
    B.check(catalogue.hexdigest() == '4bd17845e2d931c28633c6b7733bc9846fd92817756e0721da198fb9e3cf7148',
            'Catalogue ordering changed')
    return dict(census), records, catalogue.hexdigest()


def main():
    controls = audit_injections()
    census, records, catalogue_hash = build_coverage()
    coverage = {'scope':'Whole-core ambient balanced-pair certificates, for every seven-thin exterior.',
                'criterion':'max(delta_Q(u,v),delta_Q(v,u))+epsilon_u+epsilon_v <= 1',
                'degree_valid_catalogue_sha256':catalogue_hash,
                'covered_cores':len(records), 'certificates':records}
    data = gzip.compress((json.dumps(coverage,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    target = ROOT/'c4_one_defect_analytic_coverage_20261019.json.gz'
    target.write_bytes(data)
    result = {'status':'PASS', 'scope':'Exact finite controls and analytic coverage, not a complete c=4 exclusion.',
              'injection_controls':controls, 'coverage_census':census,
              'coverage_file':target.name, 'coverage_sha256':hashlib.sha256(data).hexdigest(),
              'degree_valid_catalogue_sha256':catalogue_hash,
              'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'dependency_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest()}
    (ROOT/'one_defect_coverage_controls_20261019.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
