#!/usr/bin/env python3
"""Exact adjacent-grade hole-difference compiler on authenticated native input."""
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
COUNT = 0


def need(ok, message):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(message)


def plus(a, b):
    return ((((a >> 4)+(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def minus(a, b):
    return ((((a >> 4)-(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    erasure_path = HERE.parent/'erasure/receipt.json'
    old_path = RESULTS/'native_packet_decomposition_20261023/native/domains.json'
    electric_path = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    erasure = json.loads(erasure_path.read_text())
    old = json.loads(old_path.read_text())['sources']
    electric = json.loads(electric_path.read_text())['order8']
    X = electric['higher_layers'][-1]['lifts'][0]
    Y = electric['radical_lifts_mod4']
    need(len(old) == 28 and len(erasure['sources']) == 20, 'actual native source inventory sizes')
    exceptional = erasure['adaptive_repair']['exceptional_original_centers']
    need(len(exceptional) == 6, 'six exact central-deficiency boundaries')
    # Inherited exhaustive census: no unlisted triple has a grade8 or9 hole.
    profiles = erasure['levels'][3]['profiles']
    need(sum(row['pin_sets'] for row in profiles if row['coverage'][8] < 2048 or row['coverage'][9] < 2048) == 6,
         'all1140 audited triple profiles have exactly six central bad cases')
    need(sum(row['pin_sets'] for row in profiles) == 1140, 'complete triple universe retained')
    forbidden = defaultdict(list)
    for row in exceptional:
        holes = {int(k):set(v) for k,v in row['original_holes'].items()}
        need(set(holes) == {7,8,9,10}, 'exact adjacent-grade holes retained')
        for grade in (8,9,10):
            for q in holes[grade]:
                for previous in holes[grade-1]:
                    label = minus(q, previous)
                    need(plus(label, previous) == q, 'literal difference witness')
                    forbidden[label].append({'erased':row['erased'], 'grade':grade,
                                             'hole':q, 'previous_hole':previous})
    need(len(forbidden) == 32, 'exact native forbidden repair-label count')
    # Test every abstract label by two independent formulas on the witnessed holes.
    repaired = []
    for label in range(2048):
        direct_ok = True
        for row in exceptional:
            holes = {int(k):set(v) for k,v in row['original_holes'].items()}
            for grade in (8,9,10):
                residual = holes[grade] & {plus(label,p) for p in holes[grade-1]}
                direct_ok = direct_ok and not residual
        need(direct_ok == (label not in forbidden), 'difference criterion equals direct translated-hole test')
        if direct_ok:
            repaired.append(label)
    candidates = []
    for row in old:
        source = row['source']
        code = row['code']
        label = ((X[source] % 128) << 4) | sum((Y[j][source] % 2) << (j-1) for j in range(1,5))
        need(label == (code[0] << 4) | sum(v << j for j,v in enumerate(code[1:])),
             'original physical potentials authenticate inherited paired label')
        if source in erasure['sources']:
            continue
        records = row['exact_all_shared_clocks_target_domain']
        need(len(records) == 2, 'complete inherited pinned domain is binary')
        paths = [record['path'] for record in records]
        clocks = [i for i,(a,b) in enumerate(zip(*paths)) if a != b]
        need(len(clocks) in (3,4), 'actual old outside-source switching clocks')
        need(label not in forbidden, 'every old outside source is a universal three-pin repair')
        candidates.append({'source':source, 'label':label, 'cost':len(clocks),
                           'changed_clocks':clocks, 'paths':paths})
    need(len(candidates) == 8, 'all eight outside sources inspected')
    need([c['source'] for c in candidates if c['cost'] == 3] == [1841,2087,2835,3015,3360],
         'five actual three-clock repair alternatives')
    # Corollary bounds for all numbers of additional unary outside pins.
    budgets = []
    for cost in (3,4):
        for e in range(8):
            Rcost = [2]*17+[cost]
            W = lambda t:sum(sorted(Rcost, reverse=True)[:t])
            worst = {'entry':0,'adjacent':0,'outside':0,'helper':0}
            for k in range(8,11+e):
                r0 = max(8,k-e)
                for r in range(19):
                    if not 0 <= k-r <= e:
                        continue
                    if r < 8:
                        worst['entry'] = max(worst['entry'],W(r+8)+4*(8-r))
                    elif r > 10:
                        worst['entry'] = max(worst['entry'],W(26-r)+4*(r-10))
                for r in range(r0+1,11):
                    if 0 <= k-r < e:
                        worst['adjacent'] = max(worst['adjacent'],W(min(2*r-1,37-2*r))+4)
                d = W(2*min(r0,18-r0))
                if 0 < k-r0 < e:
                    worst['outside'] = max(worst['outside'],d+8)
                worst['helper'] = max(worst['helper'],d+15)
                need(max(worst.values()) <= 49+cost, 'all interior events fit chosen repair cost')
            need(worst['entry'] <= 46+cost and worst['adjacent'] <= 36+cost and worst['outside'] <= 42+cost,
                 'declared stage bounds include all feasible grades')
            budgets.append({'repair_cost':cost,'free_outside':e,'interior_grade':[8,10+e],
                            'exact_grade_only_stage_maxima':worst,'uniform_clock_bound':49+cost})
    result = {'status':'FINITE-EXACT; inherits audited all1140 native erasure profiles',
              'forbidden_label_count':len(forbidden), 'allowed_abstract_labels':len(repaired),
              'forbidden_labels':sorted(forbidden),
              'forbidden_witnesses':{str(k):v for k,v in sorted(forbidden.items())},
              'native_candidates':candidates,'flexible_unary_pin_budgets':budgets,
              'checks_before_receipt':COUNT,
              'inputs':{str(p.relative_to(RESULTS)):digest(p) for p in [erasure_path,old_path,electric_path]},
              'source_sha256':digest(HERE/'repair_labels.py')}
    result = json.loads(json.dumps(result))
    path = HERE/'repair_labels_receipt.json'
    if args.write:
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    else:
        need(result == json.loads(path.read_text()), 'frozen exact repair-label certificate')
    print(f'PASS {COUNT:,}; 32 forbidden labels; all8 native sources repair all1140 erased triples')


if __name__ == '__main__':
    main()
