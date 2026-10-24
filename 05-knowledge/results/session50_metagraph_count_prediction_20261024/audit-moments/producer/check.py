#!/usr/bin/env python3
"""Exact finite-moment predictions. Standard library, no external files, no asserts."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)

import predictor as core
rational, ceil, floor, mul = core.rational, core.ceil, core.floor, core.mul
solve, odd_polynomial = core.solve, core.odd_polynomial
congruence, moment_bounds, predict = core.congruence, core.moment_bounds, core.predict

def narayana(n):
    c = [1]
    for j in range(1, n+1):
        c.append(c[j-1] + (c[j-3] if j >= 3 else 0))
    return c[n]

def source_histograms():
    for pin in json.loads((ROOT/'sources/SOURCE_PINS.json').read_text()):
        data = (ROOT/'sources'/pin['snapshot']).read_bytes()
        require(hashlib.sha256(data).hexdigest() == pin['sha256'], 'historical source byte pin')
        require(len(data) == pin['bytes'], 'historical source byte count')
    small = json.loads((ROOT/'sources/small_fibre_audit.json.txt').read_text())
    histograms = {}
    for record in small['sizes']:
        n = record['n']
        hist = Counter()
        for row in record['node_records']:
            mu = row['stored_tilings']
            copies = 1 if mu % 2 else 2
            require(copies == (1 if row['self_converse'] else 2), 'small parity/converse mismatch')
            require(mu % copies == 0 and (mu//copies) % 2 == 1, 'small positive odd mass')
            require(mu == row['predicted_tilings'], 'small direct-inverse audit mismatch')
            hist[mu//copies] += copies
        histograms[n] = hist
    rows = list(csv.DictReader((ROOT/'sources/merged_nodes.tsv.txt').open(), delimiter='\t'))
    hist = Counter()
    for row in rows:
        mu = int(row['fibre'])
        copies = 1 if mu % 2 else 2
        require(copies == (1 if row['class_a'] == row['class_b'] else 2), 'n8 parity/converse mismatch')
        require(mu % copies == 0 and (mu//copies) % 2 == 1, 'n8 positive odd mass')
        hist[mu//copies] += copies
    histograms[8] = hist
    return histograms

def moment_input(histograms):
    # No number of classes, maximum mass, support, or target value is emitted.
    return {'schema': 'odd_fibre_word_moments_v1', 'orders': [
        {'n': n, 'N': 1 << ((n-1)*(n-2)//2),
         'K': [sum(c*h**(r+1) for h, c in hist.items()) for r in range(12)]}
        for n, hist in sorted(histograms.items())]}

def eval_poly(coefficients, x):
    y = F(0)
    for t in reversed(coefficients):
        y = y*x + t
    return y

def check_polynomial_bounds(prediction, hist):
    a = 3 if prediction['peeled_h1_count'] else 1
    low = list(map(F, prediction['lower_q_coefficients']))
    up = list(map(F, prediction['upper_q_coefficients']))
    # Values outside the observed mass range are deliberate hostile controls.
    points = sorted(set(hist) | {a, a+1, a+2, 1000, 10001})
    for x in points:
        if x < a:
            continue
        if low:
            q = eval_poly(low, x)
            p = 2*q-x*q*q
            require(F(1,x)-p == (1-x*q)**2/x, 'lower polynomial identity')
            require(p <= F(1,x), 'lower bound point')
        if up:
            q = eval_poly(up, x)
            p = (1-2*(x-a)*q+x*(x-a)*q*q)/a
            require(p-F(1,x) == (x-a)*(1-x*q)**2/(a*x), 'upper polynomial identity')
            require(p >= F(1,x), 'upper bound point')

def exact_control_moments():
    # Two actual odd-fibre partitions, not a claim of native tournament realizability.
    left, right = {1: 15, 5: 3}, {3: 10}
    def k(hist, j): return sum(c*h**(j+1) for h, c in hist.items())
    require([k(left,j) for j in range(2)] == [k(right,j) for j in range(2)], 'hostile equal moments')
    require(sum(left.values())-sum(right.values()) == 8, 'hostile count gap')
    require(k(left,2) != k(right,2), 'hostile next moment must differ')
    # Total mass32 is possible by adjoining two singleton classes to both.
    left[1] += 2
    right[1] = 2
    require(k(left,0) == k(right,0) == 32, 'hostile power-of-two universe')
    return {'left_counts_by_odd_mass': left, 'right_counts_by_odd_mass': right,
            'K0_K1': [k(left,0),k(left,1)],
            'class_counts': [sum(left.values()),sum(right.values())],
            'count_difference': 8, 'scope': 'abstract odd-fibre partitions only'}

def hand_certificate(moments):
    q = [F(465000,50000000), F(-1250,50000000), F(1,50000000)]
    square = mul(q, q)
    coefficients = [F(1)] + [F(0)]*6
    for i, value in enumerate(q):
        coefficients[i] += 2*value
        coefficients[i+1] -= 2*value
    for i, value in enumerate(square):
        coefficients[i+1] -= value
        coefficients[i+2] += value
    upper = sum(a*b for a,b in zip(coefficients,moments))
    require(upper == F(496641366973194539,15625000000000), 'simple hand upper bound')
    require(upper < 32768, 'simple upper below modulus')
    congr = congruence(moments)
    require(congr['modulus'] == 32768 and 0 < congr['residue'] <= upper, 'hand congruence isolation')
    for x in (1,2,3,5,100,661,1000,10001):
        v = eval_poly(q,x)
        require(eval_poly(coefficients,x)-F(1,x) == (x-1)*(1-x*v)**2/x, 'hand global square identity')
    return {'q_coefficients':list(map(rational,q)),
            'upper_polynomial_coefficients':list(map(rational,coefficients)),
            'exact_upper':rational(upper), 'modular_moments':[x % 32768 for x in moments],
            'congruence':congr, 'prediction':congr['residue']}

def run():
    histograms = source_histograms()
    expected_input = moment_input(histograms)
    supplied = json.loads((ROOT/'moment_inputs.json').read_text())
    require(supplied == expected_input, 'frozen moment input does not match authenticated source')
    reports = []
    for row in supplied['orders']:
        n, moments = row['n'], row['K']
        hist = histograms[n]
        require(moments[0] == row['N'], 'fixed-path cube cardinality')
        require(hist[1] == narayana(n), 'independent rigid h1 recurrence control')
        variants = []
        for peeled in (0, narayana(n)):
            trials = [predict(moments[:k], peeled) for k in range(1, 13)]
            # Validation occurs AFTER blind prediction. Target is never an input.
            actual = sum(hist.values())
            for result in trials:
                require(F(result['lower']) <= actual <= F(result['upper']), 'held-out count not bracketed')
                require(actual % result['congruence']['modulus'] == result['congruence']['residue'], 'held-out congruence')
                if result['prediction'] is not None:
                    require(result['prediction'] == actual, 'blind prediction failed held-out validation')
                check_polynomial_bounds(result, hist)
            first = next(t['moment_count'] for t in trials if t['prediction'] is not None)
            variants.append({'peeled_h1_count': peeled, 'first_success_in_this_hierarchy': first, 'trials': trials})
        reports.append({'n': n, 'N': row['N'], 'held_out_class_count': actual,
                        'rigid_h1_count': narayana(n), 'variants': variants})
    # Algebraic fixed divisor, exact valuation, and Newton half-step denominator controls.
    for k in range(1, 25):
        q = odd_polynomial(k)
        e = 2*k-k.bit_count()
        modulus = 1 << e
        for t in range(-5, k+10):
            h = 2*t+1
            value = eval_poly(q, h)
            require(value == (2**k)*math.prod(range(t-k+1,t+1)), 'falling-factorial identity')
            require(value % modulus == 0, 'fixed divisor control')
        endpoint = int(eval_poly(q, 2*k+1))
        require(endpoint % (2*modulus) == modulus, 'fixed divisor exact at next odd point')
        half_binomial = F((-1)**k * math.comb(2*k,k), 4**k)
        require(half_binomial.denominator == modulus, 'Newton half-binomial exact denominator')
    modular_controls = []
    for row in supplied['orders']:
        m = (row['n']-1)*(row['n']-2)//2
        k = next(k for k in range(1, m+2) if 2*k-k.bit_count() >= m)
        modulus = 1 << (2*k-k.bit_count())
        answer = core.predict_modular(row['N'], [x % modulus for x in row['K'][:k]])
        require(answer['prediction'] == sum(histograms[row['n']].values()), 'universal modular prediction')
        modular_controls.append({'n': row['n'], 'result': answer})
    residue_input = json.loads((ROOT/'twelve_residues.json').read_text())
    n8 = next(row for row in supplied['orders'] if row['n'] == 8)
    require(residue_input == [x % (1 << 22) for x in n8['K']], 'twelve-residue input')
    simple_hand = hand_certificate(n8['K'][:8])
    require(simple_hand == json.loads((ROOT/'hand_certificate.json').read_text()), 'frozen hand certificate')
    hostile = exact_control_moments()
    # The theorem does not allow replacing h by the merged mass mu.
    wrong = {2:1}
    bad = congruence([sum(c*h for h,c in wrong.items())])
    require(sum(wrong.values()) % bad['modulus'] != bad['residue'], 'even-mass hostile control must fail')
    general_thresholds = []
    for n in range(2, 25):
        m = (n-1)*(n-2)//2
        k = next(k for k in range(1, m+2) if 2*k-k.bit_count() >= m)
        general_thresholds.append({'n':n, 'free_bits':m, 'moment_count':k, 'modulus_exponent':2*k-k.bit_count()})
    return {'status':'PASS', 'checks':CHECKS+core.CHECKS, 'orders':reports,
            'hostile_equal_low_moments':hostile, 'even_mass_hostile_control':'correctly rejected',
            'universal_modulus_only_thresholds':general_thresholds, 'modular_controls':modular_controls, 'simple_hand_certificate':simple_hand}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='normal.json')
    parser.add_argument('--write-input', action='store_true')
    args = parser.parse_args()
    if args.write_input:
        (ROOT/'moment_inputs.json').write_text(json.dumps(moment_input(source_histograms()),indent=2)+'\n')
    result = run()
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'], 'checks':result['checks'],
       'orders':[{'n':r['n'],'count':r['held_out_class_count'],
        'first_moments':[v['first_success_in_this_hierarchy'] for v in r['variants']]} for r in result['orders']]},sort_keys=True))
