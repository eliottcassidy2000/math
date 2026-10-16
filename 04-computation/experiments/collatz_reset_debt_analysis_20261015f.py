#!/usr/bin/env python3
"""Reanalyse frozen native rows: discounted debt and pairs of actual resets."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'05-knowledge/results'


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def ceilroot(n, d):
    a = math.isqrt(n//d)
    return a if a*a*d == n else a+1


def initial_norms(k):
    X, b = 1 << k, 1 << (k//2)
    # L0 is affine plus a period-three error; sum each bin in O(1).
    error3 = [F(4*n-3*((2*n+3)//6)+(4*n+1)//6-(2*n+1)//6)-F(10*n,3)
              for n in range(3)]
    out = []
    for a in range(X, 2*X, b):
        total = F(5*b*(2*a+b-1),3)
        total += (b//3)*sum(error3)+sum(error3[(a+j) % 3] for j in range(b % 3))
        need(total.denominator == 1, 'integral halving bin norm')
        out.append(int(total))
    return out


def analyse(rec):
    k, b = rec['k'], rec['bin_width']
    X, rho = 1 << k, F(7,10)
    initial = initial_norms(k)
    need(sum(v*v for v in initial) == rec['initial_envelope'][0], 'same initial source norms')
    history = [initial]+[v['row_norms'] for v in rec['stages']]
    padding = k+8-(len(history)-1)
    qfinal = sum(v*v for v in history[-1])
    need(qfinal < 49*b*X*X, 'identity padding has actual residual below 7X')
    tail = ceilroot(qfinal, b)
    debt = F(0)
    total = F(ceilroot(rec['initial_envelope'][0], b))
    failures, profiles = [], []
    for j, (old, new) in enumerate(zip(history, history[1:]), 1):
        r2 = sum(max(0, 10*v-7*u)**2 for u, v in zip(old, new))
        need([r2, 9*b] == rec['stages'][j-1]['residual_square'], 'frozen survivor residual')
        rupper = ceilroot(r2, 9*b)
        debt = rho*debt+(1-rho)*rupper
        total = rho*total+(1-rho)*rupper
        if r2 > 441*b*X*X:
            failures.append(dict(stage=j, residual_square_numerator=r2,
                                 denominator=9*b, seven_squared_threshold=441*b*X*X))
        # Exact cumulative source-weighted mass of weak-cancellation rows.
        profile = {}
        for delta in [F(1,10), F(1,5), F(1,4), F(3,10)]:
            mass = sum(u*u for u, v in zip(old,new) if F(u-v) <= delta*u)
            profile[str(delta)] = [mass, b*X*X]
        profiles.append(dict(stage=j, weighted_lower_tail_over_X2=profile))
    for _ in range(padding):
        debt = rho*debt+(1-rho)*tail
        total = rho*total+(1-rho)*tail
    paired, previous = [], 0
    for end in range(2, len(history)+1, 2):
        end = min(end, len(history)-1)
        old, new = history[previous], history[end]
        r2 = sum(max(0, 10*v-7*u)**2 for u,v in zip(old,new))
        paired.append(dict(start=previous, end=end, residual_square=[r2,9*b],
                           below_seven=r2 <= 441*b*X*X,
                           residual_over_X=math.sqrt(r2/(9*b))/X))
        previous = end
    if previous < len(history)-1:
        old, new = history[previous], history[-1]
        r2 = sum(max(0,10*v-7*u)**2 for u,v in zip(old,new))
        paired.append(dict(start=previous, end=len(history)-1, residual_square=[r2,9*b],
                           below_seven=r2 <= 441*b*X*X,
                           residual_over_X=math.sqrt(r2/(9*b))/X))
    return dict(k=k, first_reset_stages=len(history)-1, identity_padding=padding,
                failures=failures, max_R_over_X=max(v['residual_over_X'] for v in rec['stages']),
                final_Q_over_X=math.sqrt(qfinal/b)/X,
                weighted_debt=[debt.numerator,debt.denominator],
                weighted_debt_within_seven=debt <= 7*X*(1-rho**(k+8)),
                padded_bound=[total.numerator,total.denominator],
                padded_bound_over_X=float(total/X), bound_below_7p5X=2*total < 15*X,
                paired_stages=paired, all_paired_below_seven=all(p['below_seven'] for p in paired),
                weighted_cancellation_profiles=profiles)


if __name__ == '__main__':
    inherited = json.loads((DATA/'collatz_survivor_packets_20261015e.json').read_text())
    records = [r for r in inherited['first_reset'] if 8 <= r['k'] <= 18]
    records += [json.loads((DATA/f'collatz_compact_resets_20261015f_k{k}_r1.json').read_text())
                for k in range(19,24)]
    summaries = [analyse(rec) for rec in records]
    need([r['k'] for r in summaries if r['failures']] == [23], 'first failure in declared 8..23 census')
    need(all(r['all_paired_below_seven'] for r in summaries), 'declared paired finite census')
    need(all(r['bound_below_7p5X'] for r in summaries), 'discounted bound in declared finite census')
    result = dict(status='REFUTED one-reset uniform 7X; FINITE-EXACT repairs; all-scale replacement OPEN',
                  records=summaries)
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    path = DATA/'collatz_reset_debt_analysis_20261015f.json'
    path.write_text(data)
    for row in summaries:
        print(json.dumps({k:row[k] for k in ('k','max_R_over_X','final_Q_over_X',
                         'padded_bound_over_X','weighted_debt_within_seven','all_paired_below_seven')}))
    print(hashlib.sha256(data.encode()).hexdigest())
