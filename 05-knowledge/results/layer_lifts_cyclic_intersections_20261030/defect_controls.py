#!/usr/bin/env python3
"""Exact finite controls for the all-layer defect decomposition and hostiles."""
from pathlib import Path
from math import gcd, lcm
from fractions import Fraction
import argparse, json

checks = 0
def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

def valuation(x):
    if not x:
        raise ValueError('infinite valuation is not a finite integer')
    result = 0
    while x % 5 == 0:
        x //= 5
        result += 1
    return result

def evaluate(word):
    P = Q = 1
    B = charge = scalar = debt = kappa_sum = nonone_cost = 0
    for a in word:
        P, Q, B = 3*P, 2**a*Q, 3*B+Q
        need(P != Q, 'distinct positive powers of two and three')
        charge += valuation(gcd(abs(P-Q), B))
        scalar += min(valuation(P-Q), 3)
        if a > 1:
            Z = B+Q-P
            need(Z > 0, 'positive integral defect after a non-one letter')
            h, kappa = valuation(Z), valuation(2**a-2)
            debt += max(h-kappa, 0)
            kappa_sum += kappa
            nonone_cost += a
            need(4*kappa <= a-1, 'cost of native added defect')
    ones = word.count(1)
    need(8*charge <= 5*ones+8*kappa_sum+8*debt, 'untruncated epoch bound')
    need(8*charge <= 5*sum(word)+8*debt-3*nonone_cost-2*(len(word)-ones),
         'untruncated discounted debt bound')
    if len(word)-ones <= 1:
        need(debt == 0 and 8*charge <= 5*sum(word), 'all-layer one-defect grammar')
    return dict(clock=sum(word), charge=charge, scalar_truncated_charge=scalar,
                jump_debt=debt, kappa_sum=kappa_sum, nonone_cost=nonone_cost)

def compositions(n):
    if n == 0:
        yield ()
    else:
        for a in range(1, n+1):
            for tail in compositions(n-a):
                yield (a,)+tail

def main():
    words = 0
    for A in range(1, 17):
        for w in compositions(A):
            evaluate(w)
            words += 1
    w = (3,)+(1,)*16
    scalar = evaluate(w)
    need(scalar['scalar_truncated_charge'] == 12 and scalar['charge'] == 0,
         'multiplier-only hostile')
    need(50*12 > 31*sum(w), 'scalar bound fails at clock nineteen')
    dark = []
    r = 0
    for e in range(1, 8):
        q = 5**e
        candidates = [0] if e == 1 else [r+j*5**(e-2) for j in range(5)]
        good = [t for t in candidates if (pow(9,t,q)*pow(pow(4,t,q),-1,q)+4) % q == 0]
        need(len(good) == 1, 'unique five-adic exponent lift')
        r = good[0]
        need(9*pow(32,-1,5)*pow(9*pow(4,-1,5),r,5) % 5 == 2,
             'deep defect has nonreturning multiplier phase')
        # Full integer orbit checks through height five, modular lift beyond.
        item = dict(required_height=e, r=r, clock=2*r+5, multiplier_mod5=2)
        if e <= 5:
            result = evaluate((2,)+(1,)*(2*r)+(3,))
            need(result['charge'] == 0 and result['jump_debt'] >= e, 'dark height exact control')
            item.update(result)
        dark.append(item)
    banks = []
    for moduli,n in [((3,5,7),4), ((11,19,27),5), ((1001,343,57),7)]:
        L = lcm(*moduli)
        need(L > 2**(n*(n-1)//2), 'unscaled tournament deficit height coverage')
        need(L <= 2**((n+1)*n//2), 'next order is not covered by this bank')
        banks.append(dict(moduli=moduli,lcm=L,last_full_cube_order=n))
    L=lcm(1001,343,57)
    need(L < 1001*343*57 and all(L%m==0 for m in (1001,343,57)), 'product versus lcm hostile')
    return dict(checks=checks, positive_compositions_through_clock=16, words=words,
                scalar_hostile=dict(word=w,**scalar), dark_defects=dark, residue_banks=banks,
                scope='finite checks support the written all-word proof; no universal Collatz conclusion')

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('defects.json'))
    args=parser.parse_args()
    report=json.dumps(main(),indent=2,sort_keys=True)+'\n'
    args.output.write_text(report)
    print(report,end='')
