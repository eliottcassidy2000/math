#!/usr/bin/env python3
"""Exact finite hostile/positive tests of the all-orders norm-to-CF theorem.
No floating point, number-field package, or finite-to-infinite inference.
"""
from fractions import Fraction
from math import gcd
import json
from pathlib import Path


def check(condition, message="exact check failed"):
    if not condition:
        raise AssertionError(message)


def form(k, p, q):
    return p**k - sum(p**j * q**(k-j) for j in range(k))


def form_horner(k, p, q):
    # Independent homogeneous Horner evaluation.
    acc, qpow = 1, q
    for _ in range(k):
        acc = p*acc-qpow
        qpow *= q
    return acc


def isolate(k, bits=256):
    # h_k is strictly increasing on x>0; f_k and h_k have the same sign.
    lo, hi = Fraction(3, 2), Fraction(2)
    check(form(k, lo.numerator, lo.denominator) < 0)
    check(form(k, hi.numerator, hi.denominator) > 0)
    for _ in range(bits):
        mid = (lo+hi)/2
        value = form(k, mid.numerator, mid.denominator)
        check(value != 0, 'No positive rational root: rational root theorem.')
        if value < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def certified_cf(lo, hi, maxq=None):
    """Continued-fraction terms common to an open isolating interval."""
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    result=[]
    while True:
        a = lo.numerator//lo.denominator
        b = hi.numerator//hi.denominator
        if a != b:
            break
        p, q = a*pm1+pm2, a*qm1+qm2
        result.append((a,p,q))
        if maxq is not None and q > maxq:
            break
        pm2,pm1,qm2,qm1 = pm1,p,qm1,q
        # Irrationality and a common floor make both denominators positive.
        if lo == a:
            break
        lo,hi = 1/(hi-a),1/(lo-a)
    return result


def main():
    rows=[]
    Q=150
    for k in range(3,13):
        lo,hi=isolate(k)
        cf=certified_cf(lo,hi,Q)
        convs={(p,q) for _,p,q in cf}
        units=[]
        residual=[]
        checks=0
        for q in range(1,Q+1):
            for p in range(1,3*q+1):
                m=form(k,p,q)
                check(m == form_horner(k, p, q))
                check((p - q) * m == p ** (k + 1) - 2 * p ** k * q + q ** (k + 1))
                checks += 1
                d=gcd(p,q)
                if abs(m)==1:
                    check(d == 1 and q < p <= 2 * q)
                    check((p, q) in convs)
                    units.append([p,q,m])
                if q<p<2*q and 4*abs(m)*q*q <= p**k:
                    check((p // d, q // d) in convs)
                    # Independent rational enclosure certifies Legendre's bound.
                    x=Fraction(p,q)
                    error_upper=max(abs(x-lo),abs(x-hi))
                    check(error_upper < Fraction(1, 2 * q * q))
                    residual.append([p,q,m])
        rows.append({'k':k,'pairs_checked':checks,'unit_solutions':units,
          'residual_gate_count':len(residual), 'residual_gate_rows':residual,
          'certified_cf_through_first_q_above_150':cf,
          'root_interval':[str(lo),str(hi)]})
    # Hostiles: sign/root selection and unnecessary norm=1.
    check(form(4, -1, 1) == 1)
    lo,hi=isolate(3)
    cf=certified_cf(lo,hi,100)
    check((11, 6) in {(p, q) for _, p, q in cf})
    check(form(3, 11, 6) == -7)
    check(form(3, 3, 2) == -11)
    check((3, 2) not in {(p, q) for _, p, q in cf})
    check(cf[2][0] == 5)
    # Find a large partial quotient whose preceding convergent is not norm one.
    lo,hi=isolate(3,4096)
    cf=certified_cf(lo,hi)
    large_nonunit=[]
    for n in range(len(cf)-1):
        a,p,q=cf[n]
        next_a=cf[n+1][0]
        m=form(3,p,q)
        if next_a>=305 and abs(m)!=1:
            large_nonunit.append({'index':n,'p':str(p),'q':str(q),
                                  'norm':str(m),'next_partial_quotient':next_a,
                                  'certified_cf_prefix_through_next':[z[0] for z in cf[:n+2]]})
            break
    out={'status':'FINITE-EXACT; no all-height solution classification',
       'universe':{'k':'3..12 inclusive','q':'1..150 inclusive',
                   'p':'1..3q inclusive','pairs_total':sum(r['pairs_checked'] for r in rows)},
       'root_isolation_bits':256,'rows':rows,
       'hostiles':{'wrong_root_unit':{'k':4,'p':-1,'q':1,'norm':1},
        'convergent_nonunit':{'k':3,'p':11,'q':6,'norm':-7},
        's_supported_nonconvergent':{'k':3,'p':3,'q':2,'norm':-11,'S':[11],'U':'-12/11'},
        'unit_norm_not_305':{'k':3,'p':2,'q':1,'norm':1,'next_partial_quotient':5},
        'large_partial_quotient_nonunit':large_nonunit,
        'large_partial_quotient_search':{'k':3,'isolation_bits':4096,
                                        'certified_terms':len(cf)}}}
    path=Path(__file__).with_name('verify_uniform_norm.output.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'pairs_total':out['universe']['pairs_total'],
       'rows':[{'k':r['k'],'unit_solutions':r['unit_solutions'],
       'residual_gate_count':r['residual_gate_count']} for r in rows],
       'hostiles':out['hostiles']},indent=2))

if __name__=='__main__': main()
