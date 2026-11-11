#!/usr/bin/env python3
"""Forecast ordinary-class path-cover means using the census error envelope.
No exact A000568 or H evaluation. mpmath numerics are not interval arithmetic.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import mpmath as mp

FORECAST = Path(__file__).resolve().parents[1] / 'forecast' / 'predict.py'
spec = importlib.util.spec_from_file_location('count_forecast', FORECAST)
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


def row(n, r, support=12):
    assert n >= 3 and 1 <= r <= n
    with mp.workdps(2 * len(str(n)) + 100):
        ln2, ln10 = mp.log(2), mp.log(10)
        terms = [f.logweight(n, p) for m in range(3, min(n, support)+1) for p in f.parts(m)]
        lb = f.logsum(terms + [f.tail_log(n, max(3,support+1))])
        if lb >= 0:
            raise ValueError('This large-scale table requires delta upper bound <1.')
        lmu = mp.loggamma(n+1) + mp.loggamma(n) - mp.loggamma(r) - mp.loggamma(n-r+1) - (n-r)*ln2
        exact_universal = r >= n-1
        if exact_universal:
            lower = upper = mp.ninf
        else:
            # Conservative B instead of B/(1+B), deliberately avoids implying
            # tiny correction has been retained in a displayed floating mean.
            lower = lb/ln10
            x = (n-r)*ln2
            le = mp.log(mp.expm1(x)) if x < 100 else x + mp.log1p(-mp.exp(-x))
            upper = (lb+le)/ln10
        return {'n':str(n),'r':str(r),'mean_predictor_log10':f.fmt(lmu/ln10,60),
                'mean_predictor_mantissa':f.fmt(mp.power(10,mp.frac(lmu/ln10)),35),
                'mean_predictor_exponent':str(int(mp.floor(lmu/ln10))),
                'downward_relative_error_log10_upper':f.fmt(lower),
                'upward_relative_error_log10_upper':f.fmt(upper),
                'bounded_event_absolute_error_log10_upper':f.fmt(lb/ln10),
                'universal_coordinate':exact_universal,
                'error_scope':'symbolic mathematical bounds evaluated numerically; not interval-certified digits'}


def controls():
    # New mixture identity on a synthetic stabilizer inventory. This does not
    # re-enumerate tournaments or recheck an inherited sequence cache.
    from fractions import Fraction
    gs=[1,1,3,9,27]
    a=len(gs); b=sum(Fraction(1,g) for g in gs); delta=Fraction(a,b)-1
    bias=[(1-Fraction(1,g))/(a-b) for g in gs]
    lab=[Fraction(1,g)/b for g in gs]
    ordinary=[Fraction(1,a)]*a
    assert all(u==l/(1+delta)+delta*d/(1+delta) for u,l,d in zip(ordinary,lab,bias))
    tv=sum(abs(u-l) for u,l in zip(ordinary,lab))/2
    assert tv<=delta/(1+delta)
    # Positive multilinear contractions preserve coordinate boxes.
    import random
    rng=random.Random(271)
    for _ in range(100):
        b1,b2=[rng.uniform(1,10) for _ in range(2)],[rng.uniform(1,10) for _ in range(2)]
        lo1,lo2=0.9,0.8; hi1,hi2=1.1,1.3
        x=[v*rng.uniform(lo1,hi1) for v in b1]; y=[v*rng.uniform(lo2,hi2) for v in b2]
        c=[[rng.uniform(0,1) for j in range(2)] for i in range(2)]
        F=lambda p,q:sum(c[i][j]*p[i]*q[j] for i in range(2) for j in range(2))
        assert lo1*lo2*F(b1,b2)<=F(x,y)<=hi1*hi2*F(b1,b2)
    return {'synthetic_mixture_identity':'PASS exact rational arithmetic',
            'positive_contraction_boxes':100,'inherited_exact_counts_rechecked':0}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',nargs='+',type=int,default=[300,1000,1000000,10**12])
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    out={'controls':controls(),'predictions':[row(n,r) for n in a.n for r in sorted(set([1,n//2,n-2,n-1,n]))]}
    s=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.write_text(s)
    else:print(s,end='')

if __name__=='__main__':main()
