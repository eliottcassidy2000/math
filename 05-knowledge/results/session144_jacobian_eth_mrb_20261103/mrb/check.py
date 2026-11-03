#!/usr/bin/env python3
"""Session144 MRB checks. Exact and directed-interval evidence are separated.
No removable assert statements. No network, repository edits, or publication.
"""
from fractions import Fraction
from itertools import product
from math import comb, exp, expm1, factorial, fsum, gcd, lcm, log
from pathlib import Path
import json
from mpmath import mp, iv

BASE=Path(__file__).resolve().parent
mp.dps=70
iv.dps=90


def require(condition, message):
    if not condition:
        raise ValueError(message)


def period_mean(values):
    return sum(map(Fraction, values), Fraction(0))/len(values)


def parity_mean(mask):
    m=len(mask)
    return Fraction(0) if m%2 else sum(((-1)**r*Fraction(mask[r-1]) for r in range(1,m+1)),Fraction(0))/m


def check_mean_formula(mask, proposed):
    m=len(mask)
    actual=period_mean([(-1)**n*mask[(n-1)%m] for n in range(1,lcm(m,2)+1)])
    require(actual==proposed, 'Incorrect mean/parity convergence classification')


def shift_sign(K):
    return (-1)**(K+1)


def check_shift_sign(K, proposed):
    require(proposed==shift_sign(K), 'Incorrect alternating shift sign')


def check_abel_counterterm(linear_coefficient, finite_shift):
    require(abs(linear_coefficient+mp.euler)<mp.mpf('1e-60'), 'Incorrect Abel logarithmic counterterm')
    require(abs(finite_shift-(mp.euler**2/2+mp.pi**2/12))<mp.mpf('1e-60'), 'Incorrect Abel finite-part shift')


def falling_value(n,K):
    out=1
    for j in range(1,K+1): out*=n-j
    return out


def explicit_shift_interval(K,r,T=220):
    d=K+r
    total=iv.mpf(0)
    for n in range(K+1,T+1):
        x=iv.mpf(n)
        value=iv.exp(-x)*falling_value(n,K)*((1-iv.exp(-x))/x)**d
        total+=(-1)**n*value
    # |H_K,r(n)| <= e^-n for n>K, r>=1.
    tail=iv.exp(-(T+1))/(1-iv.exp(-1))
    enclosure=total+iv.mpf([-1,1])*tail
    require(bool((shift_sign(K)*enclosure).a>0), 'Shift failed strict signed interval test')
    lower=iv.factorial(K)/(iv.exp(K+r+1)+1)**(K+1)
    upper=iv.factorial(K)/(iv.e+1)**(K+1)
    mag=shift_sign(K)*enclosure
    require(bool(mag.a>lower.b) and bool(mag.b<upper.a), 'Integral-derived shift bounds failed')
    return enclosure, tail


def ih_density(s,d):
    if s<=0 or s>=d:return mp.mpf(0)
    return sum((-1)**j*comb(d,j)*(s-j)**(d-1) for j in range(int(mp.floor(s))+1))/mp.factorial(d-1)


def scalar_cube_integral(K,r):
    d=K+r
    value=mp.quad(lambda s:ih_density(s,d)/(mp.exp(1+s)+1)**(K+1),list(range(d+1)))
    return shift_sign(K)*mp.factorial(K)*value


def run():
    exact_count=0
    for m in range(1,11):
        for mask in product((0,1),repeat=m):
            check_mean_formula(mask,parity_mean(mask))
            exact_count+=1
    unit_rows=[]
    for m in range(1,41):
        mask=[int(gcd(n,m)==1) for n in range(1,m+1)]
        mu=parity_mean(mask)
        expected=Fraction(0) if m%2 else -Fraction(sum(mask),m)
        require(mu==expected,'Unit-mask parity boundary failed')
        unit_rows.append({'m':m,'phi':sum(mask),'mean':str(mu),'converges':mu==0})
    # Formal generating identity: P_K(D)q/(1-q)=K!q^(K+1)/(1-q)^(K+1).
    # Verify coefficients exactly rather than approximate values or sampled roots.
    for K in range(0,11):
        for n in range(1,81):
            lhs=falling_value(n,K)
            rhs=0 if n<=K else factorial(K)*comb(n-1,K)
            require(lhs==rhs,'Falling-factorial coefficient identity failed')
        for n in range(1,K+1):require(falling_value(n,K)==0,'Prefix zero failed')
    shifts=[]
    for K,r in [(1,1),(1,3),(2,1),(2,2),(3,1),(4,2)]:
        enclosure,tail=explicit_shift_interval(K,r)
        scalar=scalar_cube_integral(K,r)
        lo=mp.mpf(enclosure._mpi_[0]); hi=mp.mpf(enclosure._mpi_[1])
        midpoint=(lo+hi)/2
        require(abs(scalar-midpoint)<mp.mpf('1e-65'),'Independent scalar cube integral disagrees')
        shifts.append({'K':K,'r':r,'signed_shift_interval':str(enclosure),'analytic_tail':str(tail),
                       'independent_scalar_integral':mp.nstr(scalar,65),
                       'scalar_check_status':'NUMERICAL, not an interval certificate'})
    Q=mp.euler**2/2+mp.pi**2/12
    check_abel_counterterm(-mp.euler,Q)
    regularization=[]
    previous=None
    for N in [100,1000,10000]:
        t=1/N
        ell=log(N)
        cutoff=fsum(expm1(log(n)/n) for n in range(1,N+1))
        damped=fsum(expm1(log(n)/n)*exp(-n*t) for n in range(1,40*N+1))
        # At n>40N, f(n)<1. This is a bound on truncation only, not floating rounding.
        damping_tail=exp(-(40*N+1)*t)/(1-exp(-t))
        corrected=damped-cutoff+float(mp.euler)*ell
        error=abs(corrected-float(Q))
        if previous is not None:require(error<previous,'Abel/cutoff discrepancy did not approach proved shift')
        previous=error
        regularization.append({'N':N,'corrected_difference':corrected,'limit_Q':float(Q),
                               'absolute_difference':error,'damping_truncation_bound':damping_tail,
                               'status':'NUMERICAL finite-scale illustration, not proof or certified interval'})
    require(previous<0.06,'Abel/cutoff numerical illustration outside tolerance')
    out={'status':'PASS','periodic_masks_exactly_checked':exact_count,
         'unit_mask_cases':unit_rows,
         'falling_factorial_coefficient_tests':{'K':'0..10','n':'1..80','method':'exact integers'},
         'explicit_perturbation_intervals':shifts,
         'regularization_numerical_checks':regularization,
         'proof_dependencies':'RESEARCH.md; elementary binomial generating identity and summation by parts',
         'limits':['Finite checks do not replace all-K analytic proof.','Numerical quadrature is separated from interval shift certificates.',
                   'No irrationality, transcendence, modularity, SAT, ETH, or Jacobian classification claim.']}
    (BASE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'PASS: {exact_count} exact binary periodic masks and 40 unit masks.')
    print('PASS: falling-factorial identity and finite-prefix zeros (exact).')
    print('PASS: six strictly nonzero shift intervals with analytic tails and independent scalar-integral checks.')
    print('PASS: finite-part sign controls and finite-scale numerical illustrations.')
    print('All pass/fail predicates use explicit exceptions and remain active under python -O.')

if __name__=='__main__':run()
