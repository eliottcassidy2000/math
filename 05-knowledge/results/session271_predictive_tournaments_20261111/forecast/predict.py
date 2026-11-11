#!/usr/bin/env python3
"""Logarithmic A000568 forecasts. Never forms or rounds the integer count.
Dependency: mpmath. Model envelopes proved in PROOF.md; numeric evaluation
is arbitrary precision, not directed interval arithmetic.
"""
import argparse, json, math
from collections import Counter
from pathlib import Path
import mpmath as mp


def parts(m, lo=3):
    if m == 0:
        yield ()
    for j in range(lo, m + 1, 2):
        for p in parts(m-j, j):
            yield (j,) + p


def core(p):
    m, r = sum(p), len(p)
    e = sum((j-1)//2 for j in p)
    e += sum(math.gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
    z = math.prod(j**v * math.factorial(v) for j,v in Counter(p).items())
    return m, r, e, z


def logweight(n, p):
    m,r,e,z = core(p)
    deficit = (n-m)*(m-r) + m*(m-1)//2 - e
    # Fixed support: never subtract two huge log-gamma values.
    return sum(mp.log(n-j) for j in range(m)) - mp.log(z) - deficit*mp.log(2)


def logsum(values):
    values = list(values)
    if not values:
        return mp.ninf
    top = max(values)
    return top + mp.log(sum(mp.exp(x-top) for x in values))


def tail_log(n, first):
    """Natural log of rigorous symbolic upper bound for omitted support >=first."""
    if first > n:
        return mp.ninf
    first = max(3, first)
    def f(m):
        return m*mp.log(n) - mp.mpf(m*(2*n-m-1))*mp.log(2)/3
    return mp.log(n-first+1) + max(f(first), f(n))


def choose_cutoff(n, digits, cap=24):
    for s in range(2, min(n,cap)+1):
        if tail_log(n,s+1) <= -digits*mp.log(10):
            return s, True
    s = min(n,cap)
    return s, tail_log(n,s+1) <= -digits*mp.log(10)


def stirling(n, terms=4):
    n = mp.mpf(n)
    main = (n+mp.mpf('0.5'))*mp.log(n)-n+mp.log(2*mp.pi)/2
    for k in range(1, terms+1):
        main += mp.bernoulli(2*k)/(2*k*(2*k-1)*n**(2*k-1))
    error = abs(mp.bernoulli(2*terms+2)/((2*terms+2)*(2*terms+1)*n**(2*terms+1)))
    return main, error


def fmt(x, d=60):
    return mp.nstr(x, d)


def forecast(n, support=12, digits=50, guard=90):
    if n < 3:
        raise ValueError('forecast domain n >= 3')
    nd = len(str(n))
    with mp.workdps(2*nd + guard):
        ln10 = mp.log(10)
        base = (mp.mpf(n)*(n-1)/2*mp.log(2)-mp.loggamma(n+1))/ln10
        terms = [(p,logweight(n,p)) for m in range(3,min(n,support)+1) for p in parts(m)]
        correction = logsum(w for _,w in terms)
        tail = tail_log(n, max(3,support+1))
        total_upper = logsum([correction,tail])
        lead = logweight(n,(3,))
        # Separate logarithms retain tiny quantities that disappear in 1+epsilon.
        mantissa = mp.power(10, mp.frac(base))
        exponent = int(mp.floor(base))
        sf,se = stirling(n)
        stbase = (mp.mpf(n)*(n-1)/2*mp.log(2)-sf)/ln10
        cutoff, success = choose_cutoff(n,digits)
        row = {
            'n':str(n), 'working_decimal_digits':mp.mp.dps,
            'identity_log10':fmt(base,2*nd+55),
            'identity_scientific_mantissa':fmt(mantissa,50),
            'identity_scientific_exponent':str(exponent),
            'interpretation':'prediction, not exact integer digits; see separate model and numeric errors',
            'log10_relative_correction_leading_3cycle':fmt(lead/ln10),
            'support_cutoff':support, 'retained_nonidentity_types':len(terms),
            'log10_retained_relative_correction':fmt(correction/ln10),
            'log10_omitted_relative_tail_bound':fmt(tail/ln10),
            'log10_identity_relative_error_upper_bound':fmt(total_upper/ln10),
            'log10_log10count_correction_upper_bound':fmt(total_upper/ln10-mp.log10(ln10)),
            'stirling_four_term_log10':fmt(stbase,2*nd+55),
            'stirling_log10_absolute_error_bound':fmt(se/ln10),
            'stirling_vs_loggamma_abs_discrepancy':fmt(abs(stbase-base)),
            'adaptive_relative_decimal_digits_requested':digits,
            'adaptive_support_cutoff':cutoff, 'adaptive_bound_meets_request':bool(success),
            'numerics':'mpmath, guard digits plus precision doubling control, not certified intervals',
            'expected_labelled_automorphisms_minus_one_log10_lower':fmt(correction/ln10),
            'expected_labelled_twin_pairs_log10':fmt((mp.log(n*(n-1)//2)-(n-2)*mp.log(2))/ln10),
            'local_growth_log10_identity_ratio':fmt(n*mp.log10(2)-mp.log10(n+1)),
            'local_growth_relative_error_upper_log10':fmt(logsum([total_upper,tail_log(n+1,3)])/ln10),
        }
        # Diagnostic convergence only: bound concerns mathematical truncation.
        with mp.workdps(2*mp.mp.dps):
            b2=(mp.mpf(n)*(n-1)/2*mp.log(2)-mp.loggamma(n+1))/mp.log(10)
            row['precision_doubling_log10_abs_discrepancy']=fmt(abs(base-b2))
        return row


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n', nargs='+', default=['300','1000','10000','1000000','1000000000','1000000000000','1'+'0'*30,'1'+'0'*100,'1'+'0'*1000])
    ap.add_argument('--support',type=int,default=12)
    ap.add_argument('--relative-digits',type=int,default=50)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    lines=[json.dumps(forecast(int(n),a.support,a.relative_digits),sort_keys=True) for n in a.n]
    out='\n'.join(lines)+'\n'
    if a.output: a.output.write_text(out)
    else: print(out,end='')

if __name__=='__main__':main()
