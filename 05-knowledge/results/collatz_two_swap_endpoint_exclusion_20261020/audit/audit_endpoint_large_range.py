#!/usr/bin/env python3
"""Independent large-range replay; log3 uses log(9/8), unlike the producer."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json

checks = 0

def need(test, reason):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(reason)


def log_interval(numerator, denominator, terms=112):
    x = Q(numerator-denominator, numerator+denominator)
    power = x
    partial = Q(0)
    for j in range(terms):
        partial += 2*power/(2*j+1)
        power *= x*x
    tail = 2*power/((2*terms+1)*(1-x*x))
    return partial, partial+tail


def convergents_until(rational, cutoff):
    partials = []
    fractions = []
    numer0,numer1,denom0,denom1 = 0,1,1,0
    while True:
        a = rational.numerator // rational.denominator
        partials.append(a)
        numer0,numer1 = numer1,a*numer1+numer0
        denom0,denom1 = denom1,a*denom1+denom0
        fractions.append((numer1,denom1))
        if denom1 >= cutoff:
            return partials, fractions
        rational = 1/(rational-a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-dir', type=Path,
                        default=Path(__file__).resolve().parent.parent/'large-range')
    original = parser.parse_args().producer_dir
    # Producer files are read-only evidence, and the default survives packet relocation.
    source = original/'receipt.json'
    production = json.loads(source.read_text())
    need(source.read_bytes() == (original/'receipt.optimized.json').read_bytes(), 'producer normal/-O differ')
    script = original/'certify_endpoint_large_range.py'
    need(sha256(script.read_bytes()).hexdigest() == production['script_sha256'], 'producer source hash')
    l2,u2 = log_interval(2,1)
    l98,u98 = log_interval(9,8)
    l3,u3 = (3*l2+l98)/2,(3*u2+u98)/2
    lower,upper = l3/u2,u3/l2
    cf_lo,pairs_lo = convergents_until(lower,10**12)
    cf_hi,pairs_hi = convergents_until(upper,10**12)
    need(cf_lo == cf_hi and pairs_lo == pairs_hi, 'independent endpoint CF ambiguity')
    need(cf_lo == production['continued_fraction_prefix'], 'producer continued fraction mismatch')
    need(pairs_lo[-1] == (8573543875303,5409303924479), 'terminal pair')
    need(pairs_lo[-1][1] >= 10**12 and pairs_lo[-2][1] < 10**12, 'cutoff coverage')
    need(u2 < Q(7,10) and u3 < Q(11,10), 'Matveev A1,A2')
    need(2**9 < 23**2, 'Matveev 2^4.5 <23')
    constant = Q(7,5)*30**5*23*Q(7,10)*Q(11,10)
    need(constant == 602494200 and constant < 10**9, 'Matveev numerical constant')
    l54,u54 = log_interval(5,4)
    cutoff_log_upper = u2+12*(3*u2+u54)
    need(cutoff_log_upper < 30, 'log(2*10^12)<30')
    need(Q(10**12,9) > 31*10**9, 'cutoff inequality')
    need(Q(1,9) > Q(10**9,10**12), 'cutoff derivative positive')
    need(Q(400,2**11) < l2, 'Legendre threshold via e^(100/9)>2^11')
    need(2**11 > 2, 'signed logarithm conversion uses |U|<1/2')
    rows=[]
    for q,p in pairs_lo:
        if not 100<=p<10**12:
            continue
        a,b = q*l2-p*u3,q*u2-p*l3
        need(a>0 or b<0,'independent gap sign')
        gap_min = a if a>0 else -b
        coarse_upper = Q(2,2**min(p//9,256))
        need(gap_min > coarse_upper,'independent absolute-gap contradiction')
        rows.append({'p':p,'q':q,'sign':'positive' if a>0 else 'negative',
                     'gap_lower_bound':str(gap_min),'coarse_upper_bound':str(coarse_upper)})
    need(len(rows)==18,'independent row count')
    need([(x['q'],x['p']) for x in production['excluded_large_convergents']] == [(x['q'],x['p']) for x in rows],
         'producer row universe mismatch')
    print(json.dumps({'status':'FINITE-EXACT independent large-range audit',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'producer_script_sha256':production['script_sha256'],
        'checks':checks,'different_log_decomposition':'log3=(3 log2+log(9/8))/2; 112 terms',
        'continued_fraction_prefix':cf_lo,'terminal_convergent':pairs_lo[-1],
        'row_count':len(rows),'rows':rows,
        'external_dependency':'Bugeaud–Mignotte–Siksek, arXiv:math/0403046, Theorem 9.4 real case, printed p.16; inspected independently.',
        'logical_cautions':['Matveev applies to U=2^q/3^p-1, not directly to the additive logarithm L.',
                            'The absolute logarithm satisfies |L|<2e^(-p/9) for both signs.',
                            'The finite residual p<100 is a separate dependency.']},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
