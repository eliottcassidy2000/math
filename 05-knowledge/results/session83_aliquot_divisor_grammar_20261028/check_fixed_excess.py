#!/usr/bin/env python3
"""Scoped fixed-excess identities and complete small inverse-fiber checks."""
import argparse
import json
from math import gcd, isqrt
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0)+1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def sigma_factored(n):
    value = 1
    for p, e in factor(n).items():
        value *= (p**(e+1)-1)//(p-1)
    return value


def sigma_pairs(n):
    value = 0
    for d in range(1, isqrt(n)+1):
        if n % d == 0:
            value += d
            if d*d != n:
                value += n//d
    return value


def proper_sum_sieve(cap):
    sums = [0]*(cap+1)
    for d in range(1, cap//2+1):
        for n in range(2*d, cap+1, d):
            sums[n] += d
    return sums


def repeated_odd_prime_exceptions(c):
    """Complete b>=2 solutions 2^a*p^b for c>=0; no source cutoff."""
    out = []
    for p in range(3, c+1, 2):
        if not trial_prime(p):
            continue
        b, R = 2, 1+p
        while R <= c+1:
            if (c+1) % R == 0:
                h = (c+1)//R
                A = p+h
                if A >= 4 and A & (A-1) == 0:
                    a = A.bit_length()-2
                    n = 2**a*p**b
                    require(sigma_factored(n)-2*n == c, 'Bad exceptional reconstruction')
                    out.append({'n':n,'a':a,'p':p,'b':b})
            R += p**b
            b += 1
    return sorted(out, key=lambda x:x['n'])


def verify():
    cap = 1_000_000
    sums = proper_sum_sieve(cap)
    for n in range(1, 10001):
        require(sums[n]+n == sigma_pairs(n), f'Divisor-pair oracle at {n}')
    c2 = [n for n in range(1, cap+1) if sums[n] == n+2]
    expected = [20,104,464,650,1952,130304,522752]
    require(c2 == expected, 'Million-source excess-two census changed')
    for n in c2:
        require(sigma_factored(n)-2*n == 2, f'Independent excess at {n}')
        if len(factor(n)) <= 2:
            a = factor(n).get(2, 0)
            p = 2**(a+1)-3
            require(a>=1 and trial_prime(p) and n == 2**a*p, 'Two-prime theorem failed')
        proper = [d for d in range(1, isqrt(n)+1) if n%d==0]
        proper += [n//d for d in proper if d*d != n]
        require(all(sigma_factored(d)<2*d for d in proper if d<n), 'Primitive abundant control')
    primes = [p for p in range(3, 102, 2) if trial_prime(p)]
    identity_cases = 0
    for a in range(1, 13):
        for p in primes:
            for b in range(1, 7):
                n=2**a*p**b
                S=(2**(a+1)-1)*sum(p**j for j in range(b+1))
                R=sum(p**j for j in range(b))
                require(S-2*n+1 == (2**(a+1)-p)*R, 'Fixed-excess factorization')
                identity_cases += 1
    pair_cases = 0
    small_primes = [p for p in range(2,44) if trial_prime(p)]
    for m in range(1,41):
        S=sigma_factored(m)
        A=S-m
        D=2*m-S
        for i,p in enumerate(small_primes):
            for q in small_primes[i+1:]:
                if gcd(m,p*q) != 1:
                    continue
                n=m*p*q
                Sn=sigma_factored(n)
                c=Sn-2*n
                t=Sn-n
                require((D*p-S)*(D*q-S)==2*m*S-D*c, 'Excess factor-pair identity')
                require((A*p+S)*(A*q+S)==A*t+m*S, 'Fiber factor-pair identity')
                pair_cases += 1
    # Vanishing coefficients destroy reverse implications, not forward identities.
    require((0*2+1)*(0*3+1)==0*7+1 and sigma_factored(6)-6 != 7,
            'A=0 inverse-fiber tautology hostile')
    require((0*5-12)*(0*7-12)==2*6*12 and sigma_factored(210)-420 != 0,
            'D=0 fixed-excess tautology hostile')
    fibers={}
    expected_fibers={16:[12,26],20:[34,361],21:[18,51,91],22:[20,38],28:[28],36:[24]}
    for t,want in expected_fibers.items():
        bound=(t-1)**2
        got=[n for n in range(1,bound+1) if sums[n]==t]
        oracle=[n for n in range(1,bound+1) if sigma_pairs(n)-n==t]
        require(got==want==oracle, f'Complete fiber {t}')
        fibers[str(t)]={'sources':got,'proved_source_bound':bound}
    exceptions={str(c):repeated_odd_prime_exceptions(c) for c in [0,1,2,3,4,8,12]}
    require(exceptions['3']==[{'n':18,'a':1,'p':3,'b':2}], 'Excess-three exception')
    require(exceptions['12']==[{'n':54,'a':1,'p':3,'b':3}], 'Excess-twelve exception')
    require(all(not exceptions[str(c)] for c in [0,1,2,4,8]), 'Unexpected repeated-prime exception')
    for M in [6,28]:
        for q in small_primes:
            if M%q:
                require(sigma_factored(M*q)-2*M*q==2*M, 'Perfect-prefix dilation')
    require(sigma_factored(18)-36==3, 'q dividing perfect prefix is a hostile')
    require(sigma_factored(12)-24==4, 'Second noncoprime hostile')
    require(sigma_factored(18)-36 != 2**2-1-9, 'Composite prime-slot hostile')
    rows=[]
    for n in [12,18,20,24,28,54,88,104,464,650,1952]:
        factors=factor(n)
        tau=1
        for e in factors.values():tau*=e+1
        rows.append({'n':n,'s':sums[n],'excess':sums[n]-n,'proper_divisor_terms':tau-1,'factors':factors})
    return {'status':'PROVED identities and FINITE-EXACT explicitly bounded census',
            'source_census_cap':cap,'excess_two_sources':c2,'direct_pair_oracle_cap':10000,
            'two_prime_identity_cases':identity_cases,'guarded_factor_pair_cases':pair_cases,
            'complete_small_fibers':fibers,'repeated_odd_prime_exceptions':exceptions,'example_rows':rows,
            'hostile_controls':['650 outside two-prime family','prime slots cannot be composite',
                               'prime dilation requires coprimality','fixed-excess does not count proper divisors']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('fixed_excess_results.json'))
    args=parser.parse_args()
    result=verify()
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('Verified fixed-excess identities, complete fibers, finite repeated-prime exceptions, and census through1,000,000.')


if __name__=='__main__':main()
