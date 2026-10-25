#!/usr/bin/env python3
"""Exact controls for the repaired mixed-order binomial valuation proof."""
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib
import json


def valuation(n,p):
    if n == 0:
        raise ValueError('zero is outside this finite universe')
    out=0
    while n % p == 0:
        n//=p
        out+=1
    return out


def main():
    gates=0
    for p in (3,5,7):
        for M in range(1,100,2):
            a=valuation(M,p)
            for mu in range(1,M+1):
                actual=valuation(comb(M,mu),p)
                lower=a-valuation(mu,p)
                if actual < lower:
                    raise RuntimeError('binomial lower bound failed')
                coefficient=Fraction((p-1)*p**mu*comb(M,mu),p**(a+1))
                value=valuation(coefficient.numerator,p)-valuation(coefficient.denominator,p)
                predicted=mu-1-valuation(mu,p)
                if value < predicted or (mu==1 and value!=0) or (mu>=2 and value<1):
                    raise RuntimeError('normalized carry valuation failed')
                gates+=2
    witness=valuation(comb(15,7),3)
    if witness != 2 or valuation(15,3)-valuation(7,3) != 1:
        raise RuntimeError('strict-inequality hostile failed')
    print(json.dumps({'status':'FINITE-EXACT','gates':gates,
                     'universe':'p=3,5,7; every odd M=1..99 and mu=1..M; normalized j=1 coefficient',
                     'hostile':{'p':3,'M':15,'mu':7,'actual_binomial_valuation':2,'incorrect_equality_rhs':1},
                     'linear_normalized_valuation':0,
                     'higher_normalized_valuations':'at least 1 throughout the declared universe; general proof in repaired source',
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
