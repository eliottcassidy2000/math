#!/usr/bin/env python3
"""Independent finite identity proof of the endpoint39 certificate.

Reconstruct actual multinomial fibres at 73 distinct parameter points;
use SymPy polynomial division/inversion and Berkowitz characteristic
polynomials rather than the producer's polynomial-matrix trace recursion.
The proved degree bound <=72 upgrades these evaluations to polynomial
identities. No producer module is imported.
"""
from pathlib import Path
from math import factorial, comb, gcd
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'05-knowledge/results/long_frontier_sep06_endpoint39.out'
t = S.Symbol('t')
GATES = 0
TRACE = hashlib.sha256()


def gate(ok, label):
    global GATES
    if ok != True:
        raise RuntimeError(label+': '+str(ok))
    GATES += 1
    TRACE.update((label+'\n').encode())


def fibre(g,m):
    facts = [factorial(i) for i in range(m+1)]
    counts = {}
    for nc in range(m+1):
        for nb in range(m-nc+1):
            na=m-nb-nc
            if -39*na+(2*g-39)*nb+(3*g-39)*nc==0:
                counts[(na,nb,nc)] = S.Rational(facts[m],facts[na]*facts[nb]*facts[nc])
    return counts


def main():
    certs = [json.loads(line.split(' ',1)[1]) for line in OUTPUT.read_text().splitlines()
             if line.startswith('CHAR_CERTIFICATE ')]
    gate([c['k'] for c in certs]==list(range(1,7)), 'complete six characteristic coefficients')
    for c in certs:
        k=c['k']
        gate(c['degree']==12*k and len(c['coefficients_ascending'])==12*k+1,
             'full certified degree k='+str(k))
        gate(S.Rational(c['content'])>0 and all(v>0 for v in c['coefficients_ascending']),
             'all shifted coefficients strictly positive k='+str(k))
    primitive=0
    # Both candidate and actual characteristic polynomials have degree <=72.
    # Equality at all 73 distinct points proves their full symbolic identities.
    for x in range(1,74):
        g=x+19
        first, second=fibre(g,g),fibre(g,2*g)
        gate(len(first)==7 and len(second)==14, 'complete literal counts x='+str(x))
        p=S.Poly(sum(weight/S.Integer(comb(g,13))*t**((nc-1)//2)
                     for (na,nb,nc),weight in first.items()),t,domain=S.QQ)
        scale=S.Integer(factorial(2*g)//factorial(2*g-26))
        q_times_t=S.Poly(sum(weight/scale*t**((nc-2)//2+1)
                             for (na,nb,nc),weight in second.items()),t,domain=S.QQ)
        gate(p.degree()==6 and p.LC()==1 and p.TC()>0, 'first normalization x='+str(x))
        gate(q_times_t.degree()==13 and q_times_t.TC()>0, 'doubled carry retained x='+str(x))
        inverse=S.invert(S.Poly(t,t,domain=S.QQ),p)
        response=(q_times_t*inverse).rem(p)
        gate((response*S.Poly(t,t,domain=S.QQ)-q_times_t).rem(p).is_zero,
             'original inverse identity x='+str(x))
        cols=[(response*S.Poly(t**j,t,domain=S.QQ)).rem(p) for j in range(6)]
        matrix=S.Matrix(6,6,lambda i,j:cols[j].nth(i))
        coeffs=matrix.charpoly().all_coeffs()[1:]
        for c,actual in zip(certs,coeffs):
            value=S.Rational(c['content'])*sum(S.Integer(v)*(x-1)**j
                      for j,v in enumerate(c['coefficients_ascending']))
            gate(actual==value, 'full-degree interpolation value x='+str(x)+' k='+str(c['k']))
        primitive+=gcd(g,39)==1
    earlier=[m for m in range(1,22) if fibre(21,m)]
    gate(earlier[0]==7, 'gcd hostile actual first return is seven')
    print('UNIVERSE x=1..73, complete literal first/doubled multinomial fibres; primitive='+str(primitive))
    print('INDEPENDENT_PATH exact Poly inversion/remainder and Berkowitz characteristic polynomial')
    print('IDENTITY_CERTIFICATE each degree<=12k<=72; 73 exact equalities prove all six parameter polynomial identities')
    print('COEFFICIENTS all 258 shifted coefficients positive; original carry retained; gcd hostile g21 has first return7')
    print('PRODUCER_OUTPUT_SHA256 '+hashlib.sha256(OUTPUT.read_bytes()).hexdigest())
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__=='__main__':
    main()
