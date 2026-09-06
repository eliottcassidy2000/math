#!/usr/bin/env python3
"""Independent h=5 characteristic certificate by rational interpolation.

No producer imports. Actual first and doubled multinomial count fibres at
x=1,...,51; characteristic coefficients by principal minors; degree-bounded
Newton interpolation into monomial coefficients. Standard library only.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import factorial, comb, lcm
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0


def require(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def raw_row(g,mass):
    result={}
    for gamma in range(mass+1):
        for beta in range(mass-gamma+1):
            alpha=mass-beta-gamma
            if 2*g*beta+3*g*gamma != 33*mass:
                continue
            shift=(gamma-(1 if mass==g else 2))//2
            require(gamma-(1 if mass==g else 2)==2*shift,'literal phase exponent')
            result[shift]=F(factorial(mass),factorial(alpha)*factorial(beta)*factorial(gamma))
    return result


def multiply(a,b):
    result=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[i+j]+=x*y
    return result


def reduce_poly(a,p):
    a=a[:]
    while len(a)>5:
        leading=a.pop()
        offset=len(a)-5
        for j in range(5):
            a[offset+j]-=leading*p[j]
    return a+[F(0)]*(5-len(a))


def determinant(a):
    if not a:
        return F(1)
    a=[row[:] for row in a]
    answer=F(1)
    for j in range(len(a)):
        i=next((i for i in range(j,len(a)) if a[i][j]),None)
        if i is None:
            return F(0)
        if i!=j:
            a[j],a[i]=a[i],a[j]
            answer=-answer
        pivot=a[j][j]
        answer*=pivot
        for i in range(j+1,len(a)):
            ratio=a[i][j]/pivot
            for k in range(j+1,len(a)):
                a[i][k]-=ratio*a[j][k]
    return answer


def characteristic(a):
    answer=[F(1)]
    for k in range(1,6):
        total=F(0)
        for selected in combinations(range(5),k):
            total+=determinant([[a[i][j] for j in selected] for i in selected])
        answer.append((-1)**k*total)
    return answer


def numeric_coefficients(x):
    g=x+16
    first=raw_row(g,g)
    second=raw_row(g,2*g)
    require(tuple(sorted(first))==tuple(range(6)),'six complete first channels')
    require(tuple(sorted(second))==tuple(range(-1,11)),'twelve complete doubled channels')
    first_scale=F(comb(g,11))
    second_scale=F(factorial(2*g),factorial(2*g-22))
    p=[first[j]/first_scale for j in range(6)]
    q={j:second[j]/second_scale for j in second}
    require(p[5]==1 and p[0]>0,'monic admissible first row')
    inverse_tau=[-p[j+1]/p[0] for j in range(5)]
    require(reduce_poly([0]+inverse_tau,p)==[1,0,0,0,0],'literal inverse of tau')
    remainder=reduce_poly([q[j] for j in range(11)],p)
    remainder=[remainder[j]+q[-1]*inverse_tau[j] for j in range(5)]
    columns=[reduce_poly([F(0)]*j+remainder,p) for j in range(5)]
    matrix=[list(row) for row in zip(*columns)]
    coefficients=characteristic(matrix)
    require(all(c>0 for c in coefficients),'numeric characteristic positive')
    return coefficients


def interpolate_consecutive(values):
    # Newton forward differences on y=0,1,..., converted to powers of y.
    differences=values[:]
    degree=len(values)-1
    result=[F(0)]*(degree+1)
    basis=[F(1)]
    for r in range(degree+1):
        for j,c in enumerate(basis):
            result[j]+=differences[0]*c
        differences=[differences[j+1]-differences[j] for j in range(len(differences)-1)]
        basis=[c/F(r+1) for c in multiply(basis,[-F(r),F(1)])]
    return result


def evaluate(coefficients,y):
    value=F(0)
    for c in reversed(coefficients):
        value=value*y+c
    return value


def main():
    samples=[numeric_coefficients(x) for x in range(1,52)]
    certificates={}
    for k in range(1,6):
        degree=10*k
        values=[sample[k] for sample in samples]
        coefficients=interpolate_consecutive(values[:degree+1])
        require(len(coefficients)==degree+1 and coefficients[-1]>0,'exact characteristic degree')
        require(all(c>0 for c in coefficients),'all shifted polynomial coefficients strictly positive')
        for y,value in enumerate(values):
            require(evaluate(coefficients,y)==value,'degree bounded polynomial independent samples')
        denominator=lcm(*(c.denominator for c in coefficients))
        numerators=[int(c*denominator) for c in coefficients]
        certificates[str(k)]={'degree':degree,'shift_x':1,'denominator':str(denominator),
                              'numerators_ascending':[str(c) for c in numerators]}
    raw=json.dumps(certificates,sort_keys=True,indent=2)+'\n'
    path=Path(__file__).resolve().parents[1]/'05-knowledge/results/second_20260906_laurent_audit.json'
    path.write_text(raw)
    print('INDEPENDENT h5 endpoint33 audit: literal multinomial fibres, quotient multiplication, principal-minor characteristic, Newton interpolation.')
    print('Universe: x=1,...,51; first masses g=x+16; doubled masses2g; all count triples retained.')
    print('First channels:6; doubled channels:12 including exponent-1.')
    print('Characteristic degree bounds independently justified:10,20,30,40,50.')
    print('All155 monomial coefficients after x=y+1 are strictly positive.')
    print('certificate_sha256:',sha256(raw.encode()).hexdigest())
    print('Exact gates:',GATES)
    print('PASS')


if __name__=='__main__':
    main()
