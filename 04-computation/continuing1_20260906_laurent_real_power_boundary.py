"""Exact half-integer boundary for the frozen actual Laurent separator.

Uses rational arithmetic and integer square roots only. The positive-branch
repair has algebraic coefficients; no rational or all-six-root nonnegative
representation is claimed.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(label)


def mul(A,B):
    vals=[a*b for a in A for b in B]
    return min(vals),max(vals)


def interval(row,I):
    out=(F(0),F(0))
    for c in reversed(row):
        v=mul(out,I)
        out=(v[0]+c,v[1]+c)
    return out


def divide(A,B):
    need(not B[0]<=0<=B[1],'Lagrange denominator excludes zero')
    return mul(A,(1/B[1],1/B[0]))


def sqrt_interval(I):
    scale=10**12
    lower=isqrt(I[0].numerator*scale*scale//I[0].denominator)
    upper=isqrt(I[1].numerator*scale*scale//I[1].denominator)+1
    lo,hi=F(lower,scale),F(upper,scale)
    need(lo>0 and lo*lo<=I[0] and hi*hi>=I[1],'literal rational square-root enclosure')
    return lo,hi


def iadd(*rows):
    return sum(row[0] for row in rows),sum(row[1] for row in rows)


def ineg(row):
    return -row[1],-row[0]


def padd(*rows):
    out={}
    for row in rows:
        for powers,c in row.items():
            out[powers]=out.get(powers,F(0))+c
    return {powers:c for powers,c in out.items() if c}


def pmul(*rows):
    out={(0,0,0,0):F(1)}
    for row in rows:
        terms=[]
        for p,a in out.items():
            for q,b in row.items():
                terms.append({tuple(i+j for i,j in zip(p,q)):a*b})
        out=padd(*terms)
    return out


def main():
    name='continuing1_20260906_laurent_cone_separator_certificate.json'
    folder=Path(__file__).resolve().parent
    candidates=(folder/name,folder.parent/'05-knowledge'/'results'/name)
    path=next((p for p in candidates if p.is_file()),None)
    need(path is not None,'frozen separator certificate located')
    raw=path.read_bytes()
    need(sha256(raw).hexdigest()=='4b1ee5770b484e4164e692fbf2934f4099800b0b85d379ce01d8afef71040cc0','frozen actual certificate hash')
    bank=json.loads(raw)['cubic']
    p=list(map(F,bank['P']))
    N=list(map(F,bank['N']))
    rows=[list(map(F,row)) for row in bank['residues']]
    Is=[tuple(map(F,I)) for I in bank['intervals']]
    derivative=[i*p[i] for i in range(1,len(p))]
    for I in Is:
        need(interval(p,(I[0],I[0]))[0]*interval(p,(I[1],I[1]))[0]<0,'original polynomial root bracket')
        need(interval(N,I)[1]<0,'retained Lagrange numerator sign')
    need(Is[0][1]<Is[1][0] and Is[1][1]<Is[2][0],'three disjoint positive roots exhaust cubic')
    square_roots=[sqrt_interval(I) for I in Is]
    for r in (0,6):
        total=[F(0),F(0)]
        for I,rootI in zip(Is,square_roots):
            need(interval(rows[r],I)[1]<0,'original window remains negative at the selected roots')
            term=mul(divide(mul(interval(rows[r],I),interval(N,I)),interval(derivative,I)),rootI)
            total=[total[0]+term[0],total[1]+term[1]]
        need(total[1]<0,'strict half-power crossing of the fixed separator')
        lower=total[0].numerator//total[0].denominator
        upper=-((-total[1].numerator)//total[1].denominator)
        print('HALF_POWER r=',r,'certified_integer_enclosure=',lower,upper)
    for i,I in enumerate(square_roots):
        print('SQRT_ROOT',i+1,','.join(map(str,I)))
    a,b,c=map(F,bank['signed_multiplier'])
    E1=iadd(*square_roots)
    E2=iadd(mul(square_roots[0],square_roots[1]),mul(square_roots[0],square_roots[2]),mul(square_roots[1],square_roots[2]))
    E3=mul(mul(square_roots[0],square_roots[1]),square_roots[2])
    coefficients=[iadd((a,a),mul((c,c),mul(E1,E3))),
                  mul((c,c),iadd(E3,ineg(mul(E1,E2)))),
                  iadd((b,b),mul((c,c),iadd(mul(E1,E1),ineg(E2))))]
    for j,I in enumerate(coefficients):
        need(I[0]>0,'strict positive algebraic amplitude coefficient')
        lower=I[0]*10**6;upper=I[1]*10**6
        print('AMPLITUDE_COEFFICIENT',j,'millionth_enclosure=',lower.numerator//lower.denominator,-((-upper.numerator)//upper.denominator))
    # Formal polynomial identity in independent symbols E1,E2,E3,z.
    const=lambda value:{(0,0,0,0):F(value)}
    variables=[{tuple(int(i==j) for i in range(4)):F(1)} for j in range(4)]
    e1,e2,e3,z=variables
    qplus=padd(pmul(z,z,z),pmul(const(-1),e1,z,z),pmul(e2,z),pmul(const(-1),e3))
    A=padd(const(a),pmul(const(c),e1,e3))
    B=pmul(const(c),padd(e3,pmul(const(-1),e1,e2)))
    C=padd(const(b),pmul(const(c),padd(pmul(e1,e1),pmul(const(-1),e2))))
    original=padd(const(a),pmul(const(b),z,z),pmul(const(c),z,z,z,z))
    repaired=padd(A,pmul(B,z),pmul(C,z,z),pmul(const(c),padd(z,e1),qplus))
    need(original==repaired,'formal positive-branch cubic reduction identity')
    for candidate in (F(1),F(-1),F(1,3),F(-1,3)):
        need(interval(p,(candidate,candidate))[0]!=0,'rational-root test for the first cubic')
    need(isqrt(3)**2!=3,'norm one-third is not a rational square')
    print('CONCLUSION the fixed separator is not nonnegative for all real exponents; theta=1/2 already fails')
    print('REPAIR T=(A+B sqrt(s)+C s)J0 on the three positive roots, with explicitly positive real algebraic A,B,C')
    print('SCOPE no nonnegative half-Laurent identity on all six roots; rational-coefficient positive-branch repair is also excluded; no all-parameter closure')
    print('PASS',GATES,'always-active gates; no floating-point roots or logarithms')


if __name__=='__main__':
    main()
