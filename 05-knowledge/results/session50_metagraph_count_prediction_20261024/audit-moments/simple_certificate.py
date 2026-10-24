#!/usr/bin/env python3
"""Self-contained short certificate: eight moments only, one fixed upper polynomial."""
import json,sys
from fractions import Fraction as F

def add(a,b):return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def certify(K):
    if not isinstance(K,list) or len(K)!=8 or not all(type(x) is int and x>0 for x in K):
        raise ValueError('exactly eight positive integer word moments required')
    q=[F(465000,50000000),F(-1250,50000000),F(1,50000000)]
    P=add(add([1],mul([2,-2],q)),mul([0,-1,1],mul(q,q)))
    err=add([1],[0]+[-x for x in q])
    if add([-1],[0]+P)!=mul([-1,1],mul(err,err)):
        raise ArithmeticError('global upper square identity failed')
    U=sum(c*x for c,x in zip(P,K))
    Q=[1]
    for j in range(len(K)):Q=mul(Q,[-(2*j+1),1])
    exponent=2*len(K)-len(K).bit_count();M=1<<exponent
    residue=(-sum(Q[j+1]*K[j] for j in range(len(K)))*pow(Q[0],-1,M))%M
    first=1+(residue-1)%M;last=U.numerator//U.denominator
    count=max(0,1+(last-first)//M)
    return {'upper':str(U),'upper_polynomial_coefficients':[str(v) for v in P],'global_square_identity_verified':True,'modulus':M,'residue':residue,'integer_interval':[1,last],'candidate_count':count,'prediction':first if count==1 else None}
if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('usage: simple_certificate.py EIGHT_MOMENTS.json')
    print(json.dumps(certify(json.load(open(sys.argv[1]))),indent=2,sort_keys=True))
