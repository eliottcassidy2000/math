#!/usr/bin/env python3
from math import isqrt
import hashlib,json


def need(x,msg):
    if not x:
        raise RuntimeError(msg)


def primal(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def main():
    records=[]; pairs=states=0
    for m in range(2,513):
        a,b=0,1
        values=[]
        while True:
            values.append(a)
            a,b=b,(a+b)%m
            if (a,b)==(0,1):
                break
            need(len(values)<=m*m,'pair period finite bound')
        P=len(values)
        need(P>=3 and values[1]==1 and values[P-2]==m-1,'positive Fibonacci +/-1')
        inv=lambda x,y:((y-x)%m,x)
        need(inv(0,1)==(1,0),'first inverse step')
        need(inv(1,0)==(m-1,1),'second inverse step')
        S=set(values)-{0}
        need(S&{(-x)%m for x in S},'Fibonacci not oriented')
        records.append((m,P,P-2))
        states+=P
        if m<=32:
            for a in range(m):
                for b in range(m):
                    need(inv(b,(a+b)%m)==(a,b),'pair inverse on full state square')
                    pairs+=1
    bound=10**6
    sieve=bytearray(b'\1')*(bound+1)
    sieve[0:2]=b'\0\0'
    for p in range(2,isqrt(bound)+1):
        if sieve[p]:
            sieve[p*p:bound+1:p]=b'\0'*((bound-p*p)//p+1)
    primes=[p for p in range(bound+1) if sieve[p]]
    witness=[]
    for m in range(2,513):
        pos=neg=None
        for p in primes:
            if pos is None and p%m==1:
                pos=p
            if neg is None and p%m==m-1:
                neg=p
            if pos is not None and neg is not None:
                break
        need(pos is not None and neg is not None,'finite prime search found both residues')
        need(primal(pos) and primal(neg),'independent primality of witnesses')
        need(pos%m==1 and neg%m==m-1,'prime witness residues')
        witness.append((m,pos,neg))
    summary={'status':'PASS','nontrivial_moduli':511,'fibonacci_pair_states':states,
             'full_pair_inverse_checks':pairs,'prime_sieve_bound':bound,
             'largest_selected_prime':max(max(p,n) for m,p,n in witness),
             'fibonacci_records_sha256':hashlib.sha256(json.dumps(records).encode()).hexdigest(),
             'prime_witnesses_sha256':hashlib.sha256(json.dumps(witness).encode()).hexdigest(),
             'examples':{'fibonacci':[r for r in records if r[0] in [2,3,5,7,21,512]],
                         'primes':[r for r in witness if r[0] in [2,3,5,7,21,512]]},
             'm1':'excluded; empty connection set is a vacuous one-vertex tournament'}
    print(json.dumps(summary,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
