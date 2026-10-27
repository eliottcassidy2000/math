#!/usr/bin/env python3
import argparse
import json
import math

def require(c, msg):
    if not c:
        raise RuntimeError(msg)

def factor(n):
    out={}
    p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1
            n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def sigma(n):
    return math.prod((p**(e+1)-1)//(p-1) for p,e in factor(n).items())

def s(n):return sigma(n)-n

def prime(n):
    return n>=2 and factor(n)=={n:1}

def crt_bank(C):
    Q=1;r=0;used=set();rows=[]
    for a in range(2,C+1):
        S,T=s(a),sigma(a)
        ell=next(p for p in range(3,1000) if prime(p) and p not in used and S*T%p)
        r+=Q*((T-r)*pow(Q,-1,ell)%ell)
        Q*=ell;used.add(ell)
        rows.append({'a':a,'S':S,'T':T,'ell':ell})
    require(math.gcd(r,Q)==1,'nonprimitive CRT target')
    require(all(r%z['ell']==z['T']%z['ell'] for z in rows),'CRT owner lost')
    return Q,r,rows

def beyond_bank(m,rows):
    for z in rows:
        S,T,ell=z['S'],z['T'],z['ell']
        require(m%ell==T%ell,'wrong target progression')
        require(m>S*ell+T,'not beyond exceptional prime context')
        if m>T and (m-T)%S==0:
            p=(m-T)//S
            require(p%ell==0 and p!=ell and not prime(p),'bank obstruction failed')

def witness(m):
    for p in range(3,(m-1)//2):
        q=m-1-p
        if prime(p) and prime(q):
            require(p!=q and s(p*q)==m,'invalid distinct witness')
            return p,q
    raise RuntimeError('finite witness search failed')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--negative-control',action='store_true')
    args=parser.parse_args()
    if args.negative_control:
        require(s(9)==1+3+3,'EQUAL_PRIME_GUARD_FAILED: s(9)=4, not 7')
    rows_out=[]
    for C in [3,5,10]:
        Q,r,rows=crt_bank(C)
        B=max(z['S']*z['ell']+z['T'] for z in rows)
        m=r
        while m<=B or not prime(m):m+=Q
        beyond_bank(m,rows)
        p,q=witness(m)
        rows_out.append({'core_bank':[2,C],'Q':Q,'r':r,'B':B,'rows':rows,
                         'prime_target':m,'distinct_prime_pair':[p,q],'source':p*q})
    Q,r,rows=crt_bank(5)
    ap_targets=[]
    for m in range(r,200001,Q):
        if prime(m):
            beyond_bank(m,rows)
            ap_targets.append({'target':m,'pair':witness(m)})
    m=1063; bound=(m-1)**2
    sums=[0]*(bound+1)
    for d in range(1,bound//2+1):
        for n in range(2*d,bound+1,d):sums[n]+=d
    fiber=[n for n in range(1,bound+1) if sums[n]==m]
    require(11561 in fiber,'explicit source absent')
    require(all(not all(e>=2 for e in factor(n).values()) for n in fiber),'squarefull source found')
    # Independent complete factorization partition: fresh-prime core or squarefull.
    factored=set()
    for a in range(2,m):
        S,T=s(a),sigma(a)
        if m>T and (m-T)%S==0:
            p=(m-T)//S
            if prime(p) and a%p:factored.add(a*p)
    v=1
    while v**3<=bound:
        if all(e==1 for e in factor(v).values()):
            for u in range(1,math.isqrt(bound//v**3)+1):
                n=u*u*v**3
                if s(n)==m:factored.add(n)
        v+=1
    require(sorted(factored)==fiber,'factorized complete inverse disagrees with divisor sieve')
    print(json.dumps({'status':'FINITE-EXACT controls, not asymptotic verification',
                      'CRT_bank_certificates':rows_out,
                      'progression_prime_targets_through_200000':ap_targets,
                      'complete_target1063_fiber':fiber,'complete_source_bound':bound,
                      'no_squarefull_source_for_1063':True,
                      'independent_guarded_core_squarefull_fiber_agrees':True,
                      'large_target_complete_fiber':'NOT CLAIMED'},indent=2,sort_keys=True))

if __name__=='__main__':main()
