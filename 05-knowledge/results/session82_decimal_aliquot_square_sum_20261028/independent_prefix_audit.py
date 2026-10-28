#!/usr/bin/env python3
"""Independent exhaustive prefix sieve / monotone polynomial inversion."""
import json, math
from sympy import isprime, primerange
from check_aliquot_fibers import need, fiber, DATA_DIR

def tables(N):
    sig=[0]*(N+1); largest=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d): sig[n]+=d
    for p in range(2,N+1):
        if largest[p]==0:
            for n in range(p,N+1,p): largest[n]=p
    return sig,largest

def inverse(t,sig,largest):
    out=set()
    # Pure prime powers. Exponent 2 is the exceptional unbounded-q branch.
    if isprime(t-1): out.add((t-1)**2)
    for q in primerange(2,math.isqrt(t)+1):
        power=q*q; geom=1+q
        while geom<t:
            geom+=power; power*=q
            if geom==t: out.add(power)
    for m in range(2,t):
        S=sig[m]; sm=S-m
        if S>=t: continue
        q,rem=divmod(t-S,sm)
        if not rem and q>largest[m] and isprime(q): out.add(m*q)
        lo=largest[m]+1
        def poly(q,a):
            return sm*q**a+S*((q**a-1)//(q-1))
        a=2
        while poly(lo,a)<=t:
            L=lo; R=math.isqrt(t//sm)+1
            while L<R:
                mid=(L+R)//2
                if poly(mid,a)<t: L=mid+1
                else: R=mid
            if poly(L,a)==t and isprime(L): out.add(m*L**a)
            a+=1
    return sorted(out)

if __name__=='__main__':
    with (DATA_DIR / 'results.json').open() as h: data=json.load(h)
    sig,largest=tables(1000999)
    for t in range(2,301): need(inverse(t,sig,largest)==fiber(t)[0],f'small mismatch {t}')
    for t,d in data.items():
        got=inverse(int(t),sig,largest)
        need(got==d['fiber'],f'large mismatch {t}')
        print('PASS independent complete prefix audit',t,len(got),flush=True)
    with (DATA_DIR / 'independent_audit.json').open('w') as f:
        json.dump({'status':'PASS','small_targets':[2,300],'large_target_counts':{t:d['count'] for t,d in data.items()},'method':'independent divisor-sum prefix sieve + monotone polynomial binary search'},f,indent=2,sort_keys=True)
