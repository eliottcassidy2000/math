#!/usr/bin/env python3
"""Divisor-triangular prime-cycle derivatives versus full plethystic recovery.
Full-index reference follows Session86 A000568_COMPOSITION.md, not a new formula.
"""
from collections import Counter
from fractions import Fraction as Q
from math import factorial,gcd
import json
from check_marked_components import oddparts,channel,mul,inv,add,require

def exponent(p):
    return sum((d-1)//2 for d in p)+sum(gcd(d,e) for i,d in enumerate(p) for e in p[:i])
def zpart(p):
    z=1
    for d,m in Counter(p).items():z*=d**m*factorial(m)
    return z

def bjet(N,d):
    return [sum((Q(2**(exponent(p)+(d-1)//2+sum(gcd(d,e) for e in p)),d*zpart(p)) for p in oddparts(n)),Q()) for n in range(N+1)]
def psi(a,d):
    return [a[n//d] if n%d==0 else 0 for n in range(len(a))]
def divisor_jets(N,D):
    A=channel(N,1,0); ai2=mul(inv(A),inv(A))
    B={d:bjet(N,d) for d in range(1,D+1,2)}; H={}
    for d in B:
        value=mul(ai2,B[d]);value[0]-=int(d==1)
        for i in H:
            if d%i==0:value=add(value,mul(H[i],psi(B[d//i],i)),-1)
        H[d]=mul(value,inv(psi(B[1],d)))
    return H

def padd(a,b,scale=1):
    out=a.copy()
    for p,c in b.items():
        out[p]=out.get(p,Q())+scale*c
        if not out[p]:del out[p]
    return out

def pmul(a,b,N):
    out={}
    for p,c in a.items():
        for q,d in b.items():
            if sum(p)+sum(q)<=N:
                key=tuple(sorted(p+q));out[key]=out.get(key,Q())+c*d
    return {p:c for p,c in out.items() if c}

def prime_reference(N):
    B={p:Q(2**exponent(p),zpart(p)) for n in range(N+1) for p in oddparts(n)}
    F=padd(B,{():1},-1); bi={():Q(1)}
    for n in range(1,N+1):
        bi.update({p:-c for p,c in pmul(F,bi,n).items() if sum(p)==n})
    K=padd(padd({():Q(1)},bi,-1),{(1,):1},-1);P={}
    cache={():{():Q(1)}}
    def sub(p):
        if p not in cache:
            d=p[-1]; shifted={tuple(i*d for i in q):c for q,c in F.items() if sum(q)*d<=N}
            cache[p]=pmul(sub(p[:-1]),shifted,N)
        return cache[p]
    for n in range(3,N+1):
        row={p:c for p,c in K.items() if sum(p)==n};P.update(row)
        for p,c in row.items():K=padd(K,sub(p),-c)
    require(not K,'prime reference residual')
    return P

def reference_eval(P,N,d):
    F=channel(N,1,0);F[0]=0;out=[Q()]*(N+1)
    for p,c in P.items():
        m=p.count(d)
        if not m:continue
        q=list(p);q.remove(d);term=[c*m]+[0]*N
        for e in q:term=mul(term,psi(F,e))
        out=add(out,term)
    return out

def main():
    # Reference total prime degree must include derivative's lost weight d.
    bound=15;D=11;H=divisor_jets(bound,D);P=prime_reference(bound);checks=0
    for d,v in H.items():
        expected=reference_eval(P,bound,d)
        for n in range(bound-d+1):
            require(v[n]==expected[n],f'divisor jet mismatch {d,n}');checks+=1
    require(H[3][0]==Q(2,3),'C3 nonidentity cycle weight')
    # The H1 jet is the same child-slot derivative extracted by the Bell compiler.
    from check_prime_slot_jets import jets
    _,h=jets(bound,1)
    require(H[1]==h[1],'p1 Bell/divisor mismatch')
    print(json.dumps({'full_index_bound':bound,'odd_cycles_through':D,'coefficient_checks':checks,
        'derivative_scope':'rational aggregate cycle weight, not class count',
        'H3_prefix':[str(v) for v in H[3][:8]],'H11_constant':str(H[11][0])},indent=2))
if __name__=='__main__':main()
