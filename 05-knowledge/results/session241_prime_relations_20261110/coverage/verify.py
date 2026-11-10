#!/usr/bin/env python3
"""Exact arithmetic and literal tournament checks; no external dependencies."""
import json, math
from pathlib import Path

def need(x):
    if not x: raise RuntimeError('verification failed')
def transition(t, c):
    D,E,U,d=t
    return (D,0,2*U,d) if c=='O' else (D-E+U,d,2*U-d,d)
def state(a,b):
    t=(1,0,2,1)
    for c in 'O'*a+'Z'+'O'*b+'Z': t=transition(t,c)
    return t

def value(r,s): return ((2**r-1)*(2**s+1)+3)//2

def inverse(N):
    if N<1 or N%2==0: return None
    M=2*N-3
    if M<1: return None
    w=bin(M)[2:]
    if '0' not in w:
        if len(w)%2: return None
        r=s=len(w)//2
    else:
        leading=len(w)-len(w.lstrip('1'))
        if leading>1:
            r=leading; s=len(w)-r
        else:
            s=len(w)-len(w.rstrip('1')); r=len(w)-s-1
    if r>=3 and s>=2 and value(r,s)==N: return (r-2,s-1)
    return None

def graph(a,b):
    # Root 0 beats initial spine vertex and every O; every Z beats root.
    letters='O'+'O'*a+'Z'+'O'*b+'Z'
    n=1+len(letters); out=[0]*n
    for j,c in enumerate(letters,1):
        if c=='O': out[0]|=1<<j
        else: out[j]|=1
        for i in range(1,j): out[i]|=1<<j
    return out

def hp(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=dp[mask][v]
            if c:
                avail=out[v]&~mask
                while avail:
                    bit=avail&-avail; avail-=bit
                    dp[mask|bit][bit.bit_length()-1]+=c
    return sum(dp[-1])

def prime(p): return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))
def order(p):
    t=1; x=2%p
    while x!=1: t+=1; x=x*2%p
    return t

def valuation(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v

def main():
    seen={}; checks=0
    for r in range(3,103):
        for s in range(2,102):
            N=value(r,s); need(state(r-2,s-1)[0]==N)
            need(inverse(N)==(r-2,s-1)); need(N not in seen); seen[N]=[r,s]
            g=math.gcd(r,s)
            need(math.gcd(2**r-1,2**s+1)==(2**g+1 if (r//g)%2==0 else 1))
            checks+=1
    literal=[]
    for a in range(1,5):
        for b in range(1,5):
            n=hp(graph(a,b)); need(n==state(a,b)[0]); literal.append([a,b,n])
    bound=1000000
    forward={value(r,s):(r-2,s-1) for r in range(3,24) for s in range(2,24) if value(r,s)<=bound}
    recovered={N:inverse(N) for N in range(1,bound+1,2) if inverse(N) is not None}
    need(forward==recovered)
    periods=0
    for p in range(3,200):
        if not prime(p): continue
        t=order(p)
        for r in range(3,43):
            for s in range(2,42):
                M=2*value(r,s)-3
                predicate=(r%t==0 or (t%2==0 and s%t==t//2))
                need((M%p==0)==predicate)
                v=valuation(2**t-1,p)+valuation(r//t,p) if r%t==0 else 0
                if t%2==0 and s%t==t//2:
                    v+=valuation(2**(t//2)+1,p)+valuation(2*s//t,p)
                need(valuation(M,p)==v); periods+=1
    seeds={(1,0,2,1)}
    for _ in range(5):
        more=set()
        for D,E,U,d in seeds:
            more.update([transition((D,E,U,d),'O'),transition((D,E,U,d),'Z'),(D-E+U,D-E,D+U,D)])
        seeds |= more
    seed_checks=0
    for seed in seeds:
        D,E,U,d=seed
        for a in range(1,7):
            for b in range(1,7):
                q=seed
                for c in 'O'*a+'Z'+'O'*b+'Z': q=transition(q,c)
                need(2*(q[0]-D)-d==(2**(a+1)*U-d)*(2**(b+1)+1))
                seed_checks+=1
    # Hostiles: scalar prime/composite cannot be read off transformed product.
    need(prime(value(3,2)) and not prime(value(3,3)))
    need(inverse(7) is None and inverse(21) is None)
    counts={str(2**k):sum(value(r,s)<=2**k for r in range(3,k+3) for s in range(2,k+3)) for k in [10,20,40,80,160]}
    for k in range(4,161):
        actual=sum(value(r,s)<=2**k for r in range(3,k+3) for s in range(2,k+3))
        need(actual==(k-4)*(k-3)//2+k//2-2)
    output={'status':'PASS','arithmetic_pairs':checks,'literal_tournaments':literal,'inverse_odd_targets_checked':bound//2,'image_below_bound':len(forward),'odd_prime_period_and_valuation_checks':periods,'fixed_seed_identity_checks':seed_checks,'dyadic_formula_bounds_checked':157,'counts':counts,'hostiles':{'missing':[7,21],'prime_H_composite_shift':[19,35],'composite_H_composite_shift':[33,63]},'scope':'Restricted two-burst family only; not a proof of Chan-Pak.'}
    Path(__file__).with_name('verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
if __name__=='__main__': main()
