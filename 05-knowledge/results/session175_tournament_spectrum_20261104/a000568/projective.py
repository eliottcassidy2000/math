#!/usr/bin/env python3
"""Fixed-target projective future-response Burnside DP, standard library only."""
from math import factorial,gcd
from collections import defaultdict

def need(b,s):
    if not b: raise ValueError(s)

def odd(x): return x >> ((x & -x).bit_length()-1)

def exact_shift(x,e):
    if e>=0:return x<<e
    need((x & ((1<<(-e))-1)) == 0, 'negative shift is exact')
    return x>>(-e)

def count(n,projective=True):
    need(type(n) is int and n>=0,'nonnegative integer order required')
    # Store Q_j=R_j-jR_1, for unprocessed j>=3; degree fixes M=n-d.
    # Weight is odd(d!)*sum w(mu)*2^(M R_1(mu)).
    # Ungauged control retains R_1 and R_j and the ordinary local weight.
    lengths=list(range(3,n+1,2))[::-1]
    of=[odd(factorial(d)) for d in range(n+1)]
    initial=(0,)*len(lengths) if projective else (0,)*(len(lengths)+1)
    states={(0,initial):1};rows=[];negative=0;cached=0
    for idx,ell in enumerate(lengths):
        future=lengths[idx+1:]; new=defaultdict(int);trans=0;cache={}
        for (d,profile),w in states.items():
            M=n-d
            if M<ell:
                # Current coordinate has already been discarded by degree.
                new[d,profile]+=w;trans+=1;continue
            c=profile[0];rest=profile[1:]
            if d not in cache:
                transfers=[]
                for k in range(M//ell+1):
                    q=d+k*ell;v=(factorial(k)&-factorial(k)).bit_length()-1
                    den=of[d]*ell**k*odd(factorial(k))
                    need(of[q]%den==0,'integral odd transfer')
                    exponent=k*(ell-1)//2+ell*k*(k-1)//2-v
                    if projective:exponent+=k*(M-k*ell)
                    transfers.append((k,q,of[q]//den,exponent))
                cache[d]=transfers;cached+=len(transfers)
            admissible=[j for j in future if j<=M]
            for k,q,t,e in cache[d]:
                if projective:
                    p=tuple(r+k*(gcd(ell,j)-j) for j,r in zip(admissible,rest) if j<=n-q)
                else:
                    p=tuple(r+k*gcd(ell,j) for j,r in zip(admissible,rest[:-1]) if j<=n-q)+(rest[-1]+k,)
                shift=e+k*c
                negative+=shift<0
                new[q,p]+=exact_shift(w*t,shift);trans+=1
        states=dict(new)
        rows.append({'length':ell,'states':len(states),'transitions':trans})
    # Complete every core with its uniquely prescribed singleton tail.
    # odd(n!) is a common odd denominator; v2(r!) is removed by the tail exponent.
    total=0
    for (d,p),w in states.items():
        r=n-d;v=(factorial(r)&-factorial(r)).bit_length()-1
        den=of[d]*odd(factorial(r));need(of[n]%den==0,'tail integral transfer')
        e=r*(r-1)//2-v
        if projective:need(not p,'empty projective profile')
        else:e+=r*p[0]
        total+=exact_shift(w*(of[n]//den),e)
    need(total%of[n]==0,'final integer count')
    return total//of[n],{'peak_states':max([1]+[r['states'] for r in rows]),'transitions':sum(r['transitions'] for r in rows),'negative_shifts':negative,'precomputed_transfers':cached,'rows':rows}

def packed_count(n,projective=True):
    need(type(n) is int and n>=0,'nonnegative integer order required')
    need(type(projective) is bool,'Boolean projective mode required')
    start=3 if projective else 1
    width=max(1,n.bit_length());digit=(1<<width)-1
    facts=[factorial(d) for d in range(n+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,twos)]
    masks=[(1<<(width*c))-1 for c in range(n//2+2)]
    states={(0,0):1};rows=[];negative=0;cached=0
    for ell in range(n if n%2 else n-1,2,-2):
        idx=(ell-start)//2
        row=sum(((j-gcd(ell,j)) if projective else gcd(ell,j))<<(width*i) for i,j in enumerate(range(start,ell,2)))
        tables={}
        for degree,_ in states:
            if degree in tables:continue
            M=n-degree;entries=[]
            for k in range(M//ell+1):
                nd=degree+ell*k
                factor,remainder=divmod(odds[nd],odds[degree]*ell**k*odds[k]);need(not remainder,'integral transfer')
                base=k*(ell-1)//2+ell*k*(k-1)//2-twos[k]
                if projective:base+=k*(M-k*ell)
                c=max(0,(min(ell-2,n-nd)-start)//2+1)
                entries.append((k,nd,base,factor,masks[c],k*row));cached+=1
            tables[degree]=entries
        new=defaultdict(int);trans=0
        for (degree,response),w in states.items():
            charge=(response>>(width*idx))&digit
            if projective:charge=-charge
            for k,nd,base,factor,mask,offset in tables[degree]:
                e=base+k*charge;negative+=e<0
                new[nd,(response+offset)&mask]+=(exact_shift(w,e)*factor if e<0 else (w*factor)<<e);trans+=1
        states=dict(new)
        rows.append({'length':ell,'states':len(states),'transitions':trans})
    total=0
    for (d,p),w in states.items():
        r=n-d;factor,rem=divmod(odds[n],odds[d]*odds[r]);need(not rem,'tail integral transfer')
        e=r*(r-1)//2-twos[r]
        if projective:need(not p,'empty projective profile')
        else:e+=r*p
        total+=exact_shift(w*factor,e)
    need(total%odds[n]==0,'final integer count')
    return total//odds[n],{'peak_states':max([1]+[r['states'] for r in rows]),'transitions':sum(r['transitions'] for r in rows),'negative_shifts':negative,'precomputed_transfers':cached,'rows':rows}
