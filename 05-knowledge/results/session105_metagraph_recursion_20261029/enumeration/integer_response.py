#!/usr/bin/env python3
"""Integer-scaled projective Burnside response recursion (stdlib only)."""
from math import factorial, gcd
from collections import defaultdict


def exact_shift(value, exponent):
    if exponent >= 0:
        return value << exponent
    quotient, remainder = divmod(value, 1 << -exponent)
    if remainder:
        raise ArithmeticError('nonintegral projective transition')
    return quotient


def count(n, tail=1, cache_profiles=True):
    """Return a_n and stage statistics; analytically finish odd lengths <= tail.

    tail=1 or 3 currently. All cycle types, including composite lengths,
    are retained. Negative binary shifts use exact checked division.
    """
    if type(n) is not int or n < 0 or tail not in (1,3):
        raise ValueError('nonnegative integer n and tail in (1,3) required')
    facts=[factorial(i) for i in range(n+1)]
    lengths=tuple(range(3,n+1,2)) # unit response identically zero
    states={(n,(0,)*len(lengths)):1}
    counts=[]; transitions=0; negative=0; template_count=0
    for p in reversed(lengths):
        if p<=tail: break
        nxt=defaultdict(int)
        templates={}
        for (r,s),value in states.items():
            sp=s[-1] if p<=r else 0
            if not cache_profiles or r not in templates:
                ts=[]
                for k in range(r//p+1):
                    rr=r-k*p
                    base=k*(p-1)//2+p*k*(k-1)//2+rr*k
                    coefficient=facts[r]//(facts[rr]*p**k*facts[k])
                    delta=tuple(k*(gcd(p,j)-j) for j in range(3,min(p-2,rr)+1,2))
                    ts.append((k,rr,base,coefficient,delta))
                templates[r]=ts
                template_count+=len(ts)
            for k,rr,base,coefficient,delta in templates[r]:
                e=base+k*sp
                ss=tuple(s[i]+d for i,d in enumerate(delta))
                nxt[(rr,ss)]+=exact_shift(value*coefficient,e)
                transitions+=1;negative+=e<0
        states=nxt;counts.append(len(states))
    numerator=0
    for (r,s),value in states.items():
        if tail==1 or r<3:
            numerator+=value << (r*(r-1)//2)
        else:
            sp=s[0]
            for k in range(r//3+1):
                rr=r-3*k
                e=k+3*k*(k-1)//2+k*sp+rr*k+rr*(rr-1)//2
                coefficient=facts[r]//(facts[rr]*3**k*facts[k])
                numerator+=exact_shift(value*coefficient,e)
                transitions+=1;negative+=e<0
    answer,remainder=divmod(numerator,facts[n])
    if remainder: raise ArithmeticError('nonintegral Burnside result')
    return answer,dict(stages=counts,peak=max(counts,default=1),total_states=sum(counts),transitions=transitions,negative_shifts=negative,terminal_states=len(states),profile_templates=template_count)


def count_packed(n, odd_local=True):
    """Degree-bucket implementation; pack nonnegative deficits -s_j.

    Descending order implies every active deficit <= n, so width bit_length(n)
    prevents carries. Default scale is odd((N-r)!) times W; odd_local=False
    selects the N!/r! scale used by count().
    """
    if type(n) is not int or n < 0:
        raise ValueError('nonnegative integer n required')
    facts=[factorial(i) for i in range(n+1)]
    twos=[i-i.bit_count() for i in range(n+1)]
    odds=[facts[i]>>twos[i] for i in range(n+1)]
    width=max(1,n.bit_length()); digit=(1<<width)-1
    masks=[(1<<(width*i))-1 for i in range((n+1)//2+1)]
    states=[{} for _ in range(n+1)];states[n]={0:1}
    counts=[];transitions=0;negative=0;templates=0
    for p in range(n if n%2 else n-1,2,-2):
        shift=width*((p-3)//2)
        row=sum((j-gcd(p,j))<<(width*((j-3)//2)) for j in range(3,p,2))
        nxt=[defaultdict(int) for _ in range(n+1)]
        for r,bucket in enumerate(states):
            if not bucket:continue
            limit=r//p
            templates+=limit+1;transitions+=len(bucket)*(limit+1)
            mask=masks[max(0,(min(p-2,r)-1)//2)]
            for response,value in bucket.items():nxt[r][response&mask]+=value
            for k in range(1,limit+1):
                rr=r-p*k
                factor=(odds[n-rr]//(odds[n-r]*p**k*odds[k]) if odd_local else facts[r]//(facts[rr]*p**k*facts[k]))
                base=k*(p-1)//2+p*k*(k-1)//2+rr*k-(twos[k] if odd_local else 0)
                offset=k*row
                mask=masks[max(0,(min(p-2,rr)-1)//2)]
                target=nxt[rr]
                for response,value in bucket.items():
                    exponent=base-k*((response>>shift)&digit)
                    target[(response+offset)&mask]+=exact_shift(value*factor,exponent)
                    negative+=exponent<0
        states=nxt;counts.append(sum(map(len,states)))
    if odd_local:
        numerator=sum(exact_shift(sum(bucket.values())*(odds[n]//(odds[n-r]*odds[r])),r*(r-1)//2-twos[r]) for r,bucket in enumerate(states))
    else:
        numerator=sum(sum(bucket.values()) << (r*(r-1)//2) for r,bucket in enumerate(states))
    answer,remainder=divmod(numerator,odds[n] if odd_local else facts[n])
    if remainder:raise ArithmeticError('nonintegral Burnside result')
    return answer,dict(stages=counts,peak=max(counts,default=1),total_states=sum(counts),transitions=transitions,negative_shifts=negative,profile_templates=templates)
