"""Exact delayed-activation versions of inherited Burnside response engines.

The response quotient and integer gauges are inherited. The changes are:
keep inactive degree buckets in place, and optionally pack half-deficits
in the fixed-target projective engine. No tournament-class generation.
"""
from math import factorial, gcd
from collections import defaultdict
from baseline105 import exact_shift

def tables(n):
    if type(n) is not int or n < 0:
        raise ValueError('nonnegative integer required')
    twos=[i-i.bit_count() for i in range(n+1)]
    return twos,[factorial(i)>>twos[i] for i in range(n+1)]

def projective(n, park=True, half=True):
    twos,odds=tables(n)
    scale=2 if half else 1
    width=max(1,(n//scale).bit_length());digit=(1<<width)-1
    masks=[(1<<(width*i))-1 for i in range((n+1)//2+1)]
    states=[defaultdict(int) for _ in range(n+1)];states[n][0]=1
    visits=transitions=parked=negative=0;stages=[]
    for p in range(n if n%2 else n-1,2,-2):
        shift=width*((p-3)//2)
        row=sum(((j-gcd(p,j))//scale)<<(width*((j-3)//2)) for j in range(3,p,2))
        stop=min(p,n+1) if park else 0
        nxt=states[:stop]+[defaultdict(int) for _ in range(stop,n+1)]
        parked+=sum(len(states[r]) for r in range(stop))
        for r in range(stop,n+1):
            bucket=states[r]
            if not bucket:continue
            visits+=len(bucket);limit=r//p
            mask=masks[max(0,(min(p-2,r)-1)//2)]
            target=nxt[r]
            for response,value in bucket.items():target[response&mask]+=value
            transitions+=len(bucket)
            for k in range(1,limit+1):
                rr=r-p*k
                factor=odds[n-rr]//(odds[n-r]*p**k*odds[k])
                base=k*(p-1)//2+p*k*(k-1)//2+rr*k-twos[k]
                offset=k*row;mask=masks[max(0,(min(p-2,rr)-1)//2)]
                target=nxt[rr]
                for response,value in bucket.items():
                    exponent=base-scale*k*((response>>shift)&digit)
                    target[(response+offset)&mask]+=exact_shift(value*factor,exponent)
                    negative+=exponent<0
                transitions+=len(bucket)
        states=nxt;stages.append(sum(map(len,states)))
    numerator=sum(exact_shift(sum(bucket.values())*(odds[n]//(odds[n-r]*odds[r])),r*(r-1)//2-twos[r]) for r,bucket in enumerate(states))
    answer,rem=divmod(numerator,odds[n])
    if rem:raise ArithmeticError('nonintegral result')
    return answer,dict(stages=stages,peak=max(stages,default=1),state_visits=visits,
                       transitions=transitions,parked_visits_avoided=parked,
                       negative_shifts=negative,digit_width=width)

def prefix(n, colors=1, park=True):
    """All coefficients 0..n with the inherited raw gcd response."""
    twos,odds=tables(n)
    if type(colors) is not int or colors<1:raise ValueError('positive integer colors')
    width=max(1,n.bit_length());digit=(1<<width)-1
    masks=[(1<<(width*i))-1 for i in range(n//2+2)]
    powers=[colors**i for i in range(n+1)]
    states=[defaultdict(int) for _ in range(n+1)];states[0][0]=1
    visits=transitions=parked=0;stages=[]
    for p in range(n if n%2 else n-1,2,-2):
        shift=width*((p-1)//2)
        row=sum(gcd(p,j)<<(width*i) for i,j in enumerate(range(1,p,2)))
        active=n-p+1 if park else n+1
        nxt=[defaultdict(int) for _ in range(active)]+states[active:]
        parked+=sum(len(states[d]) for d in range(active,n+1))
        for d in range(active):
            bucket=states[d]
            if not bucket:continue
            visits+=len(bucket);limit=(n-d)//p
            mask=masks[max(0,(min(p-2,n-d)+1)//2)]
            target=nxt[d]
            for response,value in bucket.items():target[response&mask]+=value
            transitions+=len(bucket)
            for k in range(1,limit+1):
                dd=d+p*k
                factor=odds[dd]//(odds[d]*p**k*odds[k])*powers[k]
                base=k*(p-1)//2+p*k*(k-1)//2-twos[k]
                offset=k*row;mask=masks[max(0,(min(p-2,n-dd)+1)//2)]
                target=nxt[dd]
                for response,value in bucket.items():
                    exponent=base+k*((response>>shift)&digit)
                    target[(response+offset)&mask]+=(value*factor)<<exponent
                transitions+=len(bucket)
        states=nxt;stages.append(sum(map(len,states)))
    totals=[0]*(n+1)
    for d,bucket in enumerate(states):
        if not bucket:continue
        for k in range(n-d+1):
            factor=odds[d+k]//(odds[d]*odds[k])*powers[k]
            polynomial=sum(value<<(k*response) for response,value in bucket.items())
            totals[d+k]+=(polynomial*factor)<<(k*(k-1)//2-twos[k])
    answer=[]
    for d,total in enumerate(totals):
        value,rem=divmod(total,odds[d])
        if rem:raise ArithmeticError('nonintegral prefix')
        answer.append(value)
    return answer,dict(stages=stages,peak=max(stages,default=1),state_visits=visits,
                       transitions=transitions,parked_visits_avoided=parked)
