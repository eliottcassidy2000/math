"""Session87 degree buckets, with FLINT weights and an optional batched tail.

The quotient, degree-local normalization, loop order and folded terminal
observer are inherited verbatim mathematically; only arithmetic/layout vary.
"""
from collections import defaultdict
from math import factorial,gcd
from time import perf_counter
from flint import fmpz,fmpz_poly
def prefix(N,tail='fold',backend='flint'):
    if type(N) is not int or N<0:raise ValueError('nonnegative integer order')
    Z=fmpz if backend=='flint' else int
    width=max(1,N.bit_length());digit=(1<<width)-1;masks=[(1<<(width*k))-1 for k in range(N//2+2)]
    twos=[k-k.bit_count() for k in range(N+1)];odd=[factorial(k)>>twos[k] for k in range(N+1)]
    states=[{} for _ in range(N+1)];states[0][0]=Z(1);rows=[];start=perf_counter()
    for ell in reversed(range(3,N+1,2)):
        shift=width*((ell-1)//2);row=sum(gcd(ell,j)<<(width*(j//2)) for j in range(1,ell,2));nxt=[defaultdict(Z) for _ in range(N+1)];transitions=0
        for d,bucket in enumerate(states):
            if not bucket:continue
            limit=(N-d)//ell;transitions+=len(bucket)*(limit+1);mask=masks[(min(ell-2,N-d)+1)//2];target=nxt[d]
            for r,w in bucket.items():target[r&mask]+=w
            for k in range(1,limit+1):
                nd=d+ell*k;factor,rem=divmod(odd[nd],odd[d]*ell**k*odd[k])
                if rem:raise ArithmeticError('degree transfer')
                factor=Z(factor);base=k*(ell-1)//2+ell*k*(k-1)//2-twos[k];mask=masks[(min(ell-2,N-nd)+1)//2];offset=k*row;target=nxt[nd]
                for r,w in bucket.items():target[(r+offset)&mask]+=(w*factor)<<(base+k*((r>>shift)&digit))
        states=nxt;rows.append([ell,sum(map(len,states)),transitions])
    core_seconds=perf_counter()-start;start=perf_counter();totals=[Z(0)]*(N+1)
    if tail=='fold':
        for d,bucket in enumerate(states):
            if not bucket:continue
            for k in range(N-d+1):
                factor,rem=divmod(odd[d+k],odd[d]*odd[k])
                if rem:raise ArithmeticError('tail transfer')
                polynomial=sum((w<<(k*r) for r,w in bucket.items()),Z(0));totals[d+k]+=(polynomial*factor)<<(k*(k-1)//2-twos[k])
        scales=odd
    elif tail=='convolution':
        D=fmpz(odd[N]);grouped={}
        for d,bucket in enumerate(states):
            factor=D//odd[d]
            for r,w in bucket.items():
                if r not in grouped:grouped[r]=[fmpz(0)]*(N+1)
                grouped[r][d]+=w*factor
        base=[(D<<(k*(k-1)//2-twos[k]))//odd[k] for k in range(N+1)]
        for r,U in grouped.items():
            low=next(j for j,x in enumerate(U) if x);L=N-low;V=fmpz_poly([base[k]<<(k*r) for k in range(L+1)]);product=fmpz_poly(U[low:]).mul_low(V,L+1)
            for j in range(L+1):totals[low+j]+=product[j]
        scales=[D*D]*(N+1)
    else:raise ValueError('tail mode')
    values=[]
    for a,s in zip(totals,scales):
        q,r=divmod(a,s)
        if r:raise ArithmeticError('integral final count')
        values.append(int(q))
    return values,dict(N=N,tail=tail,backend=backend,rows=rows,core_seconds=core_seconds,tail_seconds=perf_counter()-start)
