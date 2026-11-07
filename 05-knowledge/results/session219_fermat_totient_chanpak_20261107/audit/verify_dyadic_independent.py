#!/usr/bin/env python3
"""Independent exact phi-preimage check: support enumeration and bounded sieve."""
import itertools, json, math
from pathlib import Path

def prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n)+1))

def pure_fiber(a):
    # If odd p|n and phi(n)=2^a, p=1+2^j, j<=a.
    # Primality forces j a power of two. At a<=63, only F0,...,F5 occur.
    assert 0 <= a <= 63
    ps=[3,5,17,257,65537]
    out=[]
    for bits in itertools.product([0,1], repeat=5):
        selected=[p for p,b in zip(ps,bits) if b]
        w=sum((p-1).bit_length()-1 for p in selected)
        if w>a: continue
        es=[0,1] if w==a else [a-w+1]
        for e in es:
            n=(1<<e)*math.prod(selected)
            ph=(1<<max(0,e-1))*math.prod(p-1 for p in selected)
            assert ph == 1<<a
            out.append(n)
    assert len(out)==len(set(out))
    return sorted(out)

assert all(prime(p) for p in [3,5,17,257,65537,641,6700417])
assert 641*6700417==(1<<32)+1
counts=[len(pure_fiber(a)) for a in range(64)]
assert counts==list(range(2,34))+[32]*32
# Complete inverse sieve: phi(n)>=sqrt(n/2), so targets<=1024 need n<=2097152.
N=2*1024**2
phi=list(range(N+1))
for p in range(2,N+1):
    if phi[p]==p:
        for n in range(p,N+1,p): phi[n]-=phi[n]//p
for a in range(11):
    expected=[n for n in range(1,N+1) if phi[n]==1<<a]
    assert expected==pure_fiber(a),(a,expected,pure_fiber(a))
receipt={"status":"PASS","pure_dyadic_exponents":[0,63],"counts":counts,"independent_sieve_exponents":[0,10],"complete_sieve_bound":N,"F5_factors":[641,6700417],"boundary_fibers":{str(a):pure_fiber(a) for a in [0,1,30,31,32,63]}}
Path(__file__).with_name('dyadic_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='boundary_fibers'},indent=2))
