#!/usr/bin/env python3
"""Exact A000568 engines. Standard library only; no floating point counting."""
from collections import defaultdict, Counter
from math import factorial, gcd, comb
from fractions import Fraction
from itertools import combinations

def require(x, label):
    if not x: raise AssertionError(label)

def integer_response(N, modulus=None, *, binary_precision=None):
    """Factorial-scaled future-response DP, all coefficients 0..N.
    modulus, if supplied, returns a_n modulo modulus using modulus*N!.
    """
    require(N >= 0, 'nonnegative order')
    require(binary_precision is None or binary_precision >= 1, 'positive binary precision')
    require(modulus is None or binary_precision is None, 'choose one residue mode')
    require(modulus is None or modulus > 0, 'positive modulus')
    fact=[factorial(n) for n in range(N+1)]
    ring=None if modulus is None else modulus*fact[N]
    bits=None if binary_precision is None else binary_precision+N-N.bit_count()
    if bits is not None: ring=1<<bits
    lengths=list(range(1,N+1,2))[::-1]
    states={(0,(0,)*len(lengths)):1}
    stats={'peak_states':1,'transitions':0}
    for i,ell in enumerate(lengths):
        future=lengths[i+1:]
        new=defaultdict(int)
        for (d,response),w in states.items():
            live=[j for j in lengths[i:] if j<=N-d]
            old=dict(zip(live,response))
            charge=old.get(ell,0)
            for k in range((N-d)//ell+1):
                q=ell*k; nd=d+q
                coefficient=fact[nd]//(fact[d]*ell**k*fact[k])
                e=k*(ell-1)//2+ell*k*(k-1)//2+k*charge
                r=tuple(old[j]+k*gcd(ell,j) for j in future if j<=N-nd)
                if bits is not None and e >= bits: continue
                term=(w*coefficient)<<e
                key=(nd,r)
                if ring: new[key]=(new[key]+term)%ring
                else: new[key]+=term
                stats['transitions']+=1
        states={key:value for key,value in new.items() if value}
        stats['peak_states']=max(stats['peak_states'],len(states))
    result=[0]*(N+1)
    for (n,r),w in states.items():
        require(not r,'terminal response')
        if binary_precision is not None:
            v=n-n.bit_count(); odd=fact[n]>>v
            require(w%(1<<v)==0,'terminal binary divisibility')
            result[n]=((w>>v)*pow(odd,-1,1<<binary_precision))% (1<<binary_precision)
        else:
            require(w%fact[n]==0,'terminal factorial divisibility')
            result[n]=w//fact[n]
            if modulus: result[n]%=modulus
    return result,stats

def core_series(N):
    """One DFS through nonunit odd cores, integer kernel K[m,t]."""
    fact=[factorial(n) for n in range(N+1)]
    kernel=defaultdict(int); nodes=0
    def visit(parts,m,t,e,c):
        nonlocal nodes
        nodes+=1; kernel[m,t]+=c<<e
        lo=parts[-1] if parts else 3
        for r in range(lo,N-m+1,2):
            k=parts.count(r)+1
            numerator=c*fact[m+r]
            denominator=fact[m]*r*k
            require(numerator%denominator==0,'class-size transition')
            visit(parts+(r,),m+r,t+1,e+(r-1)//2+sum(gcd(r,s) for s in parts),numerator//denominator)
    visit((),0,0,0,1)
    out=[]
    for n in range(N+1):
        total=sum(comb(n,m)*w << ((n-m)*(n-m-1)//2+(n-m)*t)
                  for (m,t),w in kernel.items() if m<=n)
        require(total%fact[n]==0,'core terminal divisibility')
        out.append(total//fact[n])
    return out,{'core_nodes':nodes,'kernel_states':len(kernel)}

def partitions(n,lo=1):
    if not n: yield ()
    for r in range(lo,n+1,2):
        for tail in partitions(n-r,r): yield (r,)+tail

def direct_series(N):
    out=[]
    for n in range(N+1):
        a=Fraction(0)
        for p in partitions(n):
            e=sum((r-1)//2 for r in p)+sum(gcd(r,s) for r,s in combinations(p,2))
            z=1
            for r,k in Counter(p).items(): z*=r**k*factorial(k)
            a+=Fraction(1<<e,z)
        require(a.denominator==1,'Burnside integrality')
        out.append(a.numerator)
    return out
