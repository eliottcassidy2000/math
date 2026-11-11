"""Exact constant-gcd terminal sectors; no full-count speed claim."""
from collections import defaultdict
from math import factorial,gcd,comb
from time import perf_counter
import json

def kernel(N,S,g,charges=None,*,gauge=True):
    S=tuple(sorted(set(S))); charges=charges or {}
    assert all(l>0 and l%2 and l>=g and l%g==0 for l in S) and g>0 and g%2
    assert all(gcd(a,b)==g for i,a in enumerate(S) for b in S[i+1:] if a+b<=N)
    f=[factorial(n) for n in range(N+1)]
    states={(0,0):1}; trans=0; peak=1
    for l in S:
        new=defaultdict(int); cache={}
        for (d,k),v in states.items():
            if d not in cache:
                cache[d]=[(m,d+l*m,f[d+l*m]//(f[d]*l**m*f[m]),((l-g)*m*m+(g-1)*m)//2+charges.get(l,0)*m) for m in range((N-d)//l+1)]
            for m,nd,c,e in cache[d]:
                new[nd,k+m]+=(v*c)<<(e if gauge else e+g*(m*(m-1)//2+k*m)); trans+=1
        states=new; peak=max(peak,len(states))
    out=[0]*(N+1)
    for (n,k),v in states.items():out[n]+=v<<(g*k*(k-1)//2 if gauge else 0)
    return out,{'transitions':trans,'peak_states':peak}

def brute(N,S,charges=None):
    """Independent multiplicity walk with original gcd cross terms."""
    S=tuple(S);charges=charges or {}; f=[factorial(n) for n in range(N+1)]
    out=[0]*(N+1);leaves=0
    def visit(i,d,e,z,history):
        nonlocal leaves
        if i==len(S):
            out[d]+=(f[d]//z)<<e;leaves+=1;return
        l=S[i];r=sum(gcd(l,a)*m for a,m in history)+charges.get(l,0)
        for m in range((N-d)//l+1):
            visit(i+1,d+l*m,e+(l-1)*m//2+l*m*(m-1)//2+m*r,z*l**m*f[m],history+[(l,m)])
    visit(0,0,0,1,[])
    return out,{'multiplicity_vectors':leaves}

def run():
    rows=[]
    for N,S,g,q in [(30,(1,3,5,7),1,{}),(90,(3,9,15,21),3,{9:7,15:4}), (24,(9,15,21),3,{}),(20,(9,15),1,{})]:
        a,_=kernel(N,S,g,q);b,_=brute(N,S,q);assert a==b
    try:kernel(24,(9,15),1)
    except AssertionError:pass
    else:raise AssertionError('missed hostile active gcd edge')
    for R in (70,100,130):
        S=tuple(3*p for p in (1,3,5,7,11,13,17,19,23,29,31,37,41,43));N=3*R
        t=perf_counter();a,info=kernel(N,S,3);kt=perf_counter()-t
        t=perf_counter();b,bi=brute(N,S);bt=perf_counter()-t
        assert a==b
        rows.append(dict(N=N,support=list(S),kernel_seconds=kt,brute_seconds=bt,speedup=bt/kt,**info,**bi,exact_equal=True,terminal_bits=a[N].bit_length()))
        print(json.dumps(rows[-1]),flush=True)
    return rows
if __name__=='__main__':run()
