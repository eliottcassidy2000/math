"""Exact all-order bounded-H fibres of POINTED strong one-apex tournaments.

No graph isomorphism census is hidden here: (T,x), T-x transitive, is uniquely
represented by its apex word. Strong nontrivial words begin in 1 and end in 0.
The inverse is inherited from Session124; the converse/square-shell refinement
is proved in THEORY.md. All decisions and bounds are exact integer arithmetic.
"""
from functools import lru_cache
from math import isqrt
import argparse
import hashlib
import json
from pathlib import Path


def require(c,msg):
    if not c: raise RuntimeError(msg)

@lru_cache(maxsize=None)
def f(q):
    require(type(q) is int and q>=0,'nonnegative integer address required')
    if q<=1: return 1
    n=q//2
    return f(n) if q%2 else f(n)+n+n%2

def word(q):
    require(type(q) is int and q>0,'positive address required')
    return tuple(map(int,bin(q)[2:]))

def G(q):
    require(type(q) is int and q>=1,'q must be positive')
    return q*q+q%2+2*f(q)-1

def reflected_word(q):
    u=word(q)
    return u+tuple(1-b for b in reversed(u))

def inverse(h):
    require(type(h) is int and h>=0,'nonnegative integer target required')
    if h<3 or h%2==0: return []
    return [word(2*n) for n in range(max(1,(h-1)//2),h)
            if f(n)+n+n%2==h]

def bounded_fibres(bound):
    require(type(bound) is int and bound>=1,'positive integer bound required')
    fv=[1]*bound
    counts=[0]*(bound+1)
    for n in range(1,bound):
        if n>1:
            q=n//2
            fv[n]=fv[q] if n%2 else fv[q]+q+q%2
        h=fv[n]+n+n%2
        if h<=bound: counts[h]+=1
    return counts

def summary(bound):
    counts=bounded_fibres(bound)
    fixed={G(q):q for q in range(1,isqrt(bound)+1) if G(q)<=bound}
    for h in range(bound+1):
        require(counts[h]%2==int(h in fixed),f'parity failure H={h}')
    k=isqrt(bound)
    require(len(fixed)==k-1+int(G(k)<=bound),'square-shell cumulative failure')
    rows=[]
    for h in range(3,min(bound,101)+1,2):
        ws=inverse(h)
        require(len(ws)==counts[h],'inverse mismatch')
        rows.append({'H':h,'pointed_strong_classes':counts[h],
                     'pointed_self_converse_classes':int(h in fixed),
                     'words':[''.join(map(str,w)) for w in ws]})
    return {'status':'FINITE-EXACT under proved all-order inverse bound',
            'bound_H':bound,'tested_prefix_addresses':bound-1,
            'pointed_strong_classes_total':sum(counts),
            'distinct_H_values':sum(c>0 for c in counts),
            'odd_fibre_values':len(fixed),
            'missing_odd_H_first30':[h for h in range(3,bound+1,2) if counts[h]==0][:30],
            'maximum_fibre_size':max(counts),
            'first_H_with_maximum_fibre':counts.index(max(counts)),
            'G_first25':[G(q) for q in range(1,26)],
            'counts_sha256':hashlib.sha256(','.join(map(str,counts)).encode()).hexdigest(),
            'first_odd_H_rows':rows}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=1000000)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args();report=summary(args.bound)
    body=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(body)
    else: print(body,end='')
