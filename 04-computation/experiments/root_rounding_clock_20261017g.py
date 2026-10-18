#!/usr/bin/env python3
"""One-time source rounding and exact physical-clock checks.

Finite tests never supply a universal convergence or cancellation estimate.
The complete arrays are advanced only along explicitly computed odd edges.
"""
import argparse, hashlib, json, math
from fractions import Fraction as Q
from collections import defaultdict
from pathlib import Path
import numpy as np

def need(ok,m):
    if not ok:raise RuntimeError(m)

def primitive(x):
    s=x.copy();j=np.zeros(len(s),dtype=np.int64)
    while np.any(s%8==5):
        hit=s%8==5;s[hit]=(s[hit]-1)//4;j[hit]+=1
    return s,j

def rows(N):
    s=np.arange(1,4*N+1,2,dtype=np.int64);s=s[s%8!=5]
    exact=np.zeros(len(s),dtype=np.int64);ideal=exact.copy()
    m=s.copy();hom=s.copy();phase=s%3
    while np.any(hom<=4*N):
        for out,y in ((exact,m),(ideal,hom)):
            out += np.where(y<=2*N,np.choose(phase,[0,-3,3]),
                            np.where(y<=4*N,np.choose(phase,[-1,-1,2]),0))
        hom*=4;m=4*m+1;phase=(phase+1)%3
    return s,exact,ideal

def push(s,c):
    a=np.where(s%8==1,2,1);raw=(3*s+1)//(1<<a);t,j=primitive(raw)
    ends,inv=np.unique(t,return_inverse=True);out=np.zeros(len(ends),dtype=np.int64)
    np.add.at(out,inv,np.where(a%2,-c,c))
    return ends,out

def test_rows():
    samples=sorted(set(range(1,513))|{2**k for k in range(10,19)}|{1023,1025,4095,4097,16383,16385,65535,65537})
    maximum=0;witness=None;dyadic=[];transports=[]
    for N in samples:
        s,e,h=rows(N);err=e-h;mass=int(np.abs(err).sum())
        # Exact floor(log_4(4N)); keep the factor4 outside integer division.
        bound=4*(((4*N).bit_length()-1)//2)
        need(mass<=bound,'one-time logarithmic initial error bound')
        if mass>maximum:maximum=mass;witness=N
        if N&(N-1)==0:
            need(np.all(err[s!=1]==0) and abs(int(err[0]))<=2,'dyadic error is one known ROOT atom')
            dyadic.append(dict(N=N,root_correction=int(err[0])))
        if N in (1024,4097,16384,65537):
            terms=s[err!=0];values=err[err!=0];ledger=[]
            for d in range(65):
                need(int(np.abs(values).sum())<=mass,'error does not grow under exact transport')
                if d in (0,1,2,8,32,64):ledger.append(dict(depth=d,error_l1=int(np.abs(values).sum()),atoms=len(terms)))
                terms,values=push(terms,values)
            transports.append(dict(N=N,initial_error=mass,rows=ledger))
    return dict(universe=samples,max_initial_error=maximum,first_max_witness=witness,dyadic=dyadic,transports=transports)

def clock_check(limit=4096,depth=32):
    sources=[s for s in range(1,limit+1,2) if s%8!=5]
    checks=steps=0;digest=hashlib.sha256();examples=[]
    for source in sources:
        s=source;j=0;x=source;physical=exponent=0
        for d in range(1,depth+1):
            a=2 if s%8==1 else 1;raw=(3*s+1)//2**a
            t=raw;nextj=0
            while t%8==5:t=(t-1)//4;nextj+=1
            A=a+2*nextj
            physical+=a+2*j;exponent+=A
            # Replay the actual shortcut map, without stopping at its even ROOT cycle.
            for k in range(a+2*j):x=x//2 if x%2==0 else (3*x+1)//2;steps+=1
            need(x==raw,'actual endpoint after the complete integer clock')
            need(physical==exponent-2*nextj,'terminal braid is exactly the clock coboundary')
            cA=(3+2**a-2**A)//3
            need((3*s+cA)==2**A*t,'primitive affine branch with total scale exponent')
            digest.update(f'{source},{d},{x},{physical},{t},{nextj}\n'.encode())
            if source in (3,7,27) and d<=4:examples.append([source,d,t,nextj,exponent,physical])
            s,j=t,nextj;checks+=1
    return dict(source_count=len(sources),depth=depth,checks=checks,literal_shortcut_steps=steps,sha256=digest.hexdigest(),examples=examples)

def factorization_check():
    count=0
    for q in (3,9,27,81):
        for r in range(q):
            for a in (1,2):
                for j in range(9):
                    Bminus=(Q(r)-Q(4**j-1,3))/4**j
                    lhs=(3*Bminus+1)/2**a
                    rhs=Q(3*r+1,2**(a+2*j))
                    need(lhs==rhs,'A C is exactly the Syracuse affine branch')
                    need(Q(2 if a==1 else 1,3)*Q(3,4)*Q(1,4**j)==Q(1,2**(a+2*j)), 'geometric probability factors')
                    count+=1
    # Hostile: a single-clock zero-sum input cannot contract in one step.
    out=defaultdict(Q)
    for A in range(1,17):
        a=1 if A%2 else 2;j=(A-a)//2
        for r,sign in ((0,1),(1,-1)):
            target=(Q(3*r+1,2**a)-Q(4**j-1,3))/4**j
            rr=target.numerator*pow(target.denominator,-1,9)%9
            out[rr,A-1]+=sign*Q(1,2**A)
    mass=sum(abs(x) for x in out.values())
    need(mass==2*(1-Q(1,2**16)),'finite truncation of exact one-step norm conservation')
    return dict(affine_and_probability_checks=count,single_clock_input_l1=2,single_clock_output_l1_through_A16=str(mass),full_output_l1=2)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    result=dict(status='FINITE-EXACT checks of proved one-time rounding, clock, and positive-kernel identities',rounding=test_rows(),clock=clock_check(),factorization=factorization_check())
    Path(a.output).write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rounding'},sort_keys=True))
    print(json.dumps({k:v for k,v in result['rounding'].items() if k!='universe'},sort_keys=True))
