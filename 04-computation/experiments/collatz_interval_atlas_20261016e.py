#!/usr/bin/env python3
"""Independent row-controller certificates, literal transport, and sharp cones.

The prescribed finite atlas uses stages 6,7,8 cyclically on midpoint anchors.
No eventual ROOT sign or convergence oracle is used.
"""
import argparse
from collections import defaultdict
import copy
import hashlib
import json
import math
import numpy as np
from collatz_compact_resets_20261015f import Graph,T

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def weight(n):return 2 if n%3==2 else -1
def tokens(n):return [(n,i) for i in range(abs(weight(n)))]
def norms(X,g):
    row=defaultdict(int);total=0;out=[]
    def change(n,c):
        nonlocal total
        e=int(g.endpoint[n]);old=row[e]
        row[e]+=c*weight(n)*(-1)**int(g.clock[n])
        total+=abs(row[e])-abs(old)
    for n in range(X+1,4*X+1):change(n,1)
    for N in range(X,2*X):
        out.append(total)
        need(total>=0 and total%2==0,'integer nonnegative even residual')
        if N+1<2*X:
            change(N+1,-1)
            for n in range(4*N+1,4*N+5):change(n,1)
    need(all(abs(a-b)<=4 for a,b in zip(out,out[1:])),'common-controller Lipschitz control')
    return out

def point_certificate(N,g):
    paths={};fibres=defaultdict(lambda:[[],[]]);pairs={}
    for n in range(N+1,4*N+1):
        e,d=int(g.endpoint[n]),int(g.clock[n]);sgn=weight(n)*(-1)**d
        for a in tokens(n):paths[a]=(e,d);fibres[e][sgn<0].append(a)
    for pos,neg in fibres.values():
        for a,b in zip(pos,neg):pairs[a]=b;pairs[b]=a
    return [N,paths,pairs]

def replay(cert):
    N,paths,pairs=cert
    need(set(paths)=={a for n in range(N+1,4*N+1) for a in tokens(n)},'literal annular source atoms')
    for a,(e,d) in paths.items():
        x=a[0]
        for _ in range(d):
            need(x!=1,'no step out of stopped ROOT');x=T(x)
        need(x==e and e%2==1,'actual integer endpoint and full clock; odd terminal')
    for a,b in pairs.items():
        need(pairs[b]==a and a[0]!=b[0],'unit capacity and distinct sources')
        e,d=paths[a];f,c=paths[b]
        need(e==f and (weight(a[0])>0)*2-1 == -(((weight(b[0])>0)*2-1)*(-1)**(c-d)),
             'opposite weighted clock signs at the exact same terminal')
    return len(paths)-len(pairs)

def move(cert,direction):
    N,paths,pairs=cert
    def discard(a):
        if a in pairs:
            b=pairs.pop(a);pairs.pop(b)
        paths.pop(a)
    def rename(a,b,delta):
        e,d=paths.pop(a);need(d+delta>=0,'backward transport strips two real halvings')
        paths[b]=(e,d+delta)
        if a in pairs:
            c=pairs.pop(a);pairs[c]=b;pairs[b]=c
    if direction==1:
        m=N+1
        for a in tokens(m):rename(a,(4*m,a[1]),2)
        for n in range(4*m-3,4*m):
            x=n;d=0
            while x%2==0:x//=2;d+=1
            for a in tokens(n):paths[a]=(x,d)
        cert[0]+=1
    else:
        need(direction==-1 and N>1,'valid inverse row move');m=N
        for n in range(4*m-3,4*m):
            for a in tokens(n):discard(a)
        for a in tokens(4*m):rename(a,(m,a[1]),-2)
        cert[0]-=1

def cone_envelope(X,anchors):
    infinity=10**100;out=[infinity]*X
    for a,B,stage in anchors:out[a-X]=min(out[a-X],B)
    for i in range(1,X):out[i]=min(out[i],out[i-1]+4)
    for i in range(X-2,-1,-1):out[i]=min(out[i],out[i+1]+4)
    need(all(abs(a-b)<=4 for a,b in zip(out,out[1:])),'cone envelope is Lipschitz')
    need(all(out[a-X]<=B for a,B,stage in anchors),'cone obeys every anchor bound')
    return out

def controls():
    g=Graph(512)
    for _ in range(3):
        cols,inv=np.unique(g.endpoint[1:],return_inverse=True);g.extend(cols);g.close(cols,inv)
    cases=[]
    for anchor in [65,77,96,119,127]:
        base=point_certificate(anchor,g);B=replay(base)
        for direction in [-1,1]:
            cert=copy.deepcopy(base)
            for dist in range(1,9):
                if not 64<=anchor+direction*dist<128:break
                move(cert,direction);b=replay(cert)
                need(b<=B+4*dist,'literal transported certificate pays at most four per row')
                cases.append([anchor,direction*dist,B,b])
    # Same endpoint residues at arbitrarily many chosen finite places are not equality.
    moduli=[6,30,210,2310,30030]
    hostiles=[]
    for Q in moduli:
        a=Q+3;b=2*Q+3;n=2*a;m=b
        need(n//2==a and n%2==0 and m%2==1,'native hostile clocks one and zero')
        need(a!=b and (a-b)%Q==0,'distinct exact endpoints with identical residue signatures')
        need(weight(n)==weight(m)==-1,'opposite weighted clock signs')
        need(a<n<=4*a and a<m<=4*a,'both original sources share one physical annulus')
        hostiles.append(dict(modulus=Q,row=a,sources=[n,m],ends=[a,b],clocks=[1,0]))
    for X in range(1,65):
        anchors=[(X,2,0),(2*X-1,4,0)]
        C=cone_envelope(X,anchors)
        need(C==[min(B+4*abs(N-a) for a,B,stage in anchors) for N in range(X,2*X)],
             'distance-transform envelope equals direct formula')
    return dict(transport_cases=cases,residue_hostiles=hostiles)

def run(k):
    X=1<<k;b=1<<(k//2);g=Graph(8*X)
    positions=list(range(X+(b-1)//2,2*X,b))
    requested={a:6+i%3 for i,a in enumerate(positions)}
    anchors=[];snapshots=[];last=None
    for stage in range(9):
        if stage:
            cols,inv=np.unique(g.endpoint[1:],return_inverse=True);g.extend(cols);g.close(cols,inv)
        if stage<6:continue
        ns=norms(X,g)
        if last is not None:need(all(a>=b for a,b in zip(last,ns)),'native extension contracts every row')
        last=ns
        for a in positions:
            if requested[a]==stage:
                anchors.append([a,ns[a-X],stage])
                if k<=8:
                    cert=point_certificate(a,g)
                    need(replay(cert)==ns[a-X],'independent literal anchor matching')
        g.replay(range(1,8*X+1))
        snapshots.append(dict(stage=stage,max_clock=int(g.clock.max()),
          endpoint_sha256=hashlib.sha256(g.endpoint.tobytes()).hexdigest(),
          clock_sha256=hashlib.sha256(g.clock.tobytes()).hexdigest()))
    anchors.sort();C=cone_envelope(X,anchors)
    need(all(v<=c for v,c in zip(last,C)),'full stage-eight norms lie below the independently chosen anchor cones')
    Q2=sum(b*B*B for a,B,stage in anchors)
    S1=b*b//4
    S2=b*(b*b-1)//12 if b%2 else b*(b*b+2)//12
    box=sum(b*B*B+8*B*S1+16*S2 for a,B,stage in anchors)
    assigned=[]
    for i,N in enumerate(range(X,2*X)):
        a,B,stage=anchors[i//b];assigned.append(B+4*abs(N-a))
    need(box==sum(v*v for v in assigned),'exact midpoint-bin polynomial bound')
    need(sum(c*c for c in C)<=box,'all-anchor overlap only improves the bound')
    return dict(k=k,X=X,bin_width=b,anchors=anchors,snapshots=snapshots,
        anchor_energy=Q2,cone_energy=sum(c*c for c in C),bin_energy=box,
        full_stage_eight_energy=sum(v*v for v in last),passes_7X=sum(c*c for c in C)<=49*X*X,
        cone_sha256=hashlib.sha256(json.dumps(C).encode()).hexdigest())

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--k',type=int,nargs='+',default=[6,8,10,12,14,16])
    args=ap.parse_args();ctl=controls();runs=[run(k) for k in args.k]
    print(json.dumps(dict(status='FINITE-EXACT; uniform anchor bound remains OPEN',controls=ctl,
                         runs=runs,checks=CHECKS),separators=(',',':'),sort_keys=True))
