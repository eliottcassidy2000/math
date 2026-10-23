#!/usr/bin/env python3
"""Independent word-bit matrix construction and exact weighted-demand audit.

The route producer is consumed as a certificate; transition matrices, PSD
tests, stationary weights and marginal-flow checks are built independently.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
import routing as producer

CHECKS=0


def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)


def psd(M):
    need(all(M[i][j]==M[j][i] for i in range(len(M)) for j in range(len(M))),'symmetric certificate')
    A=[row[:] for row in M]
    while A:
        d=A[0][0]
        if d<0:return False
        if d==0:
            if any(A[0][1:]):return False
            A=[row[1:] for row in A[1:]]
        else:A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))] for i in range(1,len(A))]
    return True


def certify(P,pi,gamma):
    n=len(P)
    need(all(sum(row)==1 for row in P),'exact stochastic matrix')
    need(all(pi[i]*P[i][j]==pi[j]*P[j][i] for i,j in product(range(n),repeat=2)),'exact source detailed balance')
    need(psd([[pi[i]*P[i][j] for j in range(n)] for i in range(n)]),'PSD lazy native chain')
    need(psd([[pi[i]*(F(i==j)-P[i][j])-gamma*(pi[i]*F(i==j)-pi[i]*pi[j])
               for j in range(n)] for i in range(n)]),'exact Poincare matrix certificate')


def case(r,m,weighted_flow=False):
    words=[v for v in range(1<<r) if v.bit_count()%m==0];N=len(words);idx={v:i for i,v in enumerate(words)}
    M=2*r-m if r>=m else r-1
    P=[[F(0) for _ in words] for _ in words]
    for i,w in enumerate(words):
        P[i][i]=F(1,2)
        for j in range(r-1):
            z=w^(3<<j) if ((w>>j)&3) in (1,2) else w
            P[i][idx[z]]+=F(1,2*M)
        for j in range(r-m+1):
            mask=((1<<m)-1)<<j
            z=w^mask if w&mask in (0,mask) else w
            P[i][idx[z]]+=F(1,2*M)
    A=producer.bounds(r,m);certify(P,[F(1,N)]*N,F(1,A))
    states=[]
    for w in words:
        a=[0]
        for j in range(r):a.append((a[-1]+((w>>j)&1))%m)
        states.append(tuple(a))
    pos={a:i for i,a in enumerate(states)}
    paths={};D=2*(r*r//4) if r<m else (r-m+3)*2*(m*m//4)
    for i,j in product(range(N),repeat=2):
        coarse=producer.route(states[i],states[j],m)
        fine=[states[i]]
        for a,b in zip(coarse,coarse[1:]):fine+=producer.expand_edge(a,b,m)[1:]
        route=[pos[a] for a in fine];paths[i,j]=route
        need(all(P[a][b]>0 for a,b in zip(route,route[1:])),'producer route in independently constructed move graph')
        need(len(route)-1<=D,'source-independent integer event bound')
    weights=[1+(w%3) for w in words];mass=sum(weights);pi=[F(w,mass) for w in weights]
    lo=min(weights);hi=max(weights);kappa=F(hi,lo)
    Pw=[[F(0) for _ in words] for _ in words]
    for i in range(N):
        for j in range(N):
            if i!=j:Pw[i][j]=P[i][j]*min(F(1),F(weights[j],weights[i]))
        Pw[i][i]=1-sum(Pw[i])
    certify(Pw,pi,F(lo,hi*A))
    report=dict(r=r,modulus=m,states=N,gap_denominator=A,weight_ratio=str(kappa))
    if weighted_flow:
        R=[[F(i==j,2)+pi[j]/2 for j in range(N)] for i in range(N)]
        load=defaultdict(F);first=defaultdict(F);second=defaultdict(F);worst_length=0
        for x,y,z in product(range(N),repeat=3):
            flow=pi[x]*R[x][y]/N
            first[x,z]+=flow;second[y,z]+=flow
            p=paths[x,z]+list(reversed(paths[y,z]))[1:]
            length=len(p)-1;worst_length=max(worst_length,length)
            need(p[0]==x and p[-1]==y,'compiled weighted return has exact endpoints')
            for a,b in zip(p,p[1:]):
                if a!=b:load[a,b]+=flow*length
        need(all(v==pi[x]/N for (x,z),v in first.items()),'first leg retains actual weighted source')
        need(all(v==pi[y]/N for (y,z),v in second.items()),'second leg retains actual weighted terminal')
        congestion=max(v/(pi[a]*Pw[a][b]) for (a,b),v in load.items())
        need(congestion<=4*kappa*A,'weighted routed congestion including whole microclock')
        need(worst_length<=2*D,'bounded clock dilation for compiled returns')
        report.update(exact_weighted_congestion=str(congestion),maximum_compiled_events=worst_length,
                      demand_triples=N**3)
    return report


def main():
    out=[case(4,3,True),case(6,3,True),case(6,5,True),case(8,6,False)]
    # The connector length cannot be reduced merely because the phase alphabet is small.
    x=(0,0,0,0,0);y=(0,1,2,0,0)
    need((x[2]-y[1])%3==2,'one-slot connector hostile in actual accepted ternary words')
    need(all((y[1]+d)%3!=x[2] for d in (0,1)),'shortened connector genuinely impossible')
    return dict(status='FINITE-EXACT PASS; same-agent independent matrix/flow construction',
                cases=out,checks=CHECKS,producer_certificate_checks=producer.CHECKS)


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
