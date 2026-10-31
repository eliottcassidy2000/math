#!/usr/bin/env python3
"""Finite-generator theorem controls on certified rows, plus late-item hostiles."""
from fractions import Fraction as Q
import heapq,json
from check_rows import w,row,pats,ROWS,need

def bounds(a,thresholds):
    s=w(a);L=1<<(a-1).bit_length();D=(s-a)*(L-1)
    M=L;joins=[];patterns=[]
    for x in range(1,a):
        gap=s-w(x)-w(a-x);need(gap>0,'join strictness')
        cut=int(D//gap);M=max(M,cut);joins.append([x,str(gap),cut])
    for p in pats(a):
        c=[w(p[i]) for i in (0,2,4,1,3)];half=sum(c)/2
        m=max([Q(0),half]+c+[c[i]+c[(i+2)%5] for i in range(5)])
        if half<s:
            gap=s-m;need(gap>0,'nonsharp root gap')
            cut=int(D//gap);r={'pattern':p,'gap':str(gap),'cutoff':cut}
        else:
            need(half==s,'slope invalid')
            ys=[sum((-1)**j*c[(i-j)%5] for j in range(5))/2 for i in range(5)]
            need(min(ys)>0,'positive sharp dual')
            S=int(D//min(ys));N=max(thresholds[x] for x in p)
            cut=2*N+S+2*L-1
            r={'pattern':p,'dual':[str(y) for y in ys],'S':S,'N':N,'cutoff':cut}
        M=max(M,cut);patterns.append(r)
    return M,{'a':a,'slope':str(s),'L':L,'D':str(D),'M':M,'joins':joins,'patterns':patterns}

def paths(s,L,coins):
    p,q=s.numerator,s.denominator
    need(L in coins and q*coins[L]==p*L,'missing zero-cost L coin')
    edges=[]
    for b,n in coins.items():
        need(isinstance(b,int) and b>0,'positive integer coin length')
        need(isinstance(n,int),'integer coin order')
        e=p*b-q*n;need(e>=0,'negative deficit edge');edges.append((b,e))
    dist=[None]*L;plans=[None]*L;dist[0]=(0,0);plans[0]=[]
    heap=[(0,0,0)]
    while heap:
        cost,length,r=heapq.heappop(heap)
        if dist[r]!=(cost,length):continue
        for b,e in edges:
            rr=(r+b)%L;z=(cost+e,length+b)
            if dist[rr] is None or z<dist[rr]:
                dist[rr]=z;plans[rr]=plans[r]+[b];heapq.heappush(heap,(*z,rr))
    M=max(coins)
    for r,z in enumerate(dist):
        need(z is not None,'unreachable residue')
        pp=plans[r];need(len(pp)<=L-1,'shortest path cycle')
        need(sum(pp)==z[1] and sum(p*b-q*coins[b] for b in pp)==z[0],'path realization')
        need(z[1]%L==r and z[1]<=(L-1)*M,'path length bound')
    return dist,plans

def coin_dp(coins,upto):
    f=[0]+[-10**20]*upto
    for b in range(1,upto+1):
        f[b]=max((f[b-c]+n for c,n in coins.items() if c<=b),default=-10**20)
    return f

def main():
    thresholds={1:1};records=[]
    for a in range(2,9):
        M,rec=bounds(a,thresholds);L=rec['L'];s=w(a);coins={b:row(a,b) for b in range(1,M+1)}
        dist,plans=paths(s,L,coins);N=max(1,max(z[1] for z in dist))
        upto=max(2*M+3*L,N+L)
        f=coin_dp(coins,upto)
        need(all(f[b]==row(a,b) for b in range(1,upto+1)),'coin closure vs certified row')
        for b in range(N,upto+1):
            need(s*b-Q(dist[b%L][0],s.denominator)==f[b],'residue tail')
        failures=[b for b in range(1,N) if f[b+L]!=f[b]+s*L]
        least=1+max(failures,default=0);need(least==1,'known threshold one')
        thresholds[a]=least
        rec.update({'residue_costs_and_first_weights':dist,'paths':plans,
                    'computed_sufficient_threshold':N,'coarse_threshold_bound':max(1,(L-1)*M),
                    'least_threshold':least,'coin_replay_through':upto})
        records.append(rec)
    hostiles=[]
    for M in (3,11,101):
        coins={1:1,2:4,M:2*M};s=Q(2);L=2
        dist,plans=paths(s,L,coins);f=coin_dp(coins,2*M+2)
        old=coin_dp({1:1,2:4},2*M+2)
        need(dist[1]==(0,M),'late residue first realization')
        need(f[M]==2*M and old[M]==2*M-1,'late generator counterexample')
        need(all(f[b]==2*b-int(b%2==1 and b<M) for b in range(1,len(f))),'hostile exact formula')
        need(f[1]==1 and 2-Q(dist[1][0],1)==2,'premature residue law hostile')
        hostiles.append({'M':M,'odd_first_weight':dist[1][1],'shortest_path':plans[1],
                         'prefix_omission_loss_at_M':f[M]-old[M]})
    rejected=0
    for coins in ({1:3,2:4},{1:1,2:3}):
        try:paths(Q(2),2,coins)
        except RuntimeError:rejected+=1
    need(rejected==2,'hostile guard rejection')
    print(json.dumps({'status':'PASS','scope':'finite-generator and shortest-path controls on certified rows2..8',
        'rows':records,'late_generator_hostiles':hostiles,'invalid_coin_sets_rejected':rejected},sort_keys=True,indent=2))
if __name__=='__main__':main()
