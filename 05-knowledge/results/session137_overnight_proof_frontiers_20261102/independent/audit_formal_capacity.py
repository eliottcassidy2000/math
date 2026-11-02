"""Independent actual/formal capacity firewall; never a tournament counterexample."""
from tournament_kernel import *
from functools import lru_cache
from math import gcd
import json

def run():
    masks=[1794,1792,11,19,7,95,159,63,764,1276,508];n=len(masks);N=1<<n;full=N-1
    a=tuple(tuple((masks[i]>>j)&1 for j in range(n)) for i in range(n));h=hp(a)
    @lru_cache(None)
    def starts(s,i):
        r=s^(1<<i)
        return 1 if not r else sum(starts(r,j) for j in range(n) if r>>j&1 and a[i][j])
    @lru_cache(None)
    def ends(s,i):
        r=s^(1<<i)
        return 1 if not r else sum(ends(r,j) for j in range(n) if r>>j&1 and a[j][i])
    Q=[[0]*n for _ in range(n)]
    for s in range(1,full):
        r=full^s
        for i in range(n):
            if s>>i&1:
                e=ends(s,i)
                for j in range(n):
                    if r>>j&1:Q[i][j]+=e*starts(r,j)
    c=[[Q[i][j]+Q[j][i] for j in range(n)] for i in range(n)];d=[row[:] for row in c]
    for j in (2,3,4):d[0][j]-=498;d[j][0]-=498
    for j in (5,6,7):d[0][j]+=498;d[j][0]+=498
    def response(cap,cut):return h+sum(cap[i][j] for i in range(n) for j in range(n) if a[i][j] and cut>>i&1 and not(cut>>j&1))
    real=[response(c,s) for s in range(N)];formal=[response(d,s) for s in range(N)]
    target=241137;actual_zeros=[s for s,x in enumerate(real) if x==target];formal_zeros=[s for s,x in enumerate(formal) if x==target]
    require(h==24651 and actual_zeros==[1793] and formal_zeros==[480,736,1248,1793],'actual/formal zero sets')
    require(real[480]==239643 and formal[480]==target,'formal repair value')
    for s in (0,480,736,1248,1793,full):require(hp(append(a,s))==real[s],'direct actual child control')
    gr=gf=0
    for x,y in zip(real,formal):gr=gcd(gr,x-h);gf=gcd(gf,y-h)
    require(gr==gf==6,'response gcd')
    require(sum(c[i][j] for i,j in combinations(range(n),2))==sum(d[i][j] for i,j in combinations(range(n),2))==582822,'capacity total')
    for k in range(n+1):require(sum(real[s] for s in range(N) if s.bit_count()==k)==sum(formal[s] for s in range(N) if s.bit_count()==k),'layer sums')
    maxfloor=0
    for i,j in combinations(range(n),2):
        floor=2*hp(induced(a,[v for v in range(n) if v not in (i,j)]));maxfloor=max(maxfloor,floor)
        require(d[i][j]>=floor,'pair-deletion floor')
    require(maxfloor==4194,'max floor')
    for block in ((2,3,4),(5,6,7),(8,9,10)):
        p=list(range(n))
        for i,j in zip(block,block[1:]+block[:1]):p[i]=j
        require(induced(a,p)==a and induced(c,p)==tuple(map(tuple,c)) and induced(d,p)==tuple(map(tuple,d)),'C3 block symmetry')
    nm=[30,28,8,16,4];b=tuple(tuple(nm[i]>>j&1 for j in range(5)) for i in range(5))
    nz=[s for s in range(32) if hp(append(b,s))==9]
    require(nz==[14,22,26,29],'nonstrong actual countercontrol')
    print(json.dumps({'scope':'exact one actual order11 parent and a deliberately formal modified capacity tensor; no tournament-realizability inference',
                      'parent_H':h,'actual_zero_cuts':actual_zeros,'formal_zero_cuts':formal_zeros,
                      'response_gcd':gr,'total_capacity':582822,'max_pair_floor':maxfloor,
                      'nonstrong_zero_cuts':nz,'verdict':'PASS listed invariants survive but actual zero uniqueness does not follow from them alone'},sort_keys=True,indent=2))
if __name__=='__main__':run()
