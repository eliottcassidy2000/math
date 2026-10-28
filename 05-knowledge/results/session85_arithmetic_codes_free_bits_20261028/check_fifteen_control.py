#!/usr/bin/env python3
from collections import Counter
from math import isqrt,prod

def need(b,m):
    if not b:raise RuntimeError(m)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def sigmask(n):
    c=0
    for i,p in enumerate([2,3,5,7]):
        while n%p==0:n//=p;c^=1<<i
    need(n==1,'outside squareclass basis');return c

def main():
    P=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
    H=[sigmask(p+1) for p in P];need(all(map(prime,P)),'prime list')
    need(H==[2,0,3,1,2,9,1,4,3,7,1,11,2,6,1],'columns')
    B=[2,1,6,9];change={0:0}
    for x in range(16):
        old=0
        for i,b in enumerate(B):
            if x>>i&1:old^=b
        change[old]=x
    need(len(change)==16,'invertible basis')
    edges={2:(0,1),11:(0,1),47:(0,1),7:(0,2),17:(0,2),31:(0,2),71:(0,2),59:(0,3),5:(1,2),23:(1,2),19:(1,3),29:(2,3),13:(0,4),41:(1,4),3:(0,0)}
    def vec(v):return 0 if v==0 else 1<<(v-1)
    for p,h in zip(P,H):
        u,v=edges[p];need(change[h]==vec(u)^vec(v),'incidence identity')
    weights=Counter();count=0
    for x in range(1<<15):
        h=0;S=1;deg=0
        for i,p in enumerate(P):
            if x>>i&1:
                h^=H[i];S*=p+1;u,v=edges[p];deg^=(1<<u)^(1<<v)
        square=isqrt(S)**2==S
        need(square==(h==0)==(deg==0),'complete cycle/sigma equivalence')
        if square:
            weights[x.bit_count()]+=1;count+=1
            need((x>>P.index(13)&1)==(x>>P.index(41)&1),'two-bit check')
    need(count==2048 and weights[1]==1 and weights[2]==10,'code size/distance')
    need(weights=={0:1,1:1,2:10,3:46,4:99,5:183,6:300,7:372,8:375,9:319,10:202,11:94,12:37,13:9},'complete primal enumerator')
    dual=Counter(sum((a&h).bit_count()%2 for h in H) for a in range(16))
    need(dual=={0:1,2:1,3:1,5:1,7:3,8:3,9:3,10:3},'dual enumerator')
    six=[2,5,7,19,29,59]
    need({sigmask(p+1) for p in six}=={1,2,3,4,6,7},'K4 restriction')
    need(sigmask(90)==5 and prime(89) and not prime(9) and not prime(39),'least Fano completion')
    need({sigmask(p+1) for p in six+[89]}==set(range(1,8)),'Fano columns')
    total=0
    for h in H:total^=h
    need(total==5 and total^sigmask(90)==0,'full bank repair guard')
    for p,q in [(2,29),(7,19),(59,5)]:
        need(sigmask(p+1)^sigmask(q+1)==sigmask(90),'matching boundary')
        need(isqrt((p+1)*(q+1)*90)**2==(p+1)*(q+1)*90,'new Fano circuit')
    matchings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
    from itertools import product
    for charges in product(range(-2,3),repeat=4):
        energies=[sum((charges[i]+charges[j])**2 for i,j in M) for M in matchings]
        need(sum(energies)==sum(charges)**2+2*sum(q*q for q in charges),'three-pairing energy')
    need([sum(((1,1,-1,-1)[i]+(1,1,-1,-1)[j])**2 for i,j in M) for M in matchings]==[8,0,0],'energy not boundary')
    print('PASS supersingular-bank graph on5vertices;2048words; dual',sorted(dual.items()),'primal',sorted(weights.items()))
if __name__=='__main__':main()
