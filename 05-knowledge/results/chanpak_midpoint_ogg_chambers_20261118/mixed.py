"""Choose a literal Q-probe address for every unit mod 2*3^k*5^ell."""
from math import gcd
import argparse
from replay import dp, join, h

SEEDS=[
    ((30,60,56,48,32,1),17,224),
    ((14,60,56,48,33,1),23,296),
    ((22,44,56,49,34,1),37,376),
    ((14,20,56,50,33,3),43,424)]
CHART={(x+2*y*2**j)%30:(R,x,y,j) for R,x,y in SEEDS for j in (1,3)}

def lift16(value,p,k):
    a=0;s=1
    for m in range(2,k+1):
        a=next(a+z*s for z in range(p) if pow(16,a+z*s,p**m)==value%p**m)
        s*=p
    return a,s

def select(N,k,ell):
    if min(k,ell)<1 or gcd(N,30)!=1:raise ValueError("unit target and positive local depths required")
    R,x,y,j=CHART[N%30]
    m,n=3**k,5**ell;M=m*n
    B=2*y*2**j
    z=((N-x)*pow(B,-1,M))%M
    a,u=lift16(z,3,k);b,v=lift16(z,5,ell)
    A=a+u*((b-a)*pow(u,-1,v)%v) if v!=1 else a
    t=j+4*A
    assert (x+2*y*pow(2,t,2*M)-N)%(2*M)==0
    return R,t

def Q(F):
    K=(2,4,9,3);count=len(F)
    out=[]
    V=[(0,i) for i in range(count)]+[(i,0) for i in range(1,4)]
    for v,a in V:
        out.append(sum(1<<j for j,(w,b) in enumerate(V)
                   if (v==w==0 and F[a]>>b&1) or (v!=w and K[v]>>w&1)))
    return out

def TT(t):return [((1<<t)-1)^((1<<(i+1))-1) for i in range(t)]

def check():
    assert sorted(CHART)==[1,7,11,13,17,19,23,29]
    for R,x,y in SEEDS:
        full=(1<<len(R))-1;H=[sum(row) for row in dp(R)];H[0]=1
        assert H[full]==x and sum(H[s]*H[full^s] for s in range(full+1))==y
    for N,(R,x,y,j) in sorted(CHART.items()):
        measured=h(Q(join(R,TT(j))))
        assert measured==x+2*y*2**j
        print(N, 'mod30:',x,y,j,'small literal H',measured)
    for k,l in ((1,1),(2,2),(3,2),(4,3)):
        M=2*3**k*5**l
        for N in range(1,M,2):
            if gcd(N,30)==1:select(N,k,l)
        print("all units independently addressed mod",M)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('target',type=int,nargs='?')
    p.add_argument('k',type=int,nargs='?',default=1);p.add_argument('ell',type=int,nargs='?',default=1)
    a=p.parse_args()
    if a.target is None:check()
    else:
        rows,t=select(a.target,a.k,a.ell)
        print('ordinal seed masks',rows,'transitive pad order',t,'apply Q after padding')
