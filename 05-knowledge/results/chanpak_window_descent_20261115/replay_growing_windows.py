"""First core profiles come from literal subset enumerations; then exact ordinal glue."""
from math import gcd
from replay_matched_slide import F13,D13,F29,D29,require

def baseline(k,l):
    A, pa = F13(k,k if D13>0 else 0)
    B, pb = F29(l,l if D29>0 else 0)
    ha,_,s,_,d,_=pa
    hb,Qb,_,rb,_,f=pb
    t=Qb-hb-rb
    return hb*A+ha*(B-hb)+t*s+f*d

u,v=29*(abs(D13)//2),13*(abs(D29)//2)
require((u,v,gcd(u,v))==(145,39,1))
for k,l in [(38,144),(39,144),(41,150),(60,190)]:
    X={u*i+v*j for i in range(k+1) for j in range(l+1)}
    low=(u-1)*(v-1); high=u*k+v*l-low
    require(all(y in X for y in range(low,high+1)))
    require(low-1 not in X and high+1 not in X)
    for y in range(low,high+1):
        i=(pow(u,-1,v)*y)%v; j=(pow(v,-1,u)*y)%u
        carry=(y-u*i-v*j)//(u*v)
        a=min(carry,(k-i)//v); b=carry-a
        require(0<=a and 0<=b<=(l-j)//u and u*(i+v*a)+v*(j+u*b)==y)
    N=baseline(k,l)
    print('heights',k,l,'order',17+12*k+28*l,'baseline decimal digits',len(str(N)),
          'number of consecutive odds',high-low+1)
