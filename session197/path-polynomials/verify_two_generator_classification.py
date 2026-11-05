#!/usr/bin/env python3
from verify_unit_offset_rigidity import *
from math import gcd
from fractions import Fraction

K=(0,);W=(8,1,3,6)

def build(t,s):
    U=transitive(4+t);V=transitive(t)
    for _ in range(s):U=join(U,W)
    for _ in range(s+1):V=join(V,W)
    require(U and V,'actual nonempty modules')
    return U,V

def strong(a):
    for v in range(len(a)):
        seen=1<<v
        while True:
            nex=seen
            for u in range(len(a)):
                if seen>>u&1:nex|=a[u]
            if nex==seen:break
            seen=nex
        if seen!=(1<<len(a))-1:return False
    return True

def pure25(n):
    while n%5==0:n//=5
    return n==1

def run():
    require(profile(K)==(1,2,3),'singleton physical profile')
    require(profile(W)==(5,32,117),'H5 physical profile')
    leaf_tests=0;accepted=[]
    for r in range(3,4097):
        n,d=2**r-1,2**(r-2)-1;g=gcd(n,d)
        require(g==(3 if r%2==0 else 1),'exact gcd of odd parts')
        N,D=n//g,d//g
        require(gcd(N,D)==1,'reduced coprime odd parts')
        if r>=5:require(N>1 and D>1,'both odd parts nontrivial')
        possible=pure25(N) and pure25(D)
        require(possible==(r==4),'complete transitive leaf classification control')
        if possible:accepted.append(r)
        leaf_tests+=1
    allocations=0;hits=0
    for r,uk,uw,vk,vw in product(range(1,17),range(9),range(5),range(9),range(5)):
        if uk+uw==0 or vk+vw==0:continue
        hu,bu=5**uw,2**(uk+5*uw);hv,bv=5**vw,2**(vk+5*vw)
        B=2*(2**r-1)*bv*hu+(4-2**r)*bu*hv
        expected=r==4 and uk==vk+4 and vw==uw+1
        require((B==0)==expected,'literal allocation cancellation equivalence')
        if expected:
            t,s=vk,uw;A=2*(2**r-1)*bu*bv;C=hu*hv*(-1)*(2**r-3)
            require((A,C)==(15360*4**t*1024**s,-65*25**s),'all ray coefficients')
            hits+=1
        allocations+=1
    physical=[]
    for t,s,m in ((0,0,1),(0,0,2),(1,0,1)):
        U,V=build(t,s);a=graph(transitive(m),K,U,transitive(4),V)
        H=h_subsets(a)[-1];want=15360*4**t*1024**s*3**m-65*25**s
        require(H==want and len(a)==m+15+2*t+8*s and strong(a),'physical witness count/order/strongness')
        physical.append({'t':t,'s':s,'m':m,'n':len(a),'H':H,'rows':a})
    t,s,m=0,1,1;U,V=build(t,s);a=graph(transitive(m),K,U,transitive(4),V)
    require(profile(U)==(5,512,9477),'nontrivial U join profile')
    require(profile(V)==(25,1024,13689),'nontrivial V join profile')
    require(coeff(K,U,transitive(4),V)==(15360*1024,0,-65*25),'s1 physical module coefficients')
    require(len(a)==24 and strong(a),'larger witness construction only')
    larger={'t':t,'s':s,'m':m,'n':len(a),'H_from_proved_identity':15360*1024*3-65*25,'rows':a}
    arithmetic=0
    for t,s in product(range(41),repeat=2):
        A=15360*4**t*1024**s;C=-65*25**s
        require(gcd(A,C)==5 and C%3,'coefficient gcd')
        require(gcd(3*A+C,9*A+C)==5,'exact ray content')
        for m in range(1,12):
            H=A*3**m+C
            require(H%5==0 and H%2==1 and H>=46015*4**t*1024**s,'uniform lower bound and residue')
            arithmetic+=1
    cycle=(0,2,1,5,3,4,6)
    require(all(Q[u]>>v&1 for u,v in zip(cycle,cycle[1:]+cycle[:1])),'strong quotient cycle')
    print(json.dumps({'status':'PASS','leaf_lengths_checked':leaf_tests,'accepted_leaf_lengths':accepted,'allocation_cases':allocations,'cancelling_allocations':hits,'physical_direct_counts':physical,'larger_constructed_witness':larger,'arithmetic_controls':arithmetic,'exact_common_divisor':5},sort_keys=True,indent=2))
if __name__=='__main__':run()
