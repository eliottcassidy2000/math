#!/usr/bin/env python3
from verify_extension import *
from math import gcd,factorial
from functools import reduce

def matching(a,r):
    O=[v for v in range(len(a)) if a[r]>>v&1];I=[v for v in range(len(a)) if a[v]>>r&1];mate={}
    def aug(u,seen):
        for v in I:
            if v in seen or not a[u]>>v&1:continue
            seen.add(v)
            if v not in mate or aug(mate[v],seen):mate[v]=u;return True
        return False
    return sum(aug(u,set()) for u in O)

def run():
    context_checks=0;profile_checks=0
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask)
            for r in range(n):
                op=skeleton_operator(a,r);nu=matching(a,r)
                require(max(j for _,j in op)==nu,('forward bandwidth',a,r,nu,op))
                require(all(abs(j)<=nu for _,j in op),'absolute bandwidth')
                require(all(c>0 for (power,j),c in op.items() if j==nu),'positive leading coefficient')
                context_checks+=1
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);p=bin_polynomial(a);u=coordinates(p,n%2)
            weighted=[v*factorial(2*k+1)*(1 if n%2 else k+1) for k,v in enumerate(u)]
            short=reduce(gcd,[value(p,d) for d in range(1,(n+1)//2+1)])
            longer=reduce(gcd,[value(p,d) for d in range(1,n+8)])
            require(reduce(gcd,weighted)==short==longer,('content',a,p,u))
            for k in range(len(u)):
                g=factorial(2*k+1)*(1 if n%2 else k+1)
                for d in range(-10,11):
                    require(value(basis(k,n%2),d)%g==0,'integer basis divisor')
            profile_checks+=1
    print(json.dumps({'rooted_matching_bandwidth_checks':context_checks,'profile_content_checks':profile_checks},sort_keys=True,indent=2))
if __name__=='__main__':run()
