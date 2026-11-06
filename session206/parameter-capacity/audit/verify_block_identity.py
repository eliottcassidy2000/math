"""Independent finite regression of block identity; proofs remain analytic."""
import json

def q(m):
    ans=0
    while m:
        p=m//2
        if m&1: ans+=p+(p&1)
        m=p
    return ans

def rev(r,k):
    z=0
    for _ in range(k): z=2*z+(r&1); r//=2
    return z

checks=0
for k in range(1,11):
    for r in range(1<<k):
        for v in range(65):
            actual=q((v<<k)+r)
            expected=q(v)+q(r)+v*rev(r,k)+((r>>(k-1))&1)*(v&1)
            if actual!=expected: raise RuntimeError((k,r,v,actual,expected))
            checks+=1
print(json.dumps({'status':'PASS','checks':checks,'k_max':10,'v_max':64},indent=2))
