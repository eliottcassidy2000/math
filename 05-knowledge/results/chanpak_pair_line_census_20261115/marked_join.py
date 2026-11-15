#!/usr/bin/env python3
"""Independent literal DP checks for the marked-host composition packet."""
from collections import defaultdict
from itertools import product
from repair import rows, dp, close

def packet(a,o):
    n=len(a);f=(1<<n)-1
    h,e=dp(a); _,s=dp([f^a[i]^(1<<i) for i in range(n)])
    Q=sum(h[x]*h[f^x] for x in range(f+1))
    ends=[sum(v for i,v in enumerate(r) if o>>i&1) for r in e]
    starts=[sum(v for i,v in enumerate(r) if o>>i&1) for r in s]
    E=sum(ends[x]*h[f^x] for x in range(f+1))
    M=sum(h[x]*starts[f^x] for x in range(f+1))
    W=sum(ends[x]*starts[f^x] for x in range(f+1))
    return h[f],Q,M,E,starts[f],ends[f],W

def join(a,b):
    na=len(a);nb=len(b)
    return [x|(((1<<nb)-1)<<na) for x in a]+[x<<na for x in b]

def joins(a,b):
    h,q,m,l,s,e,w=a
    H,Q,M,L,S,E,W=b
    return h*H,q*Q,h*M+m*Q,l*H+q*L,s*H,h*E,h*W+w*H+m*L+e*S

def main(hostfile, primitive=False):
    hold={}
    atoms=list(sorted(tuple(map(int,line.split())) for line in open(hostfile)))
    for n,code in atoms:
        for o in range(1,1<<n):
            p=packet(rows(n,code),o)
            if primitive and __import__("math").gcd(*p[:2])>1:continue
            if p[:-1] not in hold:hold[p[:-1]]=p[-1],n,code,o
            else:
                a,nn,cc,oo=hold[p[:-1]]
                if a!=p[-1]:
                    print("fork:",(nn,cc,oo),(n,code,o),(a,p[-1]),p[:-1])
                    return (nn,cc,oo),(n,code,o)
    return None

if __name__=="__main__":
    import sys
    f=main(sys.argv[1])
    atoms=[(rows(n,c),o) for n,c in [(1,0),(2,0),(3,2),(4,17)] for o in range(1<<n)]
    for (a,o),(b,O) in product(atoms,repeat=2):
        assert packet(join(a,b),o|(O<<len(a)))==joins(packet(a,o),packet(b,O))
    print("literal marked-join checks:",len(atoms)**2)
    assert f
    A,oa=rows(*f[0][:2]),f[0][2]
    B,ob=rows(*f[1][:2]),f[1][2]
    for mod in (rows(1,0),rows(3,2)):
        hh,qq,*_=packet(mod,0)
        r=[]
        for a,o in [(A,oa),(B,ob)]:
            r.append(dp(close(a,o,mod))[0][-1])
        p=packet(A,oa);p2=packet(B,ob)
        assert r[1]-r[0]==(p[-1]-p2[-1])*hh
        print("literal closure:",hh,qq,r)
