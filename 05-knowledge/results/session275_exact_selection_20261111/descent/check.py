"""Exact controls for projective proper-graft cofactor 21; no coverage scan."""
import json
from itertools import permutations
from math import gcd
from pathlib import Path

B=(108,41,2,4,47,12,62)

def packet(rows):
    n=len(rows); m=1<<n
    h=[0]*m; h[0]=1; dp=[[0]*n for _ in range(m)]
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,m):
        for v in range(n):
            if s>>v&1:
                q=s^(1<<v)
                if q: dp[s][v]=sum(dp[q][u] for u in range(n) if (q>>u&1) and (rows[u]>>v&1))
        h[s]=sum(dp[s])
    return h[-1],sum(h[s]*h[(m-1)^s] for s in range(m))

def literal(rows):
    n=len(rows)
    return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))

def join(a,b):
    n=len(a);m=len(b)
    return tuple(r|(((1<<m)-1)<<n) for r in a)+tuple(r<<n for r in b)

def cycle_blowup(b):
    # a->v, v->B, B->a.
    return (2,((1<<len(b))-1)<<2)+tuple((r<<2)|1 for r in b)

def require(b,msg):
    if not b: raise RuntimeError(msg)

result={"gadget_rows":B,"gadget_packet":packet(B)}
require(packet(B)==(29,580),'wrong witness')
require(literal(B)==29,'literal B')
child=cycle_blowup(B)
require(packet(child)[0]==609,'graft count')
require(literal(child)==609,'literal child')
result['child_rows']=child
result['child_H']=609
result['content']=29
result['normalized_output']=21
# Independent exhaustive controls only through four: sufficient for proof's boundary.
result['small_packets']={}
for n in range(1,5):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    values=set()
    for bits in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: i,j=j,i
            rows[i]|=1<<j
        values.add(packet(rows))
    result['small_packets'][str(n)]=sorted(values)
require(result['small_packets']['4']==[(1,16),(3,24),(5,32)],'four-vertex packets')
# Genuine finite controls for all-length homogeneous family: ordinal products.
result['tower_controls']=[]
for m,k in [(1,0),(1,1),(1,2),(2,0)]:
    b=()
    for _ in range(m): b=join(b,B)
    b=join(b,tuple(((1<<k)-1)^((1<<(i+1))-1) for i in range(k)))
    h=29**m; q=20**m*2**k
    require(packet(b)==(h,h*q),'product packet')
    result['tower_controls'].append({'m':m,'k':k,'H':h,'Q':h*q,'normalized_graft':1+q})
print(json.dumps(result,indent=2))
