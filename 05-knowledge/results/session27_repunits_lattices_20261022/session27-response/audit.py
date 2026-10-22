from equivalence import *
from itertools import product
from collections import defaultdict
import json

def H(out):
    n=len(out)
    if not n:return 1
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for u,c in enumerate(d[s]):
            if not c:continue
            z=out[u]&~s
            while z:
                bit=z&-z;z-=bit;v=bit.bit_length()-1;d[s|bit][v]+=c
    return sum(d[-1])
def response(Q,w):
    k=len(Q);o=[sum(x<<j for j,x in enumerate(row)) for row in Q]+[0]*len(w)
    for i,a in enumerate(w):
        for u in range(k):
            if a>>u&1:o[u]|=1<<(k+i)
            else:o[k+i]|=1<<u
        for j in range(i+1,len(w)):o[k+i]|=1<<(k+j)
    return H(o)
def coefficients(Q):
    k=len(Q);h=0;s=[0]*k;e=[0]*k;c=[[0]*k for _ in range(k)]
    for p in permutations(range(k)):
        bad=[j for j in range(1,k) if not Q[p[j-1]][p[j]]]
        if not bad:h+=1;s[p[0]]+=1;e[p[-1]]+=1
        if len(bad)<=1:
            for j in range(1,k):
                if not bad or j==bad[0]:c[p[j-1]][p[j]]+=1
    return h,[e[u]-s[u]+sum(c[u]) for u in range(k)],c
count=0;groupsizes={}
for k in range(1,6):
    groups=defaultdict(list)
    for mask in range(1<<(k*(k-1)//2)):
        Q=core(k,mask);h,l,c=coefficients(Q);signature=[]
        for a in range(1<<k):
            z=h+sum(l[u] for u in range(k) if a>>u&1)-sum(c[u][v]+c[v][u] for u in range(k) for v in range(u+1,k) if (a>>u&1) and (a>>v&1))
            actual=response(Q,(a,))
            if z!=actual:raise RuntimeError(('polynomial mismatch',k,mask,a,z,actual))
            signature.append(actual);count+=1
        groups[tuple(signature)].append(mask)
    groupsizes[k]=len(groups)
print('one-base polynomial independent DP checks:',count,'groups:',groupsizes,flush=True)
checks=0
for k in range(4):
    ps,size,matrices,alpha=realization(k)
    for mask in range(1<<(k*(k-1)//2)):
        Q=core(k,mask);f=final(Q,ps,size)
        for m in range(4):
            for w in product(range(1<<k),repeat=m):
                v=alpha
                for a in w:v=advance(v,matrices[a])
                z=sum(x*f[s] for s,x in v.items())
                if z!=response(Q,w):raise RuntimeError(('automaton mismatch',k,mask,w))
                checks+=1
print('automaton independent DP checks:',checks,flush=True)
for k in range(3):
    for left in range(1<<(k*(k-1)//2)):
        for right in range(1<<(k*(k-1)//2)):
            Q=core(k,left);R=core(k,right);r=equivalent(Q,R)
            if r['equivalent']!=(left==right):raise RuntimeError(('equivalence mismatch',k,left,right,r))
            if not r['equivalent'] and response(Q,r['word'])-response(R,r['word'])!=r['difference']:raise RuntimeError('bad witness')
r=equivalent(core(3,2),core(3,5))
if r['word']!=(1,2) or r['difference']!=2:raise RuntimeError('C3 witness mismatch')
print('equivalence positive/negative controls PASS; C3:',r)
