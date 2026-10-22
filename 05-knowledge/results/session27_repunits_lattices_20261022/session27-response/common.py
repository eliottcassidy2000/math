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
