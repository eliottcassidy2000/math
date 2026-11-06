from itertools import permutations
import json

def edge(u):
    digits=bin(u-1)[2:];n=len(digits)+2;rows=[0]*n;rows[0]=2
    for j,b in enumerate(digits):
        v=j+2; rows[v]|=1
        if b=='1': rows[1]|=1<<v
        else: rows[v]|=2
        for w in range(v+1,n):rows[v]|=1<<w
    return rows

def add(a,mask):
    n=len(a);a=a[:]+[mask]
    for v in range(n):
        if not(mask>>v&1):a[v]|=1<<n
    return a

def count(a):
    n=len(a);full=(1<<n)-1;D=[[0]*n for _ in range(1<<n)]
    for v in range(n):D[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            if D[S][v]:
                choices=a[v]&~S&full
                while choices:
                    x=choices&-choices;choices-=x;w=x.bit_length()-1
                    D[S|x][w]+=D[S][v]
    return sum(D[-1])

def graph(k):
    a=edge((1<<k)+1);a[1]^=4;a[2]^=2
    return add(add(a,11),83)

def require(p,m):
    if not p:raise RuntimeError(m)
a=graph(7)
for i in range(len(a)):
    require(not (a[i]>>i&1),'loop')
    for j in range(i):require(((a[i]>>j&1)+(a[j]>>i&1))==1,'edge')
require(count(a)==775,'775 repair')
control=[{'k':k,'count':count(graph(k)),'target':6*((1<<k)+1)+1} for k in range(4,10)]
require([q['count'] for q in control]==[29,91,267,775,2259,6631],'hostile ray')
print(json.dumps({'status':'VERIFIED FINITE WITNESS ONLY','u':129,'target':775,'flip':[1,2],'birth_masks':[11,83],'rows':a,'controls':control},indent=2))
