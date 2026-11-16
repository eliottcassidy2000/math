#!/usr/bin/env python3
"""Replay score-preserving recontacts and target-selected terminal slots."""
def need(x,m):
    if not x: raise ValueError(m)

def counts(a,root=None):
    n=len(a); F=(1<<n)-1; H=[1]+[0]*F
    d=[[0]*n for _ in range(F+1)]
    for S in range(1,F+1):
        for v in range(n):
            if not S>>v&1:continue
            T=S^(1<<v)
            d[S][v]=sum(d[T][w] for w in range(n) if T>>w&1 and a[v]>>w&1) if T else 1
        H[S]=sum(d[S])
    if root is None:return H
    E=O=0
    for S in range(1<<root,F+1):
        z=d[S][root]*H[F^S]
        if S.bit_count()%2:O+=z
        else:E+=z
    return H[F],E,O,H[F^(1<<root)]

def validate(a):
    need(all(not((x>>i)&1) and x<(1<<len(a)) for i,x in enumerate(a)),"row");
    need(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(len(a)) for j in range(i)),"edge")

def tri(a,ids):
    need(all(sum(a[i]>>j&1 for j in ids)==1 for i in ids),"not a directed triangle")
    return tuple(x^sum(1<<j for j in ids if i!=j) if i in ids else x for i,x in enumerate(a))

def blow(a,slot,root,M):
    v=[(i,k) for i in range(len(a)) for k in range(len(M) if i==slot else 1)]
    return tuple(sum(1<<j for j,(y,l) in enumerate(v)
        if (a[x]>>y&1) or (x==y==slot and M[k]>>l&1)) for x,k in v),v.index((root,0))

families=[
 ((4,17,82,39,105,135,171,31),2,[(2,3,4),(2,5,6)],[(283,158,187,29),(291,162,191,29),(283,158,191,33)]),
 ((8,21,73,194,45,15,179,55),3,[(1,2,3)],[(309,166,193,27),(309,166,189,23)])
]
for A,root,path,profiles in families:
    a=A
    for i,w in enumerate(profiles):
        validate(a)
        need(counts(a,root)==w,"triangle path profile")
        need([x.bit_count() for x in a]==[x.bit_count() for x in A],"degrees changed")
        need(a[0]==1<<root,"slot exit changed")
        if i<len(path):a=tri(a,path[i])
    print(A,"through",path,"to",a)
    for M in ((0,),(2,0),(6,4,0),(2,4,1)):
        h=counts(M); F=len(h)-1; u=h[-1]; R=sum(x*h[F^j] for j,x in enumerate(h))
        for b,p in ((A,profiles[0]),(a,profiles[-1])):
            H,E,O,g=p
            need(counts(*blow(b,0,root,M))==(E*R+(H-2*E)*u,E*R//2,E*R//2+g*u,g*u),"inflation law")

top=[
 (16374,32764,32760,16369,32736,32704,16256,32512,32256,31744,30720,12288,8192,16384,6217),
 (32750,65524,65528,65522,32737,65472,65408,32512,65024,64512,63488,61440,24576,16384,32768,12433),
 (16382,32764,16280,16304,16352,32708,32652,32512,32256,31744,30720,12288,8192,16384,6173),
 (16382,32764,16280,16304,32736,16324,32652,32512,32256,31744,30720,12288,8192,16384,6189)
]
for a,p in zip(top,((8053,6479,6484,5),(24159,19437,19452,15),(24161,19437,19452,15),(24163,19437,19452,15))):
    validate(a); root=len(a)-1
    need(a[root-1]==1<<root and counts(a,root)==p,"fifth height actual profile")
print("Verified terminal-slot returns and same-score recontacts.")
