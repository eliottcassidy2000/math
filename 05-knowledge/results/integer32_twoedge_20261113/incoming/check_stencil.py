#!/usr/bin/env python3
"""Independent binary-algebra replay of incoming Session321 identities."""
from itertools import combinations
from pathlib import Path
import hashlib,json
checks=0
def need(value,label):
    global checks
    checks+=1
    if not value:raise RuntimeError(label)
def add(*ms):return tuple(xor(row) for row in zip(*ms))
def xor(xs):
    a=0
    for x in xs:a^=x
    return a
def multiply(a,b):return tuple(xor(b[k] for k in range(len(a)) if row>>k&1) for row in a)
def phi(a):return add(multiply(a,a),a,tuple(1<<i for i in range(len(a))))
def rowzero(a):return tuple(r^((r.bit_count()%2)<<i) for i,r in enumerate(a))
def matrix(n,code):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:rows[i]|=1<<j
        else:rows[j]|=1<<i
    return rowzero(rows)
def transpose(a):return tuple(sum(((a[j]>>i)&1)<<j for j in range(len(a))) for i in range(len(a)))
records=[]
for n in range(2,7):
    E=tuple(combinations(range(n),2));top=(1<<len(E))-1
    for code in range(top+1):
        p=phi(matrix(n,code));q=phi(matrix(n,top^code))
        need(p==q,'full labelled Phi descends through converse')
        need(transpose(p)==p,'symmetric Phi')
        need(all(p[i]>>i&1 for i in range(n)),'unit diagonal')
        need(all(r.bit_count()%2 for r in p),'row sums one')
    spine=sum(1<<k for k,(i,j) in enumerate(E) if j==i+1)
    free=[1<<k for k,(i,j) in enumerate(E) if j>i+1]
    L=rowzero(tuple(sum(1<<j for j in range(n) if abs(i-j)==1) for i in range(n)))
    for word in range(1<<len(free)):
        t=spine|sum(bit for k,bit in enumerate(free) if word>>k&1)
        a=matrix(n,t);b=matrix(n,t^top^spine)
        delta=add(phi(a),phi(b))
        need(delta==add(multiply(a,L),multiply(L,a),multiply(L,L),L),'native Laplacian identity')
        need(delta==transpose(delta),'symmetric difference')
        need(all(not (delta[i]>>i&1) for i in range(n)),'alternating difference')
        need(all(r.bit_count()%2==0 for r in delta),'ones in radical')
        for i,j in combinations(range(1,n-1),2):
            if j-i>=3:
                expected=((a[i]>>(j-1))^(a[i]>>(j+1))^(a[i-1]>>j)^(a[i+1]>>j))&1
                need((delta[i]>>j)&1==expected,'interior stencil')
    records.append({'n':n,'all_labelled_tournaments':top+1,'fixed_spine_words':1<<len(free)})
C=matrix(3,5);P=phi(C)
need(P==(7,7,7),'triangle named-kernel hostile')
need(C!=matrix(3,2),'orientation really changes M')
need(all((r&3).bit_count()%2==0 for r in P),'named vector110 survives both orientations')
print(json.dumps({'status':'FINITE-EXACT','checks':checks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'universe':records,'minimal_nonzero_kernel_hostile':{'n':3,'codes':[5,2],'Phi_rows':[7,7,7],'named_kernel_vector':3}},sort_keys=True,indent=2))
