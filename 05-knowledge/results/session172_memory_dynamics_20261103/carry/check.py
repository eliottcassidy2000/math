#!/usr/bin/env python3
import itertools, json
from fractions import Fraction
checks=0
def require(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)
def mul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def trace(A): return A[0][0]+A[1][1]
def inv(A):
    a,b=A[0]; c,d=A[1]; det=a*d-b*c
    return ((Fraction(d,det),Fraction(-b,det)),(Fraction(-c,det),Fraction(a,det)))
I=((1,0),(0,1)); J=((0,-1),(1,0))
words={}; count=0
for length in range(1,7):
    for w in itertools.product(range(1,5),repeat=length):
        G=I
        for e in w: G=mul(((3,1),(0,2**e)),G)
        a,b=G[0]; d=G[1][1]
        require((a,d,G[1][0])==(3**length,2**sum(w),0),'word diagonal')
        require(trace(mul(J,G))==b,'J carry')
        key=(a,d,b)
        require(key not in words,'finite word faithfulness')
        words[key]=w; count+=1
for a,d in itertools.product(range(1,5),repeat=2):
    for b in range(-3,4):
        G=((a,b),(0,d)); base=((a,0),(0,d))
        for u,v,w in itertools.product(range(-1,2),repeat=3):
            L=((u,v),(0,w)); R=((w,u),(0,v))
            require(trace(mul(mul(L,G),R))==trace(mul(mul(L,base),R)),'two-sided blind')
            require(trace(mul(mul(mul(G,L),inv(G)),R))==trace(mul(mul(mul(base,L),inv(base)),R)),'inverse blind')
        for p in itertools.product(range(-1,2),repeat=4):
            P=(p[:2],p[2:]); value=trace(mul(P,G))
            require(value==p[0]*a+p[2]*b+p[3]*d,'probe formula')
            if p[2]: require(Fraction(value-p[0]*a-p[3]*d,p[2])==b,'probe inverse')
def step(n):
    m=3*n+1; e=0
    while m%2==0: m//=2; e+=1
    return m,e
for n in range(1,2048,2):
    x,e=step(n); y,f=step(x)
    require(((e,f)==(1,3))==(n%32==19),'guard13')
    require(((e,f)==(3,1))==(n%32==29),'guard31')
    if (e,f)==(1,3): require(y*16==9*n+5,'carry13')
    if (e,f)==(3,1): require(y*16==9*n+11,'carry31')
require(step(19)==(29,1) and step(29)==(11,3) and step(11)==(17,1),'literal receipt')
# All multilinear degree <=2 Boolean polynomials on F2^4.
monomials=[m for m in range(16) if m.bit_count()<=2]
code=[]
for coeff in range(1<<len(monomials)):
    code.append(tuple(sum((coeff>>j)&1 for j,m in enumerate(monomials) if (y&m)==m)%2 for y in range(16)))
require(len(set(code))==2048,'RM injective')
banks={'bulk':[],'weight_le_one':[y for y in range(16) if y.bit_count()<=1], 'weight_le_two':monomials,'all':list(range(16))}
counts={}
for name,bank in banks.items():
    signatures=sorted(set(tuple(h[y] for y in bank) for h in code))
    counts[name]=len(signatures)
    lookup={s:j for j,s in enumerate(signatures)}
    # Point-mass class encoding and deterministic channel columns.
    for h in code:
        s=tuple(h[y] for y in bank); j=lookup[s]
        for k,y in enumerate(bank): require(signatures[j][k]==h[y],'class-channel exactness')
require(counts=={'bulk':1,'weight_le_one':32,'weight_le_two':2048,'all':2048},'RM class counts')
# Small positive mixtures: a deterministic output cannot average incompatible columns.
for bit1,bit2 in itertools.product((0,1),repeat=2):
    for t in (Fraction(1,4),Fraction(1,2),Fraction(3,4)):
        out=t*bit1+(1-t)*bit2
        require((out in (0,1))==(bit1==bit2),'extremality control')
print(json.dumps({'status':'pass','checks':checks,'word_count':count,'rm_dimension':len(monomials),'response_classes':counts,'scope':'finite transcription controls; general claims proved in NOTE.md'},sort_keys=True,indent=2))
