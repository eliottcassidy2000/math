#!/usr/bin/env python3
"""Exact finite-state Collatz certificate controls; standard library only."""
from fractions import Fraction as F
from itertools import product
import json

def check(x, message):
    if not x: raise RuntimeError(message)

def shortcut(n): return n//2 if n%2==0 else (3*n+1)//2

def affine(word):
    a,b=F(1),F(0)
    for bit in word:
        a,b=(3*a/2,(3*b+1)/2) if bit else (a/2,b/2)
    return a,b

def parity(n,k):
    out=[]
    for _ in range(k): out.append(n%2); n=shortcut(n)
    return tuple(out)

def graph(width,cap):
    vertices=list(product((0,1),repeat=width-1))
    edges=[]
    for u in vertices:
        for b in (0,1):
            if sum(u)+b<=cap:
                edges.append((u,u[1:]+(b,),F(3 if b else 1,2),F(b,2)))
    return vertices,edges

def certificate(V,E,rho):
    # h(u) >= a_e h(v)/rho; each iteration adds one possible edge.
    h={u:F(1) for u in V}
    for k in range(len(V)):
        nxt={u:h[u] for u in V}
        for u,v,a,b in E: nxt[u]=max(nxt[u],a*h[v]/rho)
        if nxt==h: return h
        h=nxt
    return None

def matmul(A,B):
    n=len(A); C=[[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k]:
                for j in range(n):
                    if B[k][j]: C[i][j]+=A[i][k]*B[k][j]
    return C

report={}
# Independent arithmetic realization: all residues modulo 2^k give all words.
for k in range(1,13):
    words={parity(n,k) for n in range(1,2**k+1)}
    check(len(words)==2**k,'parity cylinder bijection')
report['parity_bijection_depths']=list(range(1,13))
# Inherited weighted de Bruijn rank-one identity, tested independently.
for L in range(1,6):
    V=list(product((0,1),repeat=L)); idx={u:i for i,u in enumerate(V)}
    for weights in [(F(1),F(1)),(F(1,2),F(3,2)),(F(1,2),F(4,3))]:
        n=len(V); K=[[F(0) for _ in V] for _ in V]
        for u in V:
            for b in (0,1): K[idx[u]][idx[u[1:]+(b,)]]+=weights[u[0]]/2
        P=[[F(i==j) for j in range(n)] for i in range(n)]
        for _ in range(L): P=matmul(P,K)
        for i,u in enumerate(V):
            expected=F(1,2**L)
            for bit in u: expected*=weights[bit]
            check(all(x==expected for x in P[i]),'weighted rank-one identity')
report['weighted_spectral_controls']='15 exact cases; L=1..5, three weight pairs'
# Five-window at most three odd steps: rational pathwise certificate.
V,E=graph(5,3); rho=F(99,100); h=certificate(V,E,rho)
check(h is not None,'positive certificate')
check(all(a*h[v]<=rho*h[u] for u,v,a,b in E),'edge inequalities')
B=max(b*h[v] for u,v,a,b in E); core=B/(1-rho)
report['window_cap_certificate']={'vertices':len(V),'edges':len(E),'rho':str(rho),'B':str(B),'weighted_core':str(core),'h':{''.join(map(str,k)):str(v) for k,v in h.items()}}
# Independent direct five-step inequality (stronger bound than generic compiler).
words=[w for w in product((0,1),repeat=5) if sum(w)<=3]
check(len(words)==26,'five-step universe')
check(all(affine(w)[0]<=F(27,32) and affine(w)[1]<=F(19,8) for w in words),'affine five-step bound')
check(affine((0,0,1,1,1))==(F(27,32),F(19,8)),'sharp affine control')
report['direct_bound']='T^5(n) <= (27*n+76)/32 < n for n >= 16, when the five-step word has <=3 ones'
# Hostile: no three consecutive odds still allows 110 repeated, slope 9/8.
Vbad,Ebad=graph(3,2)
check(certificate(Vbad,Ebad,F(1)) is None,'hostile must reject')
check(affine((1,1,0))==(F(9,8),F(5,8)),'expanding hostile')
shadow=[]
for k in range(1,9):
    word=(1,1,0)*k; modulus=2**len(word)
    # Construct matching residue by one-bit lifting instead of enumerating modulus.
    residue=0
    for j,bit in enumerate(word):
        if parity(residue,j+1)[j]!=bit: residue+=2**j
    n=residue+modulus*100
    check(parity(n,len(word))==word,'positive expanding shadow')
    a,b=affine(word); end=n
    for _ in word: end=shortcut(end)
    check(end==a*n+b and end>n,'expansion control')
    shadow.append({'blocks':k,'start':n,'end':end})
report['expanding_shadows']=shadow
print(json.dumps(report,indent=2,sort_keys=True))
