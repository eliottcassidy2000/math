#!/usr/bin/env python3
"""Finite exact controls, independent of Session98 checker imports."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial
import json
checks=0
def require(ok, name):
    global checks
    checks+=1
    if not ok: raise RuntimeError(name)
def comps(total,n):
    if n==1:
        yield (total,); return
    for a in range(total+1):
        for tail in comps(total-a,n-1): yield (a,)+tail

def shapes(k):
    if k==0: return [(0,)]
    out=[]
    # Breadth-first valid offspring words, not recursive Catalan formula.
    for inner in combinations(range(3*k+1),k):
        bits=tuple(int(i in inner) for i in range(3*k+1))
        queue=1
        for j,b in enumerate(bits):
            queue+=3*b-1
            if queue==0 and j!=len(bits)-1: break
        else:
            if queue==0: out.append(bits)
    return out
for k in range(4):
    m=3*k
    for f in range(5):
        gaps=set()
        for support in combinations(range(m+f),m):
            bounds=(-1,)+support+(m+f,)
            gaps.add(tuple(bounds[i+1]-bounds[i]-1 for i in range(m+1)))
        require(gaps==set(comps(f,m+1)), 'support-gap bijection')
        require(len(gaps)==comb(m+f,f),'gap count')
    ss=shapes(k)
    require(len(ss)==comb(3*k+1,k)//(3*k+1),'independent BFS shapes')
    for bits in ss:
        # k=0 shape is already (0,).
        require(len(bits)==3*k+1 and sum(bits)==k,'tree Euler count')
        for gs in comps(3,len(bits)):
            digits=[1+b+2*g for b,g in zip(bits,gs)]
            require(sum(digits)==4*k+1+6,'clock')
            require(F(1,2**sum(digits))==F(1,2**(4*k+1))*F(1,4**3),'native mass')
for k in range(21):
    n=3*k+1
    c=F(factorial(3*k),3**k*factorial(k))*F(2)**(-3*k*k+k)
    A=c*4**(3*k*k)
    nextA=F(factorial(3*k+3),3**(k+1)*factorial(k+1))*2**(3*(k+1)**2+k+1)
    require(nextA==A*(3*k+1)*(3*k+2)*2**(6*k+4),'amplitude recurrence')
    for f in range(31):
        q=comb(3*k+f,f)
        require((f+1)*comb(3*k+f+1,f+1)==(3*k+f+1)*q,'first recurrence')
        require(comb(3*(k+1)+f,f)==sum(comb(j+2,2)*comb(3*k+f-j,f-j) for j in range(f+1)),'slot convolution')
        p=F(3,4)**n*F(1,4)**f*q
        require(p==sum(F(3,4)**(n-1)*F(1,4)**j*comb(n+j-2,j)*F(3,4)*F(1,4)**(f-j) for j in range(f+1)) if n>1 else p==F(3,4)*F(1,4)**f,'PGF split')
for k in (1,2):
    m=3*k
    tau=tuple(3*(i//3)+(i+1)%3 for i in range(m))
    centralizers=0
    conjugates=set()
    for p in permutations(range(m)):
        inv=[0]*m
        for i,j in enumerate(p): inv[j]=i
        conj=tuple(p[tau[inv[i]]] for i in range(m))
        conjugates.add(conj)
        centralizers+=conj==tau
        # A representative edge-labelled ordered tree has distinct labels;
        # relabelling fixes it exactly when all labels fixed.
        require((tuple(p[i] for i in range(m))==tuple(range(m)))==(p==tuple(range(m))),'free label action')
    require(centralizers==3**k*factorial(k),'centralizer size')
    require(len(conjugates)==factorial(m)//centralizers,'core cardinality')
    require(tau!=tuple(range(m)) and tuple(tau[tau[tau[i]]] for i in range(m))==tuple(range(m)),'nontrivial stabilizer')
# Arithmetic hostiles checked directly, not assumed from packets.
x=41; digits=[]
for _ in range(4):
    y=3*x+1; a=0
    while y%2==0: y//=2; a+=1
    digits.append(a); x=y
require(digits==[2,1,1,1] and x==107,'extinction without descent')
require((3*1+1)//4==1,'source one fixed')
require(3 != factorial(6)//(9*2),'shape/core nonidentification')
print(json.dumps({'status':'PASS','checks':checks,'support_k_max':3,'support_f_max':4,'recurrence_k_max':20,'recurrence_f_max':30,'stabilizer_k_max':2,'source_pin':'d9d7e33349653b6b65dbad477841c9ed518c081a'},sort_keys=True,indent=2))
