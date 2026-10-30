#!/usr/bin/env python3
"""Independent audit; stdlib only, no producer imports, no removable assertions."""
from itertools import permutations, combinations, product
from math import isqrt,gcd,prod
from collections import Counter
import json

def need(x,msg):
    if not x: raise RuntimeError(msg)
def sigma(n):
    return sum(d+n//d if d*d!=n else d for d in range(1,isqrt(n)+1) if n%d==0)
def sc(n):
    out=set(); p=2
    while p*p<=n:
        while n%p==0:
            n//=p
            if p in out: out.remove(p)
            else: out.add(p)
        p+=1
    if n>1: out.add(n)
    return frozenset(out)
def span(cols):
    s={frozenset()}
    for c in cols: s|={x^c for x in s}
    return s
def dot(x,y): return (x&y).bit_count()%2
def word(z,cols): return sum(dot(z,c)<<i for i,c in enumerate(cols))
def permword(w,p): return sum(((w>>i)&1)<<p[i] for i in range(len(p)))
B=(2,5,7,19,29,59,89)
M=tuple(sum(1<<i for i,q in enumerate((2,3,5)) if q in sc(sigma(p))) for p in B)
need(M==(2,3,1,4,7,6,5),'bank')
C={word(z,M) for z in range(8)}
K={w for w in range(128) if all(dot(w,c)==0 for c in C)}
C0={2*c for c in C}; K0={2*k+b for k in K for b in (0,1)}
need((len(K),min(w.bit_count() for w in K if w))==(16,3),'Hamming')
need((len(C),min(w.bit_count() for w in C if w))==(8,4),'simplex')
need(Counter(w.bit_count() for w in K0)=={0:1,1:1,3:7,4:14,5:7,7:1,8:1},'loop enumerator')
need(len(K0)==32 and len(C0)==8,'loop dimension')
A=C0|{w^255 for w in C0}
need(Counter(w.bit_count() for w in A)=={0:1,4:14,8:1},'RM enumerator')
need(all(dot(a,b)==0 for a in A for b in A),'RM duality')
need({w>>1 for w in A}==K,'puncture')
# Full groups counted via basis generators; linearity makes this sufficient.
gens=[word(1<<i,M) for i in range(3)]
auto7=[p for p in permutations(range(7)) if all(permword(w,p) in C for w in gens)]
need(len(auto7)==168,'GL full permutation group')
count0=countA=0
for p in permutations(range(8)):
    if all(permword(w,p) in C0 for w in [2*g for g in gens]): count0+=1
    if all(permword(w,p) in A for w in [2*g for g in gens]+[255]): countA+=1
need((count0,countA)==(168,1344),'GL versus AGL')
lines={frozenset(t) for t in combinations(range(7),3) if M[t[0]]^M[t[1]]^M[t[2]]==0}
flags={(i,L) for L in lines for i in L}
point=B.index(19); stab=[p for p in auto7 if p[point]==point]
need((len(lines),len(flags),len(stab))==(7,21,24),'Fano counts')
avoid=sorted((L for L in lines if point not in L),key=lambda L:sorted(L))
actions={tuple(avoid.index(frozenset(p[i] for i in L)) for L in avoid) for p in stab}
need(len(actions)==24,'faithful S4')
f=next(iter(flags))
need(len({(p[f[0]],frozenset(p[i] for i in f[1])) for p in auto7})==21,'flag orbit')
need(len(span([sc(n) for n in (3,5,19,4,6,20)]))==16,'joint rank 4')
need(sigma(5*19*29)==3600,'root example')
need([sigma(n) for n in (3,81,343)]==[4,121,400],'three prime powers')
need(sigma(243)==364 and sigma(9)==13,'same-port hostiles')
for u in range(1,100):
    if gcd(u,21)==1:
        need(all(sc(sigma(3**a*7**b*u))==sc(sigma(u)) for a in (0,1,4) for b in (0,3)),'neutral six states')
for p in (2,3,5,7,19):
    R=lambda a:sum(p**j for j in range(a+1))
    for a,b,c in product(range(5),repeat=3):
        need((sc(R(a))^sc(R(b)))^(sc(R(b))^sc(R(c)))==sc(R(a))^sc(R(c)),'telescoping')
        for u in range(1,20):
            if gcd(p,u)==1: need(sigma(p**b*u)*R(a)==sigma(p**a*u)*R(b),'exponent update')
Q=lambda n:sc(sigma(n))
D=lambda m,n: Q(m*n)^Q(m)^Q(n)
for l,m,n in product(range(1,16),repeat=3): need(D(l,m)^D(l*m,n)==D(m,n)^D(l,m*n),'cocycle')
for m,n in product(range(1,60),repeat=2):
    if gcd(m,n)==1:
        need(sigma(m*n)==sigma(m)*sigma(n),'coprime multiplication')
        s=lambda x:sigma(x)-x
        need(s(m*n)==s(m)*s(n)+m*s(n)+n*s(m),'aliquot defect')
need(sigma(15)-15==9 and sigma(81)-81==40,'aliquot hostiles')
scan=[]
for p in range(2,101):
    if all(p%d for d in range(2,isqrt(p)+1)):
        for a in range(1,13):
            s=sum(p**j for j in range(a+1))
            if isqrt(s)**2==s: scan.append([p,a])
need(scan==[[3,1],[3,4],[7,3]],'bounded scan, not theorem proof')
# Arithmetic parity completion: independent divisor-list construction and maps.
cols=(0,)+M
bank=(3,)+B
def syn(w):
    s=0
    for i,x in enumerate(cols):
        if w>>i&1: s^=x
    return s
def parity(w): return w.bit_count()%2
fibers=Counter((parity(w),syn(w)) for w in range(256))
need(len(fibers)==16 and set(fibers.values())=={16},'joint fibers')
need({w for w in K0 if parity(w)==0}==A,'arithmetic parity kernel')
for w in range(256):
    divisors=[1]
    for i,p in enumerate(bank):
        if w>>i&1: divisors += [d*p for d in divisors]
    sqclass=sc(sum(divisors))
    need(sqclass==frozenset(p for j,p in enumerate((2,3,5)) if syn(w)>>j&1),'explicit divisor squareclass')
    need(len(divisors)==2**w.bit_count(),'squarefree divisor count')
checks=0
for perm in auto7:
    L={0:0, **{M[i]:M[perm[i]] for i in range(7)}}
    for t in range(8):
        portperm=tuple(cols.index(L[x]^t) for x in cols)
        for w in range(256):
            new=permword(w,portperm)
            need(parity(new)==parity(w),'affine parity')
            need(syn(new)==L[syn(w)]^(t if parity(w) else 0),'affine joint transport')
            checks+=1
need(checks==344064,'transport check count')
need(syn(1)==0 and syn(1<<bank.index(7))==1,'odd translation hostile')
# Character values require common units despite equal rational squareclass.
need(pow(sigma(3)%11,5,11)==1 and pow(sigma(81)%11,5,11)==0,'raw Legendre support hostile')
print(json.dumps({'status':'PASS','independent':True,'parity_completion_checks':checks,'joint_fibers':len(fibers),'code_groups':[len(auto7),count0,countA], 'point_stabilizer':len(stab),'fano_lines':len(lines),'flags':len(flags),'prime_power_scan':scan,'classification_basis':'Bennett-Levin 2015 Proposition 1 plus exponent-one elementary argument'},indent=2))
