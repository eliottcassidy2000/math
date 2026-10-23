#!/usr/bin/env python3
"""Exact finite checks. No floats, no assertion-dependent tests."""
import json
from collections import defaultdict
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
F=[0,1,1]
for k in range(3,26): F.append(F[-1]+F[-2])
def Z(n):
    if n<0: return None
    out=[]
    for i in range(len(F)-1,1,-1):
        if F[i]<=n: out.append(i);n-=F[i]
    require(n==0,'range');return tuple(sorted(out))
def r(n): return int(n>=0 and 2 in Z(n))
def fibre(n,retain):
    out=[Z(n)]
    if n>=1:
        out.append((-1,)+Z(n-1))
        if not r(n-1):out.append((1,)+Z(n-1))
    if retain and n>=2 and not r(n-2):out.append((-1,1)+Z(n-2))
    return sorted(out)
def legal(s,retain):
    ids=list(s) if retain else [i+1 if i<0 else i for i in s]
    # compressed indices: A=0, B=1,C=2,F3=3,...
    return all(b-a>1 for a,b in zip(ids,ids[1:]))
def val(i):return 1 if i==-1 else F[i]
slots=[-1]+list(range(1,15))
for retain in [False,True]:
    actual=defaultdict(list)
    def visit(j,s,n):
        if j==len(slots):actual[n].append(s);return
        visit(j+1,s,n)
        if legal(s+(slots[j],),retain):visit(j+1,s+(slots[j],),n+val(slots[j]))
    visit(0,(),0)
    for n in range(F[15]):require(sorted(actual[n])==fibre(n,retain),f'fibre {retain} {n}')
W={-1:'A',1:'B',2:'C'}
proper={}
for i in range(3,21):
    support=list(range(i-1,1,-2)); support.append(-1 if i%2 else 1);support=tuple(sorted(support))
    proper[i]=support
    require(len(fibre(F[i],True))-1==2,'proper retained count')
    require(len(fibre(F[i],False))-1==(1 if i%2 else 2),'proper deleted count')
    require(sum(val(j) for j in support)==F[i],'proper value')
    require(all(legal(support,t) for t in [False,True]),'proper legal')
    W[i]=''.join(W[j] for j in support)
    require(W[i].count('C')==F[i-2] and W[i].count('A')==F[i-2] and W[i].count('B')==F[i-3],'exact block counts')
    if i>=5:require(W[i]==W[i-2]+W[i-1],'block recurrence')
limit=W[20][::-1]
require(limit[:16]=='CABCACABCABCACAB','user prefix')
for n in range(1,F[20]):
    row=''.join(W[i] for i in Z(n))
    require(row[::-1]==limit[:n],f'prefix {n}')
    c='C' if min(Z(n))==2 else ('A' if min(Z(n))%2 else 'B')
    require(c==limit[n-1],f'lowest class {n}')
    require(len(fibre(n,False))==3-r(n-1),'compressed count')
    require(len(fibre(n,False))==3-int(c=='A'),'compressed color count')
    if n>=2:require(len(fibre(n,True))==3+int(limit[n-2]=='B'),'retained color count')
    if n>=2:require(len(fibre(n,True))==4-r(n-1)-r(n-2),'retained count')
require((-1,1) in fibre(2,True) and (-1,1) not in fibre(2,False),'zero hostile')
# Exact floor(k phi) via isqrt: irrational sqrt5 means no endpoint ambiguity.
from math import isqrt
def lower(k):return (k+isqrt(5*k*k))//2
def upper(k):return lower(k)+k
for k in range(1,1500):
    for label,n in [('C',upper(k)-1),('A',upper(k)),('B',2*lower(k)+k)]:
        if n<=len(limit):require(limit[n-1]==label,f'Beatty {label} {k}')
report={'status':'PASS','fibre_exhaustive_N':F[15]-1,'prefix_and_counts_N':F[20]-1,'prefix':limit[:80], 'row35_color':limit[34], 'row36_suffix':'canonical reverse row has same position35 C','zero_deletion_first_difference':{'n':2,'extra_support':['A','B']},'proper_expansion_nonuniqueness':{'F4':[['A','F3'],['B','F3']]}}
print(json.dumps(report,indent=2))
