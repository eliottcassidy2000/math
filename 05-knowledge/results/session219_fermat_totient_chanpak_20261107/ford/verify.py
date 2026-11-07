from itertools import permutations, combinations_with_replacement
from collections import Counter
from math import prod
import json

def check(ok, label):
    if not ok: raise RuntimeError(label)

def normal(es):
    out=[]
    for e in es: out.extend((2,2,3) if e==6 else (e,))
    return tuple(sorted(out))

def value(es): return prod((1<<e)-1 for e in es)

def fibers(n, lo=2):
    if n==1: return [()]
    out=[]
    e=lo
    while (1<<e)-1<=n:
        m=(1<<e)-1
        if n%m==0:
            out.extend((e,)+tail for tail in fibers(n//m,e))
        e+=1
    return out

def source(w):
    D,E,U,d=1,0,2,1
    for c in w:
        if c=='O': D,E,U,d=D,0,2*U,d
        elif c=='Z': D,E,U,d=D-E+U,d,2*U-d,d
        else: D,E,U,d=D-E+U,D-E,D+U,D
    return D,E,U,d

def F(s):return s[0]-s[1]+s[2]

def expected(a,b):
    if b==0:return a
    if b==1:return a-1 if a>=2 else None
    return a+b//2

report={}
report['star_matchings']={}
for r in range(1,9):
    ps=[p for p in permutations(range(r)) if all(i==0 or j==0 or i==j for i,j in enumerate(p))]
    check(len(ps)==r, 'star')
    report['star_matchings'][r]=len(ps)
seen={};count=0
for l in range(6):
    for es in combinations_with_replacement(range(2,19),l):
        n=value(es);nf=normal(es);count+=1
        check(n==value(nf),'normal value')
        check(n not in seen or seen[n]==nf,'normal uniqueness')
        seen[n]=nf
report['bounded_multisets']=count
report['power63_fibers']={}
for k in range(6):
    fs=fibers(63**k)
    check(len(fs)==k+1,'63 fiber')
    check({normal(s) for s in fs}=={(2,)*(2*k)+(3,)*k},'63 normal')
    report['power63_fibers'][k]=[list(s) for s in fs]
for a in range(1,21):
    for b in range(1,21):
        check(F(source('O'*a+'Z'*b))==2+((1<<(a+1))-1)*((1<<(b+1))-1),'two run')
report['physical_witnesses']={}
for h,w in [(3,''),(49,'OJJOO'),(63,'OZZOO'),(343,'OOJOZJJO')]:
    s=source(w);check(F(s)==h,'physical witness')
    report['physical_witnesses'][h]={'word':w,'profile':s,'closing_order':len(w)+3}
# Independent max-plus DP over the certified factor dictionary only.
cap=14;dp={(0,0):0}
for a in range(cap+1):
    for b in range(cap+1):
        if (a,b)==(0,0):continue
        vals=[dp[(a-u,b-v)]+1 for u,v in [(1,0),(0,2),(2,1),(0,3)] if a>=u and b>=v and (a-u,b-v) in dp]
        if vals:dp[(a,b)]=max(vals)
        check(dp.get((a,b))==expected(a,b),f'degree {a,b}')
report['degree_grid']={'max_a':cap,'max_b':cap,'checked':(cap+1)**2,'dictionary':[3,49,63,343]}
report['status']='PASS'
print(json.dumps(report,indent=2))
