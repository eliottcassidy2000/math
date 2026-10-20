#!/usr/bin/env python3
"""Portable exact checks for closure, asymptotics and the anti-automorphism-test boundary.
Uses only Python standard library; prints JSON; reads or writes no repository files.
"""
from itertools import product, combinations, permutations
import json

def require(ok,message):
    if not ok: raise RuntimeError(message)

def is_anti(n,S,p):
    r=dict(zip(p,p[::-1]))
    return all(((v-u)%n in S)==((r[u]-r[v])%n in S) for u in range(n) for v in range(u+1,n))

def directed(n,S,p): return len(set(p))==n and all((v-u)%n in S for u,v in zip(p,p[1:]))

out={'source_snapshot':'70ca1d8210dfd7db3a32a2687f74a331dfccf910'}
n=9;S={1,3,4,7};p=[3,1,2,6,4,5,0,7,8]
require(directed(n,S,p) and is_anti(n,S,p),'countercriterion failed')
sums=[(p[i]+p[-1-i])%n for i in range(n)]
require(len(set(sums))>1,'criterion hostile lacks unequal sums')
out['false_constant_sum_criterion']={'n':n,'S':sorted(S),'path':p,'pair_sums':sums,'directed':True,'positional_reversal_is_anti':True}
checks=[]
for n in range(7,24,2):
    count=0;bad=[]
    for choice in product([0,1],repeat=(n-3)//2):
        S={1}|{d if c==0 else n-d for d,c in zip(range(2,(n+1)//2),choice)}
        small=[d for d in range(2,(n+1)//2) if d in S]
        paths=[]
        for a,b in combinations(sorted(S-{1}),2):
            if a+b<=n+1:
                x,y=a-1,b-1;p=list(range(x))+list(range(x+y,n))+list(range(x,x+y))
                paths.append(((a,b),p))
        if not small:
            p=[0,1]+[v for j in range(1,(n-3)//2+1) for v in [n-2*j,n-2*j+1]]+[2]
            paths.append((None,p))
        for pair,p in paths:
            require(directed(n,S,p),'incoming path not directed')
            count+=1
            if is_anti(n,S,p):bad.append({'S':sorted(S),'pair':pair,'path':p})
    require(not bad,'incoming explicit path counterexample')
    checks.append({'n':n,'actual_arc_checks':count,'false_witnesses':bad})
out['actual_arc_retest_of_incoming_paths']=checks
# Every labelled tournament, with direct path enumeration; canonicalize only fully closed ones.
closed=[]
for n in range(3,6):
    pairs=list(combinations(range(n),2));perms=list(permutations(range(n)));canon={}
    for mask in range(1<<len(pairs)):
        A=[[False]*n for _ in range(n)]
        for k,(u,v) in enumerate(pairs):
            A[u][v]=bool(mask>>k&1);A[v][u]=not A[u][v]
        paths=[p for p in perms if all(A[u][v] for u,v in zip(p,p[1:]))]
        require(bool(paths),'Redei control failed')
        closedpaths=sum(A[p[-1]][p[0]] for p in paths)
        if closedpaths!=len(paths):continue
        codes=[sum(A[q[u]][q[v]]<<k for k,(u,v) in enumerate(pairs)) for q in perms]
        code=min(codes);aut=codes.count(mask)
        canon[code]={'canonical_code':code,'H':len(paths),'HC':closedpaths//n,'Aut':aut,
                     'scores':sorted(sum(row) for row in A)}
    closed.append({'n':n,'labelled_universe':1<<len(pairs),'fully_closed_classes':list(canon.values())})
out['fully_closed_small_order_audit']=closed
require([len(x['fully_closed_classes']) for x in closed]==[1,0,1],'small closure classification failed')
T=[1]+[1<<((n-1)*(n-2)//2) for n in range(1,41)]
a=[0]*41;NS=[0]*41
for n in range(1,41):
    NS[n]=sum(a[k]*T[n-k] for k in range(1,n));a[n]=T[n]-NS[n]
out['nonstrong_asymptotic_controls']=[{'n':n,'NS_over_T_previous':NS[n]/T[n-1],'NS':str(NS[n]),'T_previous':str(T[n-1])} for n in [8,10,16,24,40]]
# Exact inequality controls, not a substitute for the all-order proof.
for n in range(4,41):
    require((2*T[n-1]-NS[n])*(1<<(n-3)) <= (n-2)*T[n-1], 'lower asymptotic bound failed')
    require((NS[n]-2*T[n-1])*(1<<(n-3)) <= (n-3)*T[n-1], 'upper asymptotic bound failed')
    for k in range(2,n-1):
        require(T[k]*T[n-k]*(1<<((k-1)*(n-k-1)))==T[n-1], 'middle-term identity failed')
out['exact_asymptotic_inequalities']={'orders':[4,40],'all_passed':True}
# Positive/hostile control: regular C5 is pure-blue in the companion proof,
# but it is not fully closed and is not the exceptional nonregular Q5.
n=5;S={1,2};ps=list(permutations(range(n)))
paths=[p for p in ps if directed(n,S,p)]
closedpaths=sum((p[0]-p[-1])%n in S for p in paths)
require(len(paths)==15 and closedpaths==10,'regular C5 control failed')
out['regular_C5_control']={'H':len(paths),'HC':closedpaths//n,'open_paths':len(paths)-closedpaths}
print(json.dumps(out,indent=2))
