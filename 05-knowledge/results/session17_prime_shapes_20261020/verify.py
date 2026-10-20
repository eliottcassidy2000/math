#!/usr/bin/env python3
"""Exact independent controls for the oriented tetrahedral complement-pair map."""
from itertools import combinations, combinations_with_replacement, permutations
from collections import Counter
from math import prod, isqrt
import json

V=tuple(range(4)); E=tuple(combinations(V,2)); ALL=frozenset(E)
PERMS=tuple(permutations(V)); checks=Counter()
def require(q,tag):
    checks[tag]+=1
    if not q: raise RuntimeError(tag)
def parity(p): return sum(p[i]>p[j] for i in V for j in V if i<j)%2
def edge(a,b): return tuple(sorted((a,b)))
def acted(g,p): return frozenset(edge(p[a],p[b]) for a,b in g)
def cp(g): return tuple(sorted((tuple(sorted(g)),tuple(sorted(ALL-g)))))
def actcp(z,p): return cp(acted(frozenset(z[0]),p))
def degree(g): return tuple(sum(v in e for e in g) for v in V)
def shape(g):
    ds=sorted(degree(g))
    return {(1,1,1,3):'star',(0,2,2,2):'triangle_isolate',(1,1,2,2):'path'}[tuple(ds)]
def pathorders(g):
    return [p for p in PERMS if all(edge(p[i],p[i+1]) in g for i in range(3))]
def from_graph(g):
    ds=degree(g)
    if 3 in ds: return (ds.index(3),)*2
    if 0 in ds: return (ds.index(0),)*2
    p=pathorders(g)[0]
    return edge(p[0],p[3]) if parity(p)==0 else edge(p[1],p[2])
def fmap(z): return from_graph(frozenset(z[0]))
def jpair(a): return a if a[0]==a[1] else tuple(v for v in V if v not in a)
def actpair(a,p): return edge(p[a[0]],p[a[1]])
G=tuple(frozenset(c) for c in combinations(E,3))
CPS=tuple(sorted({cp(g) for g in G})); PAIRS=tuple(combinations_with_replacement(V,2))
require(Counter(shape(g) for g in G)=={'path':12,'star':4,'triangle_isolate':4},'shape_census')
require(len(CPS)==10,'complement_pair_count')
for g in G:
    require(from_graph(g)==from_graph(ALL-g),'representative_independence')
    if shape(g)=='path':
        orders=pathorders(g)
        require(len(orders)==2 and orders[0]==orders[1][::-1],'path_reversal')
        require(parity(orders[0])==parity(orders[1]),'reversal_parity')
require(set(map(fmap,CPS))==set(PAIRS),'bijection')
for p in PERMS:
    for z in CPS:
        lhs=fmap(actcp(z,p)); rhs=actpair(fmap(z),p)
        if parity(p): rhs=jpair(rhs)
        require(lhs==rhs,'full_twisted_equivariance')
fixed=[]
for p in PERMS:
    if sorted(Counter(p[i] for i in V).values()) != [1]*4: raise RuntimeError('not permutation')
    a=sum(actcp(z,p)==z for z in CPS)
    b=sum(actpair(z,p)==z for z in PAIRS)
    fixed.append({'permutation':p,'sign':1 if parity(p)==0 else -1,'complement_pairs_fixed':a,'simplex_pairs_fixed':b})
trans=(1,0,2,3)
require(sum(actcp(z,trans)==z for z in CPS)==2,'transposition_obstruction_source')
require(sum(actpair(z,trans)==z for z in PAIRS)==4,'transposition_obstruction_target')
# Stabilizer cycle structure for path-pair versus two-subset.
z=cp(frozenset(((0,1),(1,2),(2,3))))
stab=[p for p in PERMS if actcp(z,p)==z]
def compose(p,q): return tuple(p[q[i]] for i in V)
def order(p):
    q=V
    for k in range(1,25):
        q=compose(p,q)
        if q==V:return k
    raise RuntimeError('order')
require(sorted(map(order,stab))==[1,2,4,4],'cyclic_path_pair_stabilizer')
estab=[p for p in PERMS if actpair((0,1),p)==(0,1)]
require(sorted(map(order,estab))==[1,2,2,2],'klein_edge_stabilizer')
# Inherited divisor table, numeric example 2^2*3*5=60.
dtable={2:(0,1),3:(0,2),5:(0,3),6:(1,2),10:(1,3),15:(2,3),30:(0,0),4:(1,1),12:(2,2),20:(3,3)}
require(set(dtable)=={d for d in range(2,60) if 60%d==0},'divisor_table_domain')
require(set(dtable.values())==set(PAIRS),'divisor_table_range')
inverse={v:k for k,v in dtable.items()}
rows=[{'graph_pair':[list(g) for g in z],'simplex_pair':fmap(z),'divisor_of_60':inverse[fmap(z)]} for z in CPS]
# A different marked table preserves the arithmetic prime-to-extra involution.
def jgraph(g): return frozenset(tuple(v for v in V if v not in e) for e in g)
def jcp(z): return cp(jgraph(frozenset(z[0])))
btable={2:(0,1),3:(0,2),5:(0,3),4:(2,3),12:(1,3),20:(1,2),30:(0,0),15:(1,1),10:(2,2),6:(3,3)}
tau={2:4,4:2,3:12,12:3,5:20,20:5,6:6,10:10,15:15,30:30}
require(set(btable)==set(dtable),'balanced_table_domain')
require(set(btable.values())==set(PAIRS),'balanced_table_range')
for d in btable:require(btable[tau[d]]==jpair(btable[d]),'arithmetic_involution_conjugacy')
for z in CPS:require(fmap(jcp(z))==jpair(fmap(z)),'opposition_conjugacy')
jfixed=[z for z in CPS if jcp(z)==z]
jorbits={tuple(sorted((z,jcp(z)))) for z in CPS}
require(len(jfixed)==4 and len(jorbits)==7,'graph_7_plus_3_fibres')
for z in jfixed:require(shape(frozenset(z[0]))!='path','graph_fixed_types')
# Arithmetic checked both from direct divisors and exponent profiles.
def factor(n):
    out={}; d=2
    while d*d<=n:
        while n%d==0:out[d]=out.get(d,0)+1;n//=d
        d+=1
    if n>1:out[n]=out.get(n,0)+1
    return out
def counts(n):
    fs=factor(n);r=len(fs)
    f=prod(a+1 for a in fs.values())-2
    s=2**r-1-int(all(a==1 for a in fs.values()))
    u=r-int(len(fs)==1 and list(fs.values())==[1])
    div=[d for d in range(2,n) if n%d==0]
    direct=(len(div),sum(all(a==1 for a in factor(d).values()) for d in div),sum(len(factor(d))==1 and list(factor(d).values())==[1] for d in div))
    require((f,s,u)==direct,'direct_divisor_check')
    return {'n':n,'factorization':fs,'F':f,'S':s,'U':u,'defect':f-s-u}
arithmetic=[counts(n) for n in (5,25,125,625,776,576,60,63,27,75)]
require(625+125+25+1==776,'sum_776')
require((25-1)**2==625-50+1==576,'punctured_square_576')
require(25**2-25==600,'diagonal_deletion_600')
for p in range(2,101):
    require(p**4+p**3+p**2+1-(p*p-1)**2==p*p*(p+3),'polynomial_difference')
for n in range(2,1501):
    fs=factor(n);r=len(fs); exps=sorted(fs.values())
    f=prod(a+1 for a in exps)-2;s=2**r-1-int(all(a==1 for a in exps));u=r-int(exps==[1])
    require((f==s+u)==(exps in ([1],[3],[1,1,2])),'classification_window')
    if isqrt(n)**2==n:require(f!=s+u,'square_exclusion')
# No prime-label or support semantics inferred from integer sizes.
print(json.dumps({'status':'PASS','scope':{'graphs':'all 20 labelled three-edge graphs on four vertices','relabellings':'all 24 permutations and all 10 complement pairs','divisor_window':'2..1500; direct checks for 10 named inputs','polynomial_window':'integers p=2..100'},'check_counts':dict(checks),'total_checks':sum(checks.values()),'arithmetic':arithmetic,'fixed_point_rows':fixed,'map_rows':rows,'prime_fibre_preserving_table':btable},indent=2,sort_keys=True))
