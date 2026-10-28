#!/usr/bin/env python3
"""Exact prime -> canonical-tree cycle indices; exhaustive oracle through 6."""
from fractions import Fraction as Q
from collections import Counter
from itertools import permutations, combinations
from math import factorial, gcd
import json, sys
N=int(sys.argv[1]) if len(sys.argv)>1 else 9
if N<6: raise ValueError('N must be >=6')
# monomials are (power-sum partition, number of prime nodes, number of linear nodes)
ONE={((),0,0):Q(1)}; X={((1,),0,0):Q(1)}
def add(a,b,c=1):
    r=a.copy()
    for k,v in b.items():
        r[k]=r.get(k,0)+c*v
        if not r[k]: del r[k]
    return r
def mul(a,b):
    r={}
    for (p,i,j),v in a.items():
        for (q,k,l),w in b.items():
            if sum(p)+sum(q)<=N:
                key=(tuple(sorted(p+q)),i+k,j+l)
                r[key]=r.get(key,0)+v*w
    return {k:v for k,v in r.items() if v}
def mark(a,prime=0,linear=0): return {(p,i+prime,j+linear):v for (p,i,j),v in a.items()}
def adams(a,d,wrong=False):
    return {(tuple(d*x for x in p),i if wrong else d*i,j if wrong else d*j):v for (p,i,j),v in a.items() if d*sum(p)<=N}
def compose(p,a,wrong=False):
    r={}; cache={d:adams(a,d,wrong) for d in range(1,N+1)}
    for (part,i,j),v in p.items():
        t={((),i,j):v}
        for d in part: t=mul(t,cache[d])
        r=add(r,t)
    return r
def linear_closure(s,marked=True):
    a=s.copy(); term=s
    for k in range(2,N+1):
        term=mul(term,s)
        a=add(a,mark(term,linear=1) if marked else term)
    return a
def erase(a):
    r={}
    for (p,i,j),v in a.items(): r[(p,0,0)]=r.get((p,0,0),0)+v
    return r
def parts(n,start=1):
    if not n: yield ()
    for d in range(start,n+1,2):
        for p in parts(n-d,d): yield (d,)+p
def burnside(p):
    e=sum((d-1)//2 for d in p)+sum(gcd(d,f) for i,d in enumerate(p) for f in p[i+1:])
    z=1
    for d,m in Counter(p).items(): z*=d**m*factorial(m)
    return Q(2**e,z)
B={(p,0,0):burnside(p) for n in range(N+1) for p in parts(n)}
A=add(B,ONE,-1); inv=ONE; term=ONE
for k in range(1,N+1):
    term=mul(term,A); inv=add(inv,term,(-1)**k)
R=add(add(ONE,inv,-1),X,-1); P={}
for n in range(3,N+1):
    residual=add(R,compose(P,A),-1)
    P=add(P,{k:v for k,v in residual.items() if sum(k[0])==n})
# Forward fixed point, input only P (B is only used as an oracle below).
a={}
for iteration in range(N+1):
    s=add(X,mark(compose(P,a),prime=1)); new=linear_closure(s)
    if new==a: break
    a=new
else: raise RuntimeError('Forward iteration did not stabilize')
if erase(a)!=A: raise RuntimeError('Forward full cycle index differs from Burnside')
if any(2*b+l>sum(p)-1 or v<0 for (p,b,l),v in a.items()): raise RuntimeError('Tree support/positivity')
# Depth = maximum number of prime nodes on a root-to-leaf path.
depths=[]; ah=linear_closure(X)
for h in range((N-1)//2+1):
    depths.append(ah)
    ah=linear_closure(add(X,mark(compose(P,ah),prime=1)))
if depths[-1]!=a: raise RuntimeError('Depth saturation')

def cycle_type(p):
    seen=set(); lengths=[]
    for v in range(len(p)):
        if v in seen: continue
        w=v; length=0
        while w not in seen: seen.add(w); length+=1; w=p[w]
        lengths.append(length)
    return tuple(sorted(lengths))
def oracle(n):
    pairs=list(combinations(range(n),2)); index={p:i for i,p in enumerate(pairs)}
    transforms=[]
    for p in permutations(range(n)):
        transforms.append((p,[(index[min(p[i],p[j]),max(p[i],p[j])],p[i]>p[j]) for i,j in pairs]))
    seen=set(); bank={}; hist=Counter()
    for code in range(1<<len(pairs)):
        if code in seen: continue
        autos=[]
        for p,table in transforms:
            image=sum((((code>>k&1)^flip)<<dest) for k,(dest,flip) in enumerate(table))
            seen.add(image)
            if image==code: autos.append(p)
        def arc(i,j): return bool(code>>index[min(i,j),max(i,j)]&1)^(i>j)
        def tree(V):
            if len(V)==1: return 0,0,0
            mods=[]
            for size in range(1,len(V)+1):
                for tup in combinations(V,size):
                    M=frozenset(tup)
                    if all(len({arc(v,m) for m in M})==1 for v in V if v not in M): mods.append(M)
            strong=[M for M in mods if not any(M&D and M-D and D-M for D in mods)]
            proper=[M for M in strong if len(M)<len(V)]
            blocks=[M for M in proper if not any(M<D for D in proper)]
            reps=[next(iter(M)) for M in blocks]
            islinear=not any(arc(i,j) and arc(j,k) and arc(k,i) for i,j,k in permutations(reps,3))
            sub=[tree(tuple(sorted(M))) for M in blocks]
            return (sum(b for b,l,h in sub)+(not islinear),sum(l for b,l,h in sub)+islinear,max(h for b,l,h in sub)+(not islinear))
        b,l,h=tree(tuple(range(n))); hist[b,l,h]+=1
        for p in autos:
            key=(cycle_type(p),b,l); bank[key]=bank.get(key,0)+Q(1,len(autos))
    if len(seen)!=1<<len(pairs): raise RuntimeError('Incomplete exhaustive census')
    return bank,hist
oracle_hist={}
for n in range(1,7):
    bank,hist=oracle(n); oracle_hist[n]={','.join(map(str,k)):v for k,v in sorted(hist.items())}
    if bank!={k:v for k,v in a.items() if sum(k[0])==n}: raise RuntimeError(('Marked cycle index mismatch',n))
    for h,dh in enumerate(depths):
        expected=sum(v for (b,l,d),v in hist.items() if d<=h)
        actual=sum(v for (p,b,l),v in dh.items() if sum(p)==n)
        if actual!=expected: raise RuntimeError(('Depth census mismatch',n,h))
# Deliberately wrong weight substitution fails at the 3-cycle-of-3-cycles.
c3={((1,1,1),0,0):Q(1,3),((3,),0,0):Q(2,3)}
child=mark(c3,prime=1)
correct=mark(compose(c3,child),prime=1)
wrong=mark(compose(c3,child,wrong=True),prime=1)
if N>=9 and (sum(v for (p,b,l),v in correct.items() if b==4)!=1 or wrong==correct): raise RuntimeError('Hostile control failure')
def histogram(poly):
    r=Counter()
    for (p,b,l),v in poly.items(): r[sum(p),b,l]+=v
    if any(v.denominator!=1 for v in r.values()): raise RuntimeError('Nonintegral isotype count')
    return {','.join(map(str,k)):int(v) for k,v in sorted(r.items())}
report={'N':N,'checks':'PASS: full unmarked cycle index; full marked exhaustive indices through6; depth exhaustive through6; positivity/tree bound; integral histograms; bad-Adams hostile at9','tree_histogram_n_b_l':histogram(a),'exhaustive_n_b_l_depth':oracle_hist,'depth_cumulative_counts':{h:[int(sum(v for (p,b,l),v in dh.items() if sum(p)==n)) for n in range(N+1)] for h,dh in enumerate(depths)},'hostile_C3_C3_correct_prime_node_weights':{str(b):str(sum(v for (p,i,l),v in correct.items() if i==b)) for b in sorted({b for p,b,l in correct})},'hostile_C3_C3_wrong_prime_node_weights':{str(b):str(sum(v for (p,i,l),v in wrong.items() if i==b)) for b in sorted({b for p,b,l in wrong})}}
if '--full' in sys.argv:
    report['full_marked_cycle_index']=[{'partition':list(p),'prime_nodes':b,'linear_nodes':l,'coefficient':str(v)} for (p,b,l),v in sorted(a.items(),key=lambda kv:(sum(kv[0][0]),kv[0]))]
print(json.dumps(report,indent=2))
