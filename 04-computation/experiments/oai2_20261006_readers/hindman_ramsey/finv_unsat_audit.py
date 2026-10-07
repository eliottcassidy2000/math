#!/usr/bin/env python3
"""Independent audit of a dumped CEGAR UNSAT certificate for the fully gap-determined game Q_inv(n,t).
Checks: (1) every negative clause is a realizable triangle {d1,d2,d1+d2} inside [t]^n;
(2) every positive clause equals the gap-class set of an explicitly listed binary subgrid of [t]^n
(leaf structure re-verified from scratch); (3) the CNF is UNSAT under other SAT solvers.
Usage: finv_unsat_audit.py n t file.cnf file.leaves"""
import sys, itertools
from pysat.formula import CNF
from pysat.solvers import Solver
n=int(sys.argv[1]); t=int(sys.argv[2]); cnf=CNF(from_file=sys.argv[3]); lvf=sys.argv[4]
# variable convention: lex-positive vectors with |coords|<=t-1, in itertools.product order, numbered from 1
vecs=[d for d in itertools.product(range(-(t-1),t),repeat=n) if next((x for x in d if x!=0),0)>0]
num={d:i+1 for i,d in enumerate(vecs)}
assert len(vecs)==((2*t-1)**n-1)//2 and cnf.nv<=len(vecs)
def cls(x,y):
    d=tuple(b-a for a,b in zip(x,y))
    if next((v for v in d if v!=0),0)<0: d=tuple(-v for v in d)
    return num[d]
def is_binary_subgrid(L):
    if len(L)!=2**n or len(set(L))!=2**n: return False
    if any(not(0<=c<t) for p in L for c in p): return False
    def rec(pts,depth):
        if depth==n: return len(pts)==1
        vals=sorted({p[depth] for p in pts})
        if len(vals)!=2: return False
        return all(rec([p for p in pts if p[depth]==v],depth+1) for v in vals)
    return rec(list(L),0)
leaf_sets={}
bad_leaves=0
for line in open(lvf):
    a,b=line.split('|')
    cl=frozenset(int(x) for x in a.split())
    L=[tuple(int(c) for c in tok.split(',')) for tok in b.split()]
    if not is_binary_subgrid(L): bad_leaves+=1; continue
    recomputed=frozenset(cls(L[i],L[j]) for i in range(len(L)) for j in range(i+1,len(L)))
    if recomputed!=cl: bad_leaves+=1; continue
    leaf_sets[cl]=True
ntri=npos=0; bad=0
for c in cnf.clauses:
    if all(l<0 for l in c):
        ntri+=1
        ds=[vecs[-l-1] for l in c]
        ok=False
        if len(ds)==3:
            for d1,d2,d3 in itertools.permutations(ds):
                if tuple(a+b for a,b in zip(d1,d2))==d3 and all(max(0,d1[i],d1[i]+d2[i])-min(0,d1[i],d1[i]+d2[i])<=t-1 for i in range(n)):
                    ok=True;break
        elif len(set(ds))<3:  # degenerate clause from d1==d2 (x, x+d, x+2d)
            for d1,d2,d3 in itertools.product(set(ds),repeat=3):
                if d1==d2 and tuple(2*a for a in d1)==d3 and set(ds)=={d1,d3} and all(max(0,d1[i],2*d1[i])-min(0,d1[i],2*d1[i])<=t-1 for i in range(n)):
                    ok=True;break
        if not ok: bad+=1
    elif all(l>0 for l in c):
        npos+=1
        if frozenset(c) not in leaf_sets: bad+=1
    else: bad+=1
print(f"clauses: {len(cnf.clauses)} (triangle {ntri}, subgrid {npos}); leaf records failing re-verification: {bad_leaves}; unjustified clauses: {bad}")
for name in sys.argv[5:] if len(sys.argv)>5 else ['minisat22','maplechrono','kissat404']:
    try:
        with Solver(name=name,bootstrap_with=cnf.clauses) as s:
            print(f"solver {name}: {'SAT' if s.solve() else 'UNSAT'}",flush=True)
    except Exception as e:
        print(f"solver {name}: unavailable ({e})")
