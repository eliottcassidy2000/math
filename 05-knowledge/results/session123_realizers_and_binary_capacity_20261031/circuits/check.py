#!/usr/bin/env python3
"""Finite hostile and exhaustive checks; no optimization-removable assertions."""
import itertools, json
from pathlib import Path

def require(p,msg):
    if not p: raise RuntimeError(msg)
def parity(x): return x.bit_count()&1

def basis(vs):
    piv={}
    for x in vs:
        while x:
            i=x.bit_length()-1
            if i in piv: x^=piv[i]
            else:
                piv[i]=x
                break
    return list(piv.values())

def kernel(rows,n):
    # RREF followed by explicit free-variable solutions, not enumeration.
    rows=[x for x in rows if x]; piv=[]; k=0
    for col in range(n):
        j=next((j for j in range(k,len(rows)) if rows[j]>>col&1),None)
        if j is None: continue
        rows[k],rows[j]=rows[j],rows[k]
        for j in range(len(rows)):
            if j!=k and rows[j]>>col&1: rows[j]^=rows[k]
        piv.append(col); k+=1
    free=[j for j in range(n) if j not in piv]; out=[]
    for j in free:
        x=1<<j
        for row,col in zip(rows,piv):
            if row>>j&1: x|=1<<col
        out.append(x)
    return out

def span(bs):
    s=[0]
    for b in bs: s += [x^b for x in s]
    return set(s)
def monomials(n): return [0]+[1<<i for i in range(n)]+[(1<<i)|(1<<j) for i in range(n) for j in range(i+1,n)]
def value(terms,x): return sum((x&m)==m for m in terms)&1

def algorithm(fs,n,p):
    m=n-p; rows=[]
    for ts in fs:
        pol=[0]*n
        for t in ts:
            if t.bit_count()==2:
                a,b=[i for i in range(n) if t>>i&1]
                pol[a]^=1<<b; pol[b]^=1<<a
        rows += [r>>p for r in pol]
    rb=kernel(rows,m)
    lr=[sum((value(ts,h<<p)^value(ts,0))<<j for j,h in enumerate(rb)) for ts in fs]
    kb=kernel(lr,len(rb))
    out=[]
    for k in kb:
        h=0
        for j,b in enumerate(rb):
            if k>>j&1: h^=b
        out.append(h)
    return out

def test_family(fs,n,p):
    kb=algorithm(fs,n,p); got=span(kb)
    tables=[[value(ts,x) for x in range(1<<n)] for ts in fs]
    expected={h for h in range(1<<(n-p)) if all(all(t[x]==t[x^(h<<p)] for x in range(1<<n)) for t in tables)}
    require(got==expected,('kernel mismatch',fs,n,p,got,expected))
    # Build chart with kernel directions last, check invertibility and unused axes.
    complement=[]
    for i in range(n-p):
        if len(basis(kb+complement+[1<<i]))>len(kb)+len(complement): complement.append(1<<i)
    chart=complement+kb
    require(len(basis(chart))==n-p,'chart singular')
    for h in chart[len(complement):]: require(h in expected,'chart not invariant')

# Named hostile controls, beyond scalar rank accounting.
require(len(algorithm([[3]],2,1))==0,'protected coupling bu')
require(len(algorithm([[3],[5]],3,0))==0,'shared uv uw')
require(len(algorithm([[3,5,6]],3,0))==0,'triangle derivative one')
c4=[3,9,6,12]
require(len(algorithm([c4],4,0))==2,'C4 kernel')
require(len(algorithm([c4,[4],[8]],4,0))==0,'C4 exposed coordinates')
scalar=vector=0
for n in range(1,5):
    mons=monomials(n)
    for mask in range(1<<len(mons)):
        ts=[a for i,a in enumerate(mons) if mask>>i&1]
        for p in range(n+1): test_family([ts],n,p); scalar+=1
n=3; mons=monomials(n)
funcs=[[a for i,a in enumerate(mons) if mask>>i&1] for mask in range(1<<len(mons))]
for a,b in itertools.combinations_with_replacement(funcs,2):
    for p in range(n+1): test_family([a,b],n,p); vector+=1

def forest(es,n):
    root=list(range(n))
    def find(a):
        while a!=root[a]: a=root[a]
        return a
    for a,b in es:
        a,b=find(a),find(b)
        if a==b:return False
        root[a]=b
    return True

def maxmatch(es):
    if not es:return 0
    a,b=es[0]
    return max(maxmatch(es[1:]),1+maxmatch([e for e in es[1:] if a not in e and b not in e]))

def pivots(es):
    es=list(es); result=[]
    while es:
        deg={v:sum(v in e for e in es) for e in es for v in e}
        x=next(v for v,d in deg.items() if d==1)
        y=next(b if a==x else a for a,b in es if x in (a,b))
        ns=[b if a==y else a for a,b in es if y in (a,b)]
        result.append((y,ns))
        es=[e for e in es if x not in e and y not in e]
    return result

forest_cases=assignments=0
for n in range(1,6):
    all_edges=list(itertools.combinations(range(n),2))
    for em in range(1,1<<len(all_edges)):
        es=[e for i,e in enumerate(all_edges) if em>>i&1]
        if not forest(es,n):continue
        for pm in range(1<<len(es)):
            protected=[e for i,e in enumerate(es) if pm>>i&1]
            unprotected=[e for i,e in enumerate(es) if not pm>>i&1]
            ps=pivots(unprotected); nu=maxmatch(unprotected)
            require(len(ps)==nu,'matching count')
            for delta in (0,1):
                # Literal gate list using original inputs, preserving protected outputs.
                gates=[]; outs=[]
                def gate(op,a,b):
                    gates.append((op,a,b)); return n+delta+len(gates)-1
                for a,b in protected: outs.append(gate('and',a,b))
                for y,ns in ps:
                    z=ns[0]
                    for q in ns[1:]: z=gate('xor',z,q)
                    outs.append(gate('and',y,z))
                if delta:outs.append(n)
                z=outs[0]
                for q in outs[1:]:z=gate('xor',z,q)
                require(len(gates)==len(es)+len(protected)+nu+delta-1,'gate formula')
                for x in range(1<<(n+delta)):
                    vals=[(x>>i)&1 for i in range(n+delta)]
                    for op,a,b in gates:vals.append(vals[a]&vals[b] if op=='and' else vals[a]^vals[b])
                    want=(sum(((x>>a)&1)*((x>>b)&1) for a,b in es)+(delta*((x>>n)&1)))&1
                    require(vals[z]==want,'forest value')
                    for i,(a,b) in enumerate(protected):require(vals[outs[i]]==((x>>a)&1)*((x>>b)&1),'protected output')
                    assignments+=1
                forest_cases+=1
report={'status':'PASS','scalar_quadratic_boundary_cases':scalar,'unordered_two_output_quadratic_boundary_cases':vector,'protected_forest_rewrite_cases':forest_cases,'forest_assignment_checks':assignments,'largest_quadratic_n':4,'largest_forest_n':5,'proof_scope':'finite checks supplement general proofs; no unrestricted lower bound'}
print(json.dumps(report,indent=2,sort_keys=True))
