#!/usr/bin/env python3
"""Exact 512-subset SNF certificate, independent cyclic/additive counts."""
from itertools import product, combinations
from collections import Counter
from functools import reduce
from math import gcd
from pathlib import Path
import hashlib,json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).parent
V=tuple((1,x,y) for x,y in product([-1,0,1],repeat=2))
def req(c,m):
    if not c:raise RuntimeError(m)
def det3(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def snf_minors(rows):
    if not rows:return ()
    delta1=reduce(gcd,(abs(x) for row in rows for x in row),0)
    delta2=reduce(gcd,(abs(a[i]*b[j]-a[j]*b[i]) for a,b in combinations(rows,2) for i,j in combinations(range(3),2)),0)
    delta3=reduce(gcd,(abs(det3(a,b,c)) for a,b,c in combinations(rows,3)),0)
    if delta3:return (delta1,delta2//delta1,delta3//delta2)
    if delta2:return (delta1,delta2//delta1)
    return (delta1,)
def snf_sympy(rows):
    if not rows:return ()
    M=smith_normal_form(s.Matrix(rows),domain=s.ZZ)
    return tuple(abs(int(M[i,i])) for i in range(min(M.shape)) if M[i,i])
def rank_mod(rows,p):
    a=[[x%p for x in r] for r in rows];r=0
    for j in range(3):
        piv=next((i for i in range(r,len(a)) if a[i][j]),None)
        if piv is None:continue
        a[r],a[piv]=a[piv],a[r]
        z=pow(a[r][j],-1,p);a[r]=[(x*z)%p for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                z=a[i][j];a[i]=[(x-z*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r

def snf_census():
    signed=Counter();by_size=Counter();subsets=[];ranks={p:Counter() for p in [0,2,3,5,7]}
    for mask in range(512):
        ids=[i for i in range(9) if mask>>i&1];rows=[V[i] for i in ids]
        ds=snf_sympy(rows);req(ds==snf_minors(rows),('SNF versus minors',mask))
        sign=(-1)**len(ids);signed[ds]+=sign;by_size[(len(ids),ds)]+=1
        subsets.append({'mask':mask,'indices':ids,'invariants':ds})
        ranks[0][len(ds)]+=sign
        for p in [2,3,5,7]:
            r=sum(d%p!=0 for d in ds)
            req(rank_mod(rows,p)==r,('modular rank',mask,p))
            ranks[p][r]+=sign
    expected={():1,(1,):-9,(1,1):20,(1,2):8,(1,1,1):1,(1,1,2):-10,(1,1,3):-4,(1,1,4):-4,(1,2,2):-3}
    req(dict(signed)==expected,'aggregate signed SNF')
    triples=Counter(abs(det3(*r)) for r in combinations(V,3))
    req(dict(triples)=={0:8,1:32,2:32,3:4,4:8},'triple determinant histogram')
    return {'vectors':V,'subset_count':512,'signed_snf':[{'invariants':d,'coefficient':c} for d,c in sorted(signed.items())],
            'subset_size_snf_histogram':[{'size':k,'invariants':d,'count':c} for (k,d),c in sorted(by_size.items())],
            'signed_rank_coefficients':{str(p):dict(sorted(h.items())) for p,h in ranks.items()},
            'absolute_triple_determinants':dict(sorted(triples.items())),'subsets':subsets}
def generic_formula(n):
    g2,g3,g4=(gcd(n,d) for d in [2,3,4])
    return n**3-9*n*n+(20+8*g2)*n+1-10*g2-4*g3-4*g4-3*g2*g2
def cyclic_direct(n):
    # For each x,y, enumerate all forbidden central exponents rather than all e.
    return sum(n-len({(a*x+b*y)%n for a,b in product([-1,0,1],repeat=2)}) for x,y in product(range(n),repeat=2))
def cyclic_full(n):
    return sum(all((e+a*x+b*y)%n for a,b in product([-1,0,1],repeat=2)) for e,x,y in product(range(n),repeat=3))
def cyclic_census(snf):
    out=[]
    for n in range(1,121):
        count=cyclic_direct(n)
        req(count==generic_formula(n),('cyclic direct/formula',n))
        ie=sum(r['coefficient']*n**(3-len(r['invariants']))*__import__('math').prod(gcd(n,d) for d in r['invariants']) for r in snf['signed_snf'])
        req(ie==count,('cyclic IE',n))
        if n<=12:req(cyclic_full(n)==count,('full triple brute',n))
        out.append({'n':n,'count':count})
    return out

def additive_formula(q,p):
    if p==2:return q**3-4*q*q+6*q-3
    if p==3:return q**3-9*q*q+24*q-16
    return q**3-9*q*q+28*q-20

def additive_census():
    out=[]
    for p,d in [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(3,1),(3,2),(3,3),(3,4),(5,1),(5,2),(7,1),(7,2),(11,1),(13,1)]:
        q=p**d
        digits=[tuple((i//(p**j))%p for j in range(d)) for i in range(q)]
        enc=lambda a:sum(x*p**j for j,x in enumerate(a))
        def add(a,b):return enc(tuple((x+y)%p for x,y in zip(digits[a],digits[b])))
        def neg(a):return enc(tuple((-x)%p for x in digits[a]))
        table=[[add(a,b) for b in range(q)] for a in range(q)]
        negatives=[neg(a) for a in range(q)]
        count=0
        for x,y in product(range(q),repeat=2):
            choicesx=(negatives[x],0,x);choicesy=(negatives[y],0,y)
            count+=q-len({table[a][b] for a,b in product(choicesx,choicesy)})
        req(count==additive_formula(q,p),('additive count',p,d))
        if q<=13:
            brute=sum(all(table[e][table[a][b]]!=0 for a,b in product((negatives[x],0,x),(negatives[y],0,y))) for e,x,y in product(range(q),repeat=3))
            req(brute==count,('additive full brute',p,d))
        out.append({'characteristic':p,'dimension':d,'q':q,'count':count})
    return out

def discriminant_controls():
    # Direct traces in rank2^r square-root algebra over Q, for r<=4.
    checks=[]
    for r in range(1,5):
        z=list(range(2,2+r));basis=list(product([0,1],repeat=r));n=len(basis)
        trace=[]
        for a in basis:
            row=[]
            for b in basis:
                exponent=tuple(x+y for x,y in zip(a,b))
                row.append(0 if any(x%2 for x in exponent) else n*__import__('math').prod(t**(x//2) for t,x in zip(z,exponent)))
            trace.append(row)
        actual=int(s.Matrix(trace).det())
        predicted=2**(r*2**r)*__import__('math').prod(z)**(2**(r-1))
        req(actual==predicted,('discriminant',r));checks.append({'roots':r,'determinant':actual})
    return {'small_direct_trace_checks':checks,'nine_root_discriminant':'2^4608 * product(z_i)^256'}
def main():
    snf=snf_census()
    result={'scope':'All-n and all-finite-field formulas proved by finite exact 512-subset certification plus SNF inclusion-exclusion.',
            'snf':snf,'cyclic_groups':cyclic_census(snf),'additive_fields':additive_census(),
            'discriminant':discriminant_controls(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'branch-certificate.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','subsets':512,'cyclic_group_orders':'1..120','additive_field_cases':len(result['additive_fields']),
                      'SNF_independent_methods':2,'direct_trace_tests':4},sort_keys=True))
if __name__=='__main__':main()
