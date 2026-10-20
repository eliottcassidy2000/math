#!/usr/bin/env python3
"""Independent exact two-swap label/cut/form audit; SymPy only for norms."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import gcd
import json
import sympy as sp

Z=sp.Symbol('Z')

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def ceildiv(n,d):
    return -(-n//d)

def swap_occurrences(positions,q):
    occupied={x%q:j for j,x in enumerate(positions)}
    for edge in range(q):
        i,j=edge,(edge+1)%q
        if (i in occupied)==(j in occupied):
            continue
        if i in occupied:
            yield edge,occupied[i],1
        else:
            yield edge,occupied[j],-1

def kind_of(d):
    changed=[j for j,x in enumerate(d) if x]
    if not changed:return 'undo'
    if len(changed)==1:
        check(abs(d[changed[0]])==2,'one-label displacement')
        return 'same_plus2' if d[changed[0]]==2 else 'same_minus2'
    check(len(changed)==2 and all(abs(d[j])==1 for j in changed),'two-label displacement')
    vals=sorted(d[j] for j in changed)
    return {(-1,-1):'distinct_minus',(-1,1):'distinct_mixed',(1,1):'distinct_plus'}[tuple(vals)]

def forms(p,residues,d):
    raw=1
    for a,move in zip(residues,d):
        raw+=(Z-1)*(sp.Rational(2)**move-1)*Z**a
    ell=max(0,-min(d))
    cleared=sp.expand(2**ell*raw)
    changed=[(residues[j],x) for j,x in enumerate(d) if x]
    kind=kind_of(d)
    if kind=='undo':canonical=sp.Integer(1)
    elif kind=='same_plus2':
        a=changed[0][0];canonical=1+3*(Z-1)*Z**a
    elif kind=='same_minus2':
        a=changed[0][0];canonical=2*Z**(p-a)-3*Z+3
    elif kind=='distinct_plus':
        a,b=sorted(x[0] for x in changed);canonical=1+(Z-1)*(Z**a+Z**b)
    elif kind=='distinct_minus':
        a,b=sorted(x[0] for x in changed);canonical=Z**(p-a)-(Z-1)*(1+Z**(b-a))
    else:
        a=next(a for a,x in changed if x==1);b=next(a for a,x in changed if x==-1)
        canonical=Z**(p-b)*(1+(Z-1)*Z**a)-Z+1
    return sp.expand(cleared),sp.expand(canonical)

@lru_cache(None)
def free_norm(p,coefficients):
    poly=sum(c*Z**i for i,c in enumerate(coefficients))
    value=abs(int(sp.resultant(Z**p-2,poly,Z)))
    check(value>0,'zero norm')
    while value%2==0:value//=2
    while value%3==0:value//=3
    return value

def coefficient_key(poly):
    poly=sp.Poly(poly,Z)
    return tuple(int(poly.nth(i)) for i in range(poly.degree()+1))

def critical_legal(p,r,residues,d):
    kind=kind_of(d)
    changed=[(residues[j],x) for j,x in enumerate(d) if x]
    if kind=='undo':return True
    if kind.startswith('same_'):return False
    if kind=='distinct_plus':
        a,b=sorted(x[0] for x in changed)
        return b<r or (b>=r and b-a==r)
    if kind=='distinct_minus':
        a,b=sorted(x[0] for x in changed)
        return a>=p-r or (a<p-r and b-a==r)
    a=next(a for a,x in changed if x==1);b=next(a for a,x in changed if x==-1)
    return a<r and b>=p-r and a!=b and b!=a+p-r

def main():
    counts=Counter(); norm_types=set();cut_counts=Counter(); class_examples={}
    # Every ordered actual two-swap path, not just endpoints or necklaces.
    for p in range(2,11):
      for q in range(p+1,71):
        if gcd(p,q)!=1:continue
        k,r=divmod(q,p)
        source=[j*k+ceildiv(j*r,p)-1 for j in range(p)]
        residues=[(-j*r)%p for j in range(p)]
        eta=[ceildiv(j*r,p) for j in range(p)]
        for edge1,j1,e1 in swap_occurrences(source,q):
          middle=source.copy();middle[j1]+=e1
          for edge2,j2,e2 in swap_occurrences(middle,q):
            d=[0]*p;d[j1]+=e1;d[j2]+=e2
            target=[s+x for s,x in zip(source,d)]
            kind=kind_of(d);counts[kind]+=1
            if q-1 in (edge1,edge2):cut_counts[kind]+=1
            class_examples.setdefault(kind,{'p':p,'q':q,'edges':[edge1,edge2],'moves':d})
            gaps=[source[j]-source[j-1] if j else source[0]+q-source[-1] for j in range(p)]
            check(all(gaps[j]+d[j]-d[j-1]>=1 for j in range(p)),'cyclic label order')
            if kind.startswith('same_'):
                j=next(j for j,x in enumerate(d) if x)
                destination_gap=gaps[(j+1)%p] if d[j]>0 else gaps[j]
                check(destination_gap>=3,'same-label source gap')
            if k==1:
                check(not kind.startswith('same_'),'same-label move in k=1')
            if q==(3**p).bit_length():
                check(critical_legal(p,r,residues,d),'critical residue classifier')
            offsets=[(target[j]-target[0])%q for j in range(p)]
            check(offsets==sorted(offsets),'moving-anchor order')
            B=sum(3**(p-1-j)*2**offsets[j] for j in range(p))
            exponent=[eta[j]+d[j]-d[0] for j in range(p)]
            h=max(0,-min(exponent))
            Pvalue=sum(3**(p-1-j)*2**(h+exponent[j]+k*j) for j in range(p))
            check(Pvalue==2**h*B,'moving-anchor carry polynomial')
            # Exact algebraic reduction in Q[Z]/(Z^p-2), once per labelled type.
            key=(p,r,tuple(d))
            if key not in norm_types and p<=7 and q<=2*p+1:
                norm_types.add(key)
                cleared,canonical=forms(p,residues,d)
                m1=free_norm(p,coefficient_key(cleared));m2=free_norm(p,coefficient_key(canonical))
                check(m1==m2,'cleared/canonical 2,3-free norm mismatch')
                check(m2%gcd(abs(2**q-3**p),B)==0,'content norm gate')
                P=sum(3**(p-1-j)*2**(h+exponent[j])*Z**j for j in range(p))
                ordinary=abs(int(sp.resultant(2**r*Z**p-3**p,P,Z)))
                check(ordinary>0,'ordinary zero norm')
                sparse_raw=abs(int(sp.resultant(Z**p-2,cleared,Z)))
                ell=max(0,-min(d))
                exponent_two=r*(p-1)+p*(h-d[0]-ell)
                check(Fraction(ordinary)==Fraction(2)**exponent_two*3**(p*(p-1))*sparse_raw,
                    'exact moving-anchor resultant identity')
                while ordinary%2==0:ordinary//=2
                while ordinary%3==0:ordinary//=3
                check(ordinary==m2,'ordinary/sparse norm mismatch')
    # Converse of critical classifier: every predicted final labelled vector
    # is attained by at least one of its two orders, up to p=60.
    converse_counts=Counter()
    for p in range(3,61):
        q=(3**p).bit_length();k,r=divmod(q,p)
        if gcd(p,q)!=1:continue
        source=[j+ceildiv(j*r,p)-1 for j in range(p)]
        residues=[(-j*r)%p for j in range(p)]
        actual=set()
        for edge1,j1,e1 in swap_occurrences(source,q):
            middle=source.copy();middle[j1]+=e1
            for edge2,j2,e2 in swap_occurrences(middle,q):
                d=[0]*p;d[j1]+=e1;d[j2]+=e2;actual.add(tuple(d))
        for j in range(p):
          for l in range(j+1,p):
            for e in (-1,1):
              for f in (-1,1):
                d=[0]*p;d[j]=e;d[l]=f
                expected=critical_legal(p,r,residues,d)
                check(expected==(tuple(d) in actual),'critical converse failed')
                converse_counts['candidates']+=1
                converse_counts['legal']+=expected
    mixed_poly=Z**3+2*Z**2-3*Z+1
    mixed_norm=abs(int(sp.resultant(Z**4-2,mixed_poly,Z)))
    check(mixed_norm==343 and free_norm(4,coefficient_key(mixed_poly))==343,'mixed p4 barrier')
    mixed={'p':4,'q':7,'r':3,'plus_residue':2,'minus_residue':1,
        'ordered_edges':[0,3],'words':['0101011','1001011','1000111'],
        'norm':mixed_norm,'free_norm':343,'three_to_p':81,'B':331,'D':47,'g':1}
    valuation_poly=3-Z**12+Z**5-Z**4
    valuation_norm=abs(int(sp.resultant(Z**13-2,valuation_poly,Z)))
    check(valuation_norm==3**5*61*71,'three-adic stripping control')
    barriers=[]
    for p in (28,29,31):
        q=(3**p).bit_length();r=q-p;n=p-r;a=p-1;b=n-1
        poly=3-Z**(p-1)+Z**n-Z**(n-1)
        R=abs(int(sp.resultant(Z**p-2,poly,Z)));M=free_norm(p,coefficient_key(poly))
        check(M>3**p,'claimed barrier failed')
        label=(-pow(r,-1,p)*a)%p
        positions=[j+ceildiv(j*r,p)-1 for j in range(p)]
        leader=(label+1)%p
        edges=[positions[leader]%q,positions[label]%q]
        positions[leader]+=1;positions[label]+=1
        offsets=sorted(x%q for x in positions)
        B=sum(2**s*3**(p-1-j) for j,s in enumerate(offsets));D=2**q-3**p
        barriers.append({'p':p,'q':q,'r':r,'follower_residue':a,'leader_residue':b,
            'follower_label':label,'ordered_edges':edges,'norm':R,'free_norm':M,
            'three_to_p':3**p,'excess':M-3**p,'D':D,'g':gcd(D,B)})
    print(json.dumps({'status':'FINITE-EXACT','path_universe':'2<=p<=10; p<q<=70; gcd(p,q)=1; every ordered actual cyclic two-swap path',
        'path_counts':dict(sorted(counts.items())),'total_paths':sum(counts.values()),
        'original_cut_path_counts':dict(sorted(cut_counts.items())),
        'symbolic_type_count':len(norm_types),'critical_converse_universe':'3<=p<=60, q=bit_length(3^p), gcd(p,q)=1, all two-distinct-label signed vectors',
        'critical_converse_counts':dict(converse_counts),'norm_barriers':barriers,'mixed_small_barrier':mixed,'three_adic_control':{'p':13,'r':8,'norm':valuation_norm,'v3':5},
        'class_examples':class_examples},indent=2,sort_keys=True))

if __name__=='__main__':main()
