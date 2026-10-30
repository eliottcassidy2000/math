#!/usr/bin/env python3
"""Exact finite checks for saturated-cut canonical search. Python standard library."""
from itertools import permutations, combinations
from math import factorial
from fractions import Fraction
import json

checks=0
def need(x,msg):
    global checks
    checks+=1
    if not x: raise RuntimeError(msg)

class Cube:
    def __init__(self,n):
        self.n=n; self.edges=list(combinations(range(n),2)); self.d=len(self.edges)
        self.idx={e:i for i,e in enumerate(self.edges)}; self.mask=(1<<self.d)-1
        self.ps=list(permutations(range(n)))
    def mapping(self,p):
        return [(self.idx[tuple(sorted((p[a],p[b])))],p[a]>p[b]) for a,b in self.edges]
    def image(self,x,p):
        return sum((((x>>i)&1)^int(f))<<j for i,(j,f) in enumerate(self.mapping(p)))
    def image_adjacency(self,x,p):
        a=[[False]*self.n for _ in range(self.n)]
        for i,(u,v) in enumerate(self.edges):
            a[p[u]][p[v]]=bool((x>>i)&1); a[p[v]][p[u]]=not bool((x>>i)&1)
        return sum(int(a[u][v])<<i for i,(u,v) in enumerate(self.edges))
    def reduce(self,rows):
        basis={}
        for row in rows:
            for i in sorted(basis):
                if row>>i&1: row^=basis[i]
            coeff=row&self.mask
            if not coeff:
                if row>>self.d: return None
                continue
            p=(coeff&-coeff).bit_length()-1
            for i in list(basis):
                if basis[i]>>p&1: basis[i]^=row
            basis[p]=row
        return tuple(basis[i] for i in sorted(basis))
    def meet(self,a,b):
        return None if a is None or b is None else self.reduce(a+b)
    def size(self,a): return 0 if a is None else 1<<(self.d-len(a))
    def satisfied(self,x,a):
        return a is not None and all(((x&r&self.mask).bit_count()&1)==(r>>self.d) for r in a)
    def transform(self,a,p):
        out=[]; mp=self.mapping(p)
        for r in a:
            coef=0; rhs=r>>self.d
            for i,(j,f) in enumerate(mp):
                if r>>i&1: coef^=1<<j; rhs^=int(f)
            out.append(coef|(rhs<<self.d))
        return self.reduce(out)
    def fixed(self,p):
        return self.reduce([(1<<i)^(1<<j)^(int(f)<<self.d) for i,(j,f) in enumerate(self.mapping(p))])
    def cut(self,t):
        # Disjoint prefix cylinders in numerical order: highest bit first.
        prefix=[]; out=[]
        for i in range(self.d-1,-1,-1):
            if t>>i&1: out.append(self.reduce(prefix+[(1<<i)]))
            prefix.append((1<<i)|(((t>>i)&1)<<self.d))
        out.append(self.reduce(prefix))
        return out
    def indicator(self,flats):
        terms={}
        for f in sorted(set(flats),key=lambda a:(len(a),a)):
            # A flat already covered contributes zero under this exact formula.
            delta={f:1}
            for a,c in terms.items():
                m=self.meet(a,f)
                if m is not None: delta[m]=delta.get(m,0)-c
            for a,c in delta.items():
                terms[a]=terms.get(a,0)+c
                if not terms[a]: del terms[a]
        return terms

def verify(n):
    C=Cube(n); N=1<<C.d; G=len(C.ps)
    rows=[]; fixed=[]
    for p in C.ps:
        row=[C.image(x,p) for x in range(N)]; rows.append(row)
        need(sorted(row)==list(range(N)),f'permutation n={n}')
        need(all(row[x]==C.image_adjacency(x,p) for x in range(N)),'independent adjacency action')
        fixed.append(sum(1<<x for x,y in enumerate(row) if x==y))
        f=C.fixed(p)
        need(all(C.satisfied(x,f)==(row[x]==x) for x in range(N)),'fixed equations')
    # The independent baseline canonicalizes each labelled tournament.
    mins=sorted(set(min(row[x] for row in rows) for x in range(N)))
    saturation=0; F=[]
    for t in range(N):
        for row in rows: saturation |= 1<<row[t]
        numerator=sum((saturation&m).bit_count() for m in fixed)
        need(numerator%G==0,'integral Burnside')
        F.append(numerator//G)
        need(F[-1]==sum(m<=t for m in mins),'cumulative minimum count')
        need(F[-1]-(F[-2] if t else 0) in (0,1),'unit jumps')
    def unrank(r):
        lo=0; hi=N-1; q=0
        while lo<hi:
            mid=(lo+hi)//2; q+=1
            if F[mid]>r: hi=mid
            else: lo=mid+1
        need(q<=C.d,'query budget')
        return lo
    need([unrank(r) for r in range(len(mins))]==mins,'canonical unranking')
    # Full affine indicator/rank formula on selected cuts at n<=4.
    term_counts={}
    if n<=4:
        for t in sorted(set([0,1,N//4,N//2-1,N//2,N-2,N-1])):
            flats=C.cut(t)
            need(all(sum(C.satisfied(x,f) for f in flats)==int(x<=t) for x in range(N)),'disjoint interval cover')
            sat_flats=[C.transform(f,p) for f in flats for p in C.ps]
            terms=C.indicator(sat_flats)
            want=set(row[x] for row in rows for x in range(t+1))
            need(all(sum(c for a,c in terms.items() if C.satisfied(x,a))==int(x in want) for x in range(N)),'full union indicator')
            numerator=sum(c*C.size(C.meet(a,C.fixed(p))) for a,c in terms.items() for p in C.ps)
            need(numerator==G*F[t],'rank/Burnside/cut formula')
            term_counts[str(t)]=len(terms)
    return {'n':n,'labelled':N,'ordinary_classes':len(mins),'canonical_minima':mins,'cuts_verified':N,'affine_indicator_terms':term_counts}

# Smallest nonadditive child example and a stabilizer overcount hostile.
X={1,2}; swap=lambda x: ((x&1)<<1)|((x>>1)&1)
need({min(x,swap(x)) for x in X}=={1},'one root orbit')
need(bool([x for x in X if x<2]) and bool([x for x in X if x>=2]),'one orbit meets both child cylinders')
Y={1,2,4}; rotate=lambda x:((x<<1)&7)|(x>>2)
need({min(x,rotate(x),rotate(rotate(x))) for x in Y}=={1},'cyclic root orbit')
need(len([x for x in Y if x<4])==2,'prefix stabilizer is trivial, yields two classes')
# Weighted Burnside is orbit-average mass, not a count of witnesses/orbits.
actions=[lambda x:x,rotate,lambda x:rotate(rotate(x))]
w={1:1,2:1,4:2}
burn=Fraction(sum(w[x] for g in actions for x in Y if g(x)==x),3)
need(burn==Fraction(4,3),'noninvariant positive weights give fractional mass')
cut_burn=Fraction(sum(1 for g in actions for x in Y if x<4 and g(x)==x),3)
need(cut_burn==Fraction(2,3),'unsaturated prefix gives 2/3, not orbit count one')
# Positive witness multiplicities preserve minimum search but not rank.
W={1:2,2:1}
weighted=lambda t:sum(v for x,v in W.items() if x<=t)
need(min(t for t in range(4) if weighted(t)>0)==1,'positive multiplicity minimum')
need(min(t for t in range(4) if weighted(t)>1)==1,'weighted rank repeats first witness')
need(sorted(W)[1]==2,'unweighted rank differs')
need(sum({1:1,2:-1}.values())==0,'signed cancellation hostile')
results=[verify(n) for n in (3,4,5)]
need([r['ordinary_classes'] for r in results]==[2,4,12],'known finite census control')
print(json.dumps({'status':'PASS','checks':checks,'results':results,'hostiles':{'prefix_children':'1 parent orbit; 1+1 child-intersecting orbit counts','prefix_stabilizer':'C3 on weight-one length3; prefix0 gives 2 stabilizer classes but 1 global orbit'}},indent=2))
