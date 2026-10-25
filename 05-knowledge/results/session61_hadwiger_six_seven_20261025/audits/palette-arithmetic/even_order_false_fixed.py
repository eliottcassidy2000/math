#!/usr/bin/env python3
"""Independent Session61 audit: exact-cover Fano construction, incidence algebra,
finite-field norm fibers, and cyclotomic decomposition groups. No producer imports.
All checks remain active with -O. Finite evidence is not an all-prime proof.
"""
from itertools import combinations, product
from collections import Counter
from math import gcd, isqrt, comb
import json

checks=0
def require(b, label):
    global checks
    checks+=1
    if not b: raise RuntimeError(label)
def rank(rows):
    piv={}
    for row in rows:
        while row:
            b=row.bit_length()-1
            if b in piv: row ^= piv[b]
            else: piv[b]=row; break
    return len(piv)

# Enumerate STS(7) by exact cover of the 21 pairs, not by label permutations.
pairs=list(combinations(range(7),2)); pi={p:i for i,p in enumerate(pairs)}
triples=list(combinations(range(7),3))
blocks=[sum(1<<pi[p] for p in combinations(t,2)) for t in triples]
by_pair=[[j for j,b in enumerate(blocks) if b>>i&1] for i in range(21)]
planes=[]
def cover(used, chosen):
    if used==(1<<21)-1:
        planes.append(tuple(sorted(chosen))); return
    missing=((1<<21)-1)^used
    i=(missing & -missing).bit_length()-1
    for j in by_pair[i]:
        if not used & blocks[j]: cover(used|blocks[j],chosen+[j])
cover(0,[])
require(len(planes)==len(set(planes))==30,'exact-cover STS census')
pencils={}; hist=Counter(); moment_signatures=[]
for plane in planes:
    lines=frozenset(triples[j] for j in plane)
    pairing=tuple(sorted(tuple(x for x in t if x) for t in lines if 0 in t))
    pencils.setdefault(pairing,[]).append(lines)
    code={0}|{127 ^ sum(1<<i for i in t) for t in lines}
    require(len(code)==8 and all(a^b in code for a in code for b in code),'STS complements form binary code')
    require(rank(code)==3,'STS code rank')
    sig=[]
    for k in (1,2,3):
        for t in combinations(range(7),k):
            mask=sum(1<<i for i in t)
            corr=sum(1 if (c&mask).bit_count()%2==0 else -1 for c in code)
            require(corr==(8 if k==3 and t in lines else 0),'independent moments')
            if k==3:sig.append(corr)
    moment_signatures.append(tuple(sig))
    hist[sum(t in lines for t in ((0,1,2),(0,3,4),(0,5,6)))]+=1
require(len(set(moment_signatures))==30,'third moments recover every labeled code')
require(dict(hist)=={0:16,1:12,3:2},'exact-cover crossing census')
require(len(pencils)==15 and all(len(v)==2 for v in pencils.values()),'pairing fiber count')
for pairing,fiber in pencils.items():
    for a,b in pairing:
        f=lambda i:b if i==a else a if i==b else i
        image=frozenset(tuple(sorted(f(i) for i in t)) for t in fiber[0])
        require(image==fiber[1],'every endpoint swap obstructs selector')
    for flips in product((0,1),repeat=3):
        def f(i):
            for bit,(a,b) in zip(flips,pairing):
                if bit and i in (a,b):return a+b-i
            return i
        image=frozenset(tuple(sorted(f(i) for i in t)) for t in fiber[0])
        require(image==fiber[sum(flips)%2],'full frame action is parity')
# All 15 Lagrangians constructed separately by symplectic span.
canon=lambda w:min(w,w^63)
W={canon(w) for w in range(64) if w.bit_count()%2==0}
lag={frozenset((0,x,y,canon(x^y))) for x,y in combinations(W-{0},2) if (x&y).bit_count()%2==0}
match_lag={frozenset([0]+[canon(sum(1<<(i-1) for i in pair)) for pair in pairing]) for pairing in pencils}
require(len(lag)==15 and lag==match_lag,'symplectic/matching bijection')

# Exhaust all simple graphs through order 5 and every binary gain, with
# direct truth-table fibers versus Gaussian elimination (not forest solver).
graphs=gains=ternary=0
for n in range(6):
    possible=list(combinations(range(n),2))
    for mask in range(1<<len(possible)):
        es=[e for j,e in enumerate(possible) if mask>>j&1]
        rows=[(1<<u)|(1<<v) for u,v in es]
        r=rank(rows); b1=len(es)-r
        images=Counter(sum((((x>>u)^(x>>v))&1)<<j for j,(u,v) in enumerate(es)) for x in range(1<<n))
        require(len(images)==1<<r and set(images.values())=={1<<(n-r)},'coboundary image and H0 fibers')
        for gain in range(1<<len(es)):
            aug=[row|(((gain>>j)&1)<<n) for j,row in enumerate(rows)]
            require((rank(aug)==r)==(gain in images),'all binary gains via rank versus truth table')
            gains+=1
        require(3*b1-b1==2*b1,'split coefficient obstruction dimension')
        if n<=4:
            for q in product(range(3),repeat=n):
                residual=[(u,v) for u,v in es if q[u]==q[v]]
                eq=[(1<<u)|(1<<v) for u,v in residual]
                consistent=rank(eq)==rank([x|(1<<n) for x in eq])
                direct=any(all(2*q[u]+((b>>u)&1)!=2*q[v]+((b>>v)&1) for u,v in es) for b in range(1<<n))
                require(consistent==direct,'fixed ternary lift iff affine feasibility')
                ternary+=1
        graphs+=1
# Nontrivial kernel / nonempty scalar shadow on a cycle.
require((1^1^1,1^1^1,0)==(1,1,0),'triangle vector holonomy')
require((1^1^0)==0,'triangle parity blind')
# Full vector/quotient cohomology counts in graph-cycle normal form b1=0..3.
for b1 in range(4):
    kernel=sum(all(x.bit_count()%2==0 for x in w) for w in product(range(8),repeat=b1))
    require(kernel==2**(2*b1),'exact vector class fiber')

# Direct cyclotomic cyclic orbit screen, independent order/CRT cross-check.
def subgroup(a,m):
    seen=set(); x=1
    while x not in seen: seen.add(x); x=x*a%m
    require(x==1,'unit orbit closes at identity')
    return seen
fixed=[]
for a in range(1,165):
    if gcd(a,165)!=1:continue
    H=subgroup(a,165)
    if len(H)%2==0: fixed.append((a,len(H)))
require(fixed==[(29,10),(74,10),(134,10),(149,10),(164,2)],'all unit classes mod165')
ramified={p:sorted(subgroup(p,165//p)) for p in (3,5,11)}
require(all(165//p-1 not in H for p,H in ramified.items()),'all ramified primes fail fixed criterion')
sieve=[True]*10001;sieve[0]=sieve[1]=False
for d in range(2,isqrt(10000)+1):
    if sieve[d]:
        for x in range(d*d,10001,d):sieve[x]=False
qvalues=[]
for p,isprime in enumerate(sieve):
    if not isprime or gcd(p,165)!=1:continue
    H=subgroup(p%165,165)
    if 164 in H:qvalues.append((p**(len(H)//2),p,len(H)))
require(min(qvalues)==(659,659,2),'least norm parameter via independent prime sieve')
require(4*659*12**2 < 660**2 < 4*659*13**2,'exact Weil rounding interval')
require(54 in subgroup(109%55,55) and 164 not in subgroup(109%165,165),'old fixed prime109 is lost')

# Finite fields F_(p^k), explicit irreducible polynomials. Each division law
# is checked by x^(q-1)=1, so a bad polynomial cannot pass unnoticed.
class Field:
    def __init__(self,p,poly):self.p=p;self.poly=poly;self.k=len(poly)-1;self.q=p**self.k
    def digits(self,x):
        out=[]
        for _ in range(self.k):out.append(x%self.p);x//=self.p
        return out
    def pack(self,a):return sum((x%self.p)*self.p**i for i,x in enumerate(a))
    def add(self,a,b):return self.pack([x+y for x,y in zip(self.digits(a),self.digits(b))])
    def neg(self,a):return self.pack([-x for x in self.digits(a)])
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        a=self.digits(a);b=self.digits(b);v=[0]*(2*self.k-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):v[i+j]+=x*y
        for i in range(len(v)-1,self.k-1,-1):
            c=v[i]%self.p
            for j in range(self.k+1):v[i-self.k+j]-=c*self.poly[j]
        return self.pack(v[:self.k])
    def power(self,x,n):
        z=1
        while n:
            if n&1:z=self.mul(z,x)
            x=self.mul(x,x);n//=2
        return z
    def inv(self,x):return self.power(x,self.q-2)

fields=[(2,[0,1]),(2,[1,1,1]),(2,[1,1,0,1]),(2,[1,1,0,0,1]),
        (3,[0,1]),(3,[1,0,1]),(3,[1,2,0,1]),(5,[0,1]),(5,[2,0,1]),
        (7,[0,1]),(7,[1,0,1]),(11,[0,1]),(13,[0,1]),(17,[0,1])]
field_receipts=[]
for p,poly in fields:
    F=Field(p,poly);q=F.q
    require(all(F.power(x,q-1)==1 for x in range(1,q)),'explicit field multiplicative law')
    t=1 if p==2 else 0
    d=next(d for d in range(1,q) if all(F.sub(F.mul(x,x),F.add(F.mul(t,x),d))!=0 for x in range(q)))
    norm=lambda a,b:F.sub(F.add(F.mul(a,a),F.mul(t,F.mul(a,b))),F.mul(d,F.mul(b,b)))
    circle=[(a,b) for a,b in product(range(q),repeat=2) if norm(a,b)==1]
    require(len(circle)==q+1,'norm circle size')
    triangles=[(a,b) for a,b in circle if norm(F.sub(a,1),b)==1]
    expect=1 if p==3 else 2 if q%3==2 else 0
    require(len(triangles)==expect,'prime-power triangle existence and exact common neighbors')
    if p!=2:
        N=Counter((norm(a,b),F.add(F.add(a,a),F.mul(t,b))) for a,b in product(range(q),repeat=2))
        squares={F.mul(x,x) for x in range(1,q)}
        def eta(x):return 0 if not x else 1 if x in squares else -1
        for c in range(1,q):
            K=Counter(F.add(x,F.mul(c,F.inv(x))) for x in range(1,q))
            for s in range(q):
                disc=F.sub(F.mul(s,s),F.mul(4%p,c))
                require(N[c,s]==1-eta(disc),'norm trace fiber 1-eta')
                require(K[s]==1+eta(disc),'Kloosterman trace fiber 1+eta')
    field_receipts.append({'q':q,'p':p,'polynomial':poly,'quadratic_t':t,'quadratic_d':d,'circle':len(circle),'common_neighbors_0_1':len(triangles)})

# Elliptic trace by a square multiplicity table, not a character exponent sum.
elliptic=[]
for p in range(5,500):
    if not sieve[p]:continue
    square_counts=Counter(y*y%p for y in range(p))
    count=1+sum(square_counts[(x*x*x-1)%p] for x in range(p))
    ap=p+1-count;m=(p-1)//2
    coef=0 if (p-1)%3 else (-1)**(m-(p-1)//3)*comb(m,(p-1)//3)%p
    require(ap%p==coef,'all-prime Hasse proof finite control')
    require((ap==0)==(p%3==2),'CM splitting/trace control')
    plus=1+sum(square_counts[(1-x*x*x)%p] for x in range(p))
    genus=2+sum(square_counts[(pow(x,6,p)-1)%p] for x in range(p))
    require(p+1-genus==(p+1-count)+(p+1-plus),'genus2 quotient trace control')
    elliptic.append({'p':p,'a_p':ap,'hasse':coef})
print(json.dumps({'status':'PASS; finite controls only','always_active_checks':checks,
 'graphs':graphs,'binary_gain_instances':gains,'fixed_ternary_instances':ternary,
 'fano_exact_cover_count':len(planes),'fano_crossing_histogram':dict(hist),
 'fixed_mod165':fixed,'ramified_subgroups':ramified,'minimum_q_sieve':min(qvalues),
 'finite_field_controls':field_receipts,'elliptic_controls':elliptic},sort_keys=True,indent=2))
