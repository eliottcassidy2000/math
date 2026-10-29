#!/usr/bin/env python3
"""Exact finite controls for translated-pole memory; standard library only."""
from fractions import Fraction as F
from itertools import product
from math import comb, prod, gcd
import json

checks = 0

def eq(x, y, label):
    global checks
    checks += 1
    if x != y:
        raise RuntimeError((label, x, y))

def merged(source, n, a=F(1)):
    out = {}
    for t, c in source.items():
        for j in range(n):
            u = t + j*a
            out[u] = out.get(u,F(0)) + c
    return {t:c for t,c in out.items() if c}

def moments(source, count):
    return [sum((c*t**r for t,c in source.items()),F(0)) for r in range(count)]

def solve(matrix, rhs):
    a = [list(map(F,row))+[F(v)] for row,v in zip(matrix,rhs)]
    for j in range(len(a)):
        pivot = next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:
            raise RuntimeError('singular system')
        a[j],a[pivot]=a[pivot],a[j]
        d=a[j][j];a[j]=[v/d for v in a[j]]
        for i in range(len(a)):
            if i != j:
                d=a[i][j];a[i]=[x-d*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]

def det(matrix):
    a=[list(map(F,row)) for row in matrix];out=F(1)
    for j in range(len(a)):
        pivot=next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:return F(0)
        if pivot != j:a[j],a[pivot]=a[pivot],a[j];out=-out
        d=a[j][j];out*=d
        for i in range(j+1,len(a)):
            q=a[i][j]/d
            for k in range(j+1,len(a)):a[i][k]-=q*a[j][k]
    return out

def poly_roots(sites):
    p=[F(1)]
    for t in sites:
        q=[F(0)]*(len(p)+1)
        for j,c in enumerate(p):q[j]-=t*c;q[j+1]+=c
        p=q
    return p

def audit_hankel(source):
    s=len(source)
    m=moments(source,3*s+4)
    if not s:
        eq(m,[F(0)]*4,'zero source');return
    h=[[m[i+j] for j in range(s)] for i in range(s)]
    sites=sorted(source)
    expected=prod(source.values())*prod((sites[j]-sites[i])**2 for i in range(s) for j in range(i+1,s))
    eq(det(h),expected,'Vandermonde determinant with signed weights')
    q=poly_roots(sites)
    fitted=solve(h,[-m[s+i] for i in range(s)])+[F(1)]
    eq(fitted,q,'recover exact recurrence from 2s moments')
    for r in range(len(m)-s):eq(sum(q[j]*m[r+j] for j in range(s+1)),0,'global recurrence')
    # Independent interpolation of original weights after site recovery.
    v=[[t**r for t in sites] for r in range(s)]
    eq(solve(v,m[:s]),[source[t] for t in sites],'recover weights')

words=0
for d in range(5):
    for c in product((-1,0,1), repeat=d+1):
        if not c[0] or not c[-1]:continue
        words+=1
        source={F(i):F(x) for i,x in enumerate(c) if x}
        total=sum(c);prefix=[sum(c[:j+1]) for j in range(d)]
        for n in range(1,11):
            out=merged(source,n)
            # Direct coefficient convolution has exact total n*sum(c).
            eq(sum(out.values()),n*total,'source weight conserved with multiplicity')
            if n>d:
                count=(n-d)*bool(total)+sum(bool(x)+bool(total-x) for x in prefix)
                eq(len(out),count,'signed boundary plateau count')
            first=moments(source,7);now=moments(out,7);later=moments(merged(source,n+1),7)
            for r in range(7):
                eq(later[r]-now[r],sum(F(comb(r,l))*n**(r-l)*first[l] for l in range(r+1)),'triangular time recurrence')
        if d<=2:
            for n in range(1,5):audit_hankel(merged(source,n))

# Multiple unrelated translation orbits, nonunit translation, internal gaps.
examples=[({F(0):F(1),F(1):F(-1),F(1,2):F(2)},F(1)),
          ({F(-3):F(2),F(1):F(-2),F(0):F(1)},F(2)),
          ({F(0):F(1),F(2,3):F(-1),F(1,5):F(-2)},F(1,3))]
for src,a in examples:
    # Partition without floating-point equality or modulo approximations.
    groups=[]
    for t,c in sorted(src.items()):
        for g in groups:
            offset=(t-g[0])/a
            if offset.denominator==1:g[1][int(offset)]=c;break
        else:groups.append((t,{0:c}))
    for n in range(1,10):
        out=merged(src,n,a)
        if n>max(max(g)-min(g) for _,g in groups):
            prediction=0
            for _,g in groups:
                lo,hi=min(g),max(g);d=hi-lo
                c=[g.get(j,F(0)) for j in range(lo,hi+1)]
                total=sum(c);prefix=[sum(c[:j+1]) for j in range(d)]
                prediction+=(n-d)*bool(total)+sum(bool(v)+bool(total-v) for v in prefix)
            eq(len(out),prediction,'multiple orbit affine rank')
        first=moments(src,6);now=moments(out,6);later=moments(merged(src,n+1,a),6)
        for r in range(6):eq(later[r]-now[r],sum(comb(r,l)*(a*n)**(r-l)*first[l] for l in range(r+1)),'scaled clock recurrence')
    for n in range(1,3):audit_hankel(merged(src,n,a))

# Recurrence valid from r=0 must count the zero-site initial spike.
eq(moments({F(0):F(3)},6),[3,0,0,0,0,0],'zero site boundary')
audit_hankel({F(0):F(3)})

# Identical two-moment time traces, bounded versus growing full memory.
for n in range(1,31):
    integrable=merged({F(0):F(1),F(1):F(-1)},n)
    nonintegrable=merged({F(0):F(2),F(1,2):F(-2)},n)
    eq(moments(integrable,2),moments(nonintegrable,2),'same two query histories')
    eq(len(integrable),2,'bounded source rank')
    eq(len(nonintegrable),2*n,'unbounded source rank')

# All-time bounded-jet counterexamples, with each original site in a distinct Z-coset.
for R in range(9):
    L=R+1
    src={F(j,R+3):F((-1)**j*comb(L,j)) for j in range(L+1)}
    expected=F((-1)**L*prod(range(1,L+1)),(R+3)**L)
    eq(moments(src,L+1)[L],expected,'first visible finite-difference moment')
    for n in range(1,13):
        out=merged(src,n);m=moments(out,L+1)
        eq(m[:L],[F(0)]*L,'all-time invisible bounded jet')
        eq(m[L],n*expected,'first visible moment grows linearly')
        eq(len(out),(R+2)*n,'hidden full rank')

# Sharp 2B-moment uniqueness threshold via Lagrange weights.
for B in range(1,9):
    sites=list(map(F,range(2*B)))
    weights={t:1/prod(t-u for u in sites if u!=t) for t in sites}
    m=moments(weights,2*B)
    eq(m,[F(0)]*(2*B-1)+[F(1)],'sharp Lagrange witness')
    left={t:weights[t] for t in sites[:B]};right={t:-weights[t] for t in sites[B:]}
    eq(moments(left,2*B-1),moments(right,2*B-1),'two B-site sources share short prefix')
    if B<=4:audit_hankel(left);audit_hankel(right)

# Frobenius/modular hostile: sites, not only coefficients, can collapse.
for p in (2,3,5,7,11):
    source={F(0):F(1),F(p):F(-1)}
    eq([int(x)%p for x in moments(source,12)],[0]*12,'modular loss of pole separation')
    eq(det([[moments(source,3)[i+j] for j in range(2)] for i in range(2)]),-p*p,'integral Hankel debt')

# All cyclic source vectors over Z/M for M=2..5: direct first return vs formula.
cyclic_sources=0
for modulus in range(2,6):
    for source in product(range(modulus),repeat=modulus):
        cyclic_sources+=1
        d=next(d for d in range(1,modulus+1) if modulus%d==0 and
               all(source[(j-d)%modulus]==source[j] for j in range(modulus)))
        sd=sum(source[:d])%modulus
        expected=d*modulus//gcd(modulus,sd)
        state=[0]*modulus
        observed=None
        for n in range(1,modulus*modulus+1):
            state=[(state[j]+source[(j-(n-1))%modulus])%modulus for j in range(modulus)]
            if not any(state):observed=n;break
        eq(observed,expected,'exact finite source clock')
        eq(all(source[(j-expected)%modulus]==source[j] for j in range(modulus)),True,'returned source phase')
        # Direct d-step recurrence at a separate later time.
        n=modulus+1
        c_n=[sum(source[(j-k)%modulus] for k in range(n))%modulus for j in range(modulus)]
        c_nd=[sum(source[(j-k)%modulus] for k in range(n+d))%modulus for j in range(modulus)]
        eq(c_nd,[(v+sd)%modulus for v in c_n],'finite source block recurrence')

print(json.dumps({'status':'PASS','checks':checks,'signed_words':words,
                  'signed_word_spans':'0..4','iteration_times':'1..10',
                  'all_time_hostile_cutoffs':'R=0..8, n=1..12',
                  'cyclic_sources':cyclic_sources,'cyclic_moduli':'2..5',
                  'recurrence_convention':'all r>=0; zero site retained'},sort_keys=True))
