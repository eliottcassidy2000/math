#!/usr/bin/env python3
"""Exact mixed outward gap-one controls and smooth-numerator lattice resultants.

No floating-point arithmetic, symbolic package, or imported Collatz implementation.
The finite run checks the theorem, not a proposed all-length integer exclusion.
"""
from hashlib import sha256
from math import comb, gcd
from pathlib import Path
import json

checks = 0


def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ArithmeticError(message)


def independent(v, w):
    return any(v[i]*w[j] != v[j]*w[i]
               for i, j in [(0, 1), (0, 2), (1, 2)])


def short_pair(p, r, u):
    first = None
    for height in range(1, p+1):
        for b in range(-height, height+1):
            for c in range(-height, height+1):
                for m in range(-height, height+1):
                    v = (b, c, m)
                    if max(map(abs, v)) != height:
                        continue
                    if (m-b*r-c*u) % p:
                        continue
                    # Exactly one of v,-v is retained.
                    if next(x for x in v if x) < 0:
                        continue
                    if first is None:
                        first = v
                    elif independent(first, v):
                        h1 = max(map(abs, first))
                        need((h1*height)**3 <= p*p,
                             'successive-minimum product bound')
                        return first, v
    raise ArithmeticError('missing independent lattice pair')


def polynomial(p, r, u, family, vector):
    b, c, m = vector
    if c < 0:
        b, c, m = -b, -c, -m
    a, remainder = divmod(m-b*r-c*u, p)
    need(remainder == 0, 'not a lattice relation')
    height = max(map(abs, vector))
    need(abs(a) <= height, 'legal range did not bound a')
    kappa = 2
    e = a-b+kappa*c
    left = 2**max(-e, 0)*3**max(-b, 0)*1
    sign = -1 if c % 2 else 1
    right = sign*2**max(e, 0)*3**max(b, 0)*1
    # F=left*z^mplus*(z-1)^cplus-right*z^mminus*(z-1)^cminus.
    degree = max(max(m, 0)+max(c, 0), max(-m, 0)+max(-c, 0))
    f = [0]*(degree+1)
    for coefficient, shift, power in [
            (left, max(m, 0), max(c, 0)),
            (-right, max(-m, 0), max(-c, 0))]:
        for j in range(power+1):
            f[shift+j] += coefficient*comb(power, j)*(-1)**(power-j)
    while len(f) > 1 and not f[-1]:
        f.pop()
    need(any(f), 'identically zero cleared polynomial')
    need(f[0] != 0 and sum(f) != 0, 'extraneous root zero or one')
    need((len(f) == 1) == (m == 0 and c == 0),
         'constant-polynomial classification')
    need(len(f)-1 <= 2*height, 'degree bound')
    norm = sum(map(abs, f))
    base = 22
    need(norm <= base**height, 'coefficient norm bound')
    return f, dict(a=a, b=b, c=c, m=m, height=height,
                   degree=len(f)-1, coefficient_l1_norm=norm)


def determinant(matrix):
    n = len(matrix)
    if not n:
        return 1
    a = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for k in range(n-1):
        if not a[k][k]:
            pivot = next((j for j in range(k+1, n) if a[j][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                need(numerator % previous == 0, 'Bareiss nonexact division')
                a[i][j] = numerator//previous
        for i in range(k+1, n):
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def resultant(f, h):
    m, n = len(f)-1, len(h)-1
    need(m+n > 0, 'two constants cannot be an independent pair')
    rows = []
    for i in range(n):
        rows.append([0]*i+f[::-1]+[0]*(n-1-i))
    for i in range(m):
        rows.append([0]*i+h[::-1]+[0]*(m-1-i))
    return determinant(rows)


def carry(word):
    value = 0
    for i, bit in enumerate(word):
        value = (3 if bit else 1)*value+bit*2**i
    positional = sum(2**i*3**sum(word[i+1:])
                     for i, bit in enumerate(word) if bit)
    need(value == positional, 'carry methods disagree')
    return value


def move(source, edges):
    word = source[:]
    for i in edges:
        j = (i+1) % len(word)
        need(word[i] != word[j], 'nonactual adjacent swap')
        word[i], word[j] = word[j], word[i]
    return word


def egcd(a,b):
    if not b:return a,1,0
    d,x,y=egcd(b,a%b);return d,y,x-(a//b)*y


def actual_outwards(source):
    q=len(source); paths=set()
    for i in range(q):
        j=(i+1)%q
        if not (source[i] and source[j]):continue
        edges=((i-1)%q,j)
        if {edges[0],(edges[0]+1)%q}&{edges[1],(edges[1]+1)%q}:continue
        if all(source[e]!=source[(e+1)%q] for e in edges):
            paths.add(edges);paths.add(edges[::-1])
    return paths


def check_case(p,r,u,source,labels,pair=None):
    q=p+r
    edges=[(labels[u+r]-1)%q,labels[u]]
    need(not ({edges[0],(edges[0]+1)%q}&{edges[1],(edges[1]+1)%q}),'not disjoint')
    w=move(source,edges);need(w==move(source,edges[::-1]),'noncommuting')
    B=carry(w);D=2**q-3**p;g=gcd(abs(D),B)
    if abs(D)>1:
        modulus=abs(D);_,a,b=egcd(p,r)
        z=pow(2,a,modulus)*pow(3*pow(2,-1,modulus)%modulus,b,modulus)%modulus
        need(pow(z,p,modulus)==2%modulus and pow(z,r,modulus)==3*pow(2,-1,modulus)%modulus,'clock root')
        need(gcd(modulus,4+pow(z,u,modulus)*(z-1))==g,'whole gcd iff')
    vectors=short_pair(p,r,u) if pair is None else pair
    f,d1=polynomial(p,r,u,'outward',vectors[0]);h,d2=polynomial(p,r,u,'outward',vectors[1])
    R=resultant(f,h)
    need(R!=0 and R%g==0,'valuation coprimality or full gcd failed')
    need(abs(R)<=d1['coefficient_l1_norm']**d2['degree']*d2['coefficient_l1_norm']**d1['degree'],'Hadamard')
    return dict(p=p,r=r,u=u,q=q,edges=edges,word=''.join(map(str,w)),B=B,D=D,content=g,vectors=[d1,d2],polynomials=[f,h],resultant=R,clock_resultant_gcd=gcd(abs(D),abs(R)))

trace=sha256();shapes=pairs=ordered=cut=0; maximum=None;hits=[];depths=[];extra=[]
for p in range(2,61):
    for r in range(1,p):
        if gcd(p,r)!=1:continue
        shapes+=1;q=p+r
        source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
        labels={(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(p)}
        predicted=set()
        us=range(max(0,p-2*r),min(r,p-r)) if p>=3 else []
        for u in us:
            row=check_case(p,r,u,source,labels)
            need(3*r>p and u<p-r,'legal inequalities')
            predicted.add(tuple(row['edges']));predicted.add(tuple(row['edges'][::-1]))
            pairs+=1;ordered+=2;cut+=2*(q-1 in row['edges'])
            trace.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode()+b'\n')
            if maximum is None or row['content']>maximum['content']:maximum=row
            if row['content']==abs(row['D']):hits.append(row)
            if any(row['content']%(ell*ell)==0 for ell in range(5,100)) and len(depths)<4:depths.append(row)
            if row['clock_resultant_gcd']>row['content'] and len(extra)<4:extra.append(row)
        need(predicted==actual_outwards(source),'ordered complete legality')
need(not hits,'integer hostile found')
for a in (-1,1):
    for b in (-1,1):
        for c in (0,1):
            e=a-b+2*c
            need(2**(max(-e,0)+c)*3**max(-b,0)+2**max(e,0)*3**max(b,0)<=22,'joint coefficient corner')
# Illegal parameters deliberately violate the short-coefficient legality premise.
need((1-1*1-0*9)%3==0 and (0-0*1-2*9)%3==0,'illegal rows still satisfy relation congruence')
need(independent((1,0,1),(0,2,0)),'illegal hostile rows independent')
need(resultant([-3,2],[3,-8,4])==0,'illegal common-component hostile')
print(json.dumps(dict(status='FINITE-EXACT support for candidate smooth mixed outward content bound',script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,shapes=shapes,unordered_edit_pairs=pairs,ordered_actual_paths=ordered,cut_crossing_ordered_paths=cut,integer_hits=hits,maximum_content=maximum,prime_power_controls=depths,extra_clock_resultant_factors=extra,trace_sha256=trace.hexdigest()),indent=2))
