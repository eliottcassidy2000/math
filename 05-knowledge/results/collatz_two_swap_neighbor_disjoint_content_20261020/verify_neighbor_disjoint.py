#!/usr/bin/env python3
"""Exact neighbor-disjoint controls, short vectors, and Sylvester resultants.

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
    delta = 0 if family == 'plus' else 1
    e = a+2*b+delta*c
    fexp = c-b
    left = 2**max(-e, 0)*3**max(-fexp, 0)*7**c
    sign = -1 if family == 'plus' and c % 2 else 1
    right = sign*2**max(e, 0)*3**max(fexp, 0)
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
    base = 121 if family == 'plus' else 65
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
    if b==0: return a,1,0
    d,x,y=egcd(b,a%b)
    return d,y,x-(a//b)*y


def actual_paths(source):
    q=len(source)
    positions=[i for i,v in enumerate(source) if v]
    rank={v:i for i,v in enumerate(positions)}
    paths=set()
    for i in range(q):
        if source[i]==source[(i+1)%q]: continue
        direction=1 if source[i] else -1
        label1=i if direction==1 else (i+1)%q
        first={i,(i+1)%q}
        for j in range(q):
            if first & {j,(j+1)%q}: continue
            if source[j]==source[(j+1)%q]: continue
            second_direction=1 if source[j] else -1
            if direction != second_direction: continue
            label2=j if direction==1 else (j+1)%q
            if (rank[label2]-rank[label1])%len(positions) not in (1,len(positions)-1): continue
            paths.add((i,j))
    return paths


def check_case(p,r,u,family,source,labels,pair=None):
    q=p+r; n=p-r
    edges=[labels[u],labels[u+n]]
    if family=='minus': edges=[(i-1)%q for i in edges]
    need(not ({edges[0],(edges[0]+1)%q}&{edges[1],(edges[1]+1)%q}), 'edges overlap')
    word=move(source,edges)
    need(word==move(source,edges[::-1]),'disjoint orders do not commute')
    D=2**q-3**p; B=carry(word); g=gcd(abs(D),B)
    need(gcd(g,42)==1,'forbidden content prime')
    if abs(D)>1:
        _,a,b=egcd(p,n)
        modulus=abs(D)
        z=pow(2,a,modulus)*pow(4*pow(3,-1,modulus)%modulus,b,modulus)%modulus
        need(pow(z,p,modulus)==2%modulus and pow(z,n,modulus)==4*pow(3,-1,modulus)%modulus,'clock reconstruction')
        factor=pow(z,u,modulus)*(z-1)%modulus
        gate=(3+7*factor) if family=='plus' else (6-7*factor)
        need(gcd(modulus,gate)==g,'full-content iff gate')
    vectors=short_pair(p,n,u) if pair is None else pair
    need(independent(*vectors),'dependent vectors')
    f,d1=polynomial(p,n,u,family,vectors[0])
    h,d2=polynomial(p,n,u,family,vectors[1])
    R=resultant(f,h)
    need(R!=0 and R%g==0,'nonzero full-content resultant')
    need(abs(R)<=d1['coefficient_l1_norm']**d2['degree']*d2['coefficient_l1_norm']**d1['degree'],'Hadamard')
    return dict(p=p,r=r,n=n,u=u,family=family,q=q,edges=edges,word=''.join(map(str,word)),D=D,B=B,content=g,vectors=[d1,d2],polynomials=[f,h],resultant=R,clock_resultant_gcd=gcd(abs(D),abs(R)))

trace=sha256(); shapes=pairs=ordered=cut=0; hits=[]; maximum=None; extra=[]; depths=[]
for p in range(2,49):
    for r in range(1,p):
        if gcd(p,r)!=1: continue
        shapes+=1; q=p+r; n=p-r
        source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
        labels={(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(p)}
        predicted=set()
        for family,us in [('plus',range(max(0,r-n))),('minus',range(n,r))]:
            for u in us:
                row=check_case(p,r,u,family,source,labels)
                predicted.add(tuple(row['edges'])); predicted.add(tuple(row['edges'][::-1]))
                pairs+=1; ordered+=2; cut+=2*(q-1 in row['edges'])
                trace.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode()+b'\n')
                if maximum is None or row['content']>maximum['content']: maximum=row
                if row['content']==abs(row['D']): hits.append(row)
                if row['clock_resultant_gcd']>row['content'] and len(extra)<4: extra.append(row)
                if any(row['content']%(ell*ell)==0 for ell in range(5,100)) and len(depths)<4: depths.append(row)
        need(predicted==actual_paths(source),'complete ordered actual neighbor legality')
need(not hits,'integer hostile found')
# Exact coefficient-box corner controls; log convexity gives every real point.
for family,delta,c1,c2 in [('plus',0,112,9),('minus',1,56,16)]:
    for a in (-1,1):
        for b in (-1,1):
            for c in (0,1):
                e=a+2*b+delta*c; f=c-b
                need(2**(max(-e,0)+c)*3**max(-f,0)*7**c<=c1,'left corner bound')
                need(2**max(e,0)*3**max(f,0)<=c2,'right corner bound')
                need(2**(max(-e,0)+c)*3**max(-f,0)*7**c+2**max(e,0)*3**max(f,0)<= (121 if family=='plus' else 65),'joint corner bound')
controls=[]
p,r,u=5,3,0; q=p+r
source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
labels={(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(p)}
for pair in [((0,1,0),(5,1,0)),((5,0,0),(0,1,0))]:
    controls.append(check_case(p,r,u,'plus',source,labels,pair))
print(json.dumps(dict(status='FINITE-EXACT neighbor-disjoint content controls, no all-length exclusion claimed',script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,shapes=shapes,unordered_edit_pairs=pairs,ordered_actual_paths=ordered,cut_crossing_ordered_paths=cut,integer_hits=hits,maximum_content=maximum,extra_resultant_clock_factors=extra,prime_power_controls=depths,degenerate_projection_controls=controls,trace_sha256=trace.hexdigest()),indent=2)+'\n',end='')
