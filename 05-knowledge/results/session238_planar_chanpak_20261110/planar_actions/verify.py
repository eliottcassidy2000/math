#!/usr/bin/env python3
"""Exact controls, no dependencies; --inject-failure must exit nonzero."""
from fractions import Fraction as F
from itertools import product
import json
import sys

checks = 0

def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

def mat(rows):
    return tuple(tuple(F(x) for x in row) for row in rows)

I = mat([[1,0,0],[0,1,0],[0,0,1]])

def mul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))

def inv(a):
    rows = [list(a[i])+list(I[i]) for i in range(3)]
    for i in range(3):
        j = next(j for j in range(i,3) if rows[j][i])
        rows[i],rows[j] = rows[j],rows[i]
        q = rows[i][i]
        rows[i] = [x/q for x in rows[i]]
        for j in range(3):
            if j != i:
                q=rows[j][i]
                rows[j]=[x-q*y for x,y in zip(rows[j],rows[i])]
    return tuple(tuple(row[3:]) for row in rows)

def power(a,n):
    if n<0:
        a,n=inv(a),-n
    ans=I
    for _ in range(n):
        ans=mul(a,ans)
    return ans

def comm(a,b):
    return mul(mul(mul(a,b),inv(a)),inv(b))

def apply(a,v):
    return tuple(sum(row[j]*v[j] for j in range(3)) for row in a)

for x in (F(1,3),1,2,3,4,7):
    for d in (F(2,3),1,3):
        k=F(x,2)
        A=mat([[1,0,0],[0,2,0],[0,0,1]])
        C=mat([[1,k,0],[0,4,-2*d],[0,0,1]])
        D=mul(C,power(A,-2))
        K=comm(A,D)
        J=comm(D,K)
        check(D==mat([[1,k/4,0],[0,1,-2*d],[0,0,1]]),'D order')
        check(K==mat([[1,-k/8,-k*d/4],[0,1,-2*d],[0,0,1]]),'K formula')
        check(J==mat([[1,0,-3*k*d/4],[0,1,0],[0,0,1]]),'translation')
        check(comm(K,D)==inv(J),'reversed commutator hostile')
        check(mul(J,inv(J))==I,'inverse')
        if d==1:
            check(apply(power(J,-8),(1,2,1))==(1+3*x,2,1),'first gap')
            check(apply(power(J,-24),(1+x,6,1))==(1+10*x,6,1),'second gap')
            if x in (2,3):
                for j,n,v,expected in [(-8,3,(1,2,1),(1+3*x,16,1)),
                                        (-24,2,(1+x,6,1),(1+10*x,24,1))]:
                    out=apply(mul(power(J,j),power(A,n)),v)
                    check(out==expected,'cone padded gap')
                    check(out[0]<=x*out[1]/2,'inherited physical cone hostile')

# Enumerate every positive word of total vertex cost <=14.
words=[('',0,0,2)]
for word,n,m,r in words:
    for letter,cost in [('A',1),('C',2)]:
        if n+cost<=14:
            words.append((word+letter,n+cost,2*m if letter=='A' else 4*m+1,
                          2*r if letter=='A' else 4*r-2))
seen={}
for word,n,m,r in words:
    check(r==2**(n+1)-2*m,'normal form')
    check(not (m & (m<<1)),'separated bits')
    check(m==0 or m.bit_length()<=n-1,'highest bit')
    check(3*m<=2**n-1,'bound')
    check(r not in seen,'seed orbit injectivity')
    seen[r]=word
    bits=format(m,'0%db'%n) if n else ''
    check(bits==''.join('0' if c=='A' else '01' for c in word),'binary code')
    if word:
        check(r>2,'strict drift')

# Independent Hamilton-path subset DP. Root is vertex0; rows encode outgoing arcs.
def hp(rows,mask=None):
    n=len(rows)
    if mask is None: mask=(1<<n)-1
    if mask==0:return 1
    dp={}
    for v in range(n):
        if mask>>v&1:dp[(1<<v,v)]=1
    for size in range(1,mask.bit_count()):
        for (s,v),count in list(dp.items()):
            if s.bit_count()!=size:continue
            avail=rows[v]&mask&~s
            while avail:
                bit=avail&-avail;avail-=bit;w=bit.bit_length()-1
                dp[(s|bit,w)]=dp.get((s|bit,w),0)+count
    return sum(dp.get((mask,v),0) for v in range(n))

def state(rows):
    full=(1<<len(rows))-1
    R=0
    def walk(v,mask):
        nonlocal R
        R+=hp(rows,full^mask)
        avail=rows[v]&~mask
        while avail:
            bit=avail&-avail;avail-=bit
            walk(bit.bit_length()-1,mask|bit)
    walk(0,1)
    # e via paths on root-deleted tournament whose endpoint beats root.
    # Only e=0 is needed: a present sink other than root forces this.
    return hp(rows),R,hp(rows,full^1)

def append(rows,bridge=False):
    n=len(rows)
    ans=[r|(1<<n) for r in rows]+[0]
    if bridge:
        ans[0]&=~(1<<n);ans[n]|=1
    return ans
physical=0
for core in ([0],[2,4,1]):
    seed=append(core)
    h,R,d=state(seed)
    check((h,R,d)==((1,2,1) if len(core)==1 else (3,6,1)),'physical seed')
    for word,n,m,r in words:
        if n>6:continue
        rows=seed[:];hh,rr=h,R
        for letter in word:
            if letter=='A':
                rows=append(rows);rr*=2
            else:
                rows=append(append(rows,True));hh,rr=hh+rr,4*rr-2*d
        check(state(rows)==(hh,rr,d),'physical transfer '+word)
        check(hh not in (7,21),'small forbidden hostile')
        physical+=1

if '--inject-failure' in sys.argv:
    check(False,'injected failure: guard survives optimization')
print(json.dumps({'status':'PASS','checks':checks,'positive_words_weight_le14':len(words),
                  'independent_physical_instances_weight_le6':physical},sort_keys=True))
