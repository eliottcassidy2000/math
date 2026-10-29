#!/usr/bin/env python3
"""Exact finite replay; infinite conclusions are in RETURN_CONTROLLER_TRAP.md."""
from fractions import Fraction as F
import json

def need(x, msg):
    if not x:
        raise RuntimeError(msg)

def step(b,n):
    z=b*n+1
    a=(z & -z).bit_length()-1
    return z>>a,a

def ratio(b,n,y,delta):
    h=F(1,b-2)
    return (y+h)/(n+h)*F(b,2)**delta

def affine(b,word):
    c=d=0
    for a in word:
        c=b*c+2**d
        d+=a
    return b**len(word),c,d

def replay(b,n,word):
    for a in word:
        n,aa=step(b,n)
        need(aa==a,'inverse guard valuation')
    return n

# General symbolic threshold, checked by exact cross multiplication across
# a stated finite parameter grid, independent of the two specialized examples.
threshold_checks=0
for b in range(3,20,2):
    for r in range(1,21):
        delta=2**(r+2)-b**3
        if delta<=0: continue
        H=F(b*b+(b+2)*2**r,delta)
        A,C,D=affine(b,[r,1,1])
        need(H==F(C,2**D-A),'balanced block threshold')
        for n in range(1,202,2):
            y=F(b*n+1,2**r)
            rr=ratio(b,n,y,2)
            need((rr<1)==(n>H),'threshold strict direction')
            need((rr==1)==(n==H),'threshold equality')
            threshold_checks+=1

edges=0
for b,r,rho in [(3,3,F(27,28)),(5,5,F(875,884))]:
    for n in range(1,200002,2):
        y,a=step(b,n)
        branch=a%r==0
        rr=ratio(b,n,y,2 if branch else -1)
        if not branch:
            need(rr<=1,'leaf potential')
            need((rr==1)==(a==1),'leaf equality')
        elif b==3:
            need(n>=13 and rr<=rho,'3 branch guard')
        elif n==19:
            need(rr==F(125,116),'19 exception')
        elif n==83:
            need(rr==1,'83 exception')
        else:
            need(rr<=rho,'5 branch guard')
        edges+=1
need(ratio(3,13,5,2)==F(27,28),'sharp3')
need(ratio(5,147,23,2)==F(875,884),'sharp5')
need(F(335,512)<F(875,884),'high valuation bound')
need(F(25,4)*F(14,84)==F(25,24),'uncorrected shift obstruction')

# Direct native/queue replay, without using the special-input diagnostic.
# Repeated (n,Q) is an independent finite certificate of nontermination.
counts={}
for b,r in [(3,3),(5,5)]:
    stops=traps=0
    for source in range(1,20002,2):
        n=source;q=1; seen=set();visited83=False;c=d=j=0
        while q>0 and (n,q) not in seen:
            need(j<100000,'finite experiment cap exceeded')
            seen.add((n,q))
            visited83 |= n==83
            y,a=step(b,n)
            q+=2 if a%r==0 else -1
            c=b*c+2**d;d+=a;j+=1;n=y
            need(2**d*n==b**j*source+c,'affine identity')
        if q==0:
            stops+=1
            if b==5: need(not visited83,'escaped83 trap')
        else:
            traps+=1
            need(b==5 and visited83,'unclassified recurrent state')
    counts[str(b)]={'sources':stops+traps,'stops':stops,'traps':traps}

cylinders=0
for b,word,fixed in [(3,[3,1,1],F(49,5)),(5,[5,1,1],F(83))]:
    for k in range(1,41):
        ww=word*k
        A,C,D=affine(b,ww);M=2**(D+1)
        residue=((2**D-C)*pow(A,-1,M))%M
        expected=(fixed.numerator*pow(fixed.denominator,-1,M))%M
        need(residue==expected,'repeated block inverse residue')
        for t in range(4):
            n=residue+M*t
            y=replay(b,n,ww)
            need(F(y)==fixed+F(A,2**D)*(n-fixed),'fixedpoint affine')
            cylinders+=1
    # Rational orbit independently confirms exact valuation blocks.
    n=fixed
    for a in word:
        z=b*n+1
        numerator=z.numerator
        aa=(numerator & -numerator).bit_length()-1
        need(z.denominator%2==1 and aa==a,'rational cycle valuation')
        n=z/2**a
    need(n==fixed,'rational cycle return')

# Arbitrary positive initial queue controls at the two exceptional inputs.
for q0 in range(1,101):
    n=83;q=q0
    for j in range(300):
        n,a=step(5,n);q+=2 if a%5==0 else -1
        need(q>0,'83 queue collapse')
    need((n,q)==(83,q0),'83 queue return')
    n=19;q=q0;j=0
    while q:
        n,a=step(5,n);q+=2 if a%5==0 else -1;j+=1
        need(j<=q0+3,'19 escape clock')
    need(j==q0+3,'19 exact clock')

print(json.dumps({'status':'PASS','threshold_checks':threshold_checks,'edge_inputs':edges,'first_return_replay':counts,
    'repeated_block_sources':cylinders,'max_block_repetitions':40,
    'exceptional_initial_queues':100,'uses_assert':False},indent=2,sort_keys=True))
