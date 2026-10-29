#!/usr/bin/env python3
"""Exact finite controls for the source-preserving valuation/tree factor.
Standard library only. No assertions: checks survive python -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
import json

checks = 0
def need(ok, msg):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(msg)

def v2(n):
    if not n:
        raise ValueError('infinite valuation')
    n = abs(n)
    return (n & -n).bit_length() - 1

def step(x, kind):
    q = 3*x + (1 if kind == 'plus' else -1)
    a = v2(q)
    y = q // 2**a
    if kind == 'anti':
        y *= 3**(a-1)
    return y, a

def packet(word, kind):
    P, B, A = 1, 0, 0
    for a in word:
        if kind == 'anti':
            p, b = 3**a, -3**(a-1)
        else:
            p, b = 3, 1 if kind == 'plus' else -1
        P, B, A = p*P, p*B+b*2**A, A+a
    mod = 2**(A+1)
    return P, B, A, ((2**A-B)*pow(P, -1, mod)) % mod

def compositions(n):
    if n == 0:
        yield ()
    else:
        for a in range(1,n+1):
            for tail in compositions(n-a):
                yield (a,)+tail

def queue(word):
    q = 1
    for i,a in enumerate(word,1):
        q += 2 if a%2 == 0 else -1
        if q == 0:
            return i
    return None

# Independent direct arithmetic versus affine packet, signed controls included.
words = 0
for A in range(1,11):
    for word in compositions(A):
        words += 1
        for kind in ('plus','minus','anti'):
            P,B,total,r = packet(word,kind)
            for x in (r, r+2**(A+1), r-2**(A+1)):
                y=x
                observed=[]
                for _ in word:
                    y,a=step(y,kind)
                    observed.append(a)
                need(tuple(observed)==word, ('guard',kind,word,x))
                need(y*2**A==P*x+B, ('carry',kind,word,x))
            if kind=='plus':
                _,_,_,rm=packet(word,'minus')
                need((r+rm)%2**(A+1)==0, ('sign guard',word))

# Count all stopped valuation words by clock. Every word is one native cylinder.
counts=Counter()
for A in range(1,17):
    for w in compositions(A):
        if queue(w)==len(w):
            counts[len(w),A] += 1
for N in range(1,17):
    for A in range(1,17):
        expected=0
        if (N-1)%3==0:
            m=(N-1)//3
            d=A-(4*m+1)
            if d>=0 and d%2==0:
                h=d//2
                expected=comb(3*m+1,m)//(3*m+1)*comb(3*m+h,h)
        need(counts[N,A]==expected, ('tree clock count',N,A,counts[N,A],expected))

# Enumerate a whole dyadic residue universe, not just predicted representatives.
precision=12
source_counts={}
for kind in ('plus','minus','anti'):
    obs=Counter()
    for x in range(1,2**(precision+1),2):
        y=x; A=0; N=0; q=1
        while q>0:
            try:
                y,a=step(y,kind)
            except ValueError:
                break
            A += a; N += 1
            if A>precision:
                break
            q += 2 if a%2==0 else -1
            if q==0:
                obs[N,A] += 1
    for N,A in counts:
        if A<=precision:
            need(obs[N,A]==counts[N,A]*2**(precision-A),('source enumeration',kind,N,A))
    source_counts[kind]=sum(obs.values())
need(len(set(source_counts.values()))==1,'maps disagree in source mass')

# Independent b/G factorization and exact conditional total-clock moments.
for b in (0,1):
    for g in range(20):
        a=1+b+2*g
        pb=F(1,3) if b else F(2,3)
        pg=F(3,4)*F(1,4)**g
        need(F(1,2)**a==pb*pg, ('independent gap',a))
for m in range(15):
    N=3*m+1
    mean=N+m+F(2*N,3)
    need(mean==2*N-F(1,3),'mean clock')

# Clock/carry and positivity hostiles; no ROOT or halting conclusion from queue.
p12=packet((1,2),'plus'); p21=packet((2,1),'plus')
need(p12[:1]==p21[:1] and p12[2]==p21[2] and p12[1]!=p21[1] and p12[3]!=p21[3], 'order hostile')
shape_a=(4,1,1,1); shape_b=(2,3,1,1)
sa=packet(shape_a,'plus'); sb=packet(shape_b,'plus')
need(queue(shape_a)==4 and queue(shape_b)==4 and sum(shape_a)==sum(shape_b),'same shape and clock')
need(sa[1]!=sb[1] and sa[3]!=sb[3], 'ordered residual gaps needed')
growth_word=(2,1,1,1)
need(packet(growth_word,'plus')==(81,103,5,41),'growth cylinder packet')
need(queue(growth_word)==4,'growth cylinder tree stop')
y=41; trajectory=[y]
for a in growth_word:
    y,actual=step(y,'plus')
    need(actual==a,'growth cylinder trajectory valuation')
    trajectory.append(y)
need(trajectory==[41,31,47,71,107],'growth cylinder trajectory')
need(81>32 and 103>0,'every positive cylinder source grows at tree stop')
need(step(1,'plus')==(1,2),'positive fixed source')
need(step(-5,'plus')==(-7,1) and step(-7,'plus')==(-5,2),'five-clock roots')
need(queue((1,2)*40)==1 and queue((2,1)*40) is None,'root-sensitive tree')
anti_word=(2,)*6+(1,)*12
need(queue(anti_word) is None,'tree alive at machine barrier')
need(6+sum(2*a-3 for a in anti_word)==0,'exact antihydra barrier')
need(queue((1,))==1 and 6+2-3==5,'tree death is not machine halt')
# Exceptional odd rational source 1/3 has 3x-1=0, not finite first digit.
need(3*F(1,3)-1==0 and F(1,3).denominator%2==1,'punctured domain')
need((3*F(1,13)+1)==F(16,13),'golden rational fixed odd orbit')

# Independent bivariate Taylor PGF arithmetic; no moment recursion imported.
def mul(p,q):
    out=Counter()
    for (i,j),v in p.items():
        for (h,l),w in q.items():
            if i+h<=2 and j+l<=2:
                out[i+h,j+l] += v*w
    return dict(out)
def plus(p,q):
    out=Counter(p)
    for k,v in q.items():
        out[k] += v
    return dict(out)
# Taylor coefficients of D(1+z), B(1+z), obtained by rational division.
def rational_series(numer):
    out=[]
    for j in range(3):
        v=F(numer[j] if j<len(numer) else 0)
        if j>=1: v += 2*out[j-1]
        if j>=2: v += out[j-2]
        out.append(v/3)
    return {(0,j):v for j,v in enumerate(out)}
die=rational_series([2,2]); birth=rational_series([1,2,1])
H={(0,0):F(1),(1,0):F(1)}
mixed=[]
for k in range(9):
    ez=H.get((1,0),0); ed=H.get((0,1),0)
    vz=2*H.get((2,0),0)+ez-ez*ez
    vd=2*H.get((0,2),0)+ed-ed*ed
    cov=H.get((1,1),0)-ez*ed
    need(ez==1 and ed==2*k and vz==2*k,('generation moments',k))
    need(vd==F(8,3)*k*k*(k-1)+2*k,('clock variance',k,vd))
    need(cov==F(2,3)*k*(3*k-2),('clock population covariance',k,cov))
    mixed.append({'k':k,'clock_variance':str(vd),'clock_population_covariance':str(cov)})
    H=plus(die,mul(birth,mul(H,mul(H,H))))

print(json.dumps({'status':'PASS','checks':checks,'packet_words':words,'mixed_moments':mixed,
 'source_precision_after_known_odd_bit':precision,'source_mass_numerators':source_counts,
 'source_mass_denominator':2**precision,
 'stopped_words_clock_at_most_16':sum(counts.values()),
 'order_hostile':{'12':p12,'21':p21},
 'tree_joint_clock_examples':[
   {'N':n,'A':a,'word_count':counts[n,a],'mass':str(F(counts[n,a],2**a))}
   for n,a in ((1,1),(1,3),(4,5),(4,7),(7,9),(7,11))]},sort_keys=True,indent=2))
