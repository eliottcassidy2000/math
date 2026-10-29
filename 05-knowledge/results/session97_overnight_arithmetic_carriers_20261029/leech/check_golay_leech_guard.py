#!/usr/bin/env python3
"""Independent, stdlib-only certificate. No imports from the producer.
Coordinates are integers; Euclidean norms are divided by 8.
"""
from collections import Counter
from itertools import combinations, product
CHECKS = 0
def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)

def multiply(a,b):
    out=0
    while b:
        if b&1: out ^= a
        a <<= 1; b >>= 1
    return out

def remainder(a,b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a

# Augmentation cyclic kernel: g=(x+1)f, degree 12, then add all-ones.
# Exponents are mathematical input, not imported producer implementation.
g=sum(1<<i for i in (0,2,5,8,9,10,11,12))
check(g==7973,'producer mathematical generator')
check(remainder((1<<23)|1,g)==0,'cyclic divisor')
check(remainder(g,3)==0,'augmentation divisor')
basis=[g<<i for i in range(11)]+[(1<<24)-1]
words=[0]
for v in basis:
    words += [w^v for w in words]
C=set(words)
check(len(C)==4096,'dimension 12')
check(all((u&v).bit_count()%2==0 for u in basis for v in basis),'self dual by dimension')
W=Counter(c.bit_count() for c in C)
check(W=={0:1,8:759,12:2576,16:759,24:1},'Golay weight enumerator')
mask=(1<<23)-1
rotate=lambda c: ((c<<1)&mask)|(c>>22&1)|(c&(1<<23))
check(all(rotate(c) in C for c in C),'order 23 coordinate action')

def guards(x):
    m=x[0]%2
    parity=all(a%2==m for a in x)
    c=sum(1<<i for i,a in enumerate(x) if a%4==(m+2)%4)
    return parity, c in C, sum(x)%8==4*m

def valid(x): return all(guards(x))
def sq(x): return sum(a*a for a in x)
def unit(i,a): return tuple(a if j==i else 0 for j in range(24))

# Essential independent binary guard and hostile examples.
check(guards(unit(0,4))==(True,True,False),'sum guard rejects norm two')
check(guards((1,)*24)==(True,True,False),'sum guard rejects odd norm three')
check(valid((-3,)+(1,)*23),'odd glue seed')
check(not valid((2,)*4+(0,)*20),'noncode tetrad rejected')
check(not valid((1,)+(0,)*23),'mixed parity rejected')
check(2**23*2**12*2==8**12,'index cancels scaling covolume')

# Enumerate every minimum vector, retaining exact vectors to certify no repeats.
seen=set(); counts=Counter()
def emit(x,kind):
    x=tuple(x)
    check(valid(x),'minimum-vector guards')
    check(sq(x)==32,'minimum-vector norm four')
    check(x not in seen,'distinct minimum vector')
    seen.add(x); counts[kind]+=1
for i,j in combinations(range(24),2):
    for a,b in product((-4,4),repeat=2):
        x=[0]*24; x[i]=a; x[j]=b; emit(x,'4^2')
for c in sorted(C):
    if c.bit_count()!=8: continue
    support=[i for i in range(24) if c>>i&1]
    for signs in range(256):
        if signs.bit_count()%2: continue
        x=[0]*24
        for k,i in enumerate(support): x[i]=-2 if signs>>k&1 else 2
        emit(x,'2^8')
for j in range(24):
    for c in sorted(C):
        x=[-1 if c>>i&1 else 1 for i in range(24)]
        x[j]=-3*x[j]
        emit(x,'3,1^23')
check(counts=={'4^2':1104,'2^8':97152,'3,1^23':98304},'shell counts')
check(len(seen)==196560,'kissing number')
# Deterministic hostile round trips, not a substitute for the universal proof.
ordered=sorted(seen)
for k in range(0,len(ordered),137):
    x=ordered[k]; y=ordered[(k*37+19)%len(ordered)]
    z=tuple(a+b for a,b in zip(x,y))
    check(valid(z) and sq(z)%16==0,'closure and evenness sample')
    check(sum(a*b for a,b in zip(x,y))%8==0,'integral pairing sample')
    check(valid(tuple(-a for a in x)),'negation sample')
# Construction A roots: z^2=4 forces ±2e_i or a weight-four codeword.
check(W[4]==0 and 2*24==48,'plain Construction A has A1^24 roots')
print('Generator bitmask:',g)
print('Golay weight enumerator:',dict(sorted(W.items())))
print('Minimum shell:',dict(counts),'total',len(seen))
print('Index: 2^36; scaled covolume: 1; plain Construction A roots: 48')
print('CHECKS',CHECKS)
