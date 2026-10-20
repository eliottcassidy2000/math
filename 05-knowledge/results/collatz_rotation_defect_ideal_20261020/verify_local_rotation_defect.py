#!/usr/bin/env python3
"""Exact, self-contained controls for LOCAL_ROTATION_DEFECT_IDEAL.md."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import gcd
import json

counts=Counter()
digest=sha256()

def require(ok, message):
    counts['explicit_checks']+=1
    if not ok:
        raise RuntimeError(message)

@lru_cache(maxsize=400000)
def carry(w, base=3):
    b=0; modulus=1
    for bit in w:
        if bit=='1': b=base*b+modulus
        modulus*=2
    return b

def carry_sum(w,base=3):
    return sum((1<<i)*base**w[i+1:].count('1') for i,c in enumerate(w) if c=='1')

def decode(b,q,p,base=3):
    out=[]
    for _ in range(q):
        bit=b%2;out.append(str(bit))
        b=(b-bit*base**(p-bit))//2
        p-=bit
    return ''.join(out),b,p

def core(a,b):
    lo=0
    while a[lo]==b[lo]:lo+=1
    hi=len(a)
    while a[hi-1]==b[hi-1]:hi-=1
    return lo,hi,a[lo:hi],b[lo:hi]

def defect(a,b,base=3):
    lo,hi,u,v=core(a,b)
    require(len(u)==len(v) and u.count('1')==v.count('1'),'core shape')
    d=carry(u,base)-carry(v,base)
    require(d!=0,'injectivity at core')
    require(carry(a,base)-carry(b,base)==(1<<lo)*base**a[hi:].count('1')*d,'context factor')
    ell=len(u);h=u.count('1')
    bound=((1<<(ell-h))-1)*(base**h-(1<<h))//(base-2)
    require(0<abs(d)<=bound,'core bound')
    counts['core_edges']+=1
    return d

def check_word(w,base=3,all_pairs=False):
    q=len(w);p=w.count('1');D=(1<<q)-base**p;B=carry(w,base)
    require(B==carry_sum(w,base),'independent carry sum')
    require(decode(B,q,p,base)==(w,0,0),'carry decode')
    rs=list(dict.fromkeys(w[i:]+w[:i] for i in range(q)))
    bs=[carry(r,base) for r in rs]
    g=gcd(abs(D),B)
    require(gcd(g,2*base)==1,'content unit condition')
    for i in range(q):
        a=w[i:]+w[:i];b=w[(i+1)%q:]+w[:(i+1)%q]
        require(2*carry(b,base)==base**int(w[i])*carry(a,base)+int(w[i])*D,'rotation recurrence')
    gs=0;gp=0
    for i in range(1,len(rs)):
        gs=gcd(gs,abs(defect(rs[0],rs[i],base)))
        gp=gcd(gp,abs(defect(rs[i-1],rs[i],base)))
    require(gs==gp==g,'star/path exact local gcd')
    if all_pairs:
        ga=0
        for i in range(len(rs)):
            for j in range(i):ga=gcd(ga,abs(defect(rs[i],rs[j],base)))
        require(ga==g,'complete graph gcd')
    digest.update(f'{base}:{w}:{D}:{B}:{g}:{gs}:{gp}\n'.encode())
    counts['word_certificates']+=1
    return g

def necklaces(n):
    """FKM generation: return one least rotation and its least word period."""
    a=[0]*(n+1)
    def visit(t,p):
        if t>n:
            if n%p==0:yield ''.join(map(str,a[1:])),p
        else:
            a[t]=a[t-p]
            yield from visit(t+1,p)
            for j in range(a[t-p]+1,2):
                a[t]=j
                yield from visit(t+1,t)
    yield from visit(1,1)

all_word_counts={}
for q in range(2,13):
    seen={};count=0
    for z in range(1,(1<<q)-1):
        w=format(z,f'0{q}b');check_word(w,all_pairs=(q<=8));count+=1
        key=(w.count('1'),carry(w))
        require(key not in seen,'shape carry collision')
        seen[key]=w
    all_word_counts[str(q)]=count

necklace_counts={};nonunit_counts={}
for q in range(2,19):
    count=nonunit=0
    for w,period in necklaces(q):
        if period!=q:continue
        g=check_word(w,all_pairs=(q<=12));count+=1;nonunit+=g>1
    necklace_counts[str(q)]=count;nonunit_counts[str(q)]=nonunit

# Independent extremal scan, all shapes through 12.
extremal_shapes=0
for q in range(2,13):
    byp={}
    for z in range(1,(1<<q)-1):
        w=format(z,f'0{q}b');p=w.count('1');b=carry_sum(w)
        byp.setdefault(p,[]).append((b,w))
    for p,items in byp.items():
        low=min(items);high=max(items)
        require(low==((3**p-(1<<p)),'1'*p+'0'*(q-p)),'unique minimum')
        require(high==((1<<(q-p))*(3**p-(1<<p)),'0'*(q-p)+'1'*p),'unique maximum')
        extremal_shapes+=1

# Odd multiplier generalization; denominator 2^q-a^p is coprime to 2a.
odd_base_counts={}
for base in (5,7,9):
    n=0
    for q in range(2,9):
        for w,period in necklaces(q):
            if w.count('1') in (0,q):continue
            check_word(w,base);n+=1
    odd_base_counts[str(base)]=n

# Word-power invariance with repeated rotations retained as one vertex.
power_controls=0
for q in range(2,8):
    for w,period in necklaces(q):
        if period!=q:continue
        D=(1<<q)-3**w.count('1');B=carry(w);g=gcd(abs(D),B)
        for k in (2,3,4):
            v=w*k;Dk=(1<<len(v))-3**v.count('1');H=Dk//D
            require(Dk==D*H and H>0,'positive geometric repeat factor')
            require(carry(v)==H*B and gcd(abs(Dk),carry(v))==H*g,'power content')
            require(abs(Dk)//gcd(abs(Dk),carry(v))==abs(D)//g,'same reduced denominator')
            check_word(v);power_controls+=1

# Sparse unbalanced witness, no one pair has 2-3-smooth full difference.
w='0011011';rs=[w[i:]+w[:i] for i in range(len(w))]
d1=defect(rs[4],rs[0]);d2=defect(rs[3],rs[0])
require((d1,d2)==(-5,-63) and 25*d1-2*d2==1,'two-edge Bezout witness')
reduced=[]
for i in range(len(rs)):
    for j in range(i):
        z=abs(defect(rs[i],rs[j]))
        for p in (2,3):
            while z%p==0:z//=p
        reduced.append(z)
require(min(reduced)>1,'no individual smooth pair')

# Actual signed orbit replay, no convergence assumption.
def replay_cycle(start):
    x=start;bits=[];nodes=[]
    for _ in range(100):
        nodes.append(x);bits.append(str(x%2));x=(3*x+1)//2 if x%2 else x//2
        if x==start:return ''.join(bits),nodes
    raise RuntimeError('control did not close')
cycle_controls=[]
for start in (1,-5,-17):
    w,nodes=replay_cycle(start);D=(1<<len(w))-3**w.count('1');B=carry(w)
    require(B==start*D,'actual cycle fixed point')
    require(check_word(w)==abs(D),'integer cycle has full content')
    cycle_controls.append({'start':start,'word':w,'D':D,'B':B,'g':abs(D),'nodes':nodes})

out={
 'status':'PASS',
 'universe':{
  'all_nonconstant_binary_words_by_length':all_word_counts,
  'primitive_necklaces_by_length':necklace_counts,
  'primitive_necklaces_nonunit_content_by_length':nonunit_counts,
  'all_pairs_checked':'all binary words length <=8 and primitive necklaces length <=12, plus explicit sparse witness',
  'star_and_path_graphs':'every tested word, including primitive necklaces through length 18',
  'odd_base_necklaces_through_length_8':odd_base_counts,
  'word_power_controls':power_controls,
  'extremal_shapes_through_length_12':extremal_shapes,
 },
 'counts':dict(sorted(counts.items())),
 'certificate_digest_sha256':digest.hexdigest(),
 'sparse_unbalanced_witness':{'word':'0011011','D':47,'B':340,'g':1,'local_differences':[d1,d2],'bezout_coefficients':[25,-2],'minimum_23_stripped_pair_defect':min(reduced)},
 'nonunit_hostile':{'word':'000011','D':55,'B':80,'g':5,'reduced_denominator':11},
 'power_hostile':{'word':'0101','D':7,'B':14,'g':7,'reduced_denominator':1},
 'integer_cycle_controls':cycle_controls,
 'boundary':'Finite controls verify the formulas on the stated universes; the written proof supplies the all-word theorem. No all-positive-cycle exclusion or Collatz convergence is claimed.'
}
print(json.dumps(out,indent=2,sort_keys=True))
