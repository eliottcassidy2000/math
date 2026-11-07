#!/usr/bin/env python3
"""Literal path enumeration, no producer imports or dynamic programming."""
from itertools import product
import json

def need(p,msg):
    if not p:raise RuntimeError(msg)

def hp(out):
    n=len(out);full=(1<<n)-1
    def walk(u,used):
        if used==full:return 1
        choices=out[u]&~used;s=0
        while choices:
            b=choices&-choices;choices-=b
            s+=walk(b.bit_length()-1,used|b)
        return s
    return sum(walk(u,1<<u) for u in range(n))

def extend(out,masks):
    n=len(out);r=out[:]+[0]*len(masks)
    for j,mask in enumerate(masks):
        v=n+j
        for u in range(n):
            if mask>>u&1:r[u]|=1<<v
            else:r[v]|=1<<u
        for k in range(j+1,len(masks)):r[v]|=1<<(n+k)
    return r

def flip(out,block):
    r=out[:]
    for u in block:
        for v in block:
            if u!=v:r[u]^=1<<v
    return r

one_cases=0
for e in range(4):
    n=e+3
    pairs=[(u,v) for u in range(3,n) for v in range(u+1,n)]
    for code in range(1<<(e+len(pairs))):
        r=[2,4,1]+[0]*e
        for j in range(e):
            for u in range(3):
                if code>>j&1:r[u]|=1<<(3+j)
                else:r[3+j]|=1<<u
        for j,(u,v) in enumerate(pairs):
            if code>>(e+j)&1:r[u]|=1<<v
            else:r[v]|=1<<u
        rf=flip(r,(0,1,2))
        for b in range(1<<n):
            need(hp(extend(r,[b]))==hp(extend(rf,[b])),('module one-port',e,code,b))
            one_cases+=1

# One specific marked two-port local experiment: old->new masks (1,2).
need([hp(extend(r,[1,2])) for r in ([2,4,1],[4,1,2])]==[11,13], 'local13/11')

def words(n):
    if n==0:yield ()
    for a in (1,3):
        if n>=a:
            for w in words(n-a):yield (a,)+w

def ordinal(word,bits):
    n=sum(word);r=[0]*n;cycs=[];start=bitindex=0
    for a in word:
        for u in range(start,start+a):
            for v in range(start+a,n):r[u]|=1<<v
        if a==3:
            for u in range(3):r[start+u]|=1<<(start+(u+1)%3)
            cycs.append(tuple(range(start,start+3)))
            if bits>>bitindex&1:r=flip(r,cycs[-1])
            bitindex+=1
        if a==4:
            for u,v in ((0,1),(1,2),(2,3),(3,0),(2,0),(3,1)):
                r[start+u]|=1<<(start+v)
        start+=a
    return r,cycs

two_cases=0;cubes=0
for n in range(4,11):
    ws=[w for w in words(n) if w<w[::-1]]
    if n%2==0:
        ws += [w+(4,)+w[::-1] for w in words((n-4)//2)]
    for w in ws:
        k=w.count(3)
        if not k:continue
        signatures=set()
        for bits in range(1<<k):
            r,blocks=ordinal(w,bits);observed=[]
            extfactor=3**(k-1)*(5 if 4 in w else 1)
            for j,(u,v,t) in enumerate(blocks):
                # Old factors before the chosen C3 point to both ports.
                prefix=(1<<u)-1
                val=hp(extend(r,[prefix|(1<<u),prefix|(1<<v)]))
                expected=extfactor*(13 if bits>>j&1 else 11)
                need(val==expected,('localized probe',w,bits,j,val,expected))
                observed.append((val//2)%2)
                two_cases+=1
            signatures.add(tuple(observed))
        need(len(signatures)==1<<k,('marked cube injectivity',w))
        cubes+=1

strong4_tables=[]
for code in range(64):
    r=[0]*4;bit=0
    for u in range(4):
        for v in range(u+1,4):
            a,b=(u,v) if code>>bit&1 else (v,u)
            r[a]|=1<<b;bit+=1
    reach=[]
    for u in range(4):
        known=1<<u
        for step in range(4):
            for v in range(4):
                if known>>v&1:known|=r[v]
        reach.append(known)
    if reach==[15]*4:
        strong4_tables.append(tuple(hp(extend(r,[b])) for b in range(16)))
need(len(strong4_tables)==24 and len(set(strong4_tables))==24,'S4 full table injectivity')

print(json.dumps({'status':'PASS','one_port_all_exteriors_through':3,
                  'one_port_exact_mask_checks':one_cases,
                  'local_two_port_counts':[11,13],
                  'leaf_word_max_order':10,'flip_cubes':cubes,
                  'localized_two_port_queries':two_cases,
                  'strong_S4_literal_tables':len(strong4_tables),
                  'distinct_S4_literal_tables':len(set(strong4_tables))},indent=2))
