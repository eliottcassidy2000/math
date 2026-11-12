#!/usr/bin/env python3
"""Session283 physical superincreasing profile code, dependency-free controls."""
import hashlib
import itertools
import json
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

def require(test, message):
    if not test:
        raise RuntimeError(message)

def tournament(n, code):
    rows = [0]*n
    for k, (i,j) in enumerate(itertools.combinations(range(n),2)):
        if code >> k & 1: rows[i] |= 1<<j
        else: rows[j] |= 1<<i
    return tuple(rows)

def subset_counts(rows):
    n=len(rows); sz=1<<n
    dp=[[0]*n for _ in range(sz)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,sz):
        for v,z in enumerate(dp[mask]):
            if not z: continue
            opts=rows[v]&~mask
            while opts:
                bit=opts&-opts; opts-=bit
                dp[mask|bit][bit.bit_length()-1] += z
    return [sum(v) for v in dp]

def profile(rows):
    n=len(rows); h=subset_counts(rows)
    @lru_cache(None)
    def f(mask):
        if not mask: return (1,)+(0,)*n
        out=[0]*(n+1); anchor=mask&-mask; part=mask
        while part:
            if part&anchor:
                for k,z in enumerate(f(mask^part)[:-1]): out[k+1]+=h[part]*z
            part=(part-1)&mask
        return tuple(out)
    return tuple(factorial(k)*v for k,v in enumerate(f((1<<n)-1)))[1:]

def independent_profile(rows):
    n=len(rows); out=[0]*n
    for p in itertools.permutations(range(n)):
        b=sum(not (rows[u]>>v&1) for u,v in zip(p,p[1:]))
        for r in range(b+1,n+1): out[r-1]+=comb(n-1-b,r-1-b)
    return tuple(out)

def transitive(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def triangle(blocks):
    sizes=list(map(len,blocks)); offsets=[0,sizes[0],sizes[0]+sizes[1]]
    rows=[]
    for i,block in enumerate(blocks):
        j=(i+1)%3; targets=((1<<sizes[j])-1)<<offsets[j]
        rows.extend((r<<offsets[i])|targets for r in block)
    return tuple(rows)

def cap_weights(n,m):
    # g_s=s! S(m,s); g_s(m+1)=s(g_s(m)+g_(s-1)(m)).
    g=[1]+[0]*(n+1)
    for _ in range(m): g=[0]+[s*(g[s]+g[s-1]) for s in range(1,n+2)]
    return tuple(3*g[r]**2+2*g[r]*g[r+1]+g[r+1]**2+g[r-1]**2+2*g[r-1]*g[r] for r in range(1,n+1))

def bound(n):
    B=factorial(n)*(1<<(n-1)); target=72*(n+1)*B
    L=(target-1).bit_length()
    return B,L,n*L

def first_digit_box_cap(n):
    digits=[factorial(n)*comb(n-1,r-1) for r in range(1,n+1)]
    limit=bound(n)[2]
    for m in range(1,limit+1):
        weights=cap_weights(n,m)
        if all(weights[r]>sum(digits[s]*weights[s] for s in range(r)) for r in range(n)):
            return m, weights
    raise RuntimeError('proved cap bound not found')

def encode(F,weights): return sum(a*b for a,b in zip(F,weights))
def decode(z,weights):
    digits=[]
    for w in weights[::-1]:
        a,z=divmod(z,w); digits.append(a)
    require(z==0,'nonzero terminal remainder')
    return tuple(digits[::-1])

def attach(rows,mask):
    n=len(rows); out=list(rows)+[mask]
    for i in range(n):
        if not(mask>>i&1): out[i] |= 1<<n
    return tuple(out)

def main():
    summary={}; all_graphs=0; perm_checks=0; cap_checks=0; code_checks=0
    for n in range(1,6):
        reps={}; B,L,m=bound(n); weights=cap_weights(n,m)
        require(all(weights[r]>B*sum(weights[:r]) for r in range(n)),'superincreasing failure')
        for code in range(1<<(n*(n-1)//2)):
            rows=tournament(n,code); F=profile(rows); all_graphs+=1
            require(all(0<=v<=B for v in F),'profile digit bound')
            require(F==independent_profile(rows),'independent cut profile mismatch'); perm_checks+=1
            z=encode(F,weights); require(decode(z,weights)==F,'decoding failure'); code_checks+=1
            reps.setdefault(F,rows)
        ordered=sorted(reps,key=lambda f:f[::-1]); codes=[encode(f,weights) for f in ordered]
        require(all(a<b for a,b in zip(codes,codes[1:])),'reverse lex failure')
        for F,rows in reps.items():
            for c in range(1,3):
                literal=subset_counts(triangle([rows,transitive(c),transitive(c)]))[-1]
                require(literal==encode(F,cap_weights(n,c)),'literal cap failure'); cap_checks+=1
        summary[str(n)]={'labelled_tournaments':1<<(n*(n-1)//2),'profile_classes':len(reps),'safe_cap_m':m,'max_encoded_bits':max(codes).bit_length()}
    formal_checks=0
    for n in range(1,21):
        B,L,m=bound(n); weights=cap_weights(n,m)
        require(all(weights[r]>B*sum(weights[:r]) for r in range(n)),'safe bound failure')
        samples=[tuple([B]*n),tuple([0]*n),tuple(B if r%2 else 0 for r in range(n)),tuple((r*7357+19)%(B+1) for r in range(n))]
        for F in samples:
            require(decode(encode(F,weights),weights)==F,'formal digit failure'); formal_checks+=1
    A=(14,28,24,16,1); B=(10,28,25,16,1)
    require(profile(A)==profile(B)==(9,66,186,240,120),'H9 profile mismatch')
    menus=[]; signatures=[]
    for rows in [A,B]:
        degrees=[r.bit_count() for r in rows]; v=degrees.index(2)
        signatures.append(sorted(degrees[j] for j in range(5) if rows[v]>>j&1))
        menu={}
        for mask in range(32):
            h=subset_counts(attach(rows,mask))[-1]; menu.setdefault(str(h),[]).append(mask)
        menus.append(menu)
    require(signatures==[[1,1],[1,3]],'isomorphism separator failure')
    require('31' not in menus[0] and 4 in menus[1]['31'],'attachment loss control failure')
    digit_box_caps={}
    for n in [1,2,3,4,5,6,8,10,16,20]:
        m,weights=first_digit_box_cap(n)
        F=tuple(factorial(n)*comb(n-1,r-1) for r in range(1,n+1))
        require(decode(encode(F,weights),weights)==F,'digit-box cap decoding failure')
        digit_box_caps[str(n)]={'first_sufficient_cap':m,'universal_safe_cap':bound(n)[2]}
    receipt={'status':'finite controls pass; general theorem proved separately','universe':summary,'labelled_profiles':all_graphs,'independent_permutation_cut_profiles':perm_checks,'literal_cap_counts':cap_checks,'physical_code_roundtrips':code_checks,'formal_digit_roundtrips':formal_checks,'safe_weight_checks_orders':list(range(1,21)),'H9_profile':profile(A),'H9_attachment_menus':menus,'H9_nonisomorphism_signatures':signatures,'digit_box_caps':digit_box_caps,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
