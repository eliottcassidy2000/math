#!/usr/bin/env python3
"""Exact, standalone Hamilton-frame and C7-module inventory.
Python standard library only. All validations survive python -O.
The complete q=7 Hamilton-frame scan is the finite basis for the proved
all-order symmetric 7p inventory below 19467; no random sampling is used.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, isqrt
import hashlib
import json
from pathlib import Path

LIMIT = 19467
C_EXPECTED = [0,175,1638,6678,15120,20160,15120,5040]
EXPECTED = {4151:(10,2),5789:(10,2),14791:(12,2),16079:(12,4),17717:(12,2)}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def decode(n, code):
    r=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:r[i]|=1<<j
        else:r[j]|=1<<i
    return tuple(r)

def frame(n, code):
    r=[1<<(i+1) if i+1<n else 0 for i in range(n)]
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    for k,(i,j) in enumerate(free):
        if code>>k&1:r[i]|=1<<j
        else:r[j]|=1<<i
    return tuple(r)

def strong(r):
    n=len(r)
    a=[x|1<<i for i,x in enumerate(r)]
    for k in range(n):
        for i in range(n):
            if a[i]>>k&1:a[i]|=a[k]
    return all(x==(1<<n)-1 for x in a)

def induced(r,vs):
    return tuple(sum(1<<j for j,v in enumerate(vs) if r[u]>>v&1) for u in vs)

def ham_subsets(r):
    n=len(r); dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i,a in enumerate(dp[s]):
            if a:
                rest=r[i]&~s
                while rest:
                    bit=rest&-rest;rest-=bit;j=bit.bit_length()-1
                    dp[s|bit][j]+=a
        h[s]=sum(dp[s])
    return h

def ham_recursive(r, root=None):
    n=len(r)
    @lru_cache(None)
    def go(mask,last):
        if not mask:return 1
        rest=r[last]&mask;out=0
        while rest:
            bit=rest&-rest;rest-=bit;j=bit.bit_length()-1
            out+=go(mask^bit,j)
        return out
    full=(1<<n)-1
    if root is not None:return go(full^(1<<root),root)
    return sum(go(full^(1<<j),j) for j in range(n))

def covers(r):
    hs=ham_subsets(r);n=len(r);dp=[[0]*(n+1) for _ in hs];dp[0][0]=1
    for s in range(1,1<<n):
        lo=s&-s;t=s
        while t:
            if t&lo:
                for k in range(1,n+1):dp[s][k]+=hs[t]*dp[s^t][k-1]
            t=(t-1)&s
    return [factorial(k)*a for k,a in enumerate(dp[-1])]

def covers_literal(r):
    n=len(r);out=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum(1<<i for i in range(n-1) if not(r[p[i]]>>p[i+1]&1))
        for cuts in range(1<<(n-1)):
            if cuts&bad==bad:out[1+cuts.bit_count()]+=1
    return out

def kernel(r):
    # A_k counts walks visiting root0 k times and every other vertex once.
    # No two root symbols are consecutive. A_0 is retained as a control.
    b=len(r)-1
    @lru_cache(None)
    def go(mask,last):
        out=[0]*(mask.bit_count()+2)
        if not mask:out[0]=1
        rest=(r[last]>>1)&mask
        while rest:
            bit=rest&-rest;rest-=bit;j=bit.bit_length()
            sub=go(mask^bit,j)
            for k,v in enumerate(sub):out[k]+=v
        if last and r[last]&1:
            sub=go(mask,0)
            for k,v in enumerate(sub[:-1]):out[k+1]+=v
        return tuple(out)
    mask=(1<<b)-1;out=[0]*(b+2)
    if not b:out[0]=1
    sub=go(mask,0)
    for k,v in enumerate(sub[:-1]):out[k+1]+=v
    for j in range(1,b+1):
        sub=go(mask^(1<<(j-1)),j)
        for k,v in enumerate(sub):out[k]+=v
    return out

def kernel_literal(r):
    b=len(r)-1;out=[0]*(b+2)
    if not b:out[0]=1
    # A root may be inserted in any subset of the b+1 gaps of a fixed
    # permutation of nonroot vertices. This is independent of the DP.
    for p in permutations(range(1,b+1)):
        for mask in range(1<<(b+1)):
            w=[]
            for i in range(b+1):
                if mask>>i&1:w.append(0)
                if i<b:w.append(p[i])
            if w and all(r[u]>>v&1 for u,v in zip(w,w[1:])):
                out[mask.bit_count()]+=1
    return out

def substitute(q, m):
    n=len(m);r=list(m)+[0]*(len(q)-1)
    for i in range(n):
        for j in range(1,len(q)):
            if q[0]>>j&1:r[i]|=1<<(n+j-1)
            else:r[n+j-1]|=1<<i
    for i in range(1,len(q)):
        for j in range(1,len(q)):
            if q[i]>>j&1:r[n+i-1]|=1<<(n+j-1)
    return tuple(r)

def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def automorphisms(r):
    n=len(r)
    return sum(all(((r[i]>>j)&1)==((r[p[i]]>>p[j])&1)
                   for i in range(n) for j in range(n))
               for p in permutations(range(n)))

def canon_root(r):
    n=len(r)
    return min(tuple(sum(((r[p[i]]>>p[j])&1)<<j for j in range(n)) for i in range(n))
               for tail in permutations(range(1,n)) for p in [(0,)+tail])

def run():
    module=tuple(sum(1<<((i+j)%7) for j in (1,2,3)) for i in range(7))
    c=covers(module)
    for q in ((0,),(2,0),(0,1)):
        need(kernel(q)==kernel_literal(q),'empty and singleton exterior convention')
        need(kernel(q)[0]==ham_subsets(induced(q,list(range(1,len(q)))))[-1],'A0 convention')
    need(c==C_EXPECTED and covers_literal(module)==c,'module covers')
    orbit_table=[];regular_sets=set()
    for bits in range(8):
        steps=[j if bits>>(j-1)&1 else 7-j for j in range(1,4)]
        r=tuple(sum(1<<((i+j)%7) for j in steps) for i in range(7))
        h=ham_subsets(r)[-1];aut=automorphisms(r)
        orbit_table.append([bits,h,aut])
        if h==175:regular_sets.add(frozenset(steps))
    need(Counter((h,a) for _,h,a in orbit_table)=={(175,7):6,(189,21):2},'C7 orbit table')
    need(regular_sets=={frozenset(u*d%7 for d in (1,2,3)) for u in range(1,7)},'one ordinary175 circulant class')
    frame_data=[];target_witnesses={};min_witnesses={};kernel_controls=0
    observed=Counter();target_roots={}
    for n in range(3,8):
        total=1<<((n-1)*(n-2)//2);mn=None;strong_count=0;target=Counter();prime_frames=Counter();front=[]
        for code in range(total):
            r=frame(n,code)
            if not strong(r):continue
            strong_count+=1
            a=kernel(r);h=sum(x*y for x,y in zip(a,c))
            if a[1]>3 and prime(a[1]):prime_frames[a[1]]+=1
            positive=tuple(a[1:])
            if not any(all(x<=y for x,y in zip(v,positive)) for v,_ in front):
                front=[(v,w) for v,w in front if not all(x<=y for x,y in zip(positive,v))]
                front.append((positive,code))
            if mn is None or h<mn:mn=h;min_witnesses[n]=(code,r,a,h)
            # Independently read selected full kernels by permutation/gap words.
            if n<=5 or (n>=6 and code%(total//8)==0):
                need(a==kernel_literal(r),('kernel literal',n,code));kernel_controls+=1
            if h<LIMIT and h%7==0 and prime(h//7) and h//7>=11:
                s=ham_recursive(r,0);need(s>0,'root start')
                target[h]+=Fraction(1,s)
                target_witnesses.setdefault(h,(n,code,r,a))
                # Exact root-preserving canonicalization of successful q<=6
                # cases is an independent check of the reciprocal weighting.
                target_roots.setdefault((n,h),set()).add(canon_root(r))
        need(mn is not None,'strong universe')
        for h,v in target.items():
            need(v.denominator==1,('integer root count',n,h,v))
            need(v==len(target_roots[(n,h)]),('root orbit weight',n,h))
            observed[(h,n+6)]+=int(v)
        frame_data.append(dict(q=n,frames=total,strong_frames=strong_count,
            minimum=mn,targets={str(h):int(v) for h,v in sorted(target.items())},
            pareto_kernels=[dict(coefficients=v,frame=w) for v,w in sorted(front)],
            prime_H_strong_classes={str(h):v//h for h,v in sorted(prime_frames.items())}))
        need(all(v%h==0 for h,v in prime_frames.items()),'prime frame divisibility')
    need([x['minimum'] for x in frame_data]==[2163,4151,6489,10815,19467],'minima')
    need(dict(observed)=={(h,n):k for h,(n,k) in EXPECTED.items()},'complete strong symmetric inventory')
    # Independent labelled quotient universe throughq6, divided by (q-1)!.
    labelled=[]
    for n in range(3,7):
        target=Counter();mn=None;prime_labels=Counter()
        for code in range(1<<(n*(n-1)//2)):
            r=decode(n,code)
            if not strong(r):continue
            a=kernel(r);h=sum(x*y for x,y in zip(a,c))
            if a[1]>3 and prime(a[1]):prime_labels[a[1]]+=1
            mn=h if mn is None else min(mn,h)
            if h<LIMIT and h%7==0 and prime(h//7) and h//7>=11:target[h]+=1
        normalized={str(h):v//factorial(n-1) for h,v in sorted(target.items())}
        need(all(v%factorial(n-1)==0 for v in target.values()),'labelled quotient divisibility')
        need(normalized==frame_data[n-3]['targets'],'labelled/frame inventory')
        need(mn==frame_data[n-3]['minimum'],'labelled/frame minimum')
        need(all(v%factorial(n)==0 for v in prime_labels.values()),'prime labelled divisibility')
        need({str(h):v//factorial(n) for h,v in sorted(prime_labels.items())}==
             frame_data[n-3]['prime_H_strong_classes'],'prime labelled/frame census')
        labelled.append(dict(q=n,graphs=1<<(n*(n-1)//2),targets=normalized))
    witnesses=[]
    for h,(n,code,r,a) in sorted(target_witnesses.items()):
        actual=substitute(r,module)
        need(strong(actual),'strong lifted witness')
        need(ham_subsets(actual)[-1]==ham_recursive(actual)==h,'actual witness DP')
        witnesses.append(dict(H=h,p=h//7,order=len(actual),multiplicity=EXPECTED[h][1],
                              quotient_frame=code,quotient_rows=r,kernel=a,rows=actual))
    min_controls=[]
    for n,(code,r,a,h) in sorted(min_witnesses.items()):
        actual=substitute(r,module)
        need(ham_subsets(actual)[-1]==ham_recursive(actual)==h,'minimum actual witness')
        min_controls.append(dict(q=n,frame=code,rows=r,kernel=a,H=h))
    # All-order formula's finite hostile: identity-only weight at H4151 omits6/7
    # per symmetric class. Complete fixed C7 type with b=3 has72 graphs.
    fixed=0
    for code in range(1<<6):
        q=decode(4,code)
        if sum(x*y for x,y in zip(kernel(q),c))==4151:fixed+=6
    need(fixed==72 and Fraction(fixed,7*factorial(3))==Fraction(12,7),'C7 Burnside correction')
    # Induced root-quotient monotonicity, with root kept, has an explicit
    # coefficientwise insertion injection. Check all strong frameq<=6.
    monotone=0
    for n in range(3,7):
        for code in range(1<<((n-1)*(n-2)//2)):
            r=frame(n,code)
            if not strong(r):continue
            a=kernel(r)
            for v in range(1,n):
                rr=induced(r,[i for i in range(n) if i!=v]);b=kernel(rr)
                need(all(x<=a[k] for k,x in enumerate(b)),('kernel monotonicity',n,code,v))
                monotone+=1
    return dict(status='PASS',scope='complete symmetric strong 7p inventory for prime p>=11 and7p<19467; all-order proof in NOTE.md',
        module_rows=module,module_ordered_covers=c,circulant7_table=orbit_table,
        frame_universes=frame_data,labelled_controls=labelled,
        complete_symmetric_inventory=witnesses,minimum_controls=min_controls,
        independent_literal_kernel_controls=kernel_controls,
        coefficientwise_monotonicity_controls=monotone,
        fixed_C7_H4151=72,burnside_correction_H4151='12/7',
        hostile_identity_only_contribution_from_symmetric_H4151_classes='2/7 instead of2',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    print(json.dumps(run(),indent=2,sort_keys=True))
