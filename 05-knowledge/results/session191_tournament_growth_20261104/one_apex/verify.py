#!/usr/bin/env python3
"""Standalone exact checks for the Session 191 one-apex order-window theorem."""
import hashlib
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path

HERE=Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def word_count(w):
    h,a,p=1,0,0
    for b in w:
        if b:
            a,p=2*a+1,1
        else:
            h,a,p=h+a+p,2*a,0
    return h

def slot_count(w):
    if not w:return 1
    return (w[0]+1-w[-1]+sum(w[i]==0 and w[i+1]==1 for i in range(len(w)-1))
            +sum((1<<(i-j-1)) for j in range(len(w)) for i in range(j+1,len(w)) if w[j] and not w[i]))

def rows(w):
    n=len(w);r=[0]*(n+1)
    for i,b in enumerate(w,1):
        if b:r[0]|=1<<i
        else:r[i]|=1
        for j in range(i+1,n+1):r[i]|=1<<j
    return r

def dp_count(r):
    n=len(r);full=(1<<n)-1
    if not n:return 1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,full+1):
        for i in range(n):
            ways=dp[mask][i]
            if not ways:continue
            avail=r[i]&~mask
            while avail:
                bit=avail&-avail;avail-=bit;j=bit.bit_length()-1
                dp[mask|bit][j]+=ways
    return sum(dp[full])

def perm_count(r):
    return sum(all(r[p[i]]>>p[i+1]&1 for i in range(len(r)-1)) for p in itertools.permutations(range(len(r))))

def is_strong(r):
    if len(r)<2:return False
    full=(1<<len(r))-1
    for i in range(len(r)):
        seen=1<<i
        while True:
            nxt=seen
            for j in range(len(r)):
                if seen>>j&1:nxt|=r[j]
            if nxt==seen:break
            seen=nxt
        if seen!=full:return False
    return True

def low(n):
    if n in (2,3):return 2*n-1
    return ((1<<((n-2)//2))+1)*((1<<((n-1)//2))+1)

def high(n):
    return ((1<<(n//2))-1)*((1<<((n+1)//2))-1)+2

def low2(n):
    if n==4:return 11
    if n==5:return 17
    m=n//2
    return low(n)+(1 if n%2==0 else 3)*(1<<(m-2))

def high2(n):
    m=n//2
    return high(n)-((1<<(m-1)) if n%2==0 else 3*(1<<(m-1))-2)

def swap_min(n):
    if n==4:return 2
    if n==5:return 4
    return 2+(1<<((n-2)//2-1))+(1<<((n-1)//2-1))

def strong_words(n):
    for middle in itertools.product((0,1),repeat=n-2):
        yield (1,)+middle+(0,)

def string(w):return ''.join(map(str,w))

def bits(a):return tuple(map(int,bin(a)[2:]))
def zero_support(a):
    return tuple(1-int(b) for b in bin(a)[2:][::-1]) if a else ()
def divisors(a):
    ds=[]
    for d in range(1,math.isqrt(a)+1):
        if a%d==0:
            ds.append(d)
            if d*d!=a:ds.append(a//d)
    return sorted(ds)
def odd_part(a):
    while a%2==0:a//=2
    return a
@lru_cache(None)
def constants(a):
    m=odd_part(a)
    return tuple(sorted(set(word_count(bits(d))-d+word_count(zero_support(m//d-1)) for d in divisors(m))))

checks={"slot":0,"dp":0,"permutation":0,"strongness":0,"swaps":0,"fixed_weight":0,"first_second_endpoints":0,"strip":0,"single_run_formula":0,"dyadic_slope_rigidity":0,"single_run_classification":0,"divisor_ray_witnesses":0}
# Independent vertex computations, including nonstrong/padded words.
for n in range(11):
    for w in itertools.product((0,1),repeat=n):
        h=word_count(w);r=rows(w)
        require(h==slot_count(w),f"slot {w}");checks['slot']+=1
        require(h==dp_count(r),f"DP {w}");checks['dp']+=1
        if n<=6:
            require(h==perm_count(r),f"permutations {w}");checks['permutation']+=1
        if n>=2:
            require(is_strong(r)==bool(w[0] and not w[-1]),f"strong {w}");checks['strongness']+=1
        a=0;b=n
        while a<b and w[a]==0:a+=1
        while a<b and w[b-1]==1:b-=1
        require(h==word_count(w[a:b]),f"strip {w}");checks['strip']+=1

profiles=[];all_small={1};small_witnesses={"1":""}
for n in range(2,19):
    vals=set();by_weight={};best={};mindelta=None
    for w in strong_words(n):
        h=word_count(w);vals.add(h);o=sum(w)
        by_weight.setdefault(o,[]).append(h)
        best.setdefault(h,string(w))
        if h<=500:
            all_small.add(h);small_witnesses.setdefault(str(h),string(w))
        if n<=14:
            require(h==slot_count(w),f"strong slot {w}");checks['slot']+=1
        if n<=16:
            for i in range(1,n-2):
                if w[i:i+2]!=(0,1):continue
                new=w[:i]+(1,0)+w[i+2:]
                delta=word_count(new)-h
                u=w[:i];v=w[i+2:];A=int(string(u),2)
                expected=2*(1+A//2+sum(1<<(j-2) for j in range(2,len(v)+1) if not v[j-1]))
                require(delta==expected and delta>=2,f"swap {w} {i}")
                mindelta=delta if mindelta is None else min(mindelta,delta)
                checks['swaps']+=1
    for o,hs in sorted(by_weight.items()):
        z=n-o
        L=((1<<(o-1))+1)*((1<<(z-1))+1) if min(o,z)>=2 else (1<<(n-1))+1
        U=((1<<o)-1)*((1<<z)-1)+2
        require(min(hs)==L and max(hs)==U,f"weight {n} {o}")
        lo=(1,)+(0,)*(z-1)+(1,)*(o-1)+(0,)
        hi=(1,)*o+(0,)*z
        require(word_count(lo)==L and word_count(hi)==U,f"literal extrema {n} {o}")
        checks['fixed_weight']+=1
    vs=sorted(vals)
    require(vs[0]==low(n) and vs[-1]==high(n),f"global extrema {n}")
    rec={"chain_length":n,"tournament_order":n+1,"distinct_counts":len(vs),"minimum":vs[0],"maximum":vs[-1],"minimum_word":best[vs[0]],"maximum_word":best[vs[-1]]}
    if n>=4:
        require(vs[1]==low2(n) and vs[-2]==high2(n),f"second extrema {n}")
        rec.update(second_minimum=vs[1],second_maximum=vs[-2],second_minimum_word=best[vs[1]],second_maximum_word=best[vs[-2]])
        checks['first_second_endpoints']+=1
    if n>=4 and n<=16:
        require(mindelta==swap_min(n),f"swap minimum {n}")
        rec['minimum_internal_swap_increment']=mindelta
    profiles.append(rec)

# Exact one-run affine laws and their strong dyadic-slope rigidity.
short_words=[w for n in range(5) for w in itertools.product((0,1),repeat=n)]
for u in short_words:
    A=int(string(u),2) if u else 0
    p=u[-1] if u else 0
    for v in short_words:
        C=sum(1<<j for j,b in enumerate(v) if not b)
        eps=int(bool(v) and not v[0])
        for b in (0,1):
            if b==0:
                slope=A*(1+C);constant=word_count(u)-A+p+word_count(v)-1
            else:
                slope=(A+1)*C;constant=word_count(u)+word_count(v)-1-C+eps
            for k in range(1,7):
                w=u+(b,)*k+v
                require(word_count(w)==slope*(1<<k)+constant,"single run law")
                checks['single_run_formula']+=1
            w=u+(b,)+v
            if len(w)>=2 and w[0] and not w[-1]:
                require(slope>0 and constant in constants(slope),"single run classification")
                checks['single_run_classification']+=1
            if len(w)>=2 and w[0] and not w[-1] and slope>0 and slope&(slope-1)==0:
                require(constant==1 and (sum(w)==1 or len(w)-sum(w)==1),"dyadic-slope rigidity")
                checks['dyadic_slope_rigidity']+=1

for m in range(1,256,2):
    for d in divisors(m):
        e=m//d
        constant=word_count(bits(d))-d+word_count(zero_support(e-1))
        for t in range(5):
            for r in range(t+1):
                u=bits(d)+(0,)*r
                v=(0,)*(t-r)+zero_support(e-1)
                slope=m*(1<<t)
                require(constants(slope)==constants(m),"odd part determines constants")
                for k in range(1,5):
                    w=u+(0,)*k+v
                    require(w[0] and not w[-1],"divisor witness strong endpoints")
                    require(word_count(w)==slope*(1<<k)+constant,"divisor witness count")
                    revcomp=tuple(1-b for b in reversed(w))
                    require(word_count(revcomp)==word_count(w),"one run converse witness")
                    checks['divisor_ray_witnesses']+=1
ray_targets=[]
for a,c in ((1273,191),(5303,201)):
    require(c not in constants(a),"S190 target ray exclusion")
    ray_targets.append({"slope":a,"target_constant":c,"divisors":divisors(a),"possible_constants":list(constants(a)),"target_present":False})
require(constants(1273)==(-789,-35,37,1119),"1273 constants")
require(constants(5303)==(-4909,3045),"5303 constants")

# Exact all-order membership through 500 follows from the proved cutoff.
cutoff=next(n for n in range(2,100) if low(n)>500)
require(cutoff==11 and low(cutoff)==561,"500 cutoff")
require(13 not in all_small,"13 class hostile")
r13=[18,8,3,5,14]
require(dp_count(r13)==13 and perm_count(r13)==13,"13 actual witness")
require(is_strong(r13),"13 strong")
keep=(2,3,4)
require(any(all(r13[p[i]]>>p[j]&1 for i in range(3) for j in range(i+1,3)) for p in itertools.permutations(keep)),"13 transitive two-vertex deletion")
for remove in range(5):
    keep=tuple(i for i in range(5) if i!=remove)
    require(not any(all(r13[p[i]]>>p[j]&1 for i in range(4) for j in range(i+1,4)) for p in itertools.permutations(keep)),"13 no transitive one-deletion")

# Literal tests of formulas with arbitrary-size integers, far past the census.
large_tests=[]
for n in (19,20,31,32,63,64,99,100,1001,1002):
    o=n//2;z=n-o
    minword=(1,)+(0,)*(z-1)+(1,)*(o-1)+(0,)
    maxword=(1,)*o+(0,)*z
    require(word_count(minword)==low(n) and word_count(maxword)==high(n),"large endpoint")
    lowword=(1,)+(0,)*(z)+(1,)*(o-2)+(0,)
    require(word_count(lowword)==low2(n),"large second minimum")
    if n%2==0:topword=(1,)*(o-1)+(0,)*(z+1)
    else:topword=(1,)*(o-1)+(0,1)+(0,)*(z-1)
    require(word_count(topword)==high2(n),"large second maximum")
    large_tests.append(n)
receipt={"status":"PASS","single_run_target_certificates":ray_targets,"notation":"chain length n; tournament order n+1","checks":checks,"profiles":profiles,"one_apex_spectrum_through_500":sorted(all_small),"one_apex_missing_odds_through_500":[h for h in range(1,501,2) if h not in all_small],"one_apex_witness_words_through_500":small_witnesses,"all_order_cutoff":{"target_cap":500,"first_excluded_chain_length":cutoff,"minimum_at_cutoff":low(cutoff)},"global_13_witness":{"rows":r13,"H":13,"feedback_vertex_number":2,"deleted_to_transitive":[0,1]},"large_literal_chain_lengths":large_tests}
raw=json.dumps(receipt,sort_keys=True,indent=2)+'\n'
(HERE/'receipt.json').write_text(raw)
print('PASS: exact one-apex extrema, second extrema, adjacent-swap law, sharp swap costs, and the global H=13 hostile.')
print(json.dumps(checks,sort_keys=True))
print('Complete one-apex support through 500:',len(all_small),'values; first impossible chain length:',cutoff)
print('Receipt SHA256:',hashlib.sha256(raw.encode()).hexdigest())
