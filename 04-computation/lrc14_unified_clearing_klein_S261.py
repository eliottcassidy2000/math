# -*- coding: utf-8 -*-
# klein-2026-07-11-S261: THE UNIFIED CLEARING FORMULA on the coprime sub-family (THM-718 + opus
# auto-safe + kps two-conditions), and the SHRINK (13-runner anti-concentration -> ~4 runners).
# clearing_count_units(v,q) = phi(q) - |{+-j v_i mod q : coprime-to-q i, a unit}|  on the m=1
# window [15,27] (prime+composite) + primes {29,31}. (q=30 excluded: composite-m>=2, auto-safe fails.)
from math import gcd
import random, statistics
def danger(q):
    lo=-(-q//14); hi=(13*q)//14
    return set(range(0,lo))|set(range(hi+1,q))
def m_of(q): return -(-q//14)-1
def phi(q): return sum(1 for a in range(1,q) if gcd(a,q)==1)
def clearing_count_units(v,q):
    dg=danger(q); return sum(1 for p in range(1,q) if gcd(p,q)==1 and all((x*p)%q not in dg for x in v))
def cover_num_coprime(v,q):
    m=m_of(q); S=set()
    for x in v:
        if gcd(x,q)!=1: continue
        for j in range(1,m+1):
            r=(j*x)%q
            if gcd(r,q)==1: S.add(r); S.add((q-r)%q)
    return len(S)
def is_DC(v): return all(any(x%d==0 for x in v) for d in range(2,15))
VALID=[q for q in range(15,32) if q%14 and not (m_of(q)>=2 and any(q%d==0 for d in range(2,q)))]
if __name__=="__main__":
    random.seed(5); bad=0
    for q in VALID:
        for _ in range(5000):
            v=random.sample(range(1,300),13)
            if any(x%q==0 for x in v): continue
            if clearing_count_units(v,q)!=phi(q)-cover_num_coprime(v,q): bad+=1
    print(f"unified formula on {VALID}: {bad} failures")
    fams=[]
    while len(fams)<3000:
        v=sorted(random.sample(range(1,80),13)); g=0
        for x in v: g=gcd(g,x)
        if g==1 and is_DC(v): fams.append(v)
    mins=[min((sum(1 for x in v if gcd(x,q)==1) for q in VALID if not any(x%q==0 for x in v)),default=13) for v in fams]
    print(f"DC min-over-window #coprime: median {statistics.median(mins)}, mean {statistics.mean(mins):.1f}")
