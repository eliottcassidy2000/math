#!/usr/bin/env python3
"""Exact authored controls; exceptions remain active under python -O."""
from collections import Counter
from itertools import combinations
from math import gcd
import hashlib, json, pathlib

checks=Counter()
def need(ok, category, witness=None):
    checks[category]+=1
    if not ok: raise RuntimeError((category,witness))
def valuation(n):
    if not n: raise ValueError('zero has no finite valuation')
    n=abs(n); k=0
    while n%2==0: n//=2; k+=1
    return k
def oddmap(n,b):
    if b not in (-1,1) or n%2!=1: raise ValueError('bad signed port')
    v=valuation(3*n+b)
    return (3*n+b)//(1<<v),v
def crt_odd_dyadic(r,m,a,b):
    """Least positive source with n=r modm and exact prescribed valuation a."""
    d=1<<(a+1)
    u=(pow(3,-1,d)*((1<<a)-b))%d
    n=r+m*((u-r)*pow(m,-1,d)%d)
    return n or m*d
def orbits(perm):
    unused=set(perm); out=[]
    need(set(perm.values())==unused,'permutation')
    while unused:
        x=min(unused); start=x; cycle=[]
        while x in unused:
            unused.remove(x); cycle.append(x); x=perm[x]
        need(x==start,'cycle_closure'); out.append(cycle)
    return out
def table(mod,odd_only=False):
    states=[(n,b) for n in range(mod) if not odd_only or n%2 for b in (-1,1)]
    return {(n,b):((2*n+b)%mod,-b) for n,b in states}

# Complete, finite signed inverse-fibre and native carry universes.
for q in range(-101,102,2):
    if q%3==0: continue
    last=None
    for a in range(1,21):
        b=1 if (pow(2,a,3)*q)%3==1 else -1
        n=((1<<a)*q-b)//3
        need(oddmap(n,b)==(q,a),'inverse_fibre',(q,a))
        if last:
            pn,pb=last
            need((n,b)==(2*pn+pb,-pb),'sheet_alternation')
        last=n,b

for m in range(3,102,2):
    for b in (-1,1):
        # One counterexample at every declared modulus, not only m=5.
        ns=[crt_odd_dyadic(0,m,a,b) for a in (1,2)]
        ys=[oddmap(n,b)[0]%m for n in ns]
        need(ys[0]!=ys[1],'odd_modulus_nondescending',(m,b,ns,ys))
        for r in range(m):
            for a in range(1,9):
                n=crt_odd_dyadic(r,m,a,b)
                y,v=oddmap(n,b)
                s=((3*r+b)*pow(1<<a,-1,m))%m
                c=(3*r+b-(1<<a)*s)//m
                t=(n-r)//m
                need(v==a and n%m==r,'crt_exact_guard')
                need(3*r+b==(1<<a)*s+m*c,'signed_carry')
                need(y==s+m*((3*t+c)//(1<<a)),'height_lift')
                need((3*t+c-(1-s)*(1<<a))%(1<<(a+1))==0,'guard_carry_iff')
                # Altering the top guard bit preserves divisibility but breaks exactness.
                bad=n+m*(1<<a)
                need(valuation(3*bad+b)>a,'hostile_missing_guard_bit')

tables={}
for m,odd in [(14,True),(15,False),(210,True)]:
    perm=table(m,odd); cs=orbits(perm)
    hist=Counter(map(len,cs))
    tables[str(m)]={'states':len(perm),'cycle_histogram':dict(sorted(hist.items()))}
expected={'14':{2:1,6:2},'15':{6:1,12:2},'210':{6:7,12:14}}
for m in tables:
    need(tables[m]['cycle_histogram']==expected[m],'factor_cycle_census')
perm=table(210,True)
zmap={(n,b):(3*n+b)%315 for n,b in perm}
need(set(zmap.values())=={z for z in range(315) if z%3},'three_lift_bijection')
for p in perm:
    need(zmap[perm[p]]==2*zmap[p]%315,'doubling_intertwiner')
for cyc in orbits(perm):
    need(len(cyc)==(6 if zmap[cyc[0]]%5==0 else 12),'five_period_gate')
# CRT gluing must keep the common b; explicit inverse reconstruction.
ports={(n%14,n%15,b) for n,b in perm}
need(len(ports)==210,'crt_common_sheet')
for r14,r15,b in ports:
    n=(r14+14*((r15-r14)*pow(14,-1,15)%15))%210
    need((n,b) in perm,'crt_reconstruction')
need(14*30==2*len(ports),'hostile_independent_sheets')

# Exact anti-action, with no tournament interpretation forced at other moduli.
def ternary_sign(x):
    if not x: raise ValueError('equal vertices')
    while x%3==0: x//=3
    return 1 if x%3==1 else -1
ternary=[]
for k in range(0,5):
    m=3**k; p=table(m); cs=orbits(p)
    need(len(cs)==1 and len(cs[0])==2*m,'ternary_full_clock')
    modulus=3*m; zs=[z for z in range(modulus) if z%3]
    for x,y in combinations(zs,2):
        need(ternary_sign((2*y)%modulus-(2*x)%modulus)==-ternary_sign(y-x),'ternary_anti_action')
    ternary.append({'k':k,'vertices':2*m,'minimum_unrestricted_defect':1,
                    'native_sierpinski_defect':m*(2*m-1)})

# Independent direct orientation counts versus the inherited signed-pair formula.
def defect(n,sigma,edge):
    return sum(edge(i,j)!=edge(sigma[i],sigma[j]) for i,j in combinations(range(n),2))
def antipodal_blocks(lengths):
    sig=[]; blocks=[]; off=0
    for d in lengths:
        if d%2: raise ValueError('even cycles required')
        sig.extend(off+(i+1)%d for i in range(d))
        blocks.extend((off+i,off+i+d//2) for i in range(d//2))
        off+=d
    return sig,blocks
module_controls=[]
for lengths in [(2,),(4,),(6,),(8,),(10,),(4,2),(4,4),(6,2)]:
    sig,blocks=antipodal_blocks(lengths); n=len(sig); r=len(blocks)
    loc={v:(i,j) for i,p in enumerate(blocks) for j,v in enumerate(p)}
    pairs=list(combinations(range(r),2)); ix={p:i for i,p in enumerate(pairs)}
    pi=[loc[sig[p[0]]][0] for p in blocks]
    hist=Counter()
    for mask in range(1<<(len(pairs)+r)):
        def qe(i,j):
            return bool((mask>>ix[min(i,j),max(i,j)])&1) ^ (i>j)
        def edge(u,v):
            i,a=loc[u]; j,b=loc[v]
            if i==j: return bool((mask>>(len(pairs)+i))&1) ^ (a>b)
            return qe(i,j)
        outer=defect(r,pi,qe)
        inner=sum(edge(u,v)!=edge(sig[u],sig[v]) for u,v in blocks)
        total=defect(n,sig,edge)
        need(total==4*outer+inner,'clone_weighted_split')
        hist[total]+=1
    formula=len(lengths)+4*sum(d%4==0 for d in lengths)
    need(min(hist)==formula,'clone_minimum_formula',(lengths,hist))
    module_controls.append({'cycles':lengths,'orientations':sum(hist.values()),'minimum':min(hist)})

hostiles={
 'same_endpoint_wrong_sheet':{'plus_source':3,'minus_source':7,'endpoint':5,
     'plus_continuation':oddmap(5,1)[0],'minus_continuation':oddmap(5,-1)[0]},
 'crt_does_not_fix_valuation':[
     {'n':n,'mod14':n%14,'mod15':n%15,'sheet':1,'endpoint':oddmap(n,1)[0],'valuation':oddmap(n,1)[1]} for n in (1,211)],
 'clone_4_cycle':{'unrestricted_minimum':1,'antipodal_clone_minimum':5},
 'clone_210':{'unrestricted_minimum':21,'antipodal_clone_minimum':77,
     'quotient_even_cycles':14,'weighted_quotient_tax':56},
}
need(hostiles['same_endpoint_wrong_sheet']['plus_continuation']==1 and
     hostiles['same_endpoint_wrong_sheet']['minus_continuation']==7,'hostile_sheet_seam')
need(oddmap(1,1)==(1,2) and oddmap(211,1)==(317,1),'hostile_crt_guard')
base=pathlib.Path(__file__).resolve().parent
manifest=json.loads((base/'SOURCE_MANIFEST.json').read_text())
for row in manifest['sources']:
    data=(base/row['local']).read_bytes()
    need(hashlib.sha256(data).hexdigest()==row['sha256'] and len(data)==row['bytes'],'source_integrity')
print(json.dumps({'status':'PASS','checks':dict(sorted(checks.items())),'total_checks':sum(checks.values()),
 'tables':tables,'ternary':ternary,'clone_controls':module_controls,'hostiles':hostiles},indent=2,sort_keys=True))
