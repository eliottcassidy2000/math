#!/usr/bin/env python3
"""Independent literal audit; standard library, always-active checks."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
import hashlib, json

checks=Counter()
def need(ok, name):
    checks[name]+=1
    if not ok: raise RuntimeError(name)

def parts(n, lo=1):
    if n==0: yield ()
    for d in range(lo,n+1,2):
        for tail in parts(n-d,d): yield (d,)+tail

def exponent(p):
    return sum((d-1)//2 for d in p)+sum(gcd(a,b) for a,b in combinations(p,2))
def fixed(p): return 1<<exponent(p)
def z(p):
    out=1
    for d,m in Counter(p).items(): out*=d**m*factorial(m)
    return out

def lah(m,k):
    if m==k==0:return 1
    if not 1<=k<=m:return 0
    return factorial(m)//factorial(k)*comb(m-1,k-1)
def core_fixed(p):
    bank=list(Counter(p).items())
    ans=0
    for ks in product(*(range(1,m+1) for d,m in bank)):
        src=tuple(d for (d,m),k in zip(bank,ks) for _ in range(k))
        term=fixed(src)
        for (d,m),k in zip(bank,ks):term*=(-d)**(m-k)*lah(m,k)
        ans+=term
    return ans

def adjacency(n, code):
    rows=[0]*n
    for k,(a,b) in enumerate(combinations(range(n),2)):
        u,v=(a,b) if code>>k&1 else (b,a)
        rows[u]|=1<<v
    return rows

def hpairs(rows):
    n=len(rows); full=(1<<n)-1
    return [(a,b) for a,b in combinations(range(n),2)
            if ((rows[a]^rows[b]) & (full^(1<<a)^(1<<b)))==0]
def free(rows):
    n=len(rows);full=(1<<n)-1
    return all((rows[a]^rows[b])&(full^(1<<a)^(1<<b))
               for a,b in combinations(range(n),2))

def extension_audit():
    receipt=[]
    for n in range(1,7):
        total=0; free_count=0; min_success=1<<n; bbank=Counter()
        for code in range(1<<comb(n,2)):
            old=adjacency(n,code); pairs=hpairs(old); b=n-len(pairs)
            actual=Counter(); success=0
            for bits in range(1<<n):
                rows=[old[v] | ((1-((bits>>v)&1))<<n) for v in range(n)]+[bits]
                bp=n+1-len(hpairs(rows)); actual[bp]+=1
                if bp==n+1:success+=1
                need(bp>=b, 'core size nondecreasing')
                total+=1
            expected=Counter({b:n+b})
            for s in range(n-b+1):
                count=(1<<b)*comb(n-b,s)-(2*b if s==0 else n-b if s==1 else 0)
                if count:expected[b+1+s]+=count
            need(actual==expected,'literal extension kernel')
            alt=(1<<b)-(2*b if n==b else 1 if n==b+1 else 0)
            need(success==alt,'alternating pairfree success formula')
            if n>=2:need(success>0,'every n>=2 parent has pairfree one-vertex extension')
            min_success=min(min_success,success)
            free_count+=not pairs;bbank[b]+=1
        receipt.append(dict(n=n,parents=1<<comb(n,2),extensions=total,labelled_core_count=free_count,
                            minimum_successful_attachments=min_success,core_size_counts=dict(sorted(bbank.items()))))
    return receipt

def invariant_orientations(p):
    n=sum(p); perm=[];offset=0
    for d in p:
        perm.extend([offset+(j+1)%d for j in range(d)]);offset+=d
    unseen=set(combinations(range(n),2));orbits=[]
    while unseen:
        a,b=min(unseen);start=(a,b);cycle=[]
        while True:
            need((min(a,b),max(a,b)) in unseen,'unordered edge orbit disjointness')
            unseen.remove((min(a,b),max(a,b)));cycle.append((a,b));a,b=perm[a],perm[b]
            if (a,b)==start:break
            need((b,a)!=start,'odd type never reverses an edge')
        orbits.append(cycle)
    need(len(orbits)==exponent(p),'literal edge orbit exponent')
    count=0
    for mask in range(1<<len(orbits)):
        rows=[0]*n
        for j,orbit in enumerate(orbits):
            for a,b in orbit:
                if mask>>j&1:rows[a]|=1<<b
                else:rows[b]|=1<<a
        count+=free(rows)
    return count

def class_audit():
    out=[]
    for n in range(7):
        pairs=list(combinations(range(n),2));pairidx={p:i for i,p in enumerate(pairs)}
        perms=list(permutations(range(n)));seen=set();total=cores=0
        for code in range(1<<len(pairs)):
            if code in seen:continue
            total+=1;rows=adjacency(n,code);cores+=free(rows)
            for perm in perms:
                dest=0
                for k,(a,b) in enumerate(pairs):
                    u,v=(perm[a],perm[b]) if code>>k&1 else (perm[b],perm[a])
                    if u<v:dest|=1<<pairidx[u,v]
                seen.add(dest)
        need(len(seen)==1<<comb(n,2),'class orbits cover labelled universe')
        out.append(dict(n=n,tournaments=total,cores=cores))
    return out

def main():
    series=[];coeffs={};coeffrows=0
    for n in range(31):
        a=c=Fraction(0)
        for p in parts(n):
            f=core_fixed(p);need(isinstance(f,int) and f>=0,'fixed core nonnegative integer')
            need(f<=fixed(p),'core fixed subset bound')
            coeffs[p]=f;coeffrows+=1
            a+=Fraction(fixed(p),z(p));c+=Fraction(f,z(p))
        need(a.denominator==c.denominator==1,'ordinary tournament/core count integrality')
        series.append(dict(n=n,tournaments=int(a),cores=int(c)))
    for p,v in coeffs.items():
        bank=list(Counter(p).items()); recovered=0
        for ks in product(*(range(1,m+1) for d,m in bank)):
            src=tuple(d for (d,m),k in zip(bank,ks) for _ in range(k));term=coeffs[src]
            for (d,m),k in zip(bank,ks):term*=d**(m-k)*lah(m,k)
            recovered+=term
        need(recovered==fixed(p),'positive inverse coefficient roundtrip')
    literal=[]
    for n in range(9):
        for p in parts(n):
            if n>6 and all(d==1 for d in p):continue
            value=invariant_orientations(p)
            need(value==coeffs[p],'literal invariant core count equals signed transform')
            literal.append(dict(partition=p,tournaments=fixed(p),cores=value))
    ext=extension_audit(); classes=class_audit()
    for row in classes:
        n=row['n'];need(row['cores']==series[n]['cores'] and row['tournaments']==series[n]['tournaments'],
                       'independent class census equals cycle-index specialization')
    for n in range(2,30):
        need(series[n]['tournaments']<=(n+1)*series[n+1]['cores'],'one-point universality counting bound')
    for n in range(1,31):
        need(series[n]['tournaments']-series[n]['cores']<=(n-1)*series[n-1]['tournaments'],
             'one-pair contraction counting bound')
    out=dict(status='PASS',checks=dict(sorted(checks.items())),total_checks=sum(checks.values()),
             fixed_cycle_coefficients_checked=coeffrows,series=series,literal_fixed_cycle_census=literal,
             literal_extension_census=ext,independent_class_census=classes)
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=='__main__':main()
