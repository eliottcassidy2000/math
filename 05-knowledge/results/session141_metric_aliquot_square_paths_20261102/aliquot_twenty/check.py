#!/usr/bin/env python3
"""Exact controls; standard-library only; explicit exceptions survive -O."""
import hashlib
import json
from math import isqrt, prod
import struct

N = 200000
LISTS = {
    'first': [16,10,15,14,21,35,22,33],
    'second': [12,18,20,28,24],
    'negative': [4,9,2,8,49,121,27,169,16],
}
class Gates:
    def __init__(self): self.count=0
    def eq(self, a, b, label):
        self.count+=1
        if a!=b: raise RuntimeError(f'{label}: {a!r} != {b!r}')
    def true(self, a, label): self.eq(bool(a),True,label)
g=Gates()
# Path A: add each divisor once, including the source itself.
sigma=[0]*(N+1); tau=[0]*(N+1)
for d in range(1,N+1):
    for n in range(d,N+1,d): sigma[n]+=d; tau[n]+=1
# Path B: independent least-prime-factor sieve, then prime-factor products.
spf=list(range(N+1))
for p in range(2,isqrt(N)+1):
    if spf[p]==p:
        for n in range(p*p,N+1,p):
            if spf[n]==n: spf[n]=p
def factors(n):
    ans=[]
    while n>1:
        p=spf[n]; e=0
        while n%p==0: n//=p; e+=1
        ans.append((p,e))
    return ans
def multiplicative(n):
    fac=factors(n)
    return (prod(sum(p**i for i in range(e+1)) for p,e in fac),
            prod(e+1 for p,e in fac))
ha=hashlib.sha256(); hb=hashlib.sha256()
for n in range(1,N+1):
    sig,t=multiplicative(n)
    g.eq(sig,sigma[n],f'sigma[{n}]')
    g.eq(t,tau[n],f'tau[{n}]')
    ha.update(struct.pack('<QQ',sigma[n],tau[n]))
    hb.update(struct.pack('<QQ',sig,t))
g.eq(ha.hexdigest(),hb.hexdigest(),'independent array digests')
s=lambda n:sigma[n]-n

def proper(n):
    ds=set()
    for d in range(1,isqrt(n)+1):
        if n%d==0: ds.add(d); ds.add(n//d)
    ds.discard(n)
    return sorted(ds)
def pairs(n): return [(d,n//d) for d in range(1,isqrt(n)+1) if n%d==0]
for n in range(1,5001):
    ds=proper(n)
    g.eq(sum(ds),s(n),f'proper-sum[{n}]')
    g.eq(len(ds),tau[n]-1,f'proper-count[{n}]')
    witnesses=[(a,b) for a,b in pairs(n) if 1+a+b==s(n)]
    g.eq(bool(witnesses),tau[n]==4 or s(n)==n+2,f'compression[{n}]')

def row(n):
    return {'n':n,'factorization':factors(n),'proper_divisors':proper(n),
            's':s(n),'sigma':sigma[n],'tau':tau[n],'excess':s(n)-n,
            'prime_power':len(factors(n))==1,
            'compression_witnesses':[(a,b) for a,b in pairs(n) if 1+a+b==s(n)]}
rows={str(n):row(n) for n in sorted(set(sum(LISTS.values(),[])+[7,20,24,26,38,41]))}
g.eq(sorted(set(LISTS['first'])&set(LISTS['negative'])),[16],'literal list collision')
g.true(all(len(factors(n))==1 for n in LISTS['negative']),'negative list prime powers')
g.true(all(factors(n)[0][1]==1 and factors(n)[1][1]==1
           for n in LISTS['first'] if n!=16),'first list except16 distinct semiprimes')
expected_fibers={7:[8],15:[16,33],16:[12,26],20:[34,361],21:[18,51,91],
                 22:[20,38],24:[529],26:[46],28:[28],36:[24]}
fibers={}
for t,expected in expected_fibers.items():
    bound=(t-1)**2
    actual=[n for n in range(1,bound+1) if s(n)==t]
    g.eq(actual,expected,f'complete inverse s fiber {t}')
    fibers[str(t)]={'bound':bound,'sources':actual}
sigmafibers={}
for t,expected in {24:[14,15,23],42:[20,26,41],60:[24,38,59]}.items():
    actual=[n for n in range(1,t) if sigma[n]==t]
    g.eq(actual,expected,f'complete inverse sigma fiber {t}')
    sigmafibers[str(t)]={'exclusive_bound':t,'sources':actual}
chains={20:[20,22,14,10,8,7,1,0],26:[26,16,15,9,4,3,1,0],
        24:[24,36,55,17,1,0]}
for start,chain in chains.items():
    for a,b in zip(chain,chain[1:]):g.eq(s(a),b,f'chain {start}, edge {a}')
for n in [20,104,464,650,1952]:
    g.eq(s(n)-n,2,f'positive excess control {n}')
    g.true(tau[n]<=17,f'positive arity control {n}')
expected_low=[20,104,464,650,1952]
g.eq([n for n in range(1,N+1) if s(n)-n==2 and tau[n]<=17],expected_low,'finite low-arity census')
all_excess2=[n for n in range(1,N+1) if s(n)-n==2]
g.eq(all_excess2,expected_low+[130304],'full bounded excess2 census')
g.eq((s(130304)-130304,tau[130304]),(2,18),'sharp next-layer boundary')
g.eq((s(16),tau[16]),(15,5),'16 not compression')
g.true(tau[8]==4 and tau[27]==4,'negative-list C controls')
g.eq(len(factors(650)),3,'650 outside two-prime theorem')
g.eq(s(18),21,'composite q9 hostile actual')
g.eq((2**1-1)*9+2**2-1,12,'composite q9 hostile formal')
# All exponent patterns with product(e+1)<=17, exponents nonincreasing.
patterns={}
def extend(exponents,last,total):
    if exponents: patterns.setdefault(total,[]).append(exponents)
    for e in range(1,last+1):
        if total*(e+1)<=17: extend(exponents+[e],e,total*(e+1))
extend([],16,1)
for t in patterns: patterns[t].sort()
expected_new={13:[[12]],14:[[6,1],[13]],15:[[4,2],[14]],
              16:[[1,1,1,1],[3,1,1],[3,3],[7,1],[15]],17:[[16]]}
for t,exps in expected_new.items(): g.eq(patterns[t],exps,f'complete exponent shapes tau{t}')
g.eq(pairs(533),[(1,533),(13,41)],'positive factor pairs533')
g.true(2304<2*1155,'odd squarefree four-prime strict deficiency')
g.true(1248<2*625,'odd p>=5 cubic chart strict deficiency')
prime=lambda n:n>=2 and n<=N and spf[n]==n
params=[]
for a in range(2,8):
    q=2**(a+1)-3
    params.append({'a':a,'q':q,'prime':prime(q),'source':2**a*q})
g.eq([x['source'] for x in params if x['prime']],[20,104,464,1952],'two-prime arity parameters')
for q in range(3,1000,2):
    if prime(q):
        g.eq(s(4*q),3*q+7,f'affine chart q{q}')
        r=3*q+4
        if prime(r):g.eq(s(4*q),s(2*r),f'aliquot cross-chart q{q}')
        if (7*q+4)%3==0:
            r=(7*q+4)//3
            if prime(r):g.eq(sigma[4*q],sigma[2*r],f'sigma cross-chart q{q}')
# Independent exact F8 polynomial arithmetic, retained marking a=x.
def mul(a,b):
    c=0
    while b:
        if b&1:c^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return c
def rank(vs):
    basis={}
    for v in vs:
        while v:
            k=v.bit_length()-1
            if k in basis:v^=basis[k]
            else:basis[k]=v;break
    return len(basis)
a=2;a2=mul(a,a);a4=mul(a2,a2)
g.eq([a,a2,a4],[2,4,6],'marked F8 orbit')
g.eq(mul(a4,a4),a,'Frobenius3 identity')
g.eq(rank([a,a2,a4]),2,'Fano-line not normal basis')
b=mul(a,a2);b2=mul(b,b);b4=mul(b2,b2)
g.eq(rank([b,b2,b4]),3,'alternate generator rank3 control')
collatz=lambda n:n//2 if n%2==0 else 3*n+1
for n in [1,2,4]:g.eq(collatz(n)%7,4*n%7,f'cycle mod7 {n}')
g.true(2%7==9%7 and collatz(2)%7!=collatz(9)%7,'no global descent modulo7')
result={'status':'PASS','scope':{'source_universe':[1,N],'direct_divisor_universe':[1,5000],
 'all_size_proof':'PROOFS.md section5; computation is a finite control only'},
 'gates':g.count,'sigma_tau_digest':ha.hexdigest(),'literal_lists':LISTS,'source_rows':rows,
 'inverse_s_fibers':fibers,'inverse_sigma_fibers':sigmafibers,'aliquot_chains':chains,
 'excess2_census':all_excess2,'low_arity_excess2_census':expected_low,
 'tau13_to17_patterns':{str(t):patterns[t] for t in range(13,18)},
 'two_prime_parameters':params,'533_factor_pairs':pairs(533),
 'f8':{'marked_orbit':[a,a2,a4],'marked_rank':rank([a,a2,a4]),
 'alternate_orbit':[b,b2,b4],'alternate_rank':rank([b,b2,b4])}}
print(json.dumps(result,indent=2,sort_keys=True))
