#!/usr/bin/env python3
import json, hashlib, sys
from pathlib import Path
from math import isqrt
from aliquot_fibers import *

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def digest(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

# Independent complete brute-force divisor sieve, up to the proved common bound.
M=500; B=(M-1)**2
sig=sigma_sieve(B)
truth={m:[] for m in range(2,M+1)}
for n in range(4,B+1):
    m=sig[n]-n
    if m in truth: truth[m].append(n)
computed={m:fiber(m) for m in range(2,M+1)}
require(computed==truth,'complete individual-fiber/direct-sieve disagreement')
even=even_prefix(M)
require(even=={m:truth[m] for m in range(2,M+1,2)},'parity-prefix/direct-sieve disagreement')

# Bounds and parity over a universe strictly larger than the fiber validation.
N=500_000
sig2=sigma_sieve(N)
checks=0; guard_hostiles=[]
for n in range(4,N+1):
    m=sig2[n]-n
    if m==1: continue
    require(n<=(m-1)**2,'global inverse bound')
    r=isqrt(n)
    prime_square=(r*r==n and factors(r)=={r:1})
    if not prime_square: require(4*n<=(m-1)**2,'paired-divisor bound')
    if n%2==0 and n>=6: require(n<=2*m-6,'even-source bound')
    if n%2==0 and m%2:
        require(r*r==n or (n%2==0 and isqrt(n//2)**2==n//2),'odd target parity')
    checks+=1

# Killing missing-coordinate mutations explicitly.
require(aliquot(4)==3 and aliquot(2)*2+3==5,'coprimality hostile')
require(fiber(5)==[],'5 must stay untouchable')
require(32 in fiber(31),'squarefull u,v may share primes: 32=2^2*2^3')
require(121 in fiber(12) and 121>2*12,'linear cutoff misses odd prime square')
require(27 in fiber(13),'discarding non-square squarefull sources loses 27')
require(aliquot(2)==aliquot(3)==1 and aliquot(4)!=aliquot(6),'fiber equivalence not multiplicative congruence')
for p in [2,3,5,7,11,13,17,19]:
    require(max(fiber(p+1))==(p+1-1)**2,'sharp square-prime upper bound')

# Hyperbola family independently against factorization for a bounded named universe.
hyperbola_checks=0
for a in range(2,31):
    for m in range(2,201):
        got=two_prime_core(m,a)
        expect=[n for n in truth[m] if n%a==0 and len(factors(n//a))==2
                and all(e==1 for e in factors(n//a).values())
                and gcd(n//a,a)==1]
        require(got==expect,f'hyperbola mismatch a={a},m={m}')
        hyperbola_checks+=1

# Exact even prefix, plus explicit (distinct-prime) witnesses for every odd <=X except5.
X=100_000
prefix=even_prefix(X); spf=spf_sieve(X)
missing_odd=[]; odd_witnesses={}
for m in range(1,X+1,2):
    n=goldbach_witness(m,spf)
    if n is None: missing_odd.append(m)
    else:
        require(aliquot(n)==m,'bad Goldbach witness')
        odd_witnesses[m]=n
require(missing_odd==[5],'unexpected uncertified odd target')
U=[m for m in range(2,X+1,2) if not prefix[m]]
allU=sorted(U+[5])
require(len(allU)==13863,'published 100000 count mismatch')
# Comparison between true missing output count and naive source<=X omissions.
source_cutoff=set(sig2[n]-n for n in range(1,X+1))
false_positives=[m for m in range(2,X+1,2) if m not in source_cutoff and prefix[m]]
summary={
 'status':'FINITE-EXACT; no asymptotic density or odd-uniqueness claim',
 'fiber_validation':{'targets':[2,M],'direct_source_bound':B,'fiber_sha256':digest(computed)},
 'bound_validation':{'sources':[4,N],'composite_sources_tested':checks},
 'hyperbola_validation':{'a':[2,30],'m':[2,200],'cases':hyperbola_checks},
 'example_fibers':{str(m):fiber(m) for m in range(14) if m!=1},
 'infinite_fiber':{'1':'all primes'},
 'complete_census':{'X':X,'untouchables':len(allU),'even_untouchables':len(U),
   'odd_untouchables':[5],'fraction':len(allU)/X,'untouchables_sha256':digest(allU),
   'odd_positive_witnesses_sha256':digest(odd_witnesses),
   'naive_source_cutoff_false_positive_even_count':len(false_positives),
   'naive_source_cutoff_false_positive_even_first20':false_positives[:20]},
 'hostiles':{'omit_coprimality':'core2, p2 falsely predicts s(4)=5, actual3',
   'force_gcd_u_v_1':'misses source32=2^2*2^3 in F31',
   'linear_all_sources':'misses source121 in F12',
   'squarefull_means_square':'misses source27 in F13',
   'multiplicative_fiber_quotient':'s2=s3=1 but s4=3 != s6=6'},
 'sample_certificates':{str(m):{'fiber':fiber(m,True)[0],'enumeration':fiber(m,True)[1]}
                        for m in [2,5,52,88,96,120,124,146,31,13]}}
root=Path(__file__).resolve().parents[1]
(root/'results'/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
(root/'results'/'small_fibers.json').write_text(json.dumps(computed,indent=2)+'\n')
(root/'results'/'untouchables_100000.json').write_text(json.dumps(allU)+'\n')
print(json.dumps(summary,indent=2))
