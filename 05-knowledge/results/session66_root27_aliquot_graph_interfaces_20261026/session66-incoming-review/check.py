#!/usr/bin/env python3
"""Independent exact incoming-claim audit. No producer imports or floating proof.
Normal and python -O use the same explicit guards. Optional FFT is diagnostic.
"""
from pathlib import Path
from math import isqrt
from itertools import product
from functools import lru_cache
import hashlib, json, sys
HERE=Path(__file__).resolve().parent
CHECKS=[]
def require(ok,name):
    if not ok: raise RuntimeError(name)
    CHECKS.append(name)
def v2(n):
    if n == 0: raise ValueError('zero valuation not used for positive sources')
    return (n & -n).bit_length()-1
def orbit(z,n):
    bits=[]
    for _ in range(n):
        c=z%2; bits.append(str(c)); z=(3*z+c)//2
    return ''.join(bits),z
def Q(z,w):
    power=1
    for c in w:
        z=3*z+int(c)*power; power*=2
    return z
def lcp(a,b):
    h=0
    for x,y in zip(a,b):
        if x != y: break
        h+=1
    return h
def credit(w): return 1+2*w.count('1')-w.count('0')
def is_tree(w):
    q=1
    for i,c in enumerate(w):
        if q<=0: return False
        q+=2 if c=='1' else -1
    return q==0
@lru_cache(None)
def trees(j):
    # Independent balanced binary-word filter, not incoming recursive tree grammar.
    return tuple(''.join(w) for w in product('01',repeat=3*j+1)
                 if w.count('1')==j and is_tree(''.join(w)))
def frontier(z,j):
    if j==0: return '0'
    prefix=''; q=1; ones=0
    # A finite cap is an observation limit, never a nonhalting conclusion.
    for _ in range(10000):
        c=z%2; z=(3*z+c)//2; prefix+=str(c)
        q+=2 if c else -1; ones+=c
        if ones==j: return prefix+'0'*q
        if q==0: return None
    return None
pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
for pin in pins:
    b=(HERE/pin['snapshot']).read_bytes()
    require(hashlib.sha256(b).hexdigest()==pin['sha256'],'sealed SHA256 '+pin['snapshot'])
    require(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==pin['git_blob'],'sealed Git blob '+pin['snapshot'])
rank_cases=matching=0
for z in range(1,64,2):
    for n in range(1,10):
        kap,zn=orbit(z,n)
        for bits in product('01',repeat=n):
            w=''.join(bits); h=lcp(w,kap); val=v2(Q(z,w)); rank_cases+=1
            if h<n:
                if val!=h: raise RuntimeError(('mismatch',z,w))
            else:
                matching+=1
                if Q(z,w)!=2**n*zn or n-val!=-v2(zn): raise RuntimeError(('matching',z,w))
require(True,'piecewise mismatch / full-match rank law')
frontier_cases=0
for z in range(1,64,2):
    for j in range(6):
        w=frontier(z,j)
        if w is None: break
        require(is_tree(w),'legal completed frontier z=%s j=%s'%(z,j))
        n=len(w); r=n-v2(Q(z,w)); rows=trees(j)
        if not all(n-v2(Q(z,u))>r for u in rows if u!=w):
            raise RuntimeError(('fixed-size unique rank minimizer',z,j,w))
        frontier_cases+=1
require(True,'fixed-size canonical frontier uniquely minimizes rank')
require(is_tree('1000') and credit('1')==3,'q>0 completed frontier is still a complete tree')
require(Q(5,'1000')==432 and Q(21,'1000')==1728,'matching hostile numerators')
require(4-v2(Q(5,'1000'))==0 and 4-v2(Q(21,'1000'))==-2,'equal word and LCP but different matching ranks')
require(3-lcp('000',orbit(21,4)[0][1:])==0,'frontier q-LCP hostile gives zero rather than minus two')
require((1-v2(Q(7,'0')),4-v2(Q(7,'1000')))==(1,3),'unrestricted-size frontier comparison has no unique meaning')
# Sieve of primes, exact pair-support certificates through one million.
N=1000000
prime=bytearray(b'\1')*(N+1); prime[0:2]=b'\0\0'
for p in range(2,isqrt(N)+1):
    if prime[p]: prime[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
primes=[p for p in range(2,N+1) if prime[p]]
from array import array
witness=array('I',[0])*(N+1)
max_small=0
for e in range(8,N+1,2):
    for p in primes:
        if 2*p>=e: raise RuntimeError(('no distinct Goldbach witness',e))
        if prime[e-p]:
            witness[e]=p; max_small=max(max_small,p); break
require(all(prime[witness[e]] and prime[e-witness[e]] and 2*witness[e]<e for e in range(8,N+1,2)), 'all even 8..1e6 have exact distinct-prime witnesses')
# Exact support for odd sums: the only even prime is 2.
support=bytearray(N+1)
for e in range(N+1):
    if e%2: support[e]=int(e>=5 and prime[e-2])
    elif e>=8: support[e]=1
require(all(witness[e] for e in range(8,19563,2)),'exact small distinctness bridge 8..19562 for cited 9781 bound')
# Exact all-even-target census: even sources <=2N; odd sources are squares.
# If odd square root a is composite and p=spf(a), s(a^2)>=a^2/p>=a^(3/2),
# so a^3<=N^2. Prime roots contribute p+1 directly.
L=2*N
sig=array('Q',[0])*(L+1)
for d in range(1,L+1):
    for k in range(d,L+1,d): sig[k]+=d
hit=bytearray(N+1)
for n in range(2,L+1,2):
    m=sig[n]-n
    if m<=N and m%2==0: hit[m]=1
for p in primes:
    if p>=3 and p+1<=N: hit[p+1]=1
B=0
while (B+1)**3<=N*N: B+=1
composite_roots=0
for a in range(9,B+1,2):
    if prime[a]: continue
    x=a; sigma=1
    for p in primes:
        if p*p>x: break
        if x%p==0:
            exp=0
            while x%p==0: x//=p; exp+=1
            sigma *= (p**(2*exp+1)-1)//(p-1)
    if x>1: sigma*=1+x+x*x
    m=sigma-a*a; composite_roots+=1
    if m<=N: hit[m]=1
# Odd targets: 1,3,7 witnesses; every odd m>=9 by exact distinct pair.
for m in (1,3,7): hit[m]=1
for m in range(9,N+1,2):
    if not witness[m-1]: raise RuntimeError(('missing odd target witness',m))
    hit[m]=1
# All possible sources of 5 are <=(5-1)^2=16; none hits.
require(all(sum(d for d in range(1,n) if n%d==0)!=5 for n in range(2,17)),'5 is untouchable by complete quadratic source bound')
unt=[m for m in range(2,N+1) if not hit[m]]
require(len(unt)==150232,'independent all-target untouchable count through 1e6')
require(sum(m<=3000 for m in unt)==318,'incoming 3000-target subtotal')
require(unt[:40]==[2,5,52,88,96,120,124,146,162,188,206,210,216,238,246,248,262,268,276,288,290,292,304,306,322,324,326,336,342,372,406,408,426,430,448,472,474,498,516,518],'incoming first forty values')
require([m for m in unt if m%2]==[5],'only odd untouchable through 1e6')
expected_counts={10:2,100:5,1000:89,10000:1212,100000:13863,1000000:150232}
require(all(sum(m<=b for m in unt)==v for b,v in expected_counts.items()),'all six incoming power-of-ten census counts')
expected_residues={2:[150231,1],3:[65349,41873,43010],4:[83603,1,66628,0],5:[27675,35525,31933,33931,21168],6:[65349,0,43009,0,41873,1],7:[22234,20439,18186,26415,19465,22166,21327],12:[35017,0,17708,0,23285,1,30332,0,25301,0,18588,0],14:[22234,0,18186,0,19465,1,21327,0,20439,0,26415,0,22165,0],15:[12096,6854,7799,11090,4782,7569,17355,10682,11296,11356,8010,11316,13452,11545,5030]}
residues={}
for mod,expected in expected_residues.items():
    actual=[0]*mod
    for m in unt: actual[m%mod]+=1
    require(actual==expected,'independent incoming residue count mod '+str(mod)); residues[mod]=actual
# Small independent complete enumeration checks the even-prefix proof implementation.
M=200
small=[0]*((M-1)**2+1)
for d in range(1,len(small)):
    for k in range(2*d,len(small),d): small[k]+=d
brute=set(small[2:])
require(all(bool(hit[m])==(m in brute) for m in range(2,M+1)),'small independent complete divisor census through 200')
# Semiprime-only examples disprove the inference about coverage of EXCEPTIONS.
# 9=s(3*5), 11=s(3*7), 17=s(3*13): each has a distinct pair, not an exception.
require(all(p!=q and prime[p] and prime[q] and p+q==m-1 for m,p,q in [(9,3,5),(11,3,7),(17,3,13)]),'claimed dependency examples are nonexceptions')
report={
 'status':'PASS','checks':len(CHECKS),'source_commit':'44abcaf22b78125a11e039054f0d0ba9e4978f08',
 'rank_universe':{'odd_sources':[1,63],'binary_lengths':[1,9],'cases':rank_cases,'matching_cases':matching,'fixed_size_frontiers':frontier_cases,'internal_nodes':[0,5]},
 'hostiles':{'completed_frontier':{'prefix':'1','q':3,'word':'1000','is_tree':True},'matching_ranks':{'z5':0,'z21':-2},'z21_frontier_q_minus_LCP':0},
 'goldbach_support':{'limit':N,'even_distinct_witnesses':len(range(8,N+1,2)),'largest_first_prime':max_small,'exact_support_sha256':hashlib.sha256(support).hexdigest(),'small_cited_bridge':[8,19562]},
 'census':{'limit':N,'untouchables':len(unt),'even_untouchables':len(unt)-1,'odd_untouchables':[5],'even_source_limit':L,'composite_square_root_limit':B,'composite_square_roots_tested':composite_roots,'touchability_sha256':hashlib.sha256(hit).hexdigest(),'power_of_ten':expected_counts,'residues':residues},
 'not_replayed':['incoming 443731 non-semiprime-source target count','incoming 56267 semiprime/prime-square-only target count','incoming 2509482 composite-pair enumeration count','incoming long root27 orbit and halting-source censuses'],
 'checks_passed':CHECKS}
print(json.dumps(report,indent=2,sort_keys=True))
