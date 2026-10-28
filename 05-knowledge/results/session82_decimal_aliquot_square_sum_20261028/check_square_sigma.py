#!/usr/bin/env python3
"""Exact finite certificates, not an algorithmic proof of Ljunggren's theorem."""
from math import isqrt, gcd, prod
import json
import argparse
import importlib.util
from pathlib import Path
from hashlib import sha256

def require(test, label):
    if not test: raise RuntimeError(label)
def square(n): return n >= 0 and isqrt(n)**2 == n
def factors(n):
    out={}; d=2
    while d*d<=n:
        while n%d==0: out[d]=out.get(d,0)+1; n//=d
        d+=1
    if n>1: out[n]=out.get(n,0)+1
    return out
def sigma(n): return prod((p**(e+1)-1)//(p-1) for p,e in factors(n).items())
def primes(N):
    s=bytearray(b'\1')*(N+1); s[:2]=b'\0\0'
    for p in range(2,isqrt(N)+1):
        if s[p]: s[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
    return [p for p in range(2,N+1) if s[p]]

def dependency(lower):
    """First GF(2) relation among odd primes strictly above lower, increasing."""
    N=max(100,2*lower)
    while True:
        ps=[p for p in primes(N) if p>max(2,lower)]
        coordinates={}; pivots={}
        for j,p in enumerate(ps):
            v=0
            for q,e in factors(p+1).items():
                if e%2:
                    if q not in coordinates: coordinates[q]=len(coordinates)
                    v ^= 1 << coordinates[q]
            witness=1<<j
            while v:
                pivot=v.bit_length()-1
                if pivot not in pivots:
                    pivots[pivot]=(v,witness); break
                row,mask=pivots[pivot]; v^=row; witness^=mask
            if v==0:
                chosen=[q for i,q in enumerate(ps) if witness>>i&1]
                n=prod(chosen); sig=prod(q+1 for q in chosen)
                require(chosen and min(chosen)>lower and square(sig),'dependency')
                return {'primes':chosen,'n':n,'sigma':sig,'root':isqrt(sig),'aliquot':sig-n}
        N*=2

N=1_000_000
solutions=[x for x in range(-100,N+1) if square((x+1)*(x*x+1))]
require(solutions==[-1,0,1,7],'polynomial finite check')
ps=primes(N); pp=[]
for p in ps:
    n=p; a=1; sig=p+1
    while n<=N:
        if square(sig): pp.append({'n':n,'p':p,'a':a,'sigma':sig,'root':isqrt(sig),'aliquot':sig-n})
        a+=1; n*=p; sig+=n
pp.sort(key=lambda d:d['n'])
require([x['n'] for x in pp]==[3,81,343],'prime-power finite check')
# Direct divisor enumeration cross-checks separate from multiplicative formula.
for n in range(1,5001):
    direct=sum(d for d in range(1,isqrt(n)+1) if n%d==0 for d in ([d] if d*d==n else [d,n//d]))
    require(sigma(n)==direct,'sigma cross-check')
# All modular reductions survive independently on a substantial finite range.
for x in range(1,10001):
    require(gcd(x+1,x*x+1)==gcd(x+1,2),'gcd reduction')
    if x%2:
        require(square((x+1)*(x*x+1)) == (square((x+1)//2) and square((x*x+1)//2)),'odd reduction')
    else: require(not square((x+1)*(x*x+1)),'even exclusion')
family=[]; lower=7
for _ in range(12):
    row=dependency(lower); family.append(row); lower=max(row['primes'])
for i,a in enumerate(family):
    for b in family[:i]: require(gcd(a['n'],b['n'])==1,'coprime blocks')
    require(gcd(a['n'],343)==1 and square(400*a['sigma']),'343 extension')
    # These rows are manageable enough for full trial-factor sigma verification.
    require(sigma(a['n'])==a['sigma'],'block sigma')
orbit=[1001,343,57,23,1,0]
edges=[]
for n,t in zip(orbit,orbit[1:]):
    require(sigma(n)-n==t,'orbit')
    edges.append({'n':n,'s':t,'sigma':n+t,'square':square(n+t),'positive_simple_square_edge': t>0 and n!=t and square(n+t)})
hostile={str(n):{'sigma':sigma(n),'aliquot':sigma(n)-n,'square':square(sigma(n))} for n in [1,3,4,8,22,40,81,115,243,343,343**2]}
require(sigma(115)==144 and sigma(115)-115==29,'prime target witness')
require(sigma(243)==364 and not square(sigma(243)),'noncoprime closure failure')
require(sigma(343**2)==137257 and not square(sigma(343**2)),'power closure failure')
small_fibers={}
for t in [3,7,11,13,40,57]:
    bound=(t-1)**2; sums=[0]*(bound+1)
    for d in range(1,bound+1):
        for n in range(2*d,bound+1,d): sums[n]+=d
    fiber=[n for n in range(1,bound+1) if sums[n]==t]
    sub=[n for n in fiber if square(n+t)]
    small_fibers[t]={'bound':bound,'fiber':fiber,'square_subfiber':sub}
require(small_fibers[3]['fiber']==[4] and small_fibers[3]['square_subfiber']==[],'fiber3')
require(small_fibers[7]['fiber']==[8] and small_fibers[7]['square_subfiber']==[],'fiber7')
require(small_fibers[11]['fiber']==[21] and small_fibers[11]['square_subfiber']==[],'fiber11')
require(small_fibers[13]['fiber']==[27,35] and small_fibers[13]['square_subfiber']==[],'fiber13')
require(small_fibers[40]['fiber']==[44,74,81] and small_fibers[40]['square_subfiber']==[81],'fiber40')
require(small_fibers[57]['fiber']==[99,159,343,559,703] and small_fibers[57]['square_subfiber']==[343],'fiber57')
require(sigma(159)==216 and sigma(159)-159==57,'same observation different square-edge predicate')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent,
                    help='directory containing check_aliquot_fibers.py and its generated results.json')
args=parser.parse_args()
source_path=args.data_dir / 'results.json'
producer_path=args.data_dir / 'check_aliquot_fibers.py'
if not source_path.is_file() or not producer_path.is_file():
    raise SystemExit('Missing inverse-fiber producer/results. Run python3 check_aliquot_fibers.py, then python3 independent_prefix_audit.py in their directory first; place files adjacent or pass --data-dir DIRECTORY. Standalone small tests completed before this dependency check.')
spec=importlib.util.spec_from_file_location('verified_inverse_producer',producer_path)
producer=importlib.util.module_from_spec(spec); spec.loader.exec_module(producer)
requested_source=json.loads(source_path.read_text())
expected_counts={'1001':29,'1000001':5408,'1001000':2}
expected_exceptions={'1001':[11191],'1000001':[2660799,13011919,16390159,30001459,77519461,190124401],'1001000':[718720,1001998998001]}
require(set(requested_source)==set(expected_counts),'exact target set')
for ts,row in requested_source.items():
    require(row['count']==len(row['fiber'])==expected_counts[ts],'expected target count')
    require(row['non_semiprime_sources']==expected_exceptions[ts],'expected exception list')
    regenerated,_=producer.fiber(int(ts))
    require(row['fiber']==regenerated,'full regenerated fiber agreement')
requested={}
for ts,row in requested_source.items():
    t=int(ts); sub=[n for n in row['fiber'] if square(n+t)]
    require(sub==[], 'requested square subfiber empty')
    requested[ts]={'source_count':len(row['fiber']),'square_subfiber':sub,
                   'exceptions':[{'n':n,'sum':n+t,'floor_root':isqrt(n+t)} for n in row['non_semiprime_sources']]}
require(factors(501)=={3:1,167:1} and factors(500001)=={3:1,166667:1},'hypotenuse factors')
require(all(q%4==3 for c in [501,500001] for q in factors(c)), 'hypotenuse obstruction')
require(sigma(119)-119==25 and sigma(119)==144,'positive Pythagorean prime guard')
# c=5 has integer legs but neither orientation yields two primes h-r,h+r.
require([(r,isqrt(25-r*r)) for r in range(1,5) if square(25-r*r)]==[(3,4),(4,3)],'c5 legs')
linear_family=[]
for k in [3,6]:
    p=2*k-1; q=8*k-1; n=p*q; t=10*k-1
    require(factors(p)=={p:1} and factors(q)=={q:1},'linear family primality')
    require(sigma(n)==16*k*k and sigma(n)-n==t,'linear family identity')
    linear_family.append({'k':k,'p':p,'q':q,'n':n,'t':t,'sigma':sigma(n)})
print(json.dumps({'linear_family_controls':linear_family,'requested_square_subfibers':requested,'generated_results_sha256':sha256(source_path.read_bytes()).hexdigest(),'small_complete_fibers':small_fibers,'scope':'Finite checks supplement cited global theorem; no finite extrapolation','polynomial_range':[-100,N],'polynomial_solutions':solutions,'prime_power_bound':N,'prime_power_solutions':pp,'orbit_edges':edges,'coprime_squarefree_blocks':family,'hostile_examples':hostile,'checks':'PASS'},indent=2,sort_keys=True))
