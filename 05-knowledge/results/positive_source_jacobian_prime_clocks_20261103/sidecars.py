"""Replay inherited source/phase controls; no novelty claim for them."""
from pathlib import Path
from fractions import Fraction
import importlib.util
import json
import sys
import sympy as s
from itertools import combinations

HERE=Path(__file__).resolve().parent
prior=HERE.parent/'primegame_port_spectrum_20261103'
sys.path.insert(0,str(prior))
spec=importlib.util.spec_from_file_location('prior_spectrum_verifier',prior/'verify.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def need(x,msg):
    if not x: raise RuntimeError(msg)

def step(n,b):
    v=3*n+b;a=0
    while v%2==0:v//=2;a+=1
    return v,a

cycles=[]
for start in [227,259]:
    vals=[]; word=[];n=start
    for _ in range(5):
        vals.append(n);n,a=step(n,13);word.append(a)
    need(n==start,'cycle return')
    cycles.append(dict(values=vals,valuations=word,total=sum(vals),A=sum(word)))
a=cycles[0]['valuations'][::-1];b=cycles[1]['valuations']
need(any(a==b[j:]+b[:j] for j in range(5)),'reversal rotation')
need(all(c['total']==2499 and c['A']==8 for c in cycles),'cycle invariants')
head=[151]; valuations=[]
for _ in range(3):
    n,a=step(head[-1],1);head.append(n);valuations.append(a)
need(head==[151,227,341,1] and valuations==[1,1,10],'source splice')
need(38**2+7==1451 and 5*290-38**2==6,'norm/progression intersection')
data=json.loads((HERE/'verification.json').read_text())
for row in data['primes']['normalized_records']:
    m=row['m'];p=row['p']
    need(s.isprime(p) and int(s.primepi(p))-len(s.factorint(m))==row['eligible'],'record independent count')
    need(all(not s.isprime(q) for k in range(1,row['k']+1) for q in (k*m-1,k*m+1) if q<p),
         'record independent first hit')
need((HERE/'verification.json').read_bytes()==(HERE/'verification_optimized.json').read_bytes(),
     'normal/O byte agreement')
local_spec=importlib.util.spec_from_file_location('new_profile_checks',HERE/'verify.py')
local=importlib.util.module_from_spec(local_spec);local_spec.loader.exec_module(local)
profiles=[]
for flip in (False,True):
    a=[0]*12
    for u,v in combinations(range(12),2):
        b,i=divmod(u,4);c,j=divmod(v,4)
        arc=(c-b)%3==1 if b!=c else i==0 or (i>0 and (j-i)%3==1)
        if flip and b==c==2:arc=not arc
        a[u if arc else v]|=1<<(v if arc else u)
    profiles.append(local.profile(tuple(a)))
need(profiles[0]==profiles[1],'81/27 full unordered path covers')
print(json.dumps(dict(hostile_81_27=mod.hostile(),cycles_3n13=cycles,
                     hostile_full_path_cover_profiles=profiles,
                     source_splice=dict(values=head,valuations=valuations),
                     independent_prime_records=len(data['primes']['normalized_records']),
                     ordinary_and_optimized_identical=True),indent=2))
