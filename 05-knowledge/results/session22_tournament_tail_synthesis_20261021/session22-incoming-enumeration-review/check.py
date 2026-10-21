"""Bounded independent incoming-module review; no benchmark or upstream writes."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations, product
from math import gcd, factorial, comb
import argparse, importlib.util, json, hashlib, statistics
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo-root', required=True, type=Path, help='Path to the mathematical repository root')
ROOT=parser.parse_args().repo_root.resolve()
PACK=ROOT/'05-knowledge/results/tournament_tail_kernels_20261021b'
if not (PACK/'check.py').is_file():
 parser.error('--repo-root must contain 05-knowledge/results/tournament_tail_kernels_20261021b/check.py')
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('incoming',PACK/'check.py');local=load('local',PACK/'incoming_local.py')
checks=Counter()
def need(ok,label):
 checks[label]+=1
 if not ok:raise ValueError(label)
def parts(n,low=1):
 if n==0:yield ()
 for j in range(low,n+1,2):
  for p in parts(n-j,j):yield (j,)+p
def direct(N,core,nofixed):
 out=[F(0)]*(N+1)
 for mass in range(N+1):
  for p in parts(mass,3 if nofixed else 1):
   mult=Counter(p);z=1
   for ell,k in mult.items():z*=ell**k*factorial(k)
   e=sum((ell-1)//2 for ell in p)+sum(gcd(a,b) for a,b in combinations(p,2))
   poly={mass:F(2**e,z)}
   if core:
    for ell,k in mult.items():
     nxt={}
     for d,c in poly.items():
      for excess in range((N-d)//ell+1):
       q=d+ell*excess;nxt[q]=nxt.get(q,F(0))+c*(-1)**excess*comb(k+excess-1,excess)
     poly=nxt
   for d,c in poly.items():out[d]+=c
 return out
for N in range(25):
 for core,nofixed in product((False,True),repeat=2):
  ref=direct(N,core,nofixed)
  for fn in (m.list_poly_count,m.module_dp,local.local_dp):
   need(fn(N,core=core,no_fixed=nofixed)[0]==ref,'direct partitions all four modes')
  need(m.list_poly_count(N,core=core,no_fixed=nofixed,terminal=False)[0]==ref,'terminal fold versus full payload')
  if not nofixed:
   for switch in (1,3,5,9):need(m.hybrid_dp(N,core=core,switch=switch)[0]==ref,'hybrid exact basis change')
for mass in range(25):
 for p in parts(mass):
  t=len(p);z=1
  for ell,k in Counter(p).items():z*=ell**k*factorial(k)
  e=sum((a-1)//2 for a in p)+sum(gcd(a,b) for a,b in combinations(p,2));g=e-comb(t,2)
  need(g>=0 and (factorial(mass)*2**g)%z==0,'normalized history integer')
  for ell,k in product((1,3,5,7),range(4)):
   R=sum(gcd(a,ell) for a in p);S=R-t
   inc=k*(ell-1)//2+ell*comb(k,2)+k*R
   need(inc-comb(t+k,2)+comb(t,2)==(ell-1)*k*k//2+k*S,'normalized insertion')
   need(factorial(mass+ell*k)%(factorial(mass)*ell**k*factorial(k))==0,'local transfer integer')
need(2**comb(2,2)!=2**comb(1,2)*2**comb(1,2),'functional not multiplicative')
need(m.module_dp(3,no_fixed=True)[0][3]==F(2,3),'rational observer')
b=json.loads((PACK/'benchmarks.json').read_text());r=json.loads((PACK/'receipt.json').read_text())
for fn,p in [('check.py',PACK/'check.py'),('incoming_local.py',PACK/'incoming_local.py'),('previous_check.py',PACK.parent/'tournament_enumeration_20261021a/check.py')]:
 need(hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==b['source_hashes'][fn],'benchmark source pin')
for med in b['medians']:
 rows=[x for x in b['runs'] if x['algorithm']==med['algorithm'] and x['bound']==med['bound'] and x['core']==med['core']]
 need(len(rows)==3 and statistics.median(x['seconds'] for x in rows)==med['median_seconds'],'median exact receipt consistency')
for bound,core,name in ((200,False,'tournaments'),(100,True,'pair_free_cores')):
 rows=[x for x in b['runs'] if x['bound']==bound and x['core']==core]
 h=hashlib.sha256('\n'.join(r['series'][name]['coefficients']).encode()).hexdigest()
 need(all(x['coefficient_sha256']==h for x in rows),'all timed outputs equal saved vector')
print(json.dumps({'status':'PASS','total':sum(checks.values()),'checks':dict(checks),'source_pins':{p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in PACK.iterdir() if p.is_file()}},indent=2))
