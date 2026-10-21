#!/usr/bin/env python3
"""Independent exact controls and alternating same-host benchmark receipt."""
import argparse,hashlib,itertools,json,math,platform,statistics,time,sys
from fractions import Fraction
from pathlib import Path
from local_scale import local_dp
from baseline import packed_dp
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
checks=0

def need(p,label):
 global checks
 checks+=1
 if not p: raise ValueError(label)

def partitions(n,start=1):
 if not n: yield ()
 for j in range(start,n+1,2):
  for rest in partitions(n-j,j): yield (j,)+rest

def independent(bound,core=False,no_fixed=False):
 out=[Fraction(0)]*(bound+1)
 for mass in range(bound+1):
  for parts in partitions(mass,3 if no_fixed else 1):
   multiplicities={p:parts.count(p) for p in set(parts)}
   exponent=sum((p-1)//2 for p in parts)+sum(math.gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
   denominator=math.prod(p**k*math.factorial(k) for p,k in multiplicities.items())
   coeff=Fraction(1<<exponent,denominator)
   poly={mass:coeff}
   if core:
    for p,k in multiplicities.items():
     new={}
     for d,v in poly.items():
      for e in range((bound-d)//p+1):
       nd=d+p*e
       new[nd]=new.get(nd,0)+v*(-1)**e*math.comb(k+e-1,e)
     poly=new
   for d,v in poly.items(): out[d]+=v
 return out

def labelled_burnside(n):
 # Orbit action on unordered edges; orientation-reversing orbit has no fixed tournament.
 edges=list(itertools.combinations(range(n),2)); total=0
 for perm in itertools.permutations(range(n)):
  unseen=set(edges);cycles=0;valid=True
  while unseen:
   edge=next(iter(unseen)); cur=edge; parity=0
   while True:
    unseen.remove(cur); a,b=perm[cur[0]],perm[cur[1]]
    parity^=a>b;cur=tuple(sorted((a,b)))
    if cur==edge: break
   if parity:valid=False;break
   cycles+=1
  if valid: total+=1<<cycles
 return total//math.factorial(n)

def control():
 for n in range(27):
  for core,nofix in itertools.product([False,True],repeat=2):
   a,s=local_dp(n,core=core,no_fixed=nofix);b,t=packed_dp(n,core=core,no_fixed=nofix)
   need(a==b,'incoming coefficients')
   need(s['rows']==[{k:r[k] for k in ('cycle_length','states','transitions')} for r in t['rows']], 'same response graph')
   if n<=20:need(a==independent(n,core,nofix),'unmerged rational odd partitions')
 for n in range(8):need(local_dp(n)[0][-1]==labelled_burnside(n),'labelled edge action Burnside')
 known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
 need(local_dp(12)[0]==known,'known A000568 prefix')
 for bad in [-1,True,1.5,'3']:
  try:local_dp(bad)
  except ValueError:checks_ok=True
  else:checks_ok=False
  need(checks_ok,'invalid bound rejected')
 return checks

def run_bench(bounds,reps):
 results=[]
 for n in bounds:
  mode_runs=[False,True] if n<=100 else [False]
  for core in mode_runs:
   times={'incoming':[],'local':[]};outputs={};stats={}
   for repeat in range(reps):
    for name,fn in ([('incoming',packed_dp),('local',local_dp)] if repeat%2==0 else [('local',local_dp),('incoming',packed_dp)]):
     t=time.perf_counter();outputs[name],stats[name]=fn(n,core=core);times[name].append(time.perf_counter()-t)
    need(outputs['incoming']==outputs['local'],'benchmark exact equality')
   r={'bound':n,'core':core,'seconds':times,'median_speedup':statistics.median(times['incoming'])/statistics.median(times['local']), 'transitions':stats['local']['transitions'],'precomputed_transfers':stats['local']['precomputed_transfers'],'peak_states':stats['local']['peak_states'],'output_sha256':hashlib.sha256(json.dumps([str(x) for x in outputs['local']]).encode()).hexdigest()}
   results.append(r);print(json.dumps(r),flush=True)
 return results

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--bench',action='store_true');p.add_argument('--receipt',default='receipt.json');args=p.parse_args()
 receipt={'control_checks':control(),'python':sys.version,'platform':platform.platform(),'baseline_origin':'git HEAD:05-knowledge/results/tournament_enumeration_20261021a/check.py','hashes':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['baseline.py','local_scale.py','check.py']}}
 if args.bench:receipt['benchmarks']=run_bench([80,100,150,200],3)
 receipt['checks']=checks
 (HERE/args.receipt).write_text(json.dumps(receipt,indent=2)+'\n');print('PASS',checks,flush=True)
