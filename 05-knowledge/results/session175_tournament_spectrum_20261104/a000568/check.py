#!/usr/bin/env python3
import sys,json,time,hashlib,platform,statistics
from pathlib import Path
from math import factorial,gcd,prod
from fractions import Fraction
from collections import Counter
from projective import count,packed_count,need
from inherited_local_scale import local_dp
HERE=Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)
checks=0
def check(b,s):
 global checks
 checks+=1;need(b,s)
def parts(n,lo=1):
 if not n:yield ()
 for j in range(lo,n+1,2):
  for p in parts(n-j,j):yield (j,)+p
def direct(n):
 total=Fraction(0)
 for p in parts(n):
  e=sum((j-1)//2 for j in p)+sum(gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
  z=prod(j**k*factorial(k) for j,k in Counter(p).items())
  total+=Fraction(1<<e,z)
 return total
ref={}
for line in (HERE/'a000568-reference.txt').read_text().splitlines():
 if line.strip() and not line.startswith('#'):
  n,v=map(int,line.split());ref[n]=v
for n in range(77):
 a,s=packed_count(n);check(a==ref[n],'OEIS table')
 b,t=packed_count(n,False);check(a==b,'fixed-target ungauged')
 c,u=count(n);check(a==c,'unpacked projective')
 if n<=35:check(a==direct(n),'independent partition sum')
for bad in [-1,True,1.5,'3']:
 try:packed_count(bad)
 except ValueError:pass
 else:raise ValueError('invalid order accepted')
for n in [100,150,200]:
 a,s=packed_count(n);b,t=local_dp(n);check(a==b[-1],'inherited all-degree integer DP')
# Precise collision and hostile controls.
def response(p,j):return sum(gcd(a,j) for a in p)
def S(p,j):return j*len(p)-response(p,j)
check(sum((27,))==sum((9,9,9)) and S((27,),3)==S((9,9,9),3),'new genuine merge')
check(response((27,),1)!=response((9,9,9),1),'old profile separates witness')
check(S((11,5),3)!=S((9,7),3),'same mass cycle count is insufficient')
# Individual completion normalized cross exponents coincide exactly iff S does.
for M in range(0,26):
 for p in [(27,),(9,9,9)]:
  for k in range(M//3+1):
   r=M-3*k
   lhs=r*len(p)+k*response(p,3)-M*len(p)
   check(lhs==-k*S(p,3),'normalized completion identity')
bench=[]
for n in [80,100,150,200]:
 times={'projective':[],'ungauged_fixed':[]};stats={};vals={}
 for repeat in range(3):
  order=[('projective',True),('ungauged_fixed',False)]
  if repeat%2:order.reverse()
  for name,mode in order:
   t=time.perf_counter();vals[name],stats[name]=packed_count(n,mode);times[name].append(time.perf_counter()-t)
  check(vals['projective']==vals['ungauged_fixed'],'benchmark same answer')
 bench.append({'order':n,'seconds':times,'median_speedup':statistics.median(times['ungauged_fixed'])/statistics.median(times['projective']),'stats':stats,'answer':str(vals['projective'])})
receipt={'checks':checks,'python':sys.version,'platform':platform.platform(),'reference_url':'https://oeis.org/A000568/b000568.txt','reference_scope':'orders 0..76; inherited cached primary OEIS table','benchmarks':bench,'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'projective.py',HERE/'inherited_local_scale.py',HERE/'a000568-reference.txt',HERE/'check.py']}}
out=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'receipt.json'
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'checks':checks,'benchmarks':[{k:r[k] for k in ('order','median_speedup')} for r in bench]}))
