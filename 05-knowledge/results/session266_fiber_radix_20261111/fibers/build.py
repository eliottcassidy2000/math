"""Compile inherited exact H<=73 atom input into arbitrary-order multiplicities.
Run: python build.py [--n 30] [--h 63] [--output report.json]. No census is rerun.
"""
import json,math,argparse
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
p=Path(__file__).parent
args=argparse.ArgumentParser(description='Exact multiplicities for every order n and odd H<=73; larger H is not covered by the frozen atom input.')
args.add_argument('--n','--N',type=int,default=30,help='Largest order, nonnegative (default 30)')
args.add_argument('--h','--H',type=int,help='One odd path count from 1 through 73; omit for all')
args.add_argument('--output',default='report.json');a=args.parse_args()
if a.n<0:args.error('--n must be nonnegative')
if a.h is not None and (a.h<1 or a.h>73 or a.h%2==0):args.error('--h must be odd and between 1 and 73; this frozen inventory cannot answer other H')
targets=[a.h] if a.h is not None else range(1,74,2)
D=json.loads((p/'atoms.json').read_text());atoms=D['atoms']
@lru_cache(None)
def words(h,channel):
 if h==1:return {(0,0):Q(1)}
 out={}
 for row in atoms:
  d=row['h']
  if h%d:continue
  w=Q(row[channel]);s=row['n']
  for (k,t),v in words(h//d,channel).items():out[k+1,t+s]=out.get((k+1,t+s),Q(0))+w*v
 return out

def f(n,h,channel='ordinary',cumulative=False):
 return sum((v*math.comb(n-s+k+int(cumulative),k+int(cumulative)) for (k,s),v in words(h,channel).items() if n>=s),Q(0))
def sc(n,h):
 total=0
 for x in range(1,math.isqrt(h)+1):
  if h%(x*x):continue
  e=h//(x*x)
  centers=[(0,1),(1,1)] if e==1 else [(r['n'],r['self_converse']) for r in atoms if r['h']==e]
  for s,v in centers:
   if n>=s and (n-s)%2==0:total+=v*f((n-s)//2,x)
 return total
checks=0
out={}
for h in targets:
 rows=[]
 for n in range(a.n+1):
  ordinary=f(n,h);sym=sc(n,h);weighted=f(n,h,'weight');lab=math.factorial(n)*weighted
  if lab.denominator!=1:raise RuntimeError('nonintegral labeled count')
  if (ordinary+sym)%2:raise RuntimeError('nonintegral merged count')
  if h>1 and ordinary-2*f(n-1,h)+f(n-2,h)<0:raise RuntimeError('convexity')
  if f(n,h,cumulative=True)!=sum(f(j,h) for j in range(n+1)):raise RuntimeError('cumulative')
  rows.append({'n':n,'ordinary':int(ordinary),'self_converse':int(sym),'merged':int((ordinary+sym)/2),'labeled':int(lab),'cumulative_labeled':int(sum(math.factorial(j)*f(j,h,'weight') for j in range(1,n+1))),'cumulative_ordinary':int(f(n,h,cumulative=True)-f(0,h)),'cumulative_merged':int(sum((f(j,h)+sc(j,h))/2 for j in range(1,n+1)))})
  checks+=4
 out[h]=rows
# Invert SCC recurrence on the complete divisor-closed target range.
for h in targets:
 if h==1:continue
 for n in range(a.n+1):
  residual=f(n,h)-2*f(n-1,h)+f(n-2,h)
  for r in atoms:
   d=r['h']
   if d<h and h%d==0:residual-=r['ordinary']*(f(n-r['n'],h//d)-f(n-r['n']-1,h//d))
  expected=sum(r['ordinary'] for r in atoms if r['h']==h and r['n']==n)
  if residual!=expected:raise RuntimeError(('inverse',h,n,residual,expected))
  checks+=1
R={'scope':'All n for odd h<=73, inherited complete atom inventory; order-zero ordinary row is auxiliary empty tournament only at h=1; every cumulative field sums positive sizes 1..N and excludes it.','max_n':a.n,'active_checks':checks,'counts':out}
(p/a.output).write_text(json.dumps(R,indent=2)+'\n')
print(json.dumps({'active_checks':checks,'output':str(p/a.output),'requested_at_N':{h:rows[-1] for h,rows in out.items() if h==a.h or (a.h is None and h in (1,7,21,63))}}))
