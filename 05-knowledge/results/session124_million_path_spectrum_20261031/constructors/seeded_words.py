from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
from itertools import combinations
import json

def profile(n,code):
 out=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j) if code>>k&1 else (j,i)
  out[a]|=1<<b
 N=1<<n
 dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if not(mask>>v&1):continue
   rest=mask^(1<<v)
   if rest:dp[mask][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and out[u]>>v&1)
 H=[sum(x) for x in dp];H[0]=1
 # start through reverse relation recursion
 starts=[[0]*n for _ in range(N)]
 for v in range(n):starts[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if not(mask>>v&1):continue
   rest=mask^(1<<v)
   if rest:starts[mask][v]=sum(starts[rest][u] for u in range(n) if rest>>u&1 and out[v]>>u&1)
 full=N-1
 for cut in range(N):
  s=sum(starts[full][v] for v in range(n) if cut>>v&1)
  t=sum(dp[full][v] for v in range(n) if not(cut>>v&1))
  q=sum(H[full^Q]*sum(starts[Q][v] for v in range(n) if cut>>v&1) for Q in range(1,N))
  # F = sum ordered path covers P,Q allowed ends and starts
  F=0
  for Q in range(N):
   P=full^Q
   a=sum(dp[P][v] for v in range(n) if not(cut>>v&1)) if P else 1
   b=sum(starts[Q][v] for v in range(n) if cut>>v&1) if Q else 1
   F+=a*b
  yield (H[full],F,q,t),(n,code,cut)

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser(description='Exact all-finite-singleton-word support by state deduplication.')
 parser.add_argument('max_seed_order',type=int,nargs='?',default=5)
 parser.add_argument('ceiling',type=int,nargs='?',default=1000000)
 parser.add_argument('--output-roots',type=Path,required=True)
 args=parser.parse_args()
 m=args.max_seed_order
 if not 1<=m<=5 or args.ceiling<1:raise ValueError('Require 1<=max_seed_order<=5 and positive ceiling')
 roots={}
 for n in range(1,m+1):
  for code in range(1<<(n*(n-1)//2)):
   for state,w in profile(n,code):roots.setdefault(state,w)
 print('roots',m,len(roots))
 N=args.ceiling;seen=set();todo=list(roots);vals={}
 while todo:
  st=todo.pop()
  if st in seen:continue
  h,F,q,t=st
  if F>N:continue
  seen.add(st)
  if h<1 or h%2!=1 or F<1 or F%2!=1 or not 0<=t<=h or q<0:raise RuntimeError(('invalid_state',st))
  vals.setdefault(F,st)
  if q>N:continue
  a=(h,F,2*q+h,0)
  b=(h,F+q+h-t,2*q,h)
  if a[2]<=N:todo.append(a)
  if b[1]<=N:todo.append(b)
 print('states',len(seen),'values',len(vals),'first gaps',[k for k in range(1,N+1,2) if k not in vals and k not in(7,21)][:40])
 if m==5 and N>=103061:
  missing=[k for k in range(1,103061,2) if k not in vals and k not in(7,21)]
  if missing != [34869,85885]:raise RuntimeError(('prefix_misses',missing))
  for factors in [(3,59,197),(5,89,193)]:
   if not all(x in vals for x in factors):raise RuntimeError(('missing_factor',factors))
  print('SCC_product_completed_allowed_prefix',103059,'next_lane',103061)
 open(args.output_roots,'w').write(json.dumps([{'state':k,'witness':v}for k,v in roots.items()]))
