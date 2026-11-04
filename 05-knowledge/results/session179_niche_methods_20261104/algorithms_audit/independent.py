#!/usr/bin/env python3
from fractions import Fraction
from pathlib import Path
import json
from math import isqrt
ROOT=Path(__file__).resolve().parent.parent/'crossfield'
def require(x,m):
 if not x:raise RuntimeError(m)
r=json.loads((ROOT/'results.json').read_text())
h={int(k,2):Fraction(v) for k,v in r['window_cap_certificate']['h'].items()}
rho=Fraction(r['window_cap_certificate']['rho']);B=Fraction(r['window_cap_certificate']['B'])
E={i:[] for i in range(16)}
for i in range(16):
 for b in (0,1):
  if i.bit_count()+b<=3:
   j=((i<<1)&15)|b; E[i].append((j,b))
   require(Fraction(3 if b else 1,2)*h[j]<=rho*h[i],'edge')
   require(Fraction(b,2)*h[j]<=B,'additive')
cycles=[]
def dfs(start,at,seen,bits):
 for j,b in E[at]:
  if j==start:
   w=bits+[b];require(3**sum(w)<2**len(w),'cycle contraction');cycles.append((len(w),sum(w)))
  elif j>start and j not in seen:dfs(start,j,seen|{j},bits+[b])
for i in range(16):dfs(i,i,{i},[])
# Derive sharp finite-step endpoint independently by actual binary source lifts.
count=0;max_intercept=Fraction(0)
for residue in range(32):
 n=residue;bits=[]
 for _ in range(5):
  b=n%2;bits.append(b);n=n//2 if not b else (3*n+1)//2
 if sum(bits)<=3:
  a=Fraction(3**sum(bits),32);c=n-a*residue
  require(a<=Fraction(27,32) and c<=Fraction(19,8),'direct endpoint')
  max_intercept=max(max_intercept,c);count+=1
require(count==26 and max_intercept==Fraction(19,8),'sharp words')
# Verify each prime column with independent brute divisor sum and factor extraction.
v=json.loads((ROOT/'valuation_results.json').read_text());W={int(k):x for k,x in v['weights'].items()}
rows=0
for row in v['rows']:
 p,a=row['p'],row['a'];N=p**a
 sigma=sum(d for d in range(1,N+1) if N%d==0)
 require(sigma==row['sigma'],'independent divisor sum')
 columns={int(k):x for k,x in row['column'].items()}
 num=1;den=1
 for q,e in columns.items():
  require(q>=2 and all(q%d for d in range(2,isqrt(q)+1)),'prime coordinate')
  if e>=0:num*=q**e
  else:den*=q**(-e)
 require(Fraction(num,den)==Fraction(sigma,N),'valuation ratio')
 score=sum(W.get(q,256)*e for q,e in columns.items())
 require(score==row['dual_score'] and score>=0,'dual score')
 rows+=1
require(rows==30 and W[2]==-4,'finite separating family')
print(json.dumps({'status':'PASS','simple_cycles':len(cycles),'cycle_length_odd_counts':sorted(set(cycles)),'sharp_five_step_words':count,'independent_divisor_sum_columns':rows,'note':'No source module imported; finite checks do not prove Collatz or exclude all odd perfect numbers.'},indent=2))
