#!/usr/bin/env python3
"""Exact language-only controls; no tournament construction is assumed."""
import itertools,json,hashlib
from pathlib import Path
A='abe'
def safe(w):return 'be' not in w and 'ee' not in w
def cyclic(w):return all(w[i]+w[(i+1)%len(w)] not in ('be','ee') for i in range(len(w)))
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)]for row in A]
def power(A,n):
 R=[[int(i==j)for j in range(len(A))]for i in range(len(A))]
 while n:
  if n&1:R=mm(R,A)
  A=mm(A,A);n//=2
 return R
def trace(A):return sum(A[i][i]for i in range(len(A)))
def phi(n):return sum(__import__('math').gcd(n,k)==1 for k in range(1,n+1))
Q=[[1,1,1],[1,1,0],[1,1,0]]
B=[[0,0,0],[0,0,1],[0,0,1]]
rows=[]
for m in range(10):
 words=[''.join(w)for w in itertools.product(A,repeat=m)]
 good=sum(map(safe,words)); v=[[1,0]]; T=[[1,2],[1,1]]
 dp=sum(mm(v,power(T,m))[0]); assert good==dp
 if m>=2:assert good==2*rows[-1]['linear_safe']+rows[-2]['linear_safe']
 hits=[sum(next((j+1 for j in range(m-1)if w[j:j+2]in ('be','ee')),None)==k for w in words)for k in range(1,m)]
 formula=[sum(sum(row)for row in mm(power(Q,j-1),B))*3**(m-j-1)for j in range(1,m)]
 assert hits==formula and sum(hits)+good==3**m
 row={'length':m,'linear_safe':good,'first_bad_pair_positions':hits}
 if m:
  cwords=[w for w in words if cyclic(w)]
  rooted=trace(power(Q,m));assert rooted==len(cwords)
  necklaces={min(w[j:]+w[:j]for j in range(m))for w in cwords}
  burnside=sum(phi(d)*trace(power(Q,m//d))for d in range(1,m+1)if m%d==0)//m
  assert len(necklaces)==burnside
  row.update(cyclic_rooted=rooted,cyclic_necklaces=burnside)
 rows.append(row)
weighted_checks=[]
for x,y,z in [(1,1,1),(2,3,5),(0,2,3),(3,0,2)]:
 vals=[]
 for m in range(7):
  brute=sum(__import__('math').prod(dict(zip(A,(x,y,z)))[c]for c in w)for w in map(''.join,itertools.product(A,repeat=m))if safe(w))
  transfer=sum(mm([[1,0]],power([[x,y+z],[x,y]],m))[0]);assert brute==transfer
  if m>=2:assert brute==(x+y)*vals[-1]+x*z*vals[-2]
  vals.append(brute)
 weighted_checks.append({'weights':[x,y,z],'counts':vals})
assert safe('eab') and not safe('abe')
assert safe('eb') and not safe('be')
claimed=[1,3,7,16,36,81,182,409,919]
fit=[1,3,7]
for i in range(3,13):fit.append(2*fit[-1]+fit[-2]-fit[-3])
assert fit[:len(claimed)]==claimed
bulk=[[1,1,1],[1,1,0],[1,0,0]]
marked=[mm([[1,2,0]],power(bulk,m))[0][0]for m in range(13)]
assert marked==fit
# Independent paths with explicit duplicated initial marks.
for m in range(8):
 brute=0
 for start_mark,start_state in [('a',0),('b0',1),('b1',1)]:
  for tail in itertools.product(range(3),repeat=m):
   path=(start_state,)+tail
   brute+=path[-1]==0 and all(bulk[x][y]for x,y in zip(path,path[1:]))
 assert brute==marked[m]

out={'scope':'Language only; alphabet abe is a test assumption, not recovered insertion grammar. No matrix singularity evaluated.', 'rows':rows,'weighted_checks':weighted_checks,'candidate_recurrence':fit,'marked_boundary_walk_counts':marked,'marked_boundary_walk_bruteforce_max_edges':7,'first_free_ternary_mismatch':{'length':3,'free_safe':17,'target':16},'rotation_hostile':['eab','abe'],'reversal_hostile':['eb','be']}
path=Path(__file__).with_name('receipt.json');path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
