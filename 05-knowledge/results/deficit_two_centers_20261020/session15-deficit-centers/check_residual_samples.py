#!/usr/bin/env python3
"""Hostile sample checks of generator and C++ predicate, separate from exhaustive scan."""
import json,random,itertools
from pathlib import Path
from verify_local_configuration import strong_cuts,strong_reach,witness

def construct(row,mask):
 n=row['n'];a=[0]*n
 choices={int(k):v for k,v in row['fixed'].items()}
 choices.update({k:mask>>j&1 for j,k in enumerate(row['free'])})
 for k,bb in enumerate(row['bits']):
  for i,v in enumerate(bb[choices[k]]):a[i]|=v
 return a

def prefix(a):
 n=len(a)
 for S in itertools.combinations(range(n),4):
  if not strong_cuts(a,S):continue
  W=[v for v in range(n) if v not in S]
  if not any(sum(a[u]>>v&1 for u in S)>=3 for v in W):continue
  if strong_reach(a,W):return list(S)
 return None

def main():
 root=Path(__file__).parent;rows=json.loads((root/'residual_frame_specs.json').read_text());rng=random.Random(7344421);out=[]
 for row in rows:
  n,q=row['n'],row['q'];N=1<<len(row['free']);masks=list(range(N)) if N<=4096 else [0,N-1]+[rng.randrange(N) for _ in range(512)]
  for mask in masks:
   a=construct(row,mask)
   assert all(not (a[i]>>i&1) for i in range(n))
   assert all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(n) for j in range(i+1,n))
   assert all(a[i]>>((i+1)%n)&1 for i in range(n))
   assert all((a[i]>>j&1)==(a[(i+q)%n]>>((j+q)%n)&1) for i in range(n) for j in range(n))
   if row['anti']:assert all((a[i]>>j&1)==(a[n-1-j]>>(n-1-i)&1) for i in range(n) for j in range(n))
   if row['low_out_orbit'] is not None:assert a[row['low_out_orbit']].bit_count()==4
   assert prefix(a) is not None,(n,mask)
  out.append(dict(n=n,q=q,low=row['low_out_orbit'],tested=len(masks),frame_invariants=True,independent_prefix_predicate=True))
  print(out[-1],flush=True)
 (root/'residual_sample_checks.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
