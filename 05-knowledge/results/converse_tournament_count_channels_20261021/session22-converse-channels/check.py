#!/usr/bin/env python3
"""Independent rational-partition, permutation-action and literal census controls."""
import argparse, hashlib, itertools, json, math, platform, sys
from fractions import Fraction
from pathlib import Path
from channels import channel_dp, counts
from incoming_local_scale import local_dp
HERE=Path(__file__).resolve().parent
checks=0

def need(condition,label):
 global checks
 checks+=1
 if not condition: raise ArithmeticError(label)

def partitions(n,start=1):
 if not n: yield ()
 for j in range(start,n+1,2):
  for rest in partitions(n-j,j): yield (j,)+rest

def independent(bound,a,b,core=False,no_fixed=False):
 out=[Fraction(0)]*(bound+1)
 for n in range(bound+1):
  for parts in partitions(n,3 if no_fixed else 1):
   ks={p:parts.count(p) for p in set(parts)}
   E=sum((p-1)//2 for p in parts)+sum(math.gcd(p,q) for i,p in enumerate(parts) for q in parts[i+1:])
   coeff=Fraction(2**(a*E+b*len(parts)),math.prod(p**k*math.factorial(k) for p,k in ks.items()))
   poly={n:coeff}
   if core:
    for p,k in ks.items():
     new={}
     for d,v in poly.items():
      for e in range((bound-d)//p+1):
       nd=d+p*e
       new[nd]=new.get(nd,0)+v*(-1)**e*math.comb(k+e-1,e)
     poly=new
   for d,v in poly.items():out[d]+=v
 return out

def edge_burnside(n,anti=False):
 edges=list(itertools.combinations(range(n),2)); total=0
 for perm in itertools.permutations(range(n)):
  unseen=set(edges); orbits=0; valid=True
  while unseen:
   first=next(iter(unseen));cur=first;parity=0
   while True:
    unseen.remove(cur); u,v=perm[cur[0]],perm[cur[1]]
    parity^=(u>v)^anti;cur=tuple(sorted((u,v)))
    if cur==first:break
   if parity:valid=False;break
   orbits+=1
  if valid:total+=2**orbits
 return Fraction(total,math.factorial(n))

def census(n):
 edges=list(itertools.combinations(range(n),2));index={e:i for i,e in enumerate(edges)}
 maps=[]
 for perm in itertools.permutations(range(n)):
  maps.append([(index[tuple(sorted((perm[u],perm[v])))],perm[u]>perm[v]) for u,v in edges])
 def transform(code,mapping):
  return sum((((code>>i)&1)^swap)<<j for i,(j,swap) in enumerate(mapping))
 def strong(code):
  if not n:return False
  adj=[0]*n
  for i,(u,v) in enumerate(edges):
   if (code>>i)&1:adj[u]|=1<<v
   else:adj[v]|=1<<u
  for u in range(n):
   reached=1<<u
   while True:
    nex=reached
    for v in range(n):
     if reached>>v&1:nex|=adj[v]
    if nex==reached:break
    reached=nex
   if reached!=(1<<n)-1:return False
  return True
 unseen=set(range(1<<len(edges)));A=C=S=U=0;full=(1<<len(edges))-1
 while unseen:
  code=min(unseen);orbit={transform(code,mapping) for mapping in maps};unseen-=orbit
  sc=(code^full) in orbit;st=strong(code)
  A+=1;C+=sc;S+=st;U+=sc and st
 return {'ordinary':A,'self_converse':C,'merged':(A+C)//2,'strong':S,'strong_self_converse':U,'strong_merged':(S+U)//2}

def control():
 for n in range(21):
  for core,no_fixed in itertools.product((False,True),repeat=2):
   old,oldstats=local_dp(n,core=core,no_fixed=no_fixed)
   for a,b in [(1,0),(1,2),(2,0),(2,1),(3,2)]:
    out,stats=channel_dp(n,a,b,core=core,no_fixed=no_fixed)
    need(out==independent(n,a,b,core,no_fixed),'unmerged partitions')
    if (a,b)==(1,0):
     need(out==old and stats==oldstats,'frozen kernel exact API/telemetry match')
 rows=counts(40)
 for n in range(8):
  need(rows['ordinary'][n]==edge_burnside(n),'ordinary full permutation edge action')
  need(rows['self_converse'][n]==edge_burnside(n,True),'anti full permutation edge action')
 censuses=[]
 for n in range(7):
  row=census(n);censuses.append({'n':n,**row})
  for k,v in row.items():need(rows[k][n]==v,'literal class census '+k)
 known_sc=[1,1,1,2,2,8,12,88,176,2752,8784,279968,1492288,95458560,872687552,111698291584,1787154671104,457509297625088,13013584213369088,6662951988432581120]
 need(rows['self_converse'][:20]==known_sc,'THM-283 full published prefix')
 for bad in [-1,True,1.5,'3']:
  for kwargs in [{'bound':bad},{'bound':3,'a':bad},{'bound':3,'b':bad}]:
   try:channel_dp(**kwargs)
   except ValueError:ok=True
   else:ok=False
   need(ok,'invalid literal integer rejected')
 for kwargs in [{'bound':2,'a':0},{'bound':2,'core':1},{'bound':2,'no_fixed':0}]:
  try:channel_dp(**kwargs)
  except ValueError:ok=True
  else:ok=False
  need(ok,'invalid modes rejected')
 return {'checks':checks,'literal_census':censuses,'counts_through_40':rows}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--receipt',default='receipt.json');args=p.parse_args()
 receipt=control();receipt.update({'python':sys.version,'platform':platform.platform(),'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['channels.py','incoming_local_scale.py','check.py']}})
 (HERE/args.receipt).write_text(json.dumps(receipt,indent=2)+'\n');print('PASS',checks)
