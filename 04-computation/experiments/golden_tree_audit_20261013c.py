#!/usr/bin/env python3
"""Exact golden inverse-tree port/count/parity identities and ray control."""
import argparse,json
from functools import lru_cache

def ck(x,m):
 if not x:raise RuntimeError(m)
def G(x):
 a,b=x
 if a%2==b%2==0:return a//2,b//2
 return (b+b%2)//2,(a+b+(a+b)%2)//2
ROOTS=((0,0),(1,1),(1,2),(2,3))
def inverse(y):
 a,b=y;out=[(2*a,2*b)]
 for da,db in [(1,0),(0,1),(1,1)]:
  u,v=2*a-da,2*b-db;out.append((v-u,u))
 return out
@lru_cache(None)
def clock(x):
 if x in ROOTS:return 0
 return 1+clock(G(x))
def main():
 p=argparse.ArgumentParser();p.add_argument('--ray-census',required=True);p.add_argument('--output',required=True);args=p.parse_args()
 byroot={x:i for i,x in enumerate(ROOTS)};levels={x:0 for x in ROOTS};current=set(ROOTS);reports=[];controls=0
 for k in range(8):
  count=[0]*4;signed=[0]*4
  for x in current:
   root=byroot[x];count[root]+=1;signed[root]+=(-1)**levels[x]
   ck(clock(x)==levels[x],'strict clock depth');controls+=1
  expected=(2+3*(-4)**k)//5
  ck(count==[4**k]*4 and signed==[expected]*4,'tree count and amplitude');controls+=1
  reports.append(dict(depth=k,count_per_root=count[0],signed_per_root=signed[0],mean=signed[0]/count[0]))
  if k==7:break
  nxt=set()
  for y in current:
   inv=inverse(y);ck(len(set(inv))==4,'four distinct inverse branches');controls+=1
   for x in inv:
    ck(G(x)==y,'inverse replay');controls+=1
    if x not in levels:levels[x]=k+1;byroot[x]=byroot[y]
    else:ck(byroot[x]==byroot[y],'source port');controls+=1
    nxt.add(x)
  ck(current<=nxt,'depth nesting');current=nxt;controls+=1
 D=json.load(open(args.ray_census));ray_bins=0
 for lo,hi,count,s0,s1 in D['bins']:
  if hi>2**18:break
  vals=[sum((-1)**clock((n,r*n))for n in range(lo,hi))for r in (0,1)]
  ck(vals==[s0,s1] and count==hi-lo,'independent ray sums');controls+=count*2;ray_bins+=1
 out=dict(status='FINITE-EXACT',controls=controls,ray_bins_controlled=ray_bins,ray_source_limit_exclusive=2**18,trees=reports,formula='per root: size=4^K, signed parity=(2+3(-4)^K)/5',scope='certificate-depth balls; these are not integer-ray or Euclidean-box averages')
 with open(args.output,'w')as f:json.dump(out,f,indent=2)
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
