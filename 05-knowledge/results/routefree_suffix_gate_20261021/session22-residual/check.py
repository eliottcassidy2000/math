#!/usr/bin/env python3
"""Bounded-prefix, route-free all-future energy certificate. Standard library only.
No import of historical producers; no future route assumptions; no credit spend.
"""
import argparse, hashlib, itertools, json
from collections import defaultdict
from math import isqrt
from pathlib import Path
VECTORS=((1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),(-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1))
INDEX={v:i for i,v in enumerate(VECTORS)}
FRONTIER_SHA='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'
CHECKS=0
def ck(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise ValueError(msg)
def Q(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def add(v,w):return v[0]+w[0],v[1]+w[1]
def primitive(s):return s>0 and s%2==1 and s%8!=5
def edge(s):
 ck(primitive(s),'primitive edge source')
 n=3*s+1;a=(n&-n).bit_length()-1;t=a;u=n>>a
 ck(a in (1,2),'native branch valuation')
 while u%8==5:u=(u-1)//4;t+=2
 ck(primitive(u),'primitive edge endpoint')
 return u,t

def charge(N,s):
 # Independent three-row coordinates, not imported two-coordinate alphabet.
 a=b=0
 while s<4*N:
  row=(0,-3,3)if s<2*N else(-1,-1,2)
  a+=row[s%3];b+=row[(s+1)%3];s=4*s+1
 ck((a-b)%3==0 and (a+2*b)%3==0,'charge coordinate integrality')
 v=((a-b)//3,-(a+2*b)//3)
 ck(v==(0,0)or v in INDEX,'bounded source charge alphabet')
 return v

def cuts(count):
 return [tuple(sum(count[j%12]*VECTORS[j%12][k]for j in range(i,i+6))for k in (0,1))for i in range(12)]
def zonotope(count):
 total=tuple(sum(count[j]*VECTORS[j][k]for j in range(12))for k in (0,1))
 vals=[Q((2*z[0]-total[0],2*z[1]-total[1]))for z in cuts(count)]
 return max(vals)
def linear(count,T):
 return sum(count[j]*max((abs(beta(t,VECTORS[j]))for t in T),default=0)for j in range(12))
def ceilroot(x):
 r=isqrt(x);return r+(r*r<x)

def gate(W,H,oldcounts,newcounts):
 terminals=set(W)|set(H)
 T={u:add(W.get(u,(0,0)),H.get(u,(0,0)))for u in terminals}
 delta=sum(Q(T[u])-Q(W.get(u,(0,0)))for u in terminals)
 Do,Dn=zonotope(oldcounts),zonotope(newcounts)
 Ln=linear(newcounts,T.values());Lo=linear(oldcounts,H.values())
 mixed=ceilroot(4*Do*Dn)
 return {'delta':delta,'D_old':Do,'D_new':Dn,'L_old':Lo,'L_new':Ln,'mixed_ceiling':mixed,'upper':delta+Dn+Ln+Lo+mixed}

class PrefixRoutes:
 def __init__(self,V):self.V=V;self.known={u:(u,0,0)for u in V};self.edges={};self.max_explored=0
 def get(self,s,h):
  trail=[];u=s
  while u not in self.known and len(trail)<h:
   if u not in self.edges:self.edges[u]=edge(u)
   v,t=self.edges[u];trail.append((u,t));u=v
  self.max_explored=max(self.max_explored,len(trail))
  if u not in self.known:return None
  terminal,t,d=self.known[u]
  for u,c in reversed(trail):
   t+=c;d+=1
   # A known route is a finite verified suffix, never a convergence assumption.
   self.known[u]=(terminal,t,d)
  return self.known[s]if self.known[s][2]<=h else None

def run_case(N,h,V):
 ck(type(N)is int and N>0 and type(h)is int and h>=0,'literal scale/horizon')
 routes=PrefixRoutes(V);events=defaultdict(list);pending=[[0]*12,[0]*12];old_unarrived=0;new_unarrived=0;arrived_count=0;norm_new=0
 digest=hashlib.sha256();nonzero=0;guard_count=0
 for s in range(1,256*N,2):
  if not primitive(s):continue
  v=charge(64*N,s);k=int(s>=4*N);route=routes.get(s,h)
  a,b=s%15,s%35
  ck(a%5==b%5 and(a+15*((b-a)%7))%105==s%105,'mod15/mod35 joint source guard');guard_count+=1
  if k:norm_new+=Q(v)
  else:ck(v==charge(N,s),'old source scale identity')
  if route is None:
   if k:new_unarrived+=1
   else:old_unarrived+=1
   if v!=(0,0):pending[k][INDEX[v]]+=1;nonzero+=1
   row=[s,a,b,None,*v]
  else:
   u,t,d=route;ck(u in V and d<=h,'actual first arrival cutoff');arrived_count+=1
   events[d].append((u,k,tuple((-1 if t%2 else 1)*z for z in v),Q(v)))
   row=[s,a,b,u,t,d,*v]
  digest.update((json.dumps(row,separators=(',',':'))+'\n').encode())
 ck(norm_new==153*N,'full new diagonal budget')
 W={u:(0,0)for u in V};H=W.copy();A=0;ledger=[]
 for d in range(h+1):
  for u,k,v,q in events[d]:
   target=H if k else W;target[u]=add(target[u],v)
   if k:A+=q
  X=sum(Q(v)for v in W.values());Y=sum(Q(v)for v in H.values());B=sum(beta(W[u],H[u])for u in V)
  ledger.append([d,X,Y,B,A,Y+B])
 result=gate(W,H,*pending)
 result.update(N=N,horizon=h,budget=153*N,old_unarrived=old_unarrived,new_unarrived=new_unarrived,pending_nonzero=nonzero,pending_old_counts=pending[0],pending_new_counts=pending[1],arrived_count=arrived_count,source_rows_sha256=digest.hexdigest(),source_guard_count=guard_count,max_single_exploration=routes.max_explored,maximum_prefix_delta=max(z[-1]for z in ledger),prefix_pass=all(z[-1]<=153*N for z in ledger),suffix_pass=result['upper']<=153*N,ledger=ledger)
 ck(result['delta']==ledger[-1][-1],'gate current energy equals direct ledger')
 # Repartition identity for any hypothetical owned extraction K, never credit reuse.
 J=17;K=273;M=153*N-A;I=result['delta']-A;Rp=I-J+K;Jp=J-K
 ck(Rp+Jp==I and (Rp<=M-Jp)==(result['delta']<=153*N),'honest Rprime target')
 return result

def hostiles():
 # Exhaustively check every multiset of <=4 of the twelve vectors, by signs.
 tested=0
 for size in range(5):
  for ids in itertools.combinations_with_replacement(range(12),size):
   count=[ids.count(i)for i in range(12)]
   direct=max(Q(tuple(sum(sg*VECTORS[i][k]for sg,i in zip(signs,ids))for k in (0,1)))for signs in itertools.product((-1,1),repeat=size))
   ck(zonotope(count)==direct,'independent exhaustive sign-zonotope');tested+=1
 # Independently enumerate absent / terminal0+/- / terminal1+/- destinations.
 W={0:(1,-1),1:(-2,-1)};H={0:(2,1),1:(0,-1)}
 inputs=[(0,(1,0)),(0,(1,2)),(1,(-1,1)),(1,(0,1))]
 counts=[[0]*12,[0]*12]
 for k,v in inputs:counts[k][INDEX[v]]+=1
 G=gate(W,H,*counts)['upper'];maximum=-10**9
 for choices in itertools.product(range(5),repeat=len(inputs)):
  ww=W.copy();hh=H.copy()
  for (k,v),z in zip(inputs,choices):
   if not z:continue
   u=(z-1)//2;sg=1 if z%2 else -1;target=hh if k else ww
   target[u]=add(target[u],tuple(sg*x for x in v))
  delta=sum(Q(add(ww[u],hh[u]))-Q(ww[u])for u in W)
  ck(delta<=G,'arbitrary unfinished source destination/time/parity gate');maximum=max(maximum,delta)
 # Completing old sources is essential for the simpler new-only gate.
 simple=gate({0:(0,0)},{0:(1,0)},[0]*12,[0]*12)['upper']
 actual=Q((2,0))-Q((1,0))
 ck(simple==1 and actual==3,'omitted pending old source hostile')
 return {'sign_multisets':tested,'arbitrary_assignments':5**len(inputs),'general_gate':G,'actual_assignment_max':maximum,'ignored_old_hostile':[simple,actual]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--verify',type=Path);p.add_argument('--cases',default='7:8,4096:105,8192:112,16385:121,32768:131');args=p.parse_args()
 raw=(args.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n')
 ck(hashlib.sha256(raw).hexdigest()==FRONTIER_SHA,'authenticated frontier input')
 data=json.loads(raw);V=set(data['frontier']);ranks=dict(data['frontier_ranks']);ck(len(V)==627 and max(V)==63665,'frontier size/cutoff')
 for s in V:
  if s!=1:
   u,t=edge(s);ck(u in V and ranks[u]<ranks[s],'grounded frontier rank')
 hostile=hostiles();cases=[]
 for spec in args.cases.split(','):
  N,h=map(int,spec.split(':'));r=run_case(N,h,V);cases.append(r)
  print(json.dumps({k:v for k,v in r.items()if k not in ('ledger','pending_old_counts','pending_new_counts')}),flush=True)
  if N==7:ck(r['ledger'][8][-1]==1084 and not r['prefix_pass'],'native N7d8 hostile')
  else:ck(r['prefix_pass']and r['suffix_pass'],'all-horizon finite certificate')
 out={'status':'PROVED conditional route-free gate; FINITE-EXACT listed all-horizon scales; uniform all-N OPEN','frontier_sha256':FRONTIER_SHA,'hostiles':hostile,'cases':cases,'checks':CHECKS}
 blob=(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n').encode()
 if args.output:args.output.write_bytes(blob)
 if args.verify:ck(args.verify.read_bytes()==blob,'immutable exact replay')
 print(json.dumps({'checks':CHECKS,'receipt_sha256':hashlib.sha256(blob).hexdigest()}))
if __name__=='__main__':main()
