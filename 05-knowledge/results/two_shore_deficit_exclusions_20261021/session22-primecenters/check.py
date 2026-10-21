#!/usr/bin/env python3
"""Tiny additive/connectivity controls; no large-order center census."""
import itertools,json,pathlib,hashlib,time,datetime
ROOT=pathlib.Path(__file__).resolve().parent;checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def strong(p,D,S):
 if not S:return False
 def reach(sign):
  got={next(iter(S))};todo=list(got)
  while todo:
   x=todo.pop()
   for y in S-got:
    if (sign*(y-x))%p in D:got.add(y);todo.append(y)
  return got==S
 return reach(1) and reach(-1)
def mincut(p,D):
 V=set(range(p))
 for k in range(p-1):
  for S in itertools.combinations(range(p),k):
   if not strong(p,D,V-set(S)):return k,list(S)
 raise RuntimeError('no nonstrong cut')
t=time.monotonic();start=datetime.datetime.now(datetime.timezone.utc).isoformat();cdcounts={};rows=[]
for p in (3,5,7):
 subsets=[{i for i in range(p) if mask>>i&1} for mask in range(1,1<<p)]
 for A in subsets:
  for B in subsets:need(len({(a+b)%p for a in A for b in B})>=min(p,len(A)+len(B)-1),'prime additive inequality')
 cdcounts[str(p)]=len(subsets)**2
 for mask in range(1<<((p-1)//2)):
  D={i if mask>>(i-1)&1 else p-i for i in range(1,(p+1)//2)}
  k,S=mincut(p,D);need(k==(p-1)//2,'exact prime circulant connectivity')
  rows.append({'p':p,'connection_set':sorted(D),'minimum_cut_size':k,'cut':S})
p=9;D={1,3,4,7};S={1,4,7};X={2,5,8};Y={0,3,6};sumset={(y+d)%p for y in Y for d in D|{0}}
need(not strong(p,D,set(range(p))-S),'composite separator')
need(all((y-x)%p in D for x in X for y in Y),'literal X to Y split')
need(sumset==Y|S and len(sumset)==6<min(9,len(Y)+len(D)),'composite CD hostile')
k,cut=mincut(p,D);need(k==3<4,'composite connectivity strictly below degree')
old=set(range(11,48,2))|{12,18,20,24,28};excluded={37,41,43,47};new=old-excluded
need(len(new)==20 and max(new)==45,'updated support')
need(45+12+1==58,'safe recurrence onset')
R={'status':'FINITE-EXACT tiny controls, additive theorem cited for all primes','start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.monotonic()-t,'checks':checks,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'prime_subset_pair_counts':cdcounts,'prime_connectivity_controls':rows,'composite_hostile':{'n':9,'D':sorted(D),'S':sorted(S),'X':sorted(X),'Y':sorted(Y),'sumset':sorted(sumset),'minimum_cut_size':k},'excluded_prime_orders':sorted(excluded),'unknown_support':sorted(new),'degree_bound':45,'safe_order18_recurrence_onset':58}
(ROOT/'receipt.json').write_text(json.dumps(R,indent=2)+'\n');print(json.dumps(R,indent=2))
