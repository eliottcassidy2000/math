#!/usr/bin/env python3
"""Exact streaming replay with a single implicit cyclic sequence; Python stdlib only.

Insertion introduces just the next vertex. A2opt operation reverses one contiguous
subsequence, so it preserves one cyclic ordering and the vertex set. Only its two
boundary edges change. The tracked set of nonsquare edges is exact by induction.
Empty bad-edge set at every recorded length certifies every intermediate cycle.
"""
import argparse,gzip,hashlib,json,random,struct,time
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--trace',default='full_incremental_trace.bin.gz');p.add_argument('--initial',default='initial_cycle.txt');p.add_argument('--audit-through',type=int,default=6500);p.add_argument('--last',type=int);args=p.parse_args()
checks=0
def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def square(x):return isqrt(x)**2==x
def protected(x,y):return(x==3362 and y in(482,6242))or(y==3362 and x in(482,6242))
def edge(x,y):return(min(x,y),max(x,y))
class Sequence:
 def __init__(self,n):
  self.L=[0]*(n+1);self.R=[0]*(n+1);self.S=[1]*(n+1);self.S[0]=0;self.F=bytearray(n+1)
  rng=random.Random(991728351);self.P=[rng.getrandbits(64) for _ in range(n+1)];self.root=0
 def push(self,t):
  if self.F[t]:
   self.L[t],self.R[t]=self.R[t],self.L[t];self.F[self.L[t]]^=1;self.F[self.R[t]]^=1;self.F[t]=0
 def merge(self,a,b):
  if not a or not b:return a or b
  if (self.P[a],a)<(self.P[b],b):
   self.push(a);self.R[a]=self.merge(self.R[a],b);self.S[a]=1+self.S[self.L[a]]+self.S[self.R[a]];return a
  self.push(b);self.L[b]=self.merge(a,self.L[b]);self.S[b]=1+self.S[self.L[b]]+self.S[self.R[b]];return b
 def split(self,t,k):
  if not t:return 0,0
  self.push(t);left=self.L[t]
  if self.S[left]>=k:
   a,b=self.split(left,k);self.L[t]=b;self.S[t]=1+self.S[b]+self.S[self.R[t]];return a,t
  a,b=self.split(self.R[t],k-self.S[left]-1);self.R[t]=a;self.S[t]=1+self.S[left]+self.S[a];return t,b
 def at(self,k):
  t=self.root
  while t:
   self.push(t);left=self.L[t]
   if k<self.S[left]:t=left
   elif k==self.S[left]:return t
   else:k-=self.S[left]+1;t=self.R[t]
  raise RuntimeError('sequence index not found')
 def insert(self,k,v):
  a,b=self.split(self.root,k);self.root=self.merge(self.merge(a,v),b)
 def reverse(self,l,r):
  a,b=self.split(self.root,l);b,c=self.split(b,r-l);self.F[b]^=1;self.root=self.merge(self.merge(a,b),c)
 def values(self):
  result=[];stack=[];t=self.root
  while t or stack:
   while t:self.push(t);stack.append(t);t=self.L[t]
   t=stack.pop();result.append(t);t=self.R[t]
  return result
raw=(ROOT/args.initial).read_bytes();initial=list(map(int,raw.split()))
def validate(order,N):
 check(len(order)==N and sorted(order)==list(range(1,N+1)),'full interval coverage')
 seed=0
 for x,y in zip(order,order[1:]+order[:1]):check(square(x+y),'full square edge');seed+=protected(x,y)
 check(seed==2,'full protected seed')
file=ROOT/args.trace;f=gzip.open(file,'rb') if file.suffix=='.gz' else file.open('rb')
check(f.read(8)==b'SSINC1\0\0','trace version');first,last=struct.unpack('<II',f.read(8));check(first==len(initial)==6242 and first<=last,'declared range')
if args.last:last=min(last,args.last)
seq=Sequence(last);validate(initial,first)
for v in initial:seq.root=seq.merge(seq.root,v)
check(seq.values()==initial,'initial sequence reconstruction')
bad=set();moves=0;maxmoves=0;direct=0;audits=0;start=time.monotonic();checkpoints=[]
def replace_edges(old,new):
 for x,y in old:
  e=edge(x,y);check(not protected(x,y),'protected edge removal')
  if not square(x+y):check(e in bad,'bad-edge deletion invariant');bad.remove(e)
  else:check(e not in bad,'good-edge deletion invariant')
 for x,y in new:
  if not square(x+y):bad.add(edge(x,y))
for N in range(first+1,last+1):
 header=f.read(6);check(len(header)==6,'truncated step header');idx,count=struct.unpack('<IH',header);check(idx<N-1,'insertion index')
 a,b=seq.at(idx),seq.at((idx+1)%(N-1));replace_edges([(a,b)],[(a,N),(N,b)]);seq.insert(idx+1,N)
 check(seq.S[seq.root]==N,'new interval size');check(len(bad)<=1,'one-defect insertion');direct+=not bad
 for _ in range(count):
  chunk=f.read(8);check(len(chunk)==8,'truncated reversal');i,j=struct.unpack('<II',chunk);check(0<=i<j<N,'reversal indices')
  a,b,c,d=seq.at(i),seq.at((i+1)%N),seq.at(j),seq.at((j+1)%N)
  replace_edges([(a,b),(c,d)],[(a,c),(b,d)]);seq.reverse(i+1,j+1);check(len(bad)<=1,'one-defect repair invariant');moves+=1
 check(not bad,'every intermediate length closes');maxmoves=max(maxmoves,count)
 if N<=args.audit_through or N%8192==0:
  order=seq.values();validate(order,N);audits+=1
  if N%8192==0:
   canonical=('\n'.join(map(str,order))+'\n').encode();checkpoints.append({'N':N,'cycle_sha256':hashlib.sha256(canonical).hexdigest()})
 if N%10000==0:print(json.dumps({'verified_through':N,'checks':checks,'seconds':round(time.monotonic()-start,3)}),file=__import__('sys').stderr,flush=True)
if args.last is None:check(not f.read(1),'no trailing bytes')
order=seq.values();validate(order,last);encoded=('\n'.join(map(str,order))+'\n').encode()
if args.last is None:
 final_file=ROOT/'final_cycle.txt.gz'
 check(final_file.exists(),'stored final cycle present')
 check(encoded==gzip.decompress(final_file.read_bytes()),'independent final sequence agreement')
print(json.dumps({'status':'passed','first':first,'last':last,'cycles':last-first+1,'checks':checks,'moves':moves,'direct_insertions':direct,'max_moves':maxmoves,'full_intermediate_audits':audits,'source_sha256':hashlib.sha256(raw).hexdigest(),'trace_file_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'final_cycle_sha256':hashlib.sha256(encoded).hexdigest(),'checkpoint_cycles':checkpoints,'scope':'Every length is certified by exact permutation, single-cycle, protected-edge, and bad-edge induction; full audits additionally cover the pilot, checkpoints, and final cycle.'},indent=2,sort_keys=True))
