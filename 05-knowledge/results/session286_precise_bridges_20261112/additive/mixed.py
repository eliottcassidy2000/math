#!/usr/bin/env python3
"""Exact all-word central-residue supply control for two literal deletion cores.
Finite closure through length4 plus an all-length algebraic content bound.
No random sampling and no program-order filtering after state deduplication.
"""
from functools import lru_cache
from collections import defaultdict
from math import prod
import json

def check(p,s):
 if not p:raise RuntimeError(s)
def reverse(T):return tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T)))
def paths(T):
 n=len(T);inc=reverse(T)
 @lru_cache(None)
 def end(S,v):
  R=S^(1<<v)
  if not R:return 1
  return sum(end(R,j) for j in range(n) if (R&inc[v])>>j&1)
 return [1]+[sum(end(S,j) for j in range(n) if S>>j&1) for S in range(1,1<<n)],end
# Tuple order H,h,Q,U,V,d,e. For all nonempty deletions in this language K=0.
def packet(T):
 n=len(T);F=(1<<n)-1;D=F^1;H,E=paths(T);RH,RE=paths(reverse(T));U=V=d=e=0
 for S in range(1,F+1,2):
  sgn=1 if S.bit_count()%2 else -1
  u=RE(S,0)*RH[F^S];v=E(S,0)*H[F^S]
  U+=u;V+=v;d+=sgn*u;e+=sgn*v
 Q=sum(H[S]*H[D^S] for S in range(0,F+1,2))
 K=sum((-1)**S.bit_count()*H[S]*H[D^S] for S in range(0,F+1,2))
 check(K==0,'alternating deletion mass')
 return (H[F],H[D],Q,U,V,d,e)
def attach(D,c):return (c<<1,)+tuple((x<<1)|(1-((c>>i)&1)) for i,x in enumerate(D))
def fuse(a,b):
 H,h,Q,U,V,d,e=a;J,j,R,X,Y,f,g=b
 return (h*J+(H-h)*j+(U-h)*(Y-j)+(d-h)*(g-j),h*j,Q*R,h*X+(U-h)*R,Q*Y+(V-Q)*j,h*f,(e)*j)
def main():
 cores=[(18,8,3,5,14),(24,17,3,6,12)]
 A=sorted({packet(attach(D,c)) for D in cores for c in range(32)})
 check(len(A)==40,'full40 distinct packet alphabet')
 check({(p[1],p[2]) for p in A}=={(13,104),(15,120)},'common ratio8')
 states=set(A);rows=[]
 expected=[{13:1,15:2},{169:2,195:2,225:3},{2197:4,2535:4,2925:3,3375:4},{28561:8,32955:8,38025:6,43875:4,50625:5}]
 for k in range(1,5):
  if k>1:states={fuse(p,a) for p in states for a in A}
  groups=defaultdict(set)
  for p in states:groups[p[1:]].add(p[0])
  maxima={}
  for key,E in groups.items():
   h=key[0];base=min(E);r=len({((x-base)//2)%h for x in E})
   maxima[h]=max(maxima.get(h,0),r)
  check(maxima==expected[k-1],'finite exact residue maxima')
  check(all(v<h for h,v in maxima.items()),'no complete residue fiber')
  rows.append(dict(length=k,distinct_packets=len(states),central_fibers=len(groups),max_residues_by_deletion=maxima))
 # General formula content: valuations discard at most two occurrence factors.
 # Direct long-word controls with mixed cores, independent of finite DP pruning.
 content_controls=0
 for k in range(5,31):
  for a in range(k+1):
   types=[13]*a+[15]*(k-a)
   types=types[a//2:]+types[:a//2]
   pools={h:[p for p in A if p[1]==h] for h in (13,15)}
   word=[pools[h][(7*i+3*a+k)%len(pools[h])] for i,h in enumerate(types)]
   p=word[0]
   for z in word[1:]:p=fuse(p,z)
   g=13**max(a-2,0)*15**max(k-a-2,0)
   check(g>1 and p[0]%g==0,'all-length content control');content_controls+=1
 print(json.dumps(dict(status='PASS',literal_cores=cores,distinct_alphabet=40,finite_rows=rows,all_length_content_controls=content_controls,scope='All words in two named32-cut alphabets; no saturation theorem outside this language'),sort_keys=True,indent=2))
if __name__=='__main__':main()
