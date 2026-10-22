#!/usr/bin/env python3
"""Exact relation-cone identities; literal source tests, no convergence claim."""
import argparse, hashlib, itertools, json
from pathlib import Path
from fractions import Fraction as F
L=((1,2),(-2,-1),(1,-1)); U=((0,1),(-1,-1),(1,0))
Z=(0,0); CHECKS=0
def ck(p,m):
 global CHECKS
 CHECKS+=1
 if not p: raise ValueError(m)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def mul(n,a): return tuple(n*x for x in a)
def q(a): x,y=a; return x*x-x*y+y*y
def beta(a,b): return q(add(a,b))-q(a)-q(b)
def delta(b):w,h=b;return q(h)+beta(w,h)
def kappa(b,c):w,h=b;v,j=c;return beta(h,j)+beta(w,j)+beta(v,h)
def blockadd(b,c):return add(b[0],c[0]),add(b[1],c[1])
def blockmul(n,b):return mul(n,b[0]),mul(n,b[1])
A=(L[0],U[1]); B=(L[2],mul(-1,L[2]))
def primitive(s):return type(s) is int and s>0 and s%2==1 and s%8!=5
def edge(s):
 if not primitive(s):raise ValueError('primitive source required')
 c=2 if s%8==1 else 1; u=(3*s+1)//(2**c)
 while u%8==5:u=(u-1)//4;c+=2
 ck(primitive(u),'primitive edge target')
 return u,c

def charge(N,s):
 v=Z
 while s<4*N:
  v=add(v,(L if s<2*N else U)[s%3]);s=4*s+1
 return v

def count(lo,hi,rs):return sum((hi-1-r)//72-(lo-1-r)//72 for r in rs)
def bank(N):
 cut=(3*N+1)//2
 for tag,lo,hi,rs,coef,offset,clock in [('A',cut,2*N,(9,27,63),16,3,4),('B',N,cut,(17,35,71),64,19,6)]:
  for r in rs:
   for u in range(lo+(r-lo)%72,hi,72):yield tag,u,(coef*u+offset)//3,clock

def route(u,V,depth=500):
 c=0
 for d in range(depth+1):
  if u in V:return u,c,d
  if d==depth:return None
  u,z=edge(u);c+=z

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--out',type=Path,default=Path('verification.json'));a=p.parse_args()
 frontierpath=a.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json'
 raw=frontierpath.read_bytes();V=set(json.loads(raw)['frontier']);ck(len(V)==627 and max(V)==63665,'pinned frontier geometry')
 # Algebraic identity and sharp real maximum; integers x,y need not realize the maximum.
 cases=0
 for x,y,aa,bb in itertools.product(range(-5,6),repeat=4):
  w=add(add(mul(x,L[0]),mul(y,L[2])),add(mul(aa,L[0]),mul(bb,L[2])))
  h=add(mul(x,U[1]),mul(-y,L[2]));f=-2*x*x+3*x*y-3*y*y-3*aa*x+3*aa*y-6*bb*y
  ck(delta((w,h))==f,'arrival boundary identity')
  ck(5*f<=6*(aa*aa-aa*bb+4*bb*bb),'boundary Schur bound')
  if aa==bb==0:ck(f<=0 and (f==0)==(x==y==0),'negative-definite completed cone')
  cases+=1
 # General outside-bank identity and Schur gate, preserving all cross terms.
 outside_cases=0
 for vals in itertools.product(range(-1,2),repeat=8):
  ox,oy,hx,hy,aa,bb,x,y=vals;o=((ox,oy),(hx,hy));bbvec=add(mul(aa,L[0]),mul(bb,L[2]));c=(add(o[0],bbvec),o[1]);d=blockadd(c,blockadd(blockmul(x,A),blockmul(y,B)))
  l1=kappa(c,A);l2=kappa(c,B);base=delta(o)+beta(bbvec,o[1]);rhs=F(3*l1*l1+3*l1*l2+2*l2*l2,15)
  ck(delta(d)==base-2*x*x+3*x*y-3*y*y+l1*x+l2*y,'outside exact identity')
  ck(delta(d)<=base+rhs,'outside source-preserving Schur gate');outside_cases+=1
 # Literal native source family across every small scale and sparse very large scales.
 rows=0;guard_hash=hashlib.sha256();summaries=[]
 for N in list(range(1,1001))+[4096,8192,32768,63666,100000]:
  addresses=set();tot={'A':0,'B':0}
  for tag,u,s,c in bank(N):
   ck(u not in addresses and s not in addresses,'injective original-source ownership');addresses.update((u,s));tot[tag]+=1
   ck(edge(s)==(u,c),'actual one-edge relation and full clock')
   ck(N<=u<2*N and 4*N<=s<256*N,'original old/new carrier')
   target=A if tag=='A' else B
   ck((charge(N,u),charge(64*N,s))==target,'literal B-chain charge, not assumed palette')
   ck(charge(N,u)==charge(64*N,u),'unchanged old source')
   guard_hash.update(json.dumps([N,tag,u,s,c,u%15,u%35,s%15,s%35],separators=(',',':')).encode());rows+=1
  cut=(3*N+1)//2
  ck(tot=={'A':count(cut,2*N,(9,27,63)),'B':count(N,cut,(17,35,71))},'exact progression count')
  if N%144==0:ck(tot['A']==tot['B']==N//48,'linear all-scale count on aligned scales')
  if N>=4096:summaries.append({'N':N,'counts':tot})
 # Literal arrival records: no assumption of eventual arrival; routes are finite checks.
 route_checks=0;unresolved=0;event_hash=hashlib.sha256()
 for N in [63666,65536]:
  for tag,u,s,c in bank(N):
   ck(u not in V and s not in V,'first-arrival seam guard')
   ru=route(u,V);rs=route(s,V)
   if ru is None or rs is None:unresolved+=1;continue
   ck(rs==(ru[0],ru[1]+c,ru[2]+1),'literal first-arrival terminal clock and depth')
   # Every finite horizon checked through last arrival; retain arbitrary sign of root clock.
   for d in range(rs[2]+2):
    sig=(-1)**ru[1];w=mul(sig,charge(N,u)) if ru[2]<=d else Z;h=mul(sig,charge(64*N,s)) if rs[2]<=d else Z
    ck(delta((w,h))<=0,'single native pair has no positive partial arrival increment')
   event_hash.update(json.dumps([N,tag,u,s,ru,rs,u%15,u%35,s%15,s%35],separators=(',',':')).encode());route_checks+=1
 # Exact large N counts do not enumerate the annulus.
 large=[]
 for N in [144*64**k for k in range(20)]+[10**50+7]:
  cut=(3*N+1)//2;ca=count(cut,2*N,(9,27,63));cb=count(N,cut,(17,35,71))
  if N%144==0:ck(ca==cb==N//48,'large exact linear supply')
  large.append([N,ca,cb])
 # A wrong independent-arrival state gives +1 or +3 and is forbidden by the relation.
 ck(delta((Z,A[1]))==1 and delta((Z,B[1]))==3,'hostile child-alone states')
 # Nonpositive completed energy cannot be added to an exterior as a free credit.
 hostile=blockadd(((-10,0),(0,0)),B)
 ck(delta(hostile)>delta(((-10,0),(0,0)))+delta(B),'external interaction hostile')
 # Actual profile cone has no independent signs for parent and child.
 ck(delta((A[0],mul(-1,A[1])))==4,'wrong-sign hostile')
 out={'status':'PROVED_SCOPED_AND_FINITE_EXACT','checks':CHECKS,'boundary_algebra_cases':cases,'outside_algebra_cases':outside_cases,'literal_family_rows':rows,'finite_arrival_pair_checks':route_checks,'unresolved_at_500':unresolved,'source_and_register_hash':guard_hash.hexdigest(),'arrival_clock_register_hash':event_hash.hexdigest(),'frontier_sha256':hashlib.sha256(raw).hexdigest(),'selected_scale_counts':summaries,'large_exact_counts':large,'open':['uniform outside-plus-boundary Schur gate','uniform 153N','universal source entry'],'reported_composed_bounds':[520738,1039055,3848613],'reported_bounds_independently_recovered':False}
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
