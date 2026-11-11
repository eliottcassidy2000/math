"""Exact focused controls. Run from this directory or by full path."""
from bifan import packet,predict,curvature
from explore import H,fans,tt,flip
from itertools import combinations,permutations
from fractions import Fraction
import json
from math import isqrt

def hp(a):
 n=len(a)
 return [p for p in permutations(range(n)) if all(a[u]>>v&1 for u,v in zip(p,p[1:]))]
def delete(a,r):
 ids=[u for u in range(len(a)) if u!=r]
 return [sum(((a[u]>>v)&1)<<j for j,v in enumerate(ids)) for u in ids]
def perm_packet(a,r):
 u=delete(a,r); ids=[v for v in range(len(a)) if v!=r];A={j for j,v in enumerate(ids) if a[r]>>v&1}; k=len(u); full=(1<<k)-1
 paths={}
 for mask in range(1<<k):
  paths[mask]=[p for p in permutations([i for i in range(k) if mask>>i&1]) if all(u[x]>>y&1 for x,y in zip(p,p[1:]))]
 q=R=L=0
 for mask,pp in paths.items():
  hother=len(paths[full^mask]); q+=len(pp)*hother
  R+=sum(bool(p) and p[0] in A for p in pp)*hother
  L+=sum(bool(p) and p[-1] in A for p in pp)*hother
 h=len(paths[full]);s=sum(bool(p) and p[0] in A for p in paths[full]);e=sum(bool(p) and p[-1] in A for p in paths[full]);F=len(hp(a))
 return(q,-L-h,R-q,F+L+2*h-s-R+e)
def inverse(p,N):
 q,b,c,d=p;Delta=q*d-b*c;M=q*N-Delta;out=[]
 if M<=0:return out
 def power(v):return v>=2 and v&(v-1)==0
 for x in range(1,isqrt(M)+1):
  if M%x:continue
  for f in {x,M//x}:
   g=M//f
   if (f-c)%q or (g-b)%q:continue
   xx=(f-c)//q;yy=(g-b)//q
   if power(xx) and power(yy):out.append((xx.bit_length()-1,yy.bit_length()-1))
 return sorted(set(out))

def main():
 trials=0;floorcount=0
 for n in range(1,6):
  ed=list(combinations(range(n),2))
  for code in range(1<<len(ed)):
   a=tt(n)
   for bit,(i,j) in enumerate(ed):
    if code>>bit&1:a=flip(a,i,j)
   for r in range(n):
    p=packet(a,r);h=H(delete(a,r)) if n>1 else 1
    assert curvature(p)>=2*h;floorcount+=1
    if n>4:continue
    assert p==perm_packet(a,r)
    for m in (1,2):
     for nn in (1,2):
      f=fans(a,r,m,nn);assert H(f)==predict(p,m,nn);trials+=1
    q,b,c,d=p
    hh=[[H(fans(a,r,m,nn)) for nn in (1,2)] for m in (1,2)]
    mixed=hh[1][1]-hh[1][0]-hh[0][1]+hh[0][0]
    assert mixed==4*q
    assert Fraction(hh[0][0]*hh[1][1]-hh[0][1]*hh[1][0],mixed)==curvature(p)
    assert packet(fans(a,r,1,1),r)==(4*q,2*b,2*c,d)
 T=[96,97,59,115,35,0,52]
 p0=(168,-104,-137,109);p4=(168,-134,-107,107)
 assert len(hp(T))==19
 assert packet(T,0)==perm_packet(T,0)==p0
 assert packet(T,4)==perm_packet(T,4)==p4
 pairs=[]
 for m in (1,2,3):
  f0=fans(T,0,m,m);f4=fans(T,4,m,m)
  a,b=H(f0),H(f4);assert a-b==2
  assert(a,b)==(predict(p0,m,m),predict(p4,m,m));pairs.append([b,a])
  # Full cut switch at {0,4} followed by ten old-edge repairs.
  g=f0[:]
  for i in (0,4):
   for j in range(len(g)):
    if j not in (0,4):g=flip(g,i,j)
  for i in (0,4):
   for j in range(7):
    if j not in (0,4):g=flip(g,i,j)
  assert g==f4
 twins=[([58,36,17,6,10,28],1),([38,16,34,7,45,10],5)]
 sequences=[];two=[]
 for a,r in twins:
  assert packet(a,r)==perm_packet(a,r)
  sequences.append([H(fans(a,r,0,n)) for n in range(4)])
  two.append(H(fans(a,r,1,1)))
 assert sequences[0]==sequences[1]==[37,73,141,277]
 assert two==[215,211]
 hostile=[44,25,34,36,45,2];p=packet(hostile,2)
 assert p==perm_packet(hostile,2)==(84,-33,-60,41)
 assert H(delete(hostile,2))==9 and curvature(p)==Fraction(122,7)
 ip=packet([2,4,1],0)
 for N in range(1,500):
  expected=[(m,n) for m in range(1,8) for n in range(1,8) if predict(ip,m,n)==N]
  assert inverse(ip,N)==expected
 portal63=[16,5,17,7,74,31,175,63];portal61=flip(portal63,3,4)
 assert H(portal61)==61 and H(portal63)==63
 assert packet(portal61,4)==perm_packet(portal61,4)==(128,-18,-60,11)
 assert packet(portal63,4)==perm_packet(portal63,4)==(128,-26,-52,13)
 common=next(p for p in hp(portal61) if all(portal63[u]>>v&1 for u,v in zip(p,p[1:])))
 assert frozenset((3,4)) not in [frozenset(e) for e in zip(common,common[1:])]
 portal=[]
 for m in (1,2,3):
  aa=fans(portal61,4,m,m);bb=fans(portal63,4,m,m)
  assert flip(aa,3,4)==bb
  ha=H(aa);hb=H(bb);assert hb-ha==2
  assert ha==128*4**m-78*2**m+11;portal.append([ha,hb])
  frame=list(common)
  for v in range(8,8+m):frame.insert(frame.index(4)+1,v)
  for v in range(8+m,8+2*m):frame.insert(frame.index(4),v)
  assert all(aa[u]>>v&1 and bb[u]>>v&1 for u,v in zip(frame,frame[1:]))
 transitive_cube=[]
 for mask in range(8):
  a=[2,4,1,mask]
  for v in range(3):
   if not(mask>>v&1):a[v]|=1<<3
  q,b,c,d=packet(a,3)
  assert q==12 and b+c==-15 and d==H(a)+6
  child=H(fans(a,3,1,1));assert child==H(a)+24
  transitive_cube.append([H(a),child])
 print(json.dumps({'vertex_transitive_cube':transitive_cube,'portal_common_frame':common,'portal_lifts':portal,'small_bifan_checks' :trials,'finite_floor_rooted_inputs':floorcount,'renewable_pairs':pairs,'hidden_single_fan':sequences[0],'opposite_seam_outputs':two,'false_floor_curvature':'122/7','inverse_targets_checked':499},indent=2))
if __name__=='__main__':main()
