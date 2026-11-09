#!/usr/bin/env python3
import importlib.util,json
from pathlib import Path
spec=importlib.util.spec_from_file_location('I',Path(__file__).with_name('inherited.py'))
I=importlib.util.module_from_spec(spec);spec.loader.exec_module(I)
def induced(a,verts):
 return [sum(((a[v]>>w)&1)<<j for j,w in enumerate(verts)) for v in verts]
def reach(a,v):
 seen=todo=1<<v
 while todo:
  bit=todo&-todo;todo-=bit;i=bit.bit_length()-1
  nxt=a[i]&~seen;seen|=nxt;todo|=nxt
 return seen

def check(a):
 h,e,u,d=I.profile(a,0);delta=u-h;r=a[0].bit_count()
 I.require(delta>=r*(r-1),'degree bound')
 if r==1:
  w=a[0].bit_length()-1
  B=induced(a,list(range(1,len(a))));w-=1
  dd,q,R,unused=I.profile(B,w)
  I.require((e,h,u,delta)==(dd-q,dd-q+R,dd+R,q),'sole exit identities')
  if delta==0:
   reaches=[reach(B,j) for j in range(len(B))]
   same=[j for j in range(len(B)) if reaches[w]>>j&1 and reaches[j]>>w&1]
   P=[j for j in range(len(B)) if j in same or B[j]>>w&1]
   C=[j for j in range(len(B)) if j not in P]
   I.require(bool(C),'nonempty final suffix')
   p,ep,z,dp=I.profile(induced(B,P),P.index(w))
   cc=induced(B,C);c=I.profile(cc,0)[0];Q=I.stats(cc,0)['Q']
   I.require((h,e,u)==(p*c+z*Q,p*c,p*c+z*Q),'zero-defect physical formula')
   I.require(3*e<=h,'zero face')
   I.require((3*e==h)==(len(C)==1 and z==p),'zero face equality')
   if e==1:I.require(h>1 and ((h-1)&(h-2))==0,'e1 lacunarity')
 if r:
  L=2*h+delta-e
  I.require(2*L>=3*h+1,'old floor')
  if 2*L==3*h+1:I.require((h,e,delta)==(3,1,0),'old floor equality')
  if h>=5:
   I.require(2*L>=3*h+3,'improved floor')
   eq=2*L==3*h+3
   pred=((h,e,delta) in [(5,1,0),(9,3,0)] or (delta==1 and e==(h-1)//2 and h>=5 and ((h-1)&(h-2))==0))
   I.require(eq==pred,'improved equality classification')
   if eq:I.require(L%3==0,'equality output composite')
  return h,e,u,L
 return None

def reversed_extreme(n):
 a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 a[0]^=1<<(n-1);a[-1]|=1
 return a

def reroot(a,r):
 return induced(a,[r]+[v for v in range(len(a)) if v!=r])

def main():
 total=0;zero=0;equal=set()
 for n in range(1,7):
  for a in I.tournaments(n):
   result=check(a);total+=1
   if result:
    h,e,u,L=result
    zero+=u==h
    if h>=5 and 2*L==3*h+3:equal.add((h,e,u,L))
 families=[]
 for k in range(2,10):
  a=reroot(reversed_extreme(k+2),k)
  h,e,u,L=check(a)
  I.require((h,e,u)==(2**k+1,2**(k-1),2**k+2),'literal extremal family')
  families.append((h,e,u,L))
 return dict(status='FINITE-EXACT',rooted_labelled_universe=total,zero_defect_nonsinks=zero,improved_floor_equalities=sorted(equal),literal_infinite_family_parameters=families,all_target_coverage='OPEN')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
