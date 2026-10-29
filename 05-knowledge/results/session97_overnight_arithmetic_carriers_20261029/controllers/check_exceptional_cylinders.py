#!/usr/bin/env python3
from fractions import Fraction as F
from check_integer_termination import step, require
import json

WORDS=[(1,), (3,1,1,1),(3,3,1,1,1,1,1),(3,1,3,1,1,1,1),(3,1,1,3,1,1,1)]
def chart(w):
 c=d=0
 for a in w:c=3*c+(1<<d);d+=a
 j=len(w);mod=1<<(d+1);res=((1<<d)-c)*pow(3**j,-1,mod)%mod
 return {'word':w,'J':j,'D':d,'C':c,'residue':res,'modulus':mod}
def fast_run(n):
 q=1;w=[];states=[]
 while q:
  y,a=step(n);states.append((n,q,a));w.append(a);q+=2 if a%3==0 else -1;n=y
 return n,tuple(w),states

def main():
 charts=[chart(w) for w in WORDS]
 small=[]
 for n in range(13,125,16):
  y,w,states=fast_run(n)
  require(any(a not in (1,3) for a in w),'small arithmetic obstruction')
  small.append({'source':n,'word':w})
 require(F(3,2)*F(6,7)**3<1,'three branch descent')
 require(F(3,2)*F(117,448)<1,'high branch descent')
 counts=[0]*len(WORDS)
 for n in range(1,1000002,2):
  target,w,states=fast_run(n)
  matches=[i for i,c in enumerate(charts) if n%c['modulus']==c['residue']]
  require(len(matches)<=1,'disjoint cylinders')
  require((target>n)==bool(matches),'complete expansion classifier')
  require((target==n)==(n==1),'only equality')
  if matches:
   i=matches[0];counts[i]+=1;require(w==WORDS[i],'exact word')
  if all(a in (1,3) for a in w):
   require(all(x>=125 for x,q,a in states if a==3),'enhanced branch guard')
 # Long native expansion clock, retained without arbitrary cutoff.
 for t in range(2,301):
  n=(1<<t)-1
  for k in range(t-1):
   y,a=step(n);require(a==1 and y>n,'long expanding source');n=y
 print(json.dumps({'sources':500001,'charts':charts,'small_branch_obstructions':small,
  'expansion_counts':counts,'long_run_t_max':300},sort_keys=True))
if __name__=='__main__':main()
