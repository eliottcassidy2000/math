#!/usr/bin/env python3
"""Golden cyclic reader with exact operation-sensitive fibre tests."""
import json
from collections import defaultdict
COUNT=0
def need(p,msg):
 global COUNT
 COUNT+=1
 if not p:raise AssertionError(msg)
def mul(x,y):
 a,b=x;c,d=y
 return a*c+b*d,a*d+b*c+b*d
def legal(w,n):
 return not any((w>>i)&1 and (w>>((i+1)%n))&1 for i in range(n))
def toggle(w,i,n):
 v=w^(1<<i)
 return v if legal(v,n) else w
def words(n):
 def rec(i,w,prev):
  if i==n:
   if legal(w,n):yield w
   return
  yield from rec(i+1,w,0)
  if not prev:yield from rec(i+1,w|(1<<i),1)
 yield from rec(0,0,0)
def main():
 powers=[(1,0)]
 for _ in range(20):powers.append(mul(powers[-1],(0,1)))
 L=[2,1]
 for _ in range(2,21):L.append(L[-1]+L[-2])
 rows=[]
 for n in range(1,21):
  a,b=powers[n];a-=1;det=a*a+a*b-b*b;d=abs(det)
  def residue(w):
   x=sum(powers[i][0] for i in range(n) if (w>>i)&1)
   y=sum(powers[i][1] for i in range(n) if (w>>i)&1)
   return (((a+b)*x-b*y)%d,(-b*x+a*y)%d)
  W=list(words(n));q={w:residue(w) for w in W};fibres=defaultdict(list)
  for w in W:fibres[q[w]].append(w)
  need(len(W)==L[n],'Lucas word count')
  need(d==L[n]-1-(-1)**n,'norm count')
  need(len(fibres)==d,'reader onto')
  if n%2:
   need(all(len(f)==1 for f in fibres.values()),'odd bijection')
  else:
   alternating=[sum(1<<i for i in range(j,n,2)) for j in (0,1)]
   need(set(fibres[(0,0)])=={0,*alternating},'three zero sources')
   need(all(len(f)==1 for k,f in fibres.items() if k!=(0,0)),'other singleton fibres')
  profiles=set();fixedprofiles=[set() for _ in range(n)]
  for w in W:
   outcomes=[]
   for i in range(n):
    v=toggle(w,i,n);need(v in q,'legal toggle');need(toggle(v,i,n)==w,'toggle involution')
    outcomes.append(q[v]);fixedprofiles[i].add((q[w],q[v]))
   profiles.add((q[w],tuple(outcomes)))
  need(len(profiles)==(1 if n==2 else len(W)),'one step minimal refinement')
  for profile in fixedprofiles:need(len(profile)==(1 if n==2 else len(W)),'single fixed toggle reconstruction')
  zero_counts=None
  if n%2==0:
   zero_counts=[sum(q[toggle(w,i,n)]==(0,0) for i in range(n)) for w in [0,*alternating]]
   need(zero_counts==([2,2,2] if n==2 else [0,n//2,n//2]),'strong lumpability hostile')
  rows.append(dict(n=n,words=len(W),ring_size=d,one_step_profiles=len(profiles),zero_transition_numerators=zero_counts))
 print(json.dumps(dict(status='PASS',checks=COUNT,rows=rows),indent=2))
if __name__=='__main__':main()
