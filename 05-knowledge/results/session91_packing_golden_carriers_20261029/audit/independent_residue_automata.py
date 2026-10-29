#!/usr/bin/env python3
"""Independent reachable automata and Moore minimization for Fibonacci observers."""
from collections import deque
from math import gcd
import json
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def automaton(m,witness=False):
 states=[(x,y,p) for x in range(m) for y in range(m) for p in (0,1)]+[None];idx={s:i for i,s in enumerate(states)};dead=len(states)-1
 tr=[]
 for s in states:
  if s is None:tr.append((dead,dead));continue
  x,y,p=s
  tr.append((idx[y%m,(x+y)%m,0],dead if p else idx[(y+1)%m,(x+y+2)%m,1]))
 seen={idx[0,0,0]:''};todo=deque(seen)
 while todo:
  s=todo.popleft()
  for d,t in enumerate(tr[s]):
   if t not in seen:seen[t]=(seen[s]+str(d)) if witness else '';todo.append(t)
 ck(len(seen)==len(states),'all states reachable')
 if witness:
  for i,w in seen.items():
   if i==dead:ck('11' in w,'dead witness');continue
   f=[1,2]
   while len(f)<=len(w):f.append(f[-1]+f[-2])
   x=sum(int(d)*f[j] for j,d in enumerate(w[::-1]))%m
   y=sum(int(d)*f[j+1] for j,d in enumerate(w[::-1]))%m
   ck(states[i]==(x,y,int(w[-1]) if w else 0),'direct weighted witness')
 return states,tr

def minimize(states,tr,m,a,b,r):
 col=[int(s is not None and (a*s[0]+b*s[1]-r)%m==0) for s in states]
 while True:
  labels={};new=[]
  for i,t in enumerate(tr):
   key=(col[i],col[t[0]],col[t[1]])
   if key not in labels:labels[key]=len(labels)
   new.append(labels[key])
  if len(set(col))==len(labels):return len(labels)
  col=new

def expected(m,a,b,r):
 if r%gcd(gcd(a,b),m):return 1
 if a==b==0:return 3
 g=gcd(a,b);D=abs(a*a+a*b-b*b)
 return 2*m*m//(gcd(m,g)*gcd(m,D//g))+1

def main():
 examples=[];cases=0
 for m in range(1,10):
  st,tr=automaton(m,True)
  for a in range(m):
   for b in range(m):
    for r in range(m):
     got=minimize(st,tr,m,a,b,r);want=expected(m,a,b,r)
     ck(got==want,f'observer size m{m} a{a} b{b} r{r}: {got}!={want}');cases+=1
 for m,a,b in [(11,1,0),(11,1,4),(121,1,4),(121,1,37)]:
  st,tr=automaton(m);got=minimize(st,tr,m,a,b,0)
  ck(got==expected(m,a,b,0),'golden exceptional modulus')
  examples.append(dict(modulus=m,observer=[a,b],states=got))
 print(json.dumps(dict(status='PASS',checks=N,residue_cases=cases,examples=examples),indent=2))
if __name__=='__main__':main()
