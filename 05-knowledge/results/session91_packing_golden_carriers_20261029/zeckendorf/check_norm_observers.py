#!/usr/bin/env python3
"""Exact Zeckendorf observer carriers, Smith norm law, and joint-modulus gluing.
Standard library only. All checks remain active under python -O.
"""
from collections import deque
from itertools import product
from math import gcd, lcm
import json

CHECKS=0
def need(v,msg='check failed'):
 global CHECKS
 CHECKS+=1
 if not v: raise AssertionError(msg)

def shift(v,d,m):
 x,y=v
 return ((y+d)%m,(x+y+2*d)%m)

def read(word,m=None):
 x=y=p=0
 for d in word:
  if p and d:return None
  x,y=y+d,x+y+2*d
  p=d
 if m is not None:x%=m;y%=m
 return x,y,p

def canonical(n):
 f=[1,2]
 while f[-1]<=n:f.append(sum(f[-2:]))
 out=[]
 for a in reversed(f):
  out.append(int(a<=n))
  if a<=n:n-=a
 while out and not out[0]:out.pop(0)
 return tuple(out)

def period(m):
 # Columns of A^h. A is invertible over every Z/m.
 target=(1%m,0,0,1%m)
 q=target
 for h in range(1,100000):
  a,b,c,d=q;q=(c,d,(a+c)%m,(b+d)%m)
  if q==target:return h
 raise RuntimeError('period search cap')

def predicted(m,a,b,r):
 r%=m
 if r%gcd(gcd(abs(a),abs(b)),m):return 1
 g=gcd(abs(a),abs(b))
 if not g:return 3
 delta=a*a+a*b-b*b
 need(delta!=0,'nonzero integer norm')
 image=m*m//(gcd(m,g)*gcd(m,abs(delta)//g))
 return 2*image+1

def image_states(m,a,b):
 return {( (a*x+b*y)%m,(b*x+(a+b)*y)%m)
         for x in range(m) for y in range(m)}

def machine(m,a,b,r):
 im=image_states(m,a,b)
 states=[(u,v,p) for u,v in sorted(im) for p in (0,1)]+[None]
 ids={s:i for i,s in enumerate(states)};dead=ids[None]
 trans=[]
 for s in states:
  if s is None:trans.append((dead,dead));continue
  u,v,p=s
  row=[]
  for d in (0,1):
   if p and d:row.append(dead);continue
   z=((v+d*(a+2*b))%m,(u+v+d*(2*a+3*b))%m,d)
   need(z in ids,'closed observer image')
   row.append(ids[z])
  trans.append(tuple(row))
 accept=[s is not None and s[0]==r%m for s in states]
 reach={ids[(0,0,0)]};todo=deque(reach)
 while todo:
  for t in trans[todo.popleft()]:
   if t not in reach:reach.add(t);todo.append(t)
 need(len(reach)==len(states),'all quotient states reachable')
 return states,trans,accept

def minimize(trans,accept):
 p=list(map(int,accept))
 while True:
  types={};q=[]
  for i,row in enumerate(trans):
   key=(p[i],p[row[0]],p[row[1]])
   if key not in types:types[key]=len(types)
   q.append(types[key])
  if len(set(q))==len(set(p)):return len(set(q))
  p=q

def direct_digits():
 tested=0
 for k in range(12):
  f=[1,2]
  while len(f)<k+1:f.append(sum(f[-2:]))
  for w in product((0,1),repeat=k):
   got=read(w)
   legal=all(not (w[i] and w[i+1]) for i in range(k-1))
   need((got is not None)==legal,'legality')
   if legal:
    x=sum(d*f[i] for i,d in enumerate(reversed(w)))
    y=sum(d*f[i+1] for i,d in enumerate(reversed(w)))
    need(got[:2]==(x,y),'direct Fibonacci values')
    need(canonical(x)==tuple(w[next((i for i,d in enumerate(w) if d),len(w)):]),'canonical word')
   tested+=1
 return tested

def translations():
 rows=[]
 for m in range(1,41):
  h=period(m);L=h if m>1 else 3
  need(L>=3,'separated unit blocks')
  w1=(1,)+(0,)*(L-1);w2=(0,1)+(0,)*(L-2)
  for x,y in product(range(m),repeat=2):
   z=(x,y)
   for d in w1:z=shift(z,d,m)
   need(z==((x+1)%m,(y+1)%m),'translation u')
   z=(x,y)
   for d in w2:z=shift(z,d,m)
   need(z==(x,(y+1)%m),'translation v')
  rows.append((m,h))
 return rows

def all_observers():
 cases=0
 for m in range(1,9):
  for a,b,r in product(range(m),repeat=3):
   s,t,accept=machine(m,a,b,r)
   got=minimize(t,accept);want=predicted(m,a,b,r)
   need(got==want,f'minimality {(m,a,b,r,got,want)}')
   cases+=1
 return cases

def special_cases():
 out=[]
 for m,a,b,r in [(11,1,0,0),(11,1,4,0),(121,1,4,0),(121,1,37,0),
                  (5,1,3,0),(25,1,3,0),(8,2,2,0),(8,2,2,1)]:
  s,t,accept=machine(m,a,b,r)
  got=minimize(t,accept);need(got==predicted(m,a,b,r),'special minimality')
  out.append(dict(m=m,a=a,b=b,target=r,norm=a*a+a*b-b*b,
                  reachable_states=len(s),minimal_states=got))
 need((37*37-37-1)%121==0,'Hensel golden root')
 need((4*4-4-1)%121==11,'unlifted root defect')
 # Integer value differs within one golden-residue/legality state.
 seen={};hostile=None
 for n in range(2000):
  x,y,p=read(canonical(n));q=(5*x-2*y)%11
  key=(q,p)
  if key in seen and seen[key][0]%11!=n%11:
   hostile=(seen[key],(n,canonical(n),x,y,p));break
  seen.setdefault(key,(n,canonical(n),x,y,p))
 need(hostile is not None,'golden projection loses value residue')
 return out,hostile

def joint_moduli():
 cases=0
 for m in range(1,13):
  for n in range(1,13):
   L=lcm(m,n);g=gcd(m,n)
   actual={((x%m,y%m),(x%n,y%n)) for x in range(L) for y in range(L)}
   need(len(actual)==L*L,'CRT image')
   compatible=set()
   for x,y,u,v in product(range(m),range(m),range(n),range(n)):
    if (x-u)%g==0 and (y-v)%g==0:compatible.add(((x,y),(u,v)))
   need(actual==compatible,'joint state exact compatibility')
   need((2*m*m)*(2*n*n)==2*g*g*(2*L*L),'redundancy factor')
   cases+=1
 return cases

def prime_precision():
 rows=[]
 for p,ks,kind in [(2,4,'inert'),(3,3,'inert'),(7,3,'inert'),(5,3,'ramified'),(11,2,'split')]:
  for k in range(1,ks+1):
   m=p**k;minimum=m*m;hist={}
   for a,b in product(range(m),repeat=2):
    if a%p==0 and b%p==0:continue
    D=a*a+a*b-b*b
    J=m*m//gcd(m,abs(D))
    minimum=min(minimum,J);hist[J]=hist.get(J,0)+1
    if kind=='inert':need(J==m*m,'inert norm unit')
    if kind=='ramified':need(gcd(m,abs(D)) in (1,5),'ramified valuation one')
   expected=m*m if kind=='inert' else (m if kind=='split' else m*m//5)
   need(minimum==expected,'prime precision minimum')
   rows.append(dict(p=p,k=k,kind=kind,minimum_register_classes=minimum))
 root=4;mod=11;lifts=[]
 for k in range(1,9):
  need((root*root-root-1)%mod==0,'Hensel lift equation')
  need(root%11==4,'same root branch')
  lifts.append((k,mod,root))
  t=(-((root*root-root-1)//mod)*pow(2*root-1,-1,11))%11
  root+=mod*t;mod*=11
 return rows,lifts

def main():
 digits=direct_digits(); periods=translations(); cases=all_observers()
 specials,hostile=special_cases();joint=joint_moduli();precision,lifts=prime_precision()
 # Norm is the anti-isometry form of the golden clock.
 for a,b in product(range(-30,31),repeat=2):
  D=a*a+a*b-b*b;u,v=b,a+b
  need(u*u+u*v-v*v==-D,'golden norm sign flip')
 print(json.dumps(dict(status='PASS',checks=CHECKS,digit_words=digits,
      all_observer_cases=cases,translation_moduli=periods,joint_moduli_cases=joint,
      special_cases=specials,golden_value_hostile=hostile,prime_precision=precision,Hensel_lifts=lifts),sort_keys=True,indent=2))

if __name__=='__main__':main()
