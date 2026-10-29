#!/usr/bin/env python3
from itertools import product
from fractions import Fraction
from math import factorial
import json
COUNT=0
def need(x,msg):
 global COUNT
 COUNT+=1
 if not x:raise AssertionError(msg)
def mul(a,b):
 u,v=a;x,y=b
 return u*x+v*y,u*y+v*x+v*y
def main():
 n=9;powers=[(1,0)]
 for _ in range(n):powers.append(mul(powers[-1],(0,1)))
 need(powers[9]==(21,34),'ideal generator')
 words=[w for w in range(1<<n) if all(not(w>>i&1 and w>>((i+1)%n)&1) for i in range(n))]
 def coord(w):
  u=sum(powers[i][0] for i in range(n) if w>>i&1);v=sum(powers[i][1] for i in range(n) if w>>i&1)
  return u%2,v%2,(u+5*v)%19
 lookup={coord(w):w for w in words}
 need(len(words)==len(lookup)==76,'all CRT states')
 need(set(lookup)==set(product(range(2),range(2),range(19))),'CRT onto')
 def mirror(c):
  w=lookup[c];r=sum(((w>>i)&1)<<((-i)%n) for i in range(n));return coord(r)
 def shift(c):return ((c[0]+1)%2,c[1],(c[2]+1)%19)
 for c in lookup:
  a,b,y=c;u,v=mul((a,b),(a,b));need(mirror(c)[:2]==(u%2,v%2),'F4 Frobenius')
  need(mirror(mirror(c))==c,'reflection involution')
 for y in range(1,19):need(len({mirror((a,b,y))[2] for a,b in product(range(2),repeat=2)})==4,'four reflected coordinates')
 profiles={(c[2],mirror(c)[2],mirror(shift(c))[2]) for c in lookup}
 need(len(profiles)==76,'complete future profile')
 c=(0,0,0);orbit=[]
 while c not in orbit:orbit.append(c);c=shift(c)
 need(len(orbit)==38,'translation exact order')
 # Literal recursive permutation groups, no moment-based construction.
 def theta(v,k):
  out=0;place=1
  for _ in range(k):out+=((-v%3)%3)*place;v//=3;place*=3
  return out
 def A(k):
  if k==0:return [(0,)]
  child=A(k-1);s=3**(k-1)
  return [tuple(((i//s+r)%3)*s+gs[i//s][i%s] for i in range(3*s)) for r in range(3) for gs in product(child,repeat=3)]
 def C(k):
  if k==0:return [(0,)]
  s=3**(k-1);return [tuple(c)+tuple(s+g[i] for i in range(s))+tuple(2*s+theta(g[theta(i,k-1)],k-1) for i in range(s)) for c in C(k-1) for g in A(k-1)]
 rows=[]
 for k in range(4):
  G=C(k);N=3**k;fix=[sum(g[i]==i for i in range(N)) for g in G]
  moments=[Fraction(sum(x**r for x in fix),len(G)) for r in range(1,4)]
  target=[2*k+1,8*k*k+1,48*k**3-40*k*k+18*k+1]
  need(moments==target,'literal centralizer moments')
  need(moments[1]-moments[0]**2==4*k*(k-1),'joint excess variance')
  rows.append({'k':k,'group_order':len(G),'moments':[int(x) for x in moments]})
 for p in range(3,101,2):
  if p==5 or any(p%d==0 for d in range(2,int(p**.5)+1)):continue
  def order(a):
   x=(1,0)
   for h in range(1,p*p+1):
    u,v=mul(x,a);x=u%p,v%p
    if x==(1,0):return h
   raise AssertionError('order not found')
  h=order((1,1));d=order((p-1,p-1));expected=2*h if h%2 else (h//2 if h%4==2 else h)
  need(d==expected,'projective versus norm-one clock')
 print(json.dumps({'status':'PASS','checks':COUNT,'word_states':76,'translation_order':38,'future_profiles':len(profiles),'centralizer_moments':rows},sort_keys=True))
if __name__=='__main__':main()
