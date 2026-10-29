#!/usr/bin/env python3
from itertools import product
from collections import Counter
from math import comb
import random,json
CHECKS=0
def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise AssertionError(msg)
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(p):
 q=[0]*len(p)
 for i,j in enumerate(p):q[j]=i
 return tuple(q)
def power(p,e):
 out=tuple(range(len(p)))
 for _ in range(e):out=mul(p,out)
 return out
def digits(x,k):
 a=[0]*k
 for i in range(k-1,-1,-1):a[i]=x%3;x//=3
 return tuple(a)
def number(a):
 v=0
 for x in a:v=3*v+x
 return v
def neg(a):return tuple((-x)%3 for x in a)
def portrait(k,P):return tuple(number(tuple((d+P.get(z[:i],0))%3 for i,d in enumerate(z))) for z in (digits(v,k) for v in range(3**k)))
def rotate(k,w):return portrait(k,{w:1})
def root(p):
 out=[-1]*len(p)
 for i in range(len(p)):
  if out[i]>=0:continue
  cyc=[i];j=p[i]
  while j!=i:cyc.append(j);j=p[j]
  need(len(cyc)%2==1,'odd-order root')
  for a,v in enumerate(cyc):out[v]=cyc[(a+(len(cyc)+1)//2)%len(cyc)]
 return tuple(out)
def transport(t,s):return root(mul(t,s))
def conj(h,p):return mul(mul(h,p),inv(h))
def arc(a,b):return next(((y-x)%3==1 for x,y in zip(a,b) if x!=y),False)
def elementary(k,w):
 N=3**k;one=tuple(range(N));alpha=tuple(number(neg(digits(v,k))) for v in range(N))
 a=next(i for i,x in enumerate(w) if x);b=len(w)-a-1;q=(0,)*a+(1,)+(0,)*b
 x=rotate(k,(0,)*a);y=rotate(k,(0,)*len(w));t=mul(power(x,2),alpha);u=mul(power(y,2),alpha)
 P={}
 for i,d in enumerate(w[a+1:]):
  pref=(0,)*a+(1,)+(0,)*i
  P[pref]=d;P[neg(pref)]=(-d)%3
 h=portrait(k,P);need(mul(h,alpha)==mul(alpha,h),'endpoint constructor centralizes')
 t,u=conj(h,t),conj(h,u)
 need(mul(t,t)==one and mul(u,u)==one,'involution endpoints')
 for z in (t,u):
  for i in range(N):
   for j in range(i):need(arc(digits(i,k),digits(j,k))==arc(digits(z[j],k),digits(z[i],k)),'anti endpoint')
 steps=(transport(t,alpha),transport(u,t),transport(alpha,u))
 loop=mul(steps[2],mul(steps[1],steps[0]));desired=portrait(k,{w:1,neg(w):2})
 need(loop==desired,'individual conjugated triangle')
 return loop,tuple(reversed(tuple(inv(s) for s in steps))),steps

def main():
 rows=[];rng=random.Random(94001)
 for k in range(5):
  N=3**k;one=tuple(range(N));alpha=tuple(number(neg(digits(v,k))) for v in range(N))
  reps=[w for d in range(1,k) for w in product(range(3),repeat=d) if next((x for x in w if x),0)==1]
  reps.sort(key=lambda w:(-len(w),w));c=(N-2*k-1)//4
  need(len(reps)==c,'free portrait count')
  elementary_data={w:elementary(k,w) for w in reps}
  words=list(product(range(3),repeat=c)) if k<=3 else [tuple(rng.randrange(3) for _ in reps) for _ in range(150)]
  costs=Counter();seen=set()
  for coeff in words:
   P={};total=one;cost=0
   for w,e in zip(reps,coeff):
    if not e:continue
    P[w]=e;P[neg(w)]=(-e)%3
    loop,reverse,forward=elementary_data[w];steps=forward if e==1 else reverse
    actual=one
    for step in steps:actual=mul(step,actual);total=mul(step,total);cost+=1
    need(actual==(loop if e==1 else inv(loop)),'inverse triangle')
   expected=portrait(k,P);need(total==expected,'descending source portrait')
   need(mul(total,alpha)==mul(alpha,total),'source centralizer retained')
   need(cost==3*sum(e!=0 for e in coeff),'exact protocol cost')
   need(cost<=3*c,'linear bound');costs[cost]+=1;seen.add(total)
  if k<=3:
   need(len(seen)==3**c,'complete centralizer normal form')
   need(dict(costs)=={3*j:comb(c,j)*2**j for j in range(c+1)},'uniform binomial law')
   need(sum(cost*n for cost,n in costs.items())==2*c*3**c,'mean length')
  rows.append({'k':k,'vertices':N,'free_pairs':c,'new_bound':3*c,'old_bound':3*(1+(k-1)*N) if k else 0,'tested_corrections':len(words),'cost_counts':dict(sorted(costs.items()))})
 print(json.dumps({'status':'PASS','checks':CHECKS,'rows':rows},indent=2))
if __name__=='__main__':main()
