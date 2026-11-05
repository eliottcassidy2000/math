#!/usr/bin/env python3
"""Independent actual odd-cycle profiles and literal crossing-cycle words."""
from itertools import combinations, permutations
from functools import lru_cache
from math import comb
import json, random

def need(v,s):
 if not v:raise RuntimeError(s)
def trim(p):
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def add(*ps):
 q=[0]*max(map(len,ps))
 for p in ps:
  for i,a in enumerate(p):q[i]+=a
 return trim(q)
def neg(p):return [-a for a in p]
def mul(p,q):
 r=[0]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):r[i+j]+=a*b
 return trim(r)
def mon(n):return [0]*n+[1]
def power(p,n):
 r=[1]
 for _ in range(n):r=mul(r,p)
 return r
def shift(p,t):
 out=[0]*len(p)
 for i,a in enumerate(p):
  for j in range(i+1):out[j]+=a*comb(i,j)*t**(i-j)
 return trim(out)
def val(p,x):return sum(a*x**j for j,a in enumerate(p))
def opp(a):return tuple(row^((1<<len(a))-1)^(1<<i) for i,row in enumerate(a))
def tr(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def tour(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else (j,i);a[u]|=1<<v
 return tuple(a)
def induced(a,mask):
 vs=[v for v in range(len(a)) if mask>>v&1]
 return tuple(sum(((a[u]>>v)&1)<<i for i,v in enumerate(vs)) for u in vs)
def mirror(s,r,t):
 n=len(s);m=len(t);a=[0]*(2*n+m)
 for i in range(n):a[i]=s[i];a[n+i]=opp(s)[i]<<n
 for i in range(n):
  for j in range(n):
   if i==r and j==r:a[i]|=1<<(n+j)
   else:a[n+j]|=1<<i
 for u,row in enumerate(t):
  a[2*n+u]=row<<(2*n)|(1<<n)-1
  for i in range(n):a[n+i]|=1<<(2*n+u)
 need(all((((a[i]>>j)&1)+((a[j]>>i)&1))==1 for i,j in combinations(range(len(a)),2)),'actual tournament')
 return tuple(a)
@lru_cache(None)
def bins(a):
 n=len(a);M=1<<n;cy=[0]*M
 # A directed cycle is anchored at its least vertex. Count its directed words by DP.
 for r in range(n):
  dp={(1<<r,r):1}
  for mask in range(1<<r,M):
   if mask&((1<<r)-1) or not mask>>r&1:continue
   for u in range(r,n):
    c=dp.get((mask,u),0)
    if not c:continue
    if mask.bit_count()>=3 and mask.bit_count()%2 and a[u]>>r&1:cy[mask]+=c
    avail=a[u]&~mask&~((1<<(r+1))-1)
    while avail:
     b=avail&-avail;avail-=b;v=b.bit_length()-1
     dp[(mask|b,v)]=dp.get((mask|b,v),0)+c
 @lru_cache(None)
 def pack(mask):
  if not mask:return [1]
  b=mask&-mask;out=mul([0,1],pack(mask^b));sub=mask
  while sub:
   if sub&b and cy[sub]:out=add(out,mul([0,2*cy[sub]],pack(mask^sub)))
   sub=(sub-1)&mask
  return out
 return pack(M-1)
def paths(a,end=None):
 out=[]
 def walk(w,mask):
  if end is None or w[-1]==end:out.append(w)
  todo=a[w[-1]]&~mask
  while todo:
   b=todo&-todo;todo-=b;walk(w+(b.bit_length()-1,),mask|b)
 for r in range(len(a)):walk((r,),1<<r)
 return out
def endpoint(s,r):
 full=(1<<len(s))-1;E=[0]
 for p in paths(s,r):E=add(E,bins(induced(s,full^sum(1<<v for v in p))))
 return E
def apply(B,E,f):
 em=[a*(-1)**j for j,a in enumerate(E)]
 return add(mul(mul(B,B),f),mul([0,1],mul(mul(E,E),shift(f,1))),neg(mul([0,1],mul(mul(em,em),shift(f,-1)))))
def cross_words(s,r,t):
 a=mirror(s,r,t);n=len(s);N=len(a);actual=set()
 # Every relevant cycle contains r -> n+r; anchor there, then follow all simple words.
 def walk(w,mask):
  u=w[-1]
  if len(w)>=3 and len(w)%2 and a[u]>>r&1:actual.add(w)
  todo=a[u]&~mask&~(1<<r)
  while todo:
   b=todo&-todo;todo-=b;walk(w+(b.bit_length()-1,),mask|b)
 walk((r,n+r),(1<<r)|(1<<(n+r)))
 predicted=set();filtered=0
 ps=paths(s,r);rs=[()]+paths(t)
 for p in ps:
  for q in ps:
   for R in rs:
    w=(r,)+tuple(n+v for v in reversed(q))+tuple(2*n+v for v in R)+p[:-1]
    if len(w)%2==0:
     if len(w)==2:filtered+=1
     continue
    need(len(w)>=3 and len(set(w))==len(w),'unique marked cycle vertices')
    need(all(a[u]>>v&1 for u,v in zip(w,w[1:]+w[:1])),'actual crossing closure')
    need(w not in predicted,'unique triple reconstruction');predicted.add(w)
 need(actual==predicted,'exact crossing cycle set')
 return len(actual),filtered

def main():
 contexts=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   for r in range(n):contexts.append((tour(n,code),r))
 rng=random.Random(1971254)
 for n,count in ((4,24),(5,12)):
  for j in range(count):contexts.append((tour(n,rng.randrange(1<<(n*(n-1)//2))),j%n))
 physical=signs=cycles=twobound=0
 for s,r in contexts:
  B=bins(s);E=endpoint(s,r);e=sum(len(p)==len(s) and p[0]==r for p in paths(s))
  need(val(E,-1)==(-1)**(len(s)+1)*e,'literal endpoint sign');signs+=1
  EE=mul(E,E);U=[a if j%2==0 else 0 for j,a in enumerate(EE)];V=[a if j%2 else 0 for j,a in enumerate(EE)]
  P=add(mul(B,B),mul([0,2],V));Q=mul([0,1],U);R=mul([0,1],V)
  need(all(a>=0 for p in (P,Q,R) for a in p),'positive coefficient form')
  for t in ((),(0,),tr(2),tr(3),(2,4,1)):
   f=bins(t);actual=bins(mirror(s,r,t));want=apply(B,E,f)
   need(actual==want,'independent odd-cycle full profile')
   delta=add(shift(f,1),neg(shift(f,-1)));Delta=add(shift(f,1),shift(f,-1),[-2*a for a in f])
   need(want==add(mul(P,f),mul(Q,delta),mul(R,Delta)),'difference normal form')
   physical+=1
   if len(s)<=3:
    c,b=cross_words(s,r,t);cycles+=c;twobound+=b
 family=[]
 for k in range(1,5):
  ss=list(tr(k+2));ss[0]^=1<<(k+1);ss[-1]|=1;ss=tuple(ss)
  B=bins(ss);E=endpoint(ss,0)
  wantB=add(mon(k+2),mul([0,1],add(power([1,1],k),neg(power([-1,1],k)))))
  wantE=add(mon(k+1),power([1,1],k))
  need(B==wantB and E==wantE,'explicit side formula')
  actual=bins(mirror(ss,0,()))
  closed=add(mon(2*k+4),mul(mul([2],mon(k+2)),add(power([1,1],k+1),neg(power([-1,1],k+1)))),mul([0,1],add(power([1,1],2*k+1),power([-1,1],2*k+1),neg(mul([0,2],power([-1,0,1],k))))))
  need(actual==closed,'explicit physical quadratic profile')
  need(val(actual,2)==6*9**k+24*6**k+16*4**k-8*3**k-8*2**k+2,'B2 closed formula')
  family.append({'k':k,'profile':actual,'H':val(actual,1),'B2':val(actual,2)})
  if k==3:
   need(actual==[0,0,16,0,64,0,64,0,16,0,1],'Q3 polynomial')
   need(val(actual,1)==161 and val(actual,2)==10304 and val(actual,3)==216009,'collision exact counts')
   need(val(actual,2)==64*161 and val(actual,3)!=729*161,'collision but not full proportionality')
   need(val(apply([0,1],[1],actual),1)==10465,'literal-wrapper operator collision')
 print(json.dumps({'status':'PASS','side_contexts':len(contexts),'actual_full_profiles':physical,'endpoint_signs':signs,'literal_crossing_cycles':cycles,'excluded_singleton_two_walks':twobound,'families':family},sort_keys=True,indent=2))
if __name__=='__main__':main()
