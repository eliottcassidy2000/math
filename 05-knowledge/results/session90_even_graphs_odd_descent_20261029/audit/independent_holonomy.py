#!/usr/bin/env python3
"""Independent square-root transport audit for C3, Paley7 and C3[C3,C3,C3]."""
from itertools import permutations,product
from collections import Counter
import json
COUNT=0
def ck(x,msg):
 global COUNT
 COUNT+=1
 if not x:raise RuntimeError(msg)
def mul(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inv(p):return tuple(p.index(i) for i in range(len(p)))
def power(p,k):
 q=tuple(range(len(p)))
 while k:
  if k&1:q=mul(q,p)
  p=mul(p,p);k//=2
 return q

def inspect(name,arc,A,tau,exhaustive=False):
 n=len(tau);one=tuple(range(n));As=set(A)
 ck(len(As)==len(A), 'distinct group elements')
 ck(all(mul(a,b) in As for a in A for b in A),'closure')
 ck(all(arc(i,j)==arc(a[i],a[j]) for a in A for i in range(n) for j in range(i)),'automorphisms')
 ck(all(arc(i,j)==arc(tau[j],tau[i]) for i in range(n) for j in range(i)),'anti')
 if exhaustive:
  full={p for p in permutations(range(n)) if all(arc(i,j)==arc(p[i],p[j]) for i in range(n) for j in range(i))}
  ck(full==As,'exhaustive automorphism group')
 I=[mul(a,tau) for a in A if mul(mul(a,tau),mul(a,tau))==one]
 C={a for a in A if mul(a,tau)==mul(tau,a)}
 ck(len(I)*len(C)==len(A),'anti orbit stabilizer')
 K={mul(t,tau) for t in I}
 ck(K=={a for a in A if mul(mul(tau,a),tau)==inv(a)},'inverted locus')
 g={}
 for t,s in product(I,repeat=2):
  r=power(mul(t,s),(len(A)+1)//2)
  ck(mul(r,r)==mul(t,s),'square root')
  ck(mul(mul(r,s),inv(r))==t,'transport')
  g[t,s]=r
 H=Counter();baseH=Counter();witness=None
 for u,t,s in product(I,repeat=3):
  h=mul(mul(inv(g[u,s]),g[u,t]),g[t,s])
  ck(mul(h,s)==mul(s,h),'holonomy centralizer')
  H[h]+=1
  if s==tau:baseH[h]+=1
  if h!=one and witness is None:witness=dict(u=u,t=t,s=s,h=h)
 commute=all(mul(x,y)==mul(y,x) for x in K for y in K)
 ck((len(H)==1)==commute,'flat iff commute')
 loop_nonassociative=False
 if name=='C3[C3,C3,C3]':
  op={(x,y):power(mul(mul(x,mul(y,y)),x),(len(A)+1)//2) for x in K for y in K}
  ck(all(v in K for v in op.values()),'twisted loop closed')
  for x in K:
   ck(op[one,x]==x and op[x,one]==x,'loop identity')
   ck(len({op[x,y] for y in K})==len(K)==len({op[y,x] for y in K}),'loop latin')
  for x,y in product(K,repeat=2):
   ck(inv(op[x,y])==op[inv(x),inv(y)],'automorphic inverse')
  for x,y,z in product(K,repeat=3):
   ck(op[x,op[y,op[x,z]]]==op[op[x,op[y,x]],z],'left Bol convention')
   if op[op[x,y],z]!=op[x,op[y,z]]:loop_nonassociative=True
  ck(loop_nonassociative,'nonassociative order27 loop')
 return dict(loop_nonassociative=loop_nonassociative,name=name,vertices=n,Aut=len(A),anti_involutions=len(I),centralizer=len(C),inverted_locus=len(K),flat=commute,holonomy_distribution=sorted(H.values()),fixed_base_holonomy_distribution=sorted(baseH.values()),nonflat_witness=witness)

def main():
 rows=[]
 rows.append(inspect('C3',lambda i,j:(j-i)%3==1,[tuple((i+a)%3 for i in range(3)) for a in range(3)],tuple((-i)%3 for i in range(3)),True))
 rows.append(inspect('Paley7',lambda i,j:(j-i)%7 in {1,2,4},[tuple((a*i+b)%7 for i in range(7)) for a in (1,2,4) for b in range(7)],tuple((-i)%7 for i in range(7)),True))
 def arc(i,j):
  a,b=divmod(i,3);c,d=divmod(j,3)
  return (c-a)%3==1 if a!=c else (d-b)%3==1
 A=[tuple(3*((i//3+a)%3)+(i%3+b[i//3])%3 for i in range(9)) for a in range(3) for b in product(range(3),repeat=3)]
 tau=tuple(3*((-(i//3))%3)+(-(i%3))%3 for i in range(9))
 rows.append(inspect('C3[C3,C3,C3]',arc,A,tau,True))
 print(json.dumps(dict(status='PASS',checks=COUNT,rows=rows),indent=2))
if __name__=='__main__':main()
