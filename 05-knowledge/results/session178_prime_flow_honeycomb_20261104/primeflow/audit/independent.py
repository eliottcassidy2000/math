import json
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path

def check(x,msg):
 if not x:raise RuntimeError(msg)
def vp(n,p):
 z=0
 while n%p==0:n//=p;z+=1
 return z
def sigma(p,a):return sum(p**i for i in range(a+1))
def order(p,r):
 x=1
 for d in range(1,r):
  x=x*p%r
  if x==1:return d
 raise RuntimeError('order')
def lg(x,p):return next(k-1 for k in range(1,x+2) if p**k>x)
P=[p for p in range(3,100) if all(p%d for d in range(2,p))]
count=0
for p in P:
 for r in P:
  if p==r:continue
  d=order(p,r)
  for a in range(1,150):
   n=a+1
   predicted=0 if n%d else vp(p**d-1,r)+vp(n//d,r)-vp(p-1,r)
   check(vp(sigma(p,a),r)==predicted,('LTE',p,r,a))
   count+=1
S=[3,5,7,11,13];cert=[]
for q in [5,13]:
 E=[(p,r) for p in S for r in S if p!=r and (order(p,r)%4!=0 if p==q else order(p,r)%2!=0)]
 for p in S:
  reach={p}
  for _ in S:reach|={r for u,r in E if u in reach}
  check(reach==set(S),'initial graph')
 C=max(sum(0 if order(p,r)==1 else vp(p**order(p,r)-1,r) for p,t in E if t==r) for r in S)
 B=2
 while C+4*(lg(B+1,3)+1)>=B:B*=2
 check((C,B)==((3,16) if q==5 else (4,32)),'box')
 for A in range(B,513):check(3**(A-C)>(A+1)**4,'box inequality')
 domains={p:[a for a in range(1,B) if (a%4==1 if p==q else a%2==0)] for p in S}
 caps={r:sum(max(vp(sigma(p,a),r) for a in domains[p]) for p in S if p!=r) for r in S}
 clipped={p:[a for a in domains[p] if a<=caps[p]] for p in S}
 if q==5:check(not clipped[13] and caps[13]==1,'q5 rejection')
 else:
  check({p:max(clipped[p]) for p in S}=={3:4,5:2,7:4,11:4,13:1},'q13 bounds')
  check(all(sigma(5,a)%r for a in clipped[5] for r in S if r!=5),'q13 outdegree')
 cert.append({'q':q,'C':C,'B':B,'direct_column_maxima':caps,'remaining_exponents':clipped})
# Independent cut check by directly cancelling numerator/denominator,
# covering all nonempty cuts of six known even perfect numbers.
cuts=0
for k in [2,3,5,7,13,17]:
 q=2**k-1;N=2**(k-1)*q
 check(sigma(2,k-1)*sigma(q,1)==2*N,'perfect control')
 for p,a,epsilon in [(2,k-1,0),(q,1,1)]:
  ratio=Fraction(sigma(p,a),2**epsilon*p**a)
  O,I=ratio.numerator,ratio.denominator
  check((I<O<2*I) if not epsilon else (O<I<2*O),'cut directions')
  cuts+=1
out={'status':'PASS','independent_direct_sum_LTE_checks':count,'certificates':cert,'cut_controls':cuts,'scope':'Self-contained code, no producer imports. All guards explicit under -O.'}
print(json.dumps(out,sort_keys=True,indent=2))
