"""Literal cyclic-group audit of all-reset intersection theorem; no producer imports."""
from collections import Counter
import json
C=Counter()
def gate(x,k,ctx):
 C[k]+=1
 if not x:raise RuntimeError((k,ctx))
def v(n):
 if not n:return float('inf')
 z=0
 while n%5==0:z+=1;n//=5
 return z
def coeff(w):
 B=0;c=0
 for i,a in enumerate(w):B+=3**(len(w)-1-i)*2**c;c+=a
 return 3**len(w),2**c,B
def pair(w,m):
 P,Q,B=coeff(w);q=pow(Q,-1,m);return P*q%m,B*q%m
def mul(x,y,m):return x[0]*y[0]%m,(x[0]*y[1]+x[1])%m
def group(x,m):
 S=set();y=(1,0)
 while y not in S:S.add(y);y=mul(x,y,m)
 return S
def params(L,a):
 P,Q,B=coeff((1,)*L+(a,));t=next(t for t in (1,2,4) if pow(P*pow(Q,-1,5)%5,t,5)==1)
 return t,min(v(P-Q),v(B)) if t==1 else v(P**t-Q**t),v(2**a-2)
def g(e,t,r,h):
 if t==1:return 2*5**max(0,min(r,e)-1,e+r-h-1)
 if e<=h:return 5**max(0,min(r,e)-1)
 return 2*5**max(0,e-h-1,min(r,e)-1)
def words(n):
 yield ()
 for a in range(1,n+1):
  for w in words(n-a):yield (a,)+w
entries=list(words(4))
for L in range(16):
 for a in list(range(2,42))+[101,501,2501]:
  w=(1,)*L+(a,);t,r,h=params(L,a);gs=[g(e,t,r,h) for e in range(1,15)]
  gate(all(y%x==0 for x,y in zip(gs,gs[1:])), 'divisibility_chain',(L,a))
  distinct=[1]
  for x in gs:
   if x!=distinct[-1]:distinct.append(x)
  while distinct[-1]<=L:distinct.append(5*distinct[-1])
  radices=[y//x for x,y in zip(distinct,distinct[1:])]
  digits=[];q=L
  for b in radices:digits.append(q%b);q//=b
  digitcount=1
  for d in digits:digitcount*=d+1
  good=[j for j in range(L+1) if all(j%x<=L%x for x in gs)]
  gate(len(good)==digitcount,'mixed_radix_count',(L,a,gs))
  for j in range(L+1):
   q=j;ds=[]
   for b in radices:ds.append(q%b);q//=b
   gate((j in good)==all(x<=y for x,y in zip(ds,digits)), 'mixed_radix_predicate',(L,a,j))
  for e in range(1,5):
   m=5**e;S=group(pair(w,m),m);A=pair((1,),m);power=(1,0);actual=None
   for j in range(1,2*5**(e-1)+1):
    power=mul(A,power,m)
    if power in S:actual=j;break
   ge=g(e,t,r,h)
   gate(actual==ge,'literal_intersection_modulus',(L,a,e,actual,ge))
   gate(len(S)==t*5**max(0,e-r),'literal_affine_order',(L,a,e))
   for u in entries:
    n=sum(pair((1,)*j+u,m) in S for j in range(L+1))
    gate(n<=1+L//ge,'all_entry_layer_cap',(L,a,e,u))
    if not u:gate(n==1+L//ge,'empty_entry_layer_max',(L,a,e))
   for j in range(L+1):
    # Completion entry: 1^j u=w.
    u=w[j:]
    n=sum(pair((1,)*k+u,m) in S for k in range(L+1))
    gate(n==1+j//ge+(L-j)//ge,'completion_entry_count',(L,a,e,j))
print(json.dumps({'status':'PASS','gates':sum(C.values()),'counts':dict(C),'universe':{'L':'0..15','a':'2..41 plus101,501,2501','depth':'1..4','entry_clock_max':4,'digit_depth':14}},indent=2))
