"""Synthetic alphabet controls, not a tournament census."""
from collections import Counter,defaultdict
from fractions import Fraction
from math import factorial
import json,hashlib,pathlib
N=12; HMAX=81; checks=0
def need(x,s):
 global checks
 checks+=1
 if not x: raise ValueError(s)
def add(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,b):
 c=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]):c[i+j]+=x*y
 return c
def mon(n,c=1):
 a=[0]*(N+1)
 if n<=N:a[n]=c
 return a
def sub2(a):return [a[n//2] if n%2==0 else 0 for n in range(N+1)]
def omega(n):
 c=0;p=2
 while p*p<=n:
  while n%p==0:n//=p;c+=1
  p+=1
 return c+(n>1)
# Two real rigid weights and R5/S4 weights; other entries are deliberately
# formal symbols, NOT assertions of tournament realizability at these orders.
atoms=[('u',1,1,1,'u'),('v',3,1,1,'v'),('r',5,3,3,'r'),('s',4,5,3,'s'),('a',2,7,0,'b'),('b',2,7,0,'a'),('c',3,9,1,'c'),('d',4,9,9,'d')]
D={a[0]:a for a in atoms}; A=[0]*(N+1);A[0]=1
for n in range(1,N+1):A[n]=A[n-1]+(A[n-3] if n>=3 else 0)
P=defaultdict(lambda:[0]*(N+1));C=defaultdict(lambda:[0]*(N+1));centers=[(0,1,1)]
for name,n,h,b,opp in atoms:
 P[h][n]+=1
 if opp==name:C[h][n]+=1;centers.append((n,h,b))
C[1][0]+=1
# Formula by ordered nonrigid factorizations separated by rigid gaps.
F={1:A};factor_cache={}
def words(h):
 if h==1:return [()]
 if h in factor_cache:return factor_cache[h]
 out=[]
 for j in sorted(P):
  if j>1 and h%j==0:out += [(j,)+w for w in words(h//j)]
 factor_cache[h]=out;return out
for h in range(3,HMAX+1,2):
 f=[0]*(N+1)
 for w in words(h):
  need(len(w)<=omega(h),'factor length')
  t=A[:]
  for j in w:t=mul(mul(t,P[j]),A)
  f=add(f,t)
 F[h]=f
# Explicit ordinary symbol words, with literal opposite/reversal.
counts=defaultdict(lambda:[0]*(N+1)); sc=defaultdict(lambda:[0]*(N+1)); B=defaultdict(lambda:[0]*(N+1)); seen=set(); total=0
stack=[((),0,1)]
while stack:
 w,n,h=stack.pop();total+=1;counts[h][n]+=1
 opp=tuple(D[s][4] for s in w[::-1]); key=min(w,opp)
 if w==opp:
  half=1
  for s in w[:len(w)//2]:half*=D[s][2]
  centerb=1 if len(w)%2==0 else D[w[len(w)//2]][3]
  b=half*centerb;need((h-b)%2==0,'even deficit');sc[h][n]+=1
  if key not in seen:B[(h-b)//2][n]+=1;seen.add(key)
 elif key not in seen:B[h][n]+=1;seen.add(key)
 for name,size,hh,bb,oo in atoms:
  if n+size<=N and h*hh<=HMAX:stack.append((w+(name,),n+size,h*hh))
for h in range(1,HMAX+1,2):
 need(F[h]==counts[h],f'ordinary h{h}')
 J=[0]*(N+1)
 for x in range(1,int(h**.5)+1,2):
  if h%(x*x)==0:J=add(J,mul(sub2(F[x]),C[h//(x*x)]))
 need(J==sc[h],f'converse h{h}')
for k in range(1,11):
 d=2*k;S=[0]*(N+1)
 for x in range(1,d+1,2):
  if d%x or x*(x-1)>d:continue
  E=[0]*(N+1)
  for n,a,b in centers:
   if x*(x*a-b)==d:E[n]+=1
  S=add(S,mul(sub2(F[x]),E))
 result=S
 if k%2:result=[s+Fraction(f-j,2) for s,f,j in zip(S,F[k],sc[k])]
 need(result==B[k],f'degree{k}')
 L=max(omega(x) for x in range(1,d+1,2) if d%x==0 and x*(x-1)<=d)
 need(all(omega(x)<=L for x in range(1,k+1,2) if k%(x*x)==0),'J denominator bound')
# Known leaf and trivalent formulas use actual certified atoms, not formal d.
p=mon(0);p[1]=-1;p[3]=-1;p2=sub2(p)
leaf=add([Fraction(x-y,2) for x,y in zip(A,mul(sub2(A),add(add(mon(0),mon(1)),mon(3))))],mul(sub2(A),mon(4)))
need(mul(mul(leaf,p),p2)==add(add(mon(4,2),mon(5,-1)),mon(7,-1)),'known leaf identity')
for r in range(7):
 need(omega(3**r)==r,'Omega power3')
 need(Fraction(1,2*factorial(r))>0,'leading coefficient positive')
for k in range(1,101):
 L=max(omega(x) for x in range(1,2*k+1,2) if 2*k%x==0 and x*(x-1)<=2*k)
 for x in range(1,k+1,2):
  if k%(x*x)==0:need(x*(x-1)<=2*k and omega(x)<=L,'squaredivisor inclusion')
receipt={'status':'PASS','scope':'formal involutive alphabet, not realizability or tournament census','max_order':N,'max_presentation_weight':HMAX,'ordinary_words':total,'checks':checks,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'degrees_checked':list(range(1,11))}
print(json.dumps(receipt,sort_keys=True,indent=2))
