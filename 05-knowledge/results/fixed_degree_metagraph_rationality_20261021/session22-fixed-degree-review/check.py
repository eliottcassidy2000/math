"""Independent finite formal-alphabet controls; no tournament realizability claim."""
from collections import Counter,defaultdict
from functools import lru_cache
from math import prod
import json, hashlib
from pathlib import Path
N=10
# name,size,h,b,converse index. Distinct same-signature letters test multiplicities.
letters=[('u',1,1,1,0),('v',3,1,1,1),('a',2,3,1,2),('b',2,3,1,3),('c',4,3,3,4),('d',3,9,3,5),('e',2,5,0,7),('f',2,5,0,6),('g',4,15,1,8)]
zero=lambda:[0]*(N+1)
def mon(n):
 a=zero()
 if n<=N:a[n]=1
 return a
def add(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(N+1)]
def sub2(a):return [a[n//2] if n%2==0 else 0 for n in range(N+1)]
def power(a,r):
 out=mon(0)
 for _ in range(r):out=mul(out,a)
 return out
P=defaultdict(zero); C=defaultdict(zero); centers=[(0,1,1)]
for _,n,h,b,t in letters:
 P[h][n]+=1
for i,(_,n,h,b,t) in enumerate(letters):
 if t==i:C[h][n]+=1;centers.append((n,h,b))
C[1][0]+=1
A=zero();A[0]=1
for n in range(1,N+1):A[n]=sum(A[n-t] for t in (1,3) if n>=t)
@lru_cache(None)
def factorizations(j):
 if j==1:return ((),)
 return tuple((t,)+tail for t in P if t>1 and j%t==0 for tail in factorizations(j//t))
@lru_cache(None)
def F(j):
 out=zero()
 for fs in factorizations(j):
  term=power(A,len(fs)+1)
  for t in fs:term=mul(term,P[t])
  out=add(out,term)
 return out
@lru_cache(None)
def J(j):
 out=zero()
 for x in range(1,int(j**.5)+1,2):
  if j%(x*x)==0:out=add(out,mul(sub2(F(x)),C[j//(x*x)]))
 return out
@lru_cache(None)
def S(d):
 out=zero()
 for x in range(1,d+1,2):
  if d%x or x*(x-1)>d:continue
  E=zero()
  for n,a,b in centers:
   if x*(x*a-b)==d:E[n]+=1
  out=add(out,mul(sub2(F(x)),E))
 return out
def B(k):return [S(2*k)[n]+((F(k)[n]-J(k)[n])//2 if k%2 else 0) for n in range(N+1)]
ordinary=defaultdict(zero);sc=defaultdict(zero);deficit=defaultdict(zero);degree=defaultdict(zero)
count=0
def visit(w,n,h):
 global count
 count+=1;ordinary[h][n]+=1
 op=tuple(letters[t][4] for t in reversed(w))
 if w==op:
  x=prod(letters[t][2] for t in w[:len(w)//2]);b=x*(letters[w[len(w)//2]][3] if len(w)%2 else 1)
  sc[h][n]+=1;deficit[h-b][n]+=1;degree[(h-b)//2][n]+=1
 elif w<op:degree[h][n]+=1
 for t,(_,m,g,_,_) in enumerate(letters):
  if n+m<=N:visit(w+(t,),n+m,h*g)
visit((),0,1)
checks=0
def check(a,b,label):
 global checks
 if a!=b:raise RuntimeError((label,a,b))
 checks+=1
for j in range(1,501,2):check(F(j),ordinary[j],('F',j));check(J(j),sc[j],('J',j))
for j in range(1,501,2):check(all((a-b)%2==0 for a,b in zip(F(j),J(j))),True,('parity',j))
check(F(7),zero(),'absent h7')
check(S(14)[4]>0,True,'degree7 exists without h7 in formal alphabet')
for d in range(2,201,2):check(S(d),deficit[d],('S',d))
for k in range(1,101):check(B(k),degree[k],('B',k))
# Exact denominator controls via polynomial arithmetic in SymPy, independent of word counts.
import sympy as sp
z=sp.symbols('z');p=1-z-z**3;p2=p.subs(z,z*z)
def omega(n):
 r=0;t=2
 while t*t<=n:
  while n%t==0:r+=1;n//=t
  t+=1
 return r+(n>1)
Ps={j:sum(v[n]*z**n for n in range(N+1)) for j,v in P.items()}
@lru_cache(None)
def Fs(j):return sp.sympify(sum(prod(Ps[t] for t in fs)/p**(len(fs)+1) for fs in factorizations(j)))
def Js(j):return sum(Fs(x).subs(z,z*z)*sum(C[j//(x*x)][n]*z**n for n in range(N+1)) for x in range(1,int(j**.5)+1,2) if j%(x*x)==0)
def Ss(d):return sum(Fs(x).subs(z,z*z)*sum(z**n for n,a,b in centers if x*(x*a-b)==d) for x in range(1,d+1,2) if d%x==0 and x*(x-1)<=d)
for k in range(1,41):
 L=max(omega(x) for x in range(1,2*k+1,2) if 2*k%x==0 and x*(x-1)<=2*k)
 expr=Ss(2*k)+((Fs(k)-Js(k))/2 if k%2 else 0)
 den=p2**(L+1)*(p**(omega(k)+1) if k%2 else 1)
 check(sp.fraction(sp.cancel(den*expr))[1].has(z),False,('denominator',k))
# Known actual-tournament specializations, leaving R6 arbitrary.
R6=sp.symbols('R6'); AA=1/p;AA2=1/p2
B1=sp.cancel(z**4*AA2+(AA-(1+z+z**3)*AA2)/2)
knownB1=(2*z**4-z**5-z**7)/(p*p2)
check(sp.cancel(B1-knownB1),0,'B1')
B3=R6*AA2+z**10*AA2**2*(1+z+z**3)+(z**5*AA**2-z**5*AA2)/2
knownB3=(z**5*AA**2-z**5*AA2)/2+R6*AA2+z**10*AA2**2*(1+z+z**3)
check(sp.cancel(B3-knownB3),0,'B3')
# Independent local-pole normalization: p(z) ~ rho*(1+3rho^2)*(1-z/rho).
rho=sp.symbols('rho',positive=True)
for r in range(10):
 c=rho**(5*r)/(2*sp.factorial(r)*(rho*(1+3*rho**2))**(r+1))
 check(sp.cancel(c-rho**(4*r-1)/(2*sp.factorial(r)*(1+3*rho**2)**(r+1))),0,('constant',r))
receipt={'status':'PASS','synthetic_only':True,'max_order':N,'synthetic_letters':letters,'literal_words':count,'checks':checks,'F_J_odd_h_through':499,'S_even_deficit_through':200,'B_degree_through':100,'symbolic_denominator_k_through':40,'known_specializations':['B1','B3 with arbitrary R6'],'asymptotic_constant_r_through':9,'proof_sha256':hashlib.sha256(Path('session22-fixed-degree/PROOF.md').read_bytes()).hexdigest()}
Path('session22-fixed-degree-review/receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
