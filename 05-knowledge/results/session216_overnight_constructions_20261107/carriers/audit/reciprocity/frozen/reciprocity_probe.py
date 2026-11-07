from itertools import permutations
import sympy as sp
from inherited import *
def polynomials(g,r=0):
 n=len(g);N=1<<n;H=[0]*N;S=[0]*N;E=[0]*N
 for A in range(N):
  v=[i for i in range(n) if A>>i&1]
  for p in permutations(v):
   if all(g[p[i]]>>p[i+1]&1 for i in range(len(p)-1)):
    H[A]+=1
    if p and p[0]==r:S[A]+=1
    if p and p[-1]==r:E[A]+=1
 # nonempty path polynomial Q=P-1; powers truncated; F(x)=S*(1+Q)^(x-1)
 Q=H[:];Q[0]=0;power=[0]*N;power[0]=1;coefS=[];coefE=[]
 for k in range(n):
  coefS.append(sum(S[A]*power[(N-1)^A] for A in range(N)))
  coefE.append(sum(E[A]*power[(N-1)^A] for A in range(N)))
  nxt=[0]*N
  for A in range(N):
   B=A
   while B:nxt[A]+=Q[B]*power[A^B];B=(B-1)&A
  power=nxt
 x=sp.Symbol('x');F=0;G=0;fall=1
 for k in range(n):
  if k:fall*=x-k
  F+=coefS[k]*fall/sp.factorial(k);G+=coefE[k]*fall/sp.factorial(k)
 return sp.Poly(F,x),sp.Poly(G,x)
if __name__=='__main__':
 for n in range(1,5):
  for bits in range(1<<(n*(n-1)//2)):
   g=[0]*n;z=0
   for i in range(n):
    for j in range(i+1,n):
     if bits>>z&1:g[i]|=1<<j
     else:g[j]|=1<<i
     z+=1
   F,G=polynomials(g);x=F.gens[0]
   if sp.expand(F.as_expr().subs(x,1-x)-(-1)**(n-1)*G.as_expr())!=0:print('FAIL',g,F,G);raise SystemExit
 print('PASS')
