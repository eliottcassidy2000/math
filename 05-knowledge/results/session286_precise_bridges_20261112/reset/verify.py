#!/usr/bin/env python3
"""Exact reset/Ising transfer and literal controls. Standard library only."""
from itertools import product
from fractions import Fraction
from array import array
from math import gcd
import json
import core as v
import transfer as t
from flow import Flow
D=(18,8,3,5,14); CUTS={'A':9,'B':3,'E':18,'C':24}
TS={c:v.attach(D,x) for c,x in CUTS.items()};PS={c:v.packet(T) for c,T in TS.items()}
TABS={c:t.table(T) for c,T in TS.items()};FINAL=(t.table(TS['A'],1),t.table(TS['A'],2))
ST=list(product(range(3),range(3),range(2)));IX={x:i for i,x in enumerate(ST)};N=18

def check(q,msg):
 if not q:raise RuntimeError(msg)
def mat(c):
 A=[[0]*N for _ in range(N)]
 for j,x in enumerate(ST):
  for y,z in t.step({x:1},TABS[c]).items():A[IX[y]][j]=z
 return A
def add(A,B,a=1,b=1):return [[a*x+b*y for x,y in zip(r,s)] for r,s in zip(A,B)]
def mul(A,B):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*B)] for r in A]
def scale(A,c):return [[c*x for x in r] for r in A]
def zero(A):return all(not x for r in A for x in r)
I=[[int(i==j) for j in range(N)] for i in range(N)]
M={c:mat(c) for c in CUTS};F=M['A'];F2=mul(F,F);F3=mul(F,F2)
Z=add(mul(M['C'],M['B']),F2,b=-1);K=mul(Z,Z)
R=add(add(add(mul(mul(M['C'],M['E']),M['B']),mul(F,Z),b=-1),mul(Z,F),b=-1),F3,b=-1)
s0=[[int(s==(0,0,0))] for s in ST]
L=[[sum(t.finish({s:1},f) for f in FINAL) for s in ST]]

def spinword(bits):
 ss=[0]+list(bits)+[0,0] # final A; len(word)=len(bits)+2
 return ''.join({(0,0):'A',(0,1):'B',(1,1):'E',(1,0):'C'}[a,b] for a,b in zip(ss,ss[1:]))
def V(W):
 st={(0,0,0):1}
 for c in W[:-1]:st=t.step(st,TABS[c])
 check(W[-1]=='A','final A')
 return sum(t.finish(st,f) for f in FINAL)
def H(W):return v.word([PS[c] for c in W])['H']
def count(W):return H(W)+2*V(W)
def response(a,b):
 num=-34104960*471**a*104**b+3945600*104**a*471**b-140160*104**(a+b)
 check(num%367==0,'response integrality')
 return num//367
def coefficients(k):
 check(k>=3,'k domain')
 c=count('A'*k);linear=[2*13**(k-2)+2*response(k-i-3,i) for i in range(k-2)]
 pair={}
 for i in range(k-2):
  for j in range(i+1,k-2):
   d=j-i
   pair[i,j]=(-1658880*104**(k-4)-18*13**(k-2)) if d==1 else -55296000*104**(k-d-3)*471**(d-2)
 return c,linear,pair
def energy(bits,C):
 c,a,b=C;return c+sum(x*y for x,y in zip(a,bits))+sum(x*bits[i]*bits[j] for (i,j),x in b.items())
def mincut(C):
 c,a,b=C;n=len(a);G=Flow(n+2);S=n;T=n+1;constant=2*c
 for i in range(n):
  u=2*a[i]+sum(x for (p,q),x in b.items() if i in (p,q))
  if u>=0:G.add(i,T,u)
  else:G.add(S,i,-u);constant+=u
 for (i,j),x in b.items():G.add(i,j,-x);G.add(j,i,-x)
 value,left=G.cut(S,T);bits=[int(i in left) for i in range(n)]
 check((constant+value)%2==0,'cut parity');return (constant+value)//2,bits

def decode_target(k,target,C=None):
 C=coefficients(k) if C is None else C;c,a,b=C;n=k-2
 if k==3:
  for bit in (0,1):
   if c+a[0]*bit==target:return [bit]
  return None
 known={}
 for depth in range((n+1)//2):
  left=depth;right=n-1-depth
  residual=target-c-sum(a[i]*x for i,x in known.items())-sum(w*known[i]*known[j] for (i,j),w in b.items() if i in known and j in known)
  if residual%13**depth:return None
  lam=pow(8,depth,13)*pow(471,k-depth-5,13)%13
  digit=(residual//13**depth)*pow(lam,-1,13)%13
  if left==right:
   coeff=(5-6*(known[0]+known[n-1]))%13
   options=[(x,) for x in (0,1) if coeff*x%13==digit]
  elif depth==0:
   options=[(x,y) for x,y in product(range(2),repeat=2) if (3*x+2*y-6*x*y)%13==digit]
  else:
   l=(3-6*known[n-1])%13;r=(2-6*known[0])%13
   options=[(x,y) for x,y in product(range(2),repeat=2) if (l*x+r*y)%13==digit]
  check(len(options)<=1,'unique local13adic digit')
  if not options:return None
  known[left]=options[0][0];known[right]=options[0][-1]
 bits=[known[i] for i in range(n)]
 return bits if energy(bits,C)==target else None

def literal_H(T):
 n=len(T);full=(1<<n)-1;dp=array('Q',[0])*(n*(1<<n))
 for j in range(n):dp[(1<<j)*n+j]=1
 for S in range(1,full+1):
  rem=full^S
  for j in range(n):
   val=dp[S*n+j]
   if not val:continue
   out=T[j]&rem
   while out:
    bit=out&-out;out-=bit;z=bit.bit_length()-1;dp[(S|bit)*n+z]+=val
 return sum(dp[full*n:])
def graph(W):
 T=TS[W[0]]
 for c in W[1:]:T=v.fuse(T,TS[c])
 return T
def reroot(T,r):
 ids=[r]+[i for i in range(len(T)) if i!=r]
 return tuple(sum(1<<j for j,b in enumerate(ids) if T[a]>>b&1) for a in ids)

def main():
 # Degree-six annihilator is a full matrix identity, not fitted scalar data.
 P=I
 for lam in (0,13,13,104,104,471):P=mul(P,add(F,I,b=-lam))
 check(zero(P),'full matrix annihilator')
 check(mul(F,K)==scale(K,104) and mul(K,F)==scale(K,104),'rankone outer propagation')
 check(zero(mul(Z,K)) and zero(mul(K,Z)),'third disjoint response vanishes')
 check(mul(mul(L,K),s0)[0][0]==-27648000,'K boundary')
 check(mul(mul(L,R),s0)[0][0]==-829440,'R boundary')
 check(mul(mul(L,R),F)==scale(mul(L,R),104),'R observed right propagation')
 powers=[I]
 for _ in range(6):powers.append(mul(F,powers[-1]))
 for a in range(6):
  check(mul(mul(Z,powers[a]),Z)==scale(K,471**a),'Z middle propagation certificate')
  check(mul(mul(L,powers[a]),R)==scale(mul(L,R),104**a),'R observed left propagation certificate')
  for b in range(6):
   check(mul(mul(mul(mul(L,powers[a]),Z),powers[b]),s0)[0][0]==response(a,b),'closed response rectangle')
 # Since proposed exponentials satisfy the annihilator, certificates imply all powers.
 check(Fraction(471,104)<5 and Fraction(471,104)**2>20,'rational scale separation')
 neg=Fraction(7891200*5-68209920,367)+26
 pos=Fraction((7891200-Fraction(55296000,471))*20-(68209920+Fraction(55296000,471))-280320,367)-Fraction(3317760,104)-442
 check(neg<0 and pos>0,'uniform marginal signs by exact rational bounds')
 enumerated=0;small=[]
 for k in range(3,11):
  C=coefficients(k);vals=[]
  for bits in product(range(2),repeat=k-2):
   q=energy(bits,C);check(q==count(spinword(bits)),'all spin transfer versus quadratic');check(decode_target(k,q,C)==list(bits),'exact inverse small');vals.append(q);enumerated+=1
  lo=[int(i<k//2) for i in range(k-2)];hi=[1-x for x in lo]
  check(energy(lo,C)==min(vals) and vals.count(min(vals))==1,'unique minimum')
  check(energy(hi,C)==max(vals) and vals.count(max(vals))==1,'unique maximum')
  m,bits=mincut(C);check(m==min(vals) and bits==lo,'flow control');check(len(set(vals))==2**(k-2),'injective full small cube');check(decode_target(k,min(vals)-2,C) is None,'negative target control')
  small.append(dict(k=k,minimum=min(vals),maximum=max(vals),distinct=len(set(vals)),coefficient_gcd=gcd(*C[1],*C[2].values())))
 rows=[]
 for k in (12,20,50,100):
  C=coefficients(k);lo=[int(i<k//2) for i in range(k-2)];hi=[1-x for x in lo];m,bits=mincut(C)
  check(bits==lo and m==energy(lo,C),'large exact mincut')
  for i,a in enumerate(C[1]):
   marginal=a+sum(x for (p,q),x in C[2].items() if i in (p,q))
   check(a<0 if i<k//2 else marginal>0,'uniform coordinate monotonicity')
  for trial in range(6):
   probe=[((i*i+trial*i+trial)//(trial+1))%2 for i in range(k-2)];check(decode_target(k,energy(probe,C),C)==probe,'large exact inverse')
  rows.append(dict(k=k,min_word=spinword(lo),max_word=spinword(hi),minimum=m,maximum=energy(hi,C)))
 controls=[]
 for W in ('A','AA','AAA','BCA'):
  T=graph(W);r=1+5*(len(W)-1);G=v.fuse((2,4,1),reroot(T,r));h=literal_H(G)
  check(h==count(W),'literal reroot triangle')
  controls.append(dict(word=W,order=len(G),H=h))
 print(json.dumps(dict(status='PASS',matrix_annihilator_roots=[0,13,13,104,104,471],response_rectangle=[6,6],all_spin_checks=enumerated,literal_controls=controls,small=small,extrema=rows,decoder_residues={'outer':[0,3,2,12],'later_coefficients':[[3,2],[10,2],[3,9],[10,9]],'single':[5,12,6]},scope='Exact reset Ising/extrema/injective inverse family; no universal prescribed target coverage'),sort_keys=True,indent=2))
if __name__=='__main__':main()
