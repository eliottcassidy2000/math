#!/usr/bin/env python3
"""Independent path-list three-strand matrix, not importing reset producers."""
from itertools import product,permutations
from collections import Counter
import json
D=(18,8,3,5,14); ST=list(product(range(3),range(3),range(2))); IX={s:i for i,s in enumerate(ST)}

def must(x,m):
 if not x:raise RuntimeError(m)

def paths(S):
 return [p for p in permutations([i for i in range(5)if S>>i&1])if all(D[a]>>b&1 for a,b in zip(p,p[1:]))]
P={S:paths(S)for S in range(32)}
TR=[]
for cols in product(range(3),repeat=5):
 masks=[sum(1<<i for i,x in enumerate(cols)if x==c)for c in range(3)]
 TR.extend(product(*(P[s]for s in masks)))

def advance(s,ps,cut):
 p,q,r=ps
 return (1 if p and not(cut>>p[-1]&1)else 2 if p else s[0],s[1]or(1 if q and cut>>q[0]&1 else 2 if q else 0),int(s[2]or bool(r)))

def M(cut):
 out=[[0]*18 for _ in ST]
 for j,s in enumerate(ST):
  for ps in TR:out[IX[advance(s,ps,cut)]][j]+=1
 return out

def L(cut=9,r=0):
 out=[]
 for s in ST:
  c=0
  for ps in TR:
   state=advance(s,ps,cut)
   if state[0]==2 or state[1]==2:continue
   c+=int(bool(ps[1])and ps[1][-1]==r)+int(bool(ps[2])and ps[2][-1]==r)
  out.append(c)
 return out

def mul(A,B):
 return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def sub(A,B):return [[a-b for a,b in zip(x,y)]for x,y in zip(A,B)]
def scale(a,A):return [[a*x for x in r]for r in A]
def add(A,B):return [[a+b for a,b in zip(x,y)]for x,y in zip(A,B)]
def ident():return [[int(i==j)for j in range(18)]for i in range(18)]
def power(A,n):
 B=ident()
 for _ in range(n):B=mul(A,B)
 return B

def run():
 F=M(9);B=M(3);C=M(24);E=M(18);Z=sub(mul(C,B),mul(F,F));K=mul(Z,Z);zero=scale(0,K)
 must(len(TR)==471,'three-strand total')
 must(mul(F,K)==scale(104,K)and mul(K,F)==scale(104,K),'K eigen')
 must(mul(Z,K)==zero and mul(K,Z)==zero,'cubic nilpotency')
 identities=[]
 for m in range(19):
  must(mul(mul(Z,power(F,m)),Z)==scale(471**m,K),'ZFmZ')
  identities.append(m)
 R=sub(sub(sub(mul(mul(C,E),B),mul(F,Z)),mul(Z,F)),power(F,3))
 # R is not a full matrix eigenvector; only boundary contractions are used.
 l=L();s0=[[int(s==(0,0,0))]for s in ST]
 value=lambda X:mul(mul([l],X),s0)[0][0]
 must(value(K)==-27648000,'K scalar');must(value(R)==-829440,'R scalar')
 must(mul(mul([l],R),F)==scale(104,mul([l],R)),'LR right eigen')
 for a in range(19):must(mul(mul([l],power(F,a)),R)==scale(104**a,mul([l],R)),'left boundary eigen')
 # Full annihilator and rectangle are all-order recurrence certificates.
 ann=F
 for root in (13,13,104,104,471):ann=mul(ann,sub(F,scale(root,ident())))
 must(ann==zero,'degree-six annihilator')
 for a in range(6):
  for b in range(6):
   num=-34104960*471**a*104**b+3945600*104**a*471**b-140160*104**(a+b)
   must(num%367==0,'S integer')
   must(value(mul(mul(power(F,a),Z),power(F,b)))==num//367,'S rectangle')
 # Independent full Boolean Mobius transforms; no inferred pair-only assumption.
 spin_checks=[]; mats={(0,0):F,(0,1):B,(1,1):E,(1,0):C}
 for k in range(4,10):
  n=k-2;values=[]
  for mask in range(1<<n):
   spins=[0]+[(mask>>i)&1 for i in range(n)]+[0,0]
   vec=s0
   for i in range(k-1):vec=mul(mats[(spins[i],spins[i+1])],vec)
   j=sum(spins)-9*sum(a*b for a,b in zip(spins,spins[1:]))
   values.append(2*mul([l],vec)[0][0]+2*13**(k-2)*j)
  coef=values[:]
  for i in range(n):
   for mask in range(1<<n):
    if mask>>i&1:coef[mask]-=coef[mask^(1<<i)]
  for mask,c in enumerate(coef):
   if mask.bit_count()>2:must(c==0,'reset Boolean degree')
   elif mask.bit_count()==2:
    i,j=[i for i in range(n)if mask>>i&1];d=j-i
    want=-1658880*104**(k-4)-18*13**(k-2) if d==1 else -55296000*104**(k-d-3)*471**(d-2)
    must(c==want,'pair coefficient')
  spin_checks.append([k,1<<n])
 return dict(spin_checks=spin_checks,path_triples=len(TR),K_scalar=value(K),R_scalar=value(R),matrix_identities=identities,status='PASS',note='Nineteen consecutive exact identities independently establish recurrences via Cayley-Hamilton on the18-state matrix plus the named scalar eigenvalue; Boolean-degree tests are finite controls, all-order proof separately reviewed')
if __name__=='__main__':print(json.dumps(run(),indent=2))
