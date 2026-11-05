from fractions import Fraction as F
from collections import deque
from math import gcd
from pathlib import Path
import json
from core import constants,glue,gadget,inventories,cycles,require
P=Path(__file__).parent
I=tuple(tuple(F(i==j)for j in range(3))for i in range(3))
def mul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j]for k in range(3))for j in range(3))for i in range(3))
def power(A,k):
 out=I
 while k:
  if k&1:out=mul(out,A)
  A=mul(A,A);k//=2
 return out
def inv(A):
 W=[list(A[i])+list(I[i])for i in range(3)]
 for j in range(3):
  k=next(k for k in range(j,3)if W[k][j]);W[j],W[k]=W[k],W[j];a=W[j][j];W[j]=[x/a for x in W[j]]
  for i in range(3):
   if i!=j:
    a=W[i][j];W[i]=[x-a*y for x,y in zip(W[i],W[j])]
 return tuple(tuple(row[3:])for row in W)
def comm(A,B):return mul(mul(mul(A,B),inv(A)),inv(B))
def E(i,j,a=1):return tuple(tuple(I[x][y]+a*(x==i and y==j)for y in range(3))for x in range(3))
def M(L):
 g,q,b=constants(L);return tuple(tuple(map(F,row))for row in [[g,-(L%2),L],[1,-int(L==1),0],[q,0,b]])
def apply(A,v):return tuple(sum(A[i][j]*v[j]for j in range(3))for i in range(3))
def mod(A,m):return tuple(tuple((x.numerator*pow(x.denominator,-1,m))%m for x in row)for row in A)if m>1 else tuple((0,0,0)for _ in range(3))
A=M(1);T=E(0,1);A2=power(A,2)
H=comm(mul(M(9),inv(M(11))),mul(M(11),inv(M(13))))
require(H==E(0,1,F(3,4)),'exact central transvection')
require(mul(comm(T,mul(mul(A2,T),inv(A2))),power(T,2))==E(2,1,-4),'recover E32')
T32=E(2,1);S=mul(mul(inv(A),T32),A)
T31=mul(comm(T32,S),T32)
require(T31==E(2,0),'recover E31')
R=mul(T31,inv(T32));T12=mul(mul(mul(A,R),inv(A)),inv(T32))
require(T12==T,'recover E12')
T21=mul(mul(inv(T12),S),T12);require(T21==E(1,0),'recover E21')
T13=mul(mul(inv(T31),mul(mul(A,T32),inv(A))),T31);require(T13==E(0,2),'recover E13')
require(comm(T21,T13)==E(1,2),'recover E23')
false_ratio=mul(M(5),inv(M(7)))
require(false_ratio!=E(2,1,2),'reject false first-row omission')
for z in range(3,35,2):
 v=(z,1,z)
 for G in [E(0,1,-z),E(2,1,-z),E(0,1),E(1,0,-1)]:v=apply(G,v)
 require(v==(1,0,0),'explicit first seed reduction')
 v=(10*z-1,z,15*z)
 for G in [E(0,1,-10),E(1,0,z),E(2,0,15*z),power(mul(mul(E(0,1),E(1,0,-1)),E(0,1)),2)]:v=apply(G,v)
 require(v==(1,0,0),'explicit second seed reduction')
for r in range(42):
 G=mul(mul(mul(mul(E(1,0,17),E(0,1)),E(1,0,2)),E(0,1,-1)),E(1,0,-r))
 require(apply(G,(1,0,0))==(2*r+3,35*r+53,0),'explicit target vector map')
relative=0
for m in range(1,102,2):
 q=m;a=0
 while q%3==0 and q>1:q//=3;a+=1
 if q==1:continue
 d=3**max(a-1,0);k=d*((4*pow(3*d,-1,q))%q)
 R12=power(H,k)
 require(mod(R12,q)==mod(E(0,1),q) and mod(R12,3**a)==mod(I,3**a),'relative CRT transvection');relative+=1
# All residue outputs from actual positive words, with no inverse or formal vector as a seed.
mods=[];allw={};gens=(1,9,11,13);mat={L:tuple(tuple(int(x)for x in row)for row in M(L))for L in gens}
for m in range(1,42,2):
 seed=(3%m,1%m,3%m);seen={seed:()};queue=deque([seed]);found={};found1={}
 while queue and (len(found)<m or len(found1)<m):
  v=queue.popleft();w=seen[v];r=(53*v[0]-3*v[1]+63*v[2])%m
  if r not in found or v[0] not in found1:
   exact=(3,1,3)
   for L in w:exact=apply(M(L),exact)
   value=int(53*exact[0]-3*exact[1]+63*exact[2]);require(value%2==1 and value%m==r,'exact physical residue')
   found.setdefault(r,dict(word=list(w),H=value));found1.setdefault(v[0],dict(word=list(w),H=int(exact[0])))
  for L in gens:
   z=tuple(sum(mat[L][i][j]*v[j]for j in range(3))%m for i in range(3))
   if z not in seen:seen[z]=w+(L,);queue.append(z)
 require(len(found)==m and len(found1)==m,'every odd-modulus residue observed in both classes')
 mods.append(dict(modulus=m,C1_residues=len(found1),C2_residues=len(found),states_discovered=len(seen),max_word=max(len(v['word'])for v in found.values())))
 allw[str(m)]={kind:{str(r):bank[r]for r in sorted(bank)}for kind,bank in [('C1',found1),('C2',found)]}
# Exact seed identities and primitive contents, through large 3-power unit cycles.
seeds=0
for a in range(1,7):
 m=3**a;values0=set();values1=set();dvalues0=set();dvalues1=set()
 for j in range(1,2*3**(a-1)+1):
  z=(1<<j)+1;v=(z,1,z);w=apply(M(5),v)
  require(gcd(*map(int,w))==1,'primitive optional M5 seed')
  f0=53*v[0]-3*v[1]+63*v[2];f1=53*w[0]-3*w[1]+63*w[2]
  require(f0==113+116*(1<<j)and f1==1419+1472*(1<<j),'exact three-adic seed formulas')
  values0.add(f0%m);values1.add(f1%m);dvalues0.add(v[0]%m);dvalues1.add(int(w[0])%m);seeds+=1
 require(values0=={r for r in range(m)if r%3!=2}and values1=={r for r in range(m)if r%3!=0},'complete three-adic seed residue sets')
 require(dvalues0=={r for r in range(m)if r%3!=1}and dvalues1=={r for r in range(m)if r%3!=0},'complete C1 three-adic seed residues')
# Literal graphs for the final readout, without relying on the symbolic formula.
def converse(T):return [sum(1<<j for j in range(len(T))if T[j]>>i&1)for i in range(len(T))]
def graft(T,entries,G,r):
 n=len(T);out=list(T)+[x<<n for x in G]
 for x in range(n):
  for y in range(len(G)):
   if y==r and x in entries:out[n+y]|=1<<x
   else:out[x]|=1<<(n+y)
 return out
literal=[]
for recipe in [(),(1,),(5,),(7,),(9,),(11,),(13,)]:
 V,r=glue([0],{0},gadget(2));v=(3,1,3)
 for L in recipe:V,r=glue(V,{r},gadget(L));v=apply(M(L),v)
 sourceH=inventories(V)['H'];require(sourceH==v[0] and cycles(V)==1,'literal C1 source readout')
 V,r=glue(V,{r},gadget(5));left,b=glue([2,0],{0},[0]);T=graft(left,{0,b},converse(V),r)
 h=inventories(T)['H'];c=cycles(T);want=int(53*v[0]-3*v[1]+63*v[2])
 require(h==want and c==2,'literal final C2 readout');literal.append(dict(recipe=list(recipe),order=len(T),H=h,C=c,source_H=sourceH,source_C=1))
receipt=dict(status='PASS',relative_modulus_checks=relative,modular_orbit_controls=mods,three_adic_seed_controls=seeds,literal_controls=literal,false_ratio=[[str(x)for x in row]for row in false_ratio])
(P/'WITNESSES.json').write_text(json.dumps(allw,sort_keys=True,indent=2)+'\n')
s=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(s);print(s,end='')
