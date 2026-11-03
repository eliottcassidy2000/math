#!/usr/bin/env python3
"""Literal subset-DP checks; no external dependencies or repository imports."""
import itertools as it,json,math
from common import B,cover,evalF,ext,ordinal,paths,strong,require
C=((5,3,0),(3,7,2),(0,2,2))
def wrapper(parts,r=2,split=1):
 a=ordinal(parts);n=len(a);cut=sum(map(len,parts[:split]));zm=((1<<r)-1)<<n
 out=[v|(zm if i>=cut else 0) for i,v in enumerate(a)]
 out += [((1<<cut)-1)|sum(1<<(n+j) for j in range(i+1,r)) for i in range(r)]
 return out

def survive(a,deleted):
 remain=[i for i in range(len(a)) if i not in deleted]
 return strong([sum(1<<j for j,y in enumerate(remain) if a[x]>>y&1) for x in remain])
def bilinear(p,q):
 x=[math.factorial(i)*p[i] for i in (1,2,3)]
 y=[math.factorial(i)*q[i] for i in (1,2,3)]
 return sum(C[i][j]*x[i]*y[j] for i in range(3) for j in range(3))
# Exhaustively enumerate all run skeletons for two distinguishable feedback vertices.
counts={}
def walk(last,left,i,j):
 if not left and i and j:counts[i,j]=counts.get((i,j),0)+1
 for x in ['A','L',*sorted(left)]:
  if x=='A' and i==3 or x=='L' and j==3:continue
  ok=last is None or (last=='A' and x=='L') or (last=='L' and isinstance(x,int)) or (isinstance(last,int) and x=='A') or (isinstance(last,int) and isinstance(x,int) and last<x)
  if ok:walk(x,left-{x} if isinstance(x,int) else left,i+(x=='A'),j+(x=='L'))
walk(None,{0,1},0,0)
require(counts=={(i+1,j+1):C[i][j] for i in range(3) for j in range(3) if C[i][j]},'skeleton coefficients')
require(all(any(((B[0][i]>>j)&1)!=((B[1][perm[i]]>>perm[j])&1) for i in range(5) for j in range(5)) for perm in it.permutations(range(5))),'nonisomorphic seeds')
require(all(strong(x) for x in B),'strong seeds')
seed=[cover(x) for x in B]; query=[cover(ext(x)) for x in B]
require(seed[0]==seed[1],'equal full seed profile')
require(query[1]==[0,19,109,150,75,15,1],'query1 profile')
require(query[0]==[0,17,95,138,73,15,1],'query0 profile')
records=[];plain=None;graphs={}
for bits in it.product(range(2),repeat=2):
 for address in (None,0,1):
  parts=[B[b] for b in bits]
  if address is not None:parts[address]=ext(parts[address])
  a=wrapper(parts);p=cover(a)
  require(p[1]==bilinear(cover(parts[0]),cover(parts[1])),'Hamilton formula')
  require(all(survive(a,{v}) for v in range(len(a))),'2-strong')
  cuts=[list(S) for S in it.combinations(range(len(a)),2) if not survive(a,set(S))]
  require(cuts==[[len(a)-2,len(a)-1]],'unique minimum cut')
  if address is None:
   if plain is None:plain=p
   require(p==plain,'full-profile equality')
   require(p[1]==152757,'storage count')
   graphs[bits]=a
  else:require(p[1]==585033+60882*bits[address],'one-vertex read')
  records.append({'bits':bits,'address':address,'H':p[1],'profile':p,'cuts':cuts})
for b,a in graphs.items():
 for c,d in graphs.items():require(sum((x^y).bit_count() for x,y in zip(a,d))//2==sum(x!=y for x,y in zip(b,c)),'arc isometry')
larger=[]
for address in (None,0,1):
 bits=(0,1,0);parts=[B[b] for b in bits]
 if address is not None:parts[address]=ext(parts[address])
 h=paths(wrapper(parts))[-1]
 b=0 if address is None else bits[address]
 expected=(33*9**2-651*84**2+504*411**2) if address is None else ((689+2*b)*9**2-(3071+254*b)*84**2+(2036+200*b)*411**2 if address==0 else (561+66*b)*9-(145824+20832*b)*84+(730296+81648*b)*411)
 require(h==expected,'M3 literal Hamilton read')
 larger.append({'bits':bits,'address':address,'H':h})
res=[]
for M in range(2,31):
 h=33*9**(M-1)-651*84**(M-1)+504*411**(M-1)
 require(h%8==(5 if M==2 else 1),'mod8 clock')
 require(h%3**M==0 and h%3**(M+1)!=0,'exact 3 valuation')
 res.append({'M':M,'mod8':h%8,'mod9':h%9,'v3':M})
# Literal resilience checks beyond r=2, no full-profile calculation here.
resilience=[]
for r in (1,3,4):
 a=wrapper([B[0],B[1]],r)
 for k in range(r):require(all(survive(a,set(S)) for S in it.combinations(range(len(a)),k)),'r resilience')
 cuts=[list(S) for S in it.combinations(range(len(a)),r) if not survive(a,set(S))]
 require(cuts==[list(range(10,10+r))],'unique r cut')
 resilience.append({'r':r,'order':len(a),'unique_cut':cuts[0],'H':paths(a)[-1]})
print(json.dumps({'skeleton':[[i,j,v] for (i,j),v in sorted(counts.items())],'graphs':records,'clocks':res,'resilience':resilience,'M3_literal':larger},sort_keys=True,indent=2))
