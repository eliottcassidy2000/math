#!/usr/bin/env python3
"""Exact portable controls. No assertions; all gates survive python -O."""
import json,hashlib
from itertools import product
from pathlib import Path
from blocks import graph
from inherited import profile,require
from transfer import S,T,mm,mv,dvalue,ell,v

def add(A,B,c=1):return [[a+c*b for a,b in zip(x,y)] for x,y in zip(A,B)]
I=[[int(i==j) for j in range(6)] for i in range(6)]
A=mm(S[1],S[0]); zero=[[0]*6 for _ in range(6)]
require(mm(mm(A,add(A,I,-4)),mm(add(A,I,-1),add(A,I,-1)))==zero,'matrix annihilator')
Fgrid=[[9,27,93,351],[29,95,353,1379],[107,365,1391,5489],[413,1439,5537,21923]]
Ggrid=[[5,23,89,347],[17,83,341,1367],[59,317,1343,5441],[221,1247,5345,21731]]
for j,k in product(range(4),repeat=2):
 r=j+k+2;B=(4**r+6*r-1)//3
 require(dvalue([0,1]*r,2*j)==Fgrid[j][k]==B+4**j-2+(j==0),'F initial grid')
 require(dvalue([0,1]*r,2*j+1)==Ggrid[j][k]==B-2*4**j-2,'G initial grid')
checks=0; controls=[]
for m in range(2,8):
 for bits in product(range(2),repeat=m):
  U=3*(1+sum(b*2**(m-1-i) for i,b in enumerate(bits)))
  for p in range(m-1):
   rows=graph(bits,p); x=profile(rows,0)
   require(x[2:]==(U,3),'universal sidecar')
   require(x[0]==dvalue(bits,p),'literal vs transfer')
   if bits[-1]==1 and p<m-2:require(x[1]==0,'E boundary')
   checks+=1
for r in range(2,6):
 B=(4**r+6*r-1)//3; vals=[]
 for p in range(2*r-2):
  D=B if p==0 else B+(-2)**p-2
  x=profile(graph([0,1]*r,p),0)
  require(x==(D,0,4**r+2,3),'all scalar coordinates small literal')
  vals.append(D)
 require(len(set(vals))==2*r-2,'distinct scalar fiber')
 controls.append({'r':r,'sidecar':[0,4**r+2,3],'scalars':vals})
long_checks=0
for r in range(2,51):
 B=(4**r+6*r-1)//3
 for p in range(2*r-2):
  require(dvalue([0,1]*r,p)==(B if p==0 else B+(-2)**p-2),'larger transfer control')
  long_checks+=1
files=['PROOF.md','verify.py','blocks.py','transfer.py','inherited.py']
print(json.dumps({'status':'PASS','scope':'All-order proof uses matrix identity and recurrence; finite literal checks do not establish universal quantifiers.','universal_word_literal_checks':checks,'alternating_literal_controls':controls,'long_transfer_controls':long_checks,'source_sha256':{f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}},sort_keys=True,indent=2))
