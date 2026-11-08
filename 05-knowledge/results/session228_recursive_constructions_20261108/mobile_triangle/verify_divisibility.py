#!/usr/bin/env python3
import json,hashlib
from itertools import product
from pathlib import Path
from inherited import profile,require
from transfer import S,T,mv,ell,v

def scalar(seq):
 x=v[:]
 for kind,k in seq:x=mv(S[k] if kind=='s' else T[k],x)
 return sum(a*b for a,b in zip(ell,x))
def graph(seq):
 blocks=[([0] if typ=='s' else [2,4,1],k) for typ,k in seq]
 n=1+sum(len(a) for a,k in blocks);a=[0];off=1
 for b,k in blocks:
  for j,row in enumerate(b):
   a.append((row<<off)|sum(1<<q for q in range(off+len(b),n))|(1 if j>=k else 0))
   if j<k:a[0]|=1<<(off+j)
  off+=len(b)
 return a
alphabet=[('s',0),('s',1)]+[('t',k) for k in range(4)]
controls=0;literal=0
for length in range(1,6):
 for seq in product(alphabet,repeat=length):
  t=sum(x[0]=='t' for x in seq);D=scalar(seq)
  require(D%3**max(t-2,0)==0,'unrestricted factor')
  if length<=3:require(D==profile(graph(seq),0)[0],'literal transfer agreement');literal+=1
  controls+=1
sharp=[]
for t in range(2,41):
 seq=[('t',0)]*(t-2)+[('t',1),('t',1)]
 D=scalar(seq);require(D==41*3**(t-2),'sharp prefix family')
 if t<=3:require(profile(graph(seq),0)[0]==D,'literal sharp family')
 if t<=6:sharp.append({'t':t,'D':D})
# The stronger valuation of the separated family, independently on its explicit formula.
sep=0
for r in range(3,31):
 J=list(range(1,r-1,2));H=(4**r+6*r-1)//9
 for t in range(1,len(J)+1):
  P=J[:t];D=3**t*H+3**(t-1)*sum(4**j-2 for j in P)
  require(D%3**(t-1)==0,'separated common factor')
  if t%3:require(D%3**t!=0,'separated exact valuation')
  sep+=1
files=['DIVISIBILITY.md','verify_divisibility.py','transfer.py','inherited.py']
print(json.dumps({'status':'PASS','block_sequences_lengths1_to5':controls,'literal_profiles_lengths1_to3':literal,'sharp_family_t2_to40':39,'sharp_examples':sharp,'separated_formula_controls':sep,'sha256':{f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}},sort_keys=True,indent=2))
