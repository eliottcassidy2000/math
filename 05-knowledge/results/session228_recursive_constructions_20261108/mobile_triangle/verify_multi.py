#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
from itertools import combinations
from inherited import profile,require
from transfer import S,T,mm,mv,ell,v

def sub(A,B,c=1):return [[a-c*b for a,b in zip(x,y)] for x,y in zip(A,B)]
def gmulti(bits,positions):
 positions=set(positions);blocks=[];outs=[];i=0
 while i<len(bits):
  if i in positions:blocks.append([2,4,1]);outs.append(set(range(2*bits[i]+bits[i+1])));i+=2
  else:blocks.append([0]);outs.append({0} if bits[i] else set());i+=1
 n=1+sum(map(len,blocks));rows=[0];off=1
 for b,o in zip(blocks,outs):
  for j,row in enumerate(b):
   rows.append((row<<off)|sum(1<<q for q in range(off+len(b),n))|(1 if j not in o else 0))
   if j in o:rows[0]|=1<<(off+j)
  off+=len(b)
 return rows

def val(bits,positions):
 x=v[:];positions=set(positions);i=0
 while i<len(bits):
  if i in positions:x=mv(T[2*bits[i]+bits[i+1]],x);i+=2
  else:x=mv(S[bits[i]],x);i+=1
 return sum(a*b for a,b in zip(ell,x))
def induced(a,vs):return [sum(1<<j for j,w in enumerate(vs) if a[u]>>w&1) for u in vs]
def strong(a):
 for r in range(len(a)):
  seen=1<<r
  while True:
   new=seen
   for i in range(len(a)):
    if seen>>i&1:new|=a[i]
   if new==seen:break
   seen=new
  if seen!=(1<<len(a))-1:return False
 return True
A=mm(S[1],S[0]);F=sub(T[1],A,3);I=[[int(i==j) for j in range(6)] for i in range(6)];zero=[[0]*6 for _ in range(6)]
require(mm(mm(A,sub(A,I,4)),mm(sub(A,I),sub(A,I)))==zero,'annihilator')
P=I
for h in range(1,5):
 P=mm(P,A);require(mm(mm(F,P),F)==zero,'separated defect zero product')
for r in range(1,51):require(val([0,1]*r,[])==(4**r+6*r-1)//9,'base scalar')
literal=0;transfers=0;fibers=[]
for r in range(3,19):
 J=list(range(1,r-1,2));H=(4**r+6*r-1)//9
 for t in range(1,len(J)+1):
  ds=set();count=0
  for P in combinations(J,t):
   D=3**t*H+3**(t-1)*sum(4**j-2 for j in P);pos=[2*j for j in P]
   require(val([0,1]*r,pos)==D,'multi additive scalar');transfers+=1;ds.add(D);count+=1
   if r<=6:
    a=gmulti([0,1]*r,pos);x=profile(a,0)
    require(x==(D,0,3**(t-1)*(4**r+2),3**t),'literal multi profile')
    b=induced(a,[0]+list(range(2,len(a)-1)));y=profile(b,0)
    require(strong(b),'multi strongness')
    require(y==(D,3**t,3**(t-1)*(2**(2*r-1)+1),3**t),'literal strong profile');literal+=2
  require(len(ds)==count,'scalar injection')
  if r in (6,10,18):fibers.append({'r':r,'t':t,'cardinality':count,'minimum_D':min(ds),'maximum_D':max(ds)})
# Adjacent interactions are real, and prevent a naive unrestricted extension.
bits=[0,1,0,1];adj=profile(gmulti(bits,[0,2]),0)[0]
require(adj==41,'adjacent hostile literal')
require(3*val(bits,[0])+3*val(bits,[2])-9*val(bits,[])==39,'false additive prediction')
files=['MULTI.md','verify_multi.py','STRONG.md','MEMBERSHIP.md','PROOF.md','transfer.py','inherited.py']
print(json.dumps({'status':'PASS','literal_profiles':literal,'separated_transfer_controls':transfers,'base_controls':50,'selected_fibers':fibers,'adjacent_hostile':{'actual':41,'false_additive_prediction':39},'sha256':{f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}},sort_keys=True,indent=2))
