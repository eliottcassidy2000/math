#!/usr/bin/env python3
"""Exact integer short-period and structural-reset replay."""
import hashlib
import math
import sys
from itertools import product
from pathlib import Path
from random import Random
import numpy as np

D=62
EXPECTED={
"0.2":(3,4),"0.12":(4,4),"0.123":(5,5),"1.2":(4,4),"2.13":(5,5),
"0.0.3":(3,2,3),"0.0.13":(2,3,2),"0.0.123":(3,3,3),
"0.1.3":(3,4,3),"0.1.23":(2,3,2),"0.1.123":(3,4,3),
"0.12.3":(3,3,2),"0.3.1":(3,2,3),"0.3.12":(3,4,4),
"0.13.2":(3,4,4),"0.13.12":(2,3,3),"1.2.3":(4,3,3),"1.12.3":(4,4,3),
}

def read(path):
 raw=Path(path).read_bytes(); it=iter(map(int,raw.split()))
 assert next(it)==D and next(it)==8
 init=np.array([next(it) for _ in range(D)],dtype=object)
 fin=np.array([next(it) for _ in range(D)],dtype=object)
 M={}
 for _ in range(8):
  k=next(it);M[k]=np.array([next(it) for _ in range(D*D)],dtype=np.int64).reshape(D,D)
 assert set(M)==set(range(8)) and next(it,None) is None
 assert max(int(np.max(np.abs(A))) for A in M.values())<=5
 return init,fin,M,hashlib.sha256(raw).hexdigest()

def matrix(word,M,integer=object):
 A=np.eye(D,dtype=integer)
 for k in word:A=M[k]@A
 return A

def rank_one(M):
 nz=np.argwhere(M)
 if not len(nz):return False
 i,j=nz[0]
 return bool((M[i,j]*M==M[:,j:j+1]*M[i:i+1,:]).all())

def fmt(word):
 return '.'.join(''.join(str(j+1) for j in range(3) if k>>j&1) or '0' for k in word)

def phases(M):
 out={};bystep={};checked=0
 for k,lam in [(2,3),(3,5)]:
  for U in product(range(8),repeat=k):
   if U!=min(U[j:]+U[:j] for j in range(k)):continue
   checked+=1
   # At most 3 short matrix factors, entries bounded by 5:
   # |entry|<=62^2*5^3, safe before switching to Python integers.
   B=matrix(U,M,np.int64).astype(object)
   if np.trace(B)!=lam:continue
   P=B.copy()
   for e in range(1,8):
    if np.trace(P)!=lam**e:
     bystep[fmt(U)]=e;break
    PN=B@P
    if rank_one(P) and (PN==lam*P).all():
     es=[]
     for j in range(k):
      A=matrix(U[j:]+U[:j],M);Q=A.copy()
      for r in range(1,8):
       R=A@Q
       if rank_one(Q) and (R==lam*Q).all():es.append(r);break
       Q=R
      else:raise AssertionError(("unresolved phase",U,j))
     out[fmt(U)]=tuple(es);break
    P=PN
   else:raise AssertionError(("unresolved class", U))
 assert out==EXPECTED,(out,EXPECTED)
 assert len(bystep)==7 and checked==212
 print("checked classes",checked,"trace-power rejects passing first trace",bystep)
 for word,es in out.items():print(word,es)

def encode(vertices, arc):
 word=[]
 for j,v in enumerate(vertices):
  k=0
  for i,u in enumerate(vertices[:j]):
   if arc(v,u):
    if j-i>3:return None
    k|=1<<(j-i-1)
  word.append(k)
 return word

def physical(init,fin,M):
 rng=Random(4880);tests=0;multiple=0
 def H(word):
  v=init.copy()
  for a in word:v=M[a]@v
  return int(fin@v)
 for links in ([4]*7,[3,*([4]*6)],[3,3,*([4]*5)],[5,*([4]*5)],[4]*10):
  marks=[0]
  for d in links:marks.append(marks[-1]+d-1)
  end=marks[-1];exceptions=set(zip(marks,marks[1:]))
  def inner(a,b):return (a<b) != ((min(a,b),max(a,b)) in exceptions)
  pos=list(range(end+1))
  if links[0]==5:pos[0],pos[1]=pos[1],pos[0]
  assert encode(pos,inner) is not None
  for times in range(5):
   for _ in range(40):
    i=rng.randrange(end);pos[i],pos[i+1]=pos[i+1],pos[i]
    if encode(pos,inner) is None:pos[i],pos[i+1]=pos[i+1],pos[i]
   P=rng.randrange(5);S=rng.randrange(5)
   ids=[*range(-P,0),*pos,*range(end+1,end+1+S)]
   rank={v:k for k,v in enumerate(ids)}
   flipped={(u,v) for j,v in enumerate(ids) for u in ids[max(0,j-3):j]
       if not(0<=u<=end and 0<=v<=end) and rng.randrange(2)}
   def arc(a,b):
    if 0<=a<=end and 0<=b<=end:return inner(a,b)
    i,j=rank[a],rank[b]
    return (i<j)!=((a,b) in flipped or (b,a) in flipped)
   whole=encode(ids,arc);uword=whole[P:P+end+1]
   floor=math.prod((1<<(d-2))+1 for d in links)
   assert H(uword)==floor and rank_one(matrix(uword,M))
   safe=[t for t in marks if 8<=t<=end-8]
   t=safe[0]
   assert H(whole)==H(encode([v for v in ids if v<=t],arc))*H(encode([v for v in ids if v>=t],arc))
   e=max(0,(end+1-27+2)//3)
   assert H(whole)%5**e==0
   multiple+=e>0;tests+=1
 print("ambient structural resets",tests,"nontrivial certified 5-power",multiple)

def main():
 init,fin,M,sha=read(sys.argv[1])
 print("input SHA256",sha)
 phases(M)
 physical(init,fin,M)

if __name__=="__main__":main()
