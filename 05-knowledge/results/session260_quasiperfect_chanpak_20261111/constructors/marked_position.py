#!/usr/bin/env python3
"""Exact two-module path-word compiler; unbounded marked-position theorem."""
from collections import defaultdict
from math import comb,factorial
from pathlib import Path
import json
from module_ladders import inflate,data
HERE=Path(__file__).resolve().parent

def mul(A,B):
 C=defaultdict(int)
 for (r,s),x in A.items():
  for (u,v),y in B.items():C[r+u,s+v]+=x*y
 return dict(C)
def add(A,B,scale=1):
 for key,x in B.items():A[key]+=scale*x

def compile_core(T,O,v):
 # Fixed original vertices, including distinguished selected v.
 # P precedes v; S follows v; all three share exterior neighbourhood.
 k=len(T);P=k;S=k+1;adj=list(T)+[0,0]
 for a in range(k):
  if a==v:adj[a]|=1<<S
  elif T[a]>>v&1:adj[a]|=(1<<P)|(1<<S)
 for a in range(k):
  if a==v:adj[P]|=1<<a
  elif T[v]>>a&1:adj[P]|=1<<a;adj[S]|=1<<a
 adj[P]|=1<<S
 # All directed path words. Reusable P/S symbols never repeat consecutively;
 # P->S is the only symbol-symbol edge. Between fixed vertices there are <=2.
 paths=defaultdict(lambda:defaultdict(int));front=defaultdict(int)
 for a in range(k+2):
  mask=1<<a if a<k else 0;r=int(a==P);s=int(a==S)
  front[mask,a,r,s,int(a<k and O>>a&1)]+=1
 while front:
  nex=defaultdict(int)
  for (mask,a,r,s,first),x in front.items():
   paths[mask,first,int(a<k and O>>a&1)][r,s]+=x
   for b in range(k+2):
    if not(adj[a]>>b&1) or b<k and mask>>b&1:continue
    nex[mask|(1<<b if b<k else 0),b,r+int(b==P),s+int(b==S),first]+=x
  front=nex
  if any(r+s>2*(k+1) for _,_,r,s,_ in front):raise RuntimeError('DAG bound')
 full=(1<<k)-1;H=defaultdict(int);W=defaultdict(int);U=defaultdict(int)
 for f in (0,1):
  for e in (0,1):add(H,paths[full,f,e])
 for A in range(1<<k):
  B=full^A;fst=defaultdict(int);end=defaultdict(int);start=defaultdict(int)
  for f in (0,1):
   for e in (0,1):
    add(fst,paths[A,f,e])
    if e:add(end,paths[A,f,e])
    if f:add(start,paths[B,f,e])
  if not A:fst[0,0]+=1 # first path may be empty in U only
  add(W,mul(end,start));add(U,mul(fst,start))
 D=defaultdict(int,H);add(D,W,-1)
 return {key:{rs:x for rs,x in poly.items() if x} for key,poly in [('H',H),('D',D),('U',U)]}

def exponential(poly):
 out=defaultdict(int)
 for (r,s),x in poly.items():
  for a in range(r+1):
   for b in range(s+1):out[a,b]+=x*(-1)**(r+s-a-b)*comb(r,a)*comb(s,b)
 return {ab:x for ab,x in out.items() if x}
def ev(poly,j,k):return sum(x*a**j*b**k for (a,b),x in poly.items())

def main():
 T=(24,1,3,6,14);out=[]
 for O,v in [(6,1),(24,3)]:
  raw=compile_core(T,O,v);expo={key:exponential(poly) for key,poly in raw.items()}
  checks=0
  for j in range(5):
   for k in range(5):
    B,_=inflate(T,0,v,j+k+1);old=[a for a in range(5) if a!=v]
    sel=sum(1<<i for i,a in enumerate(old) if O>>a&1)|(1<<(4+j))
    h,w,u=data(B,sel)
    if [ev(expo[key],j,k) for key in ['H','D','U']]!=[h,h-w,u]:raise RuntimeError('literal control')
    checks+=1
  out.append(dict(rows=T,selection=O,vertex=v,exponential_coefficients={key:[[a,b,x] for (a,b),x in sorted(poly.items())] for key,poly in expo.items()},literal_controls=checks))
 receipt={'status':'PASS','scope':'Symbolic finite path-word compiler, all prefix/suffix lengths; literal controls are additional','families':out}
 (HERE/'marked_position_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
