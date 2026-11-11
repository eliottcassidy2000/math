#!/usr/bin/env python3
"""Exact skeleton certificate for a freely marked transitive module."""
from collections import defaultdict
from math import comb
from itertools import product
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'constructors'))
from module_ladders import inflate,data
HERE=Path(__file__).resolve().parent

def add(A,B,scale=1):
 for key,x in B.items():A[key]+=scale*x
def mul(A,B):
 C=defaultdict(int)
 for r,x in A.items():
  for s,y in B.items():C[tuple(a+b for a,b in zip(r,s))]+=x*y
 return C

def compile(T,v):
 # Four external vertices, two fixed selected clones; three reusable gap types.
 fixed,_=inflate(T,0,v,2);k=len(fixed);d=3
 adj=list(fixed)+[0]*d
 for a in range(k-2):
  if fixed[a]>>(k-2)&1:adj[a]|=((1<<d)-1)<<k
 for i in range(d):
  adj[k+i]|=fixed[k-2]&((1<<(k-2))-1)
  for j in range(i+1,d):adj[k+i]|=1<<(k+j)
  for j in range(2):
   if i<=j:adj[k+i]|=1<<(k-2+j)
   else:adj[k-2+j]|=1<<(k+i)
 O=(1<<(k-2))|(1<<(k-1));zero=(0,)*d
 paths=defaultdict(lambda:defaultdict(int));front=defaultdict(int)
 for a in range(k+d):front[1<<a if a<k else 0,a,tuple(int(a==k+i) for i in range(d)),int(a<k and O>>a&1)]+=1
 layers=0
 while front:
  nex=defaultdict(int);layers+=1
  if layers>k+d*(k+1):raise RuntimeError('DAG bound')
  for (mask,a,r,first),x in front.items():
   paths[mask,first,int(a<k and O>>a&1)][r]+=x
   for b in range(k+d):
    if not(adj[a]>>b&1) or b<k and mask>>b&1:continue
    rr=tuple(r[i]+int(b==k+i) for i in range(d))
    nex[mask|(1<<b if b<k else 0),b,rr,first]+=x
  front=nex
 full=(1<<k)-1;H=defaultdict(int);W=defaultdict(int);U=defaultdict(int)
 for f,e in product((0,1),repeat=2):add(H,paths[full,f,e])
 for A in range(1<<k):
  B=full^A;fst=defaultdict(int);end=defaultdict(int);start=defaultdict(int)
  for f,e in product((0,1),repeat=2):
   add(fst,paths[A,f,e])
   if e:add(end,paths[A,f,e])
   if f:add(start,paths[B,f,e])
  if not A:fst[zero]+=1
  add(W,mul(end,start));add(U,mul(fst,start))
 D=defaultdict(int,H);add(D,W,-1)
 return {key:{r:x for r,x in p.items() if x} for key,p in [('H',H),('D',D),('U',U)]}

def exp(poly):
 out=defaultdict(int)
 for r,x in poly.items():
  for a in product(*(range(t+1) for t in r)):
   c=x
   for i in range(len(r)):c*=(-1)**(r[i]-a[i])*comb(r[i],a[i])
   out[a]+=c
 return {a:x for a,x in out.items() if x}
def ev(poly,lens):
 out=0
 for r,x in poly.items():
  for a,b in zip(r,lens):x*=a**b
  out+=x
 return out

def main():
 T=(24,1,3,6,14);allout=[]
 for v in (1,2,3,4):
  raw=compile(T,v);expo={key:exp(p) for key,p in raw.items()}
  C={1:12,2:18,3:20,4:30}[v];delta={1:0,2:1,3:0,4:3}[v]
  expected={'H':{(1,1,1):-1,(2,2,2):24},
    'D':{(1,0,1):1,(1,1,1):-1,(1,2,1):1,(2,1,2):-12,(2,2,2):24,(2,3,2):-12},
    'U':{(2,2,3):2*C,(2,3,3):3*C}}
  if delta:expected['U'][0,1,1]=-delta
  if expo!=expected:raise RuntimeError(('closed formula certificate',v,expo))
  checks=0
  for a,c,b in product(range(3),repeat=3):
   B,_=inflate(T,0,v,a+c+b+2);O=(1<<(4+a))|(1<<(5+a+c));h,w,u=data(B,O)
   actual=[ev(expo[key],(a,c,b)) for key in ('H','D','U')]
   if actual!=[h,h-w,u]:raise RuntimeError((v,a,c,b,actual,(h,w,u)))
   checks+=1
  allout.append({'core':T,'vertex':v,'exponential':{key:[[list(r),x] for r,x in sorted(p.items())] for key,p in expo.items()},'literal_checks':checks})
 out={'status':'PASS','scope':'Finite symbolic path skeleton certificate valid at every nonnegative gap length','families':allout}
 (HERE/'marked_module_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
