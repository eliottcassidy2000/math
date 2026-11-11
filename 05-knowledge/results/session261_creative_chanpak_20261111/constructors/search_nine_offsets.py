#!/usr/bin/env python3
from module_tools import *
import json
HERE=Path(__file__).resolve().parent

def data_all(T):
 n=len(T);e=m.endpoints(T);h=sum(e[-1]);R=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n));s=m.endpoints(R)
 F=(1<<n)-1;K=[[0]*n for _ in range(n)];U=[0]*n
 for A in range(1<<n):
  B=F^A
  if not B:continue
  hh=sum(e[A]) if A else 1
  for j in range(n):
   if not s[B][j]:continue
   U[j]+=hh*s[B][j]
   for i in range(n):K[i][j]+=e[A][i]*s[B][j]
 us=[0]*(1<<n);ws=[0]*(1<<n)
 for O in range(1,1<<n):
  bit=O&-O;i=bit.bit_length()-1;rest=O^bit;us[O]=us[rest]+U[i];ws[O]=ws[rest]+sum(K[i][j]+K[j][i] for j in range(n) if rest>>j&1)
 return h,ws,us

def main():
 seen5={canon(tour(5,b)) for b in range(1024)};six=set()
 for A in seen5:
  for b in range(32):
   T=tuple(A[i]|((1<<5) if b>>i&1 else 0) for i in range(5))+(sum(1<<i for i in range(5) if not b>>i&1),)
   six.add(canon(T))
 best={};hits={};tested=0
 for order in (6,7):
  for A in sorted(six):
   for bits in range(1 if order==6 else 64):
    T=A if order==6 else tuple(A[i]|((1<<6) if bits>>i&1 else 0) for i in range(6))+(sum(1<<i for i in range(6) if not bits>>i&1),)
    tested+=1;h,w,u=data_all(T)
    for O in range(1,1<<order):
     D=h-w[O];q=u[O]
     if D<=1 or (D-1)&(D-2):continue
     ss=(D-1).bit_length()-1
     if ss<2 or not q:continue
     j=(q&-q).bit_length()-1
     if q>>j!=9:continue
     z={'s':ss,'rows':T,'O':O,'H':h,'W':w[O],'D':D,'U':q,'j':j};hits[ss]=hits.get(ss,0)+1
     key=(order,j,T,O)
     if ss not in best or key<best[ss][0]:best[ss]=(key,z);print(z,flush=True)
 out={'status':'PASS','scope':'Complete order6 classes and every one-vertex extension, all selected subsets; complete through order7, not an all-order exclusion','six_class_count':len(six),'presentations_tested':tested,'hit_presentations':hits,'best':{s:z for s,(key,z) in best.items()}}
 HERE.joinpath('nine_offsets_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('END',json.dumps(out),flush=True)
if __name__=='__main__':main()
