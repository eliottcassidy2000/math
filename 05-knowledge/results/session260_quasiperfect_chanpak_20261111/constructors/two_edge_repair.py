#!/usr/bin/env python3
"""Apex all-cut gap theorem controls and symbolic four-edge rank-template."""
import json
from collections import defaultdict
from math import comb
from pathlib import Path
from module_ladders import endpoints
HERE=Path(__file__).resolve().parent

def base(a):
 n=2*a+2;adj=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if j==n-1 and i<a:adj[j]|=1<<i
   else:adj[i]|=1<<j
 return adj

def cut_h(mask,n):
 h=1;weighted=prev=0
 for j in range(n):
  bit=mask>>j&1
  if not bit:h+=weighted+prev
  weighted=2*weighted+bit;prev=bit
 return h

def template_coefficients():
 # Fixed vertices (Afirst,Apenultimate,Alast,Blast,apex), with two
 # repeatable independent run symbols P=remaining A, S=B except last.
 # Their only mutual arc P->S gives a finite word language.
 k=5;P=5;S=6
 adj=[sum(1<<b for b in z) for z in [[1,2,3,4,5,6],[2,3,4,6],[4,6],[2,4],[5],[1,2,3,6],[3,4]]]
 front=defaultdict(int)
 for a in range(7):front[1<<a if a<k else 0,a,int(a==P),int(a==S)]+=1
 H=defaultdict(int)
 while front:
  nex=defaultdict(int)
  for (mask,a,r,s),x in front.items():
   if mask==31:H[r,s]+=x
   for b in range(7):
    if adj[a]>>b&1 and not(b<k and mask>>b&1):nex[mask|(1<<b if b<k else 0),b,r+int(b==P),s+int(b==S)]+=x
  front=nex
 E=defaultdict(int)
 for (r,s),x in H.items():
  for a in range(r+1):
   for b in range(s+1):E[a,b]+=x*(-1)**(r+s-a-b)*comb(r,a)*comb(s,b)
 return {ab:x for ab,x in E.items() if x}

def main():
 apex=[]
 for a in range(2,9):
  n=2*a+1;vals=[cut_h(mask,n) for mask in range(1<<n)];levels=sorted(set(vals),reverse=True)
  high=(2**a-1)*(2**(a+1)-1)+2;gap=3*2**(a-1)-2
  if levels[:2]!=[high,high-gap]:raise RuntimeError('all-apex spectrum top gap')
  apex.append(dict(a=a,cuts=1<<n,maximum=high,second_maximum=levels[1],gap=gap,target_hits=vals.count(high-2)))
 E=template_coefficients();cases=[]
 for a in range(3,8):
  T=base(a);r=2*a+1;flips=[(0,r),(a-2,r),(a-1,2*a),(a-1,r)]
  for i,j in flips:T[i]^=1<<j;T[j]^=1<<i
  actual=sum(endpoints(T)[-1]);symbolic=sum(x*b**(a-3)*c**a for (b,c),x in E.items())
  if actual!=symbolic:raise RuntimeError('template graph interface')
  target=2**a*(2**(a+1)-3)+1
  cases.append(dict(a=a,flips=flips,H=actual,target=target,excess=actual-target,rows=T if a==4 else None))
 out=dict(status='PASS',apex_cut_controls=apex,template_exponential_coefficients=[[b,c,x] for (b,c),x in sorted(E.items())],template_cases=cases)
 (HERE/'two_edge_structure_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
