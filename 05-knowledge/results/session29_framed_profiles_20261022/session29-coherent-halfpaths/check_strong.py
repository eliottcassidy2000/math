#!/usr/bin/env python3
"""Literal reachability checks for strong binary realizations; no exhaustive census."""
import json
from pathlib import Path
from check_blue import binary
from check_model import tournament,hp

def strong(T):
 n=len(T);R=[row|1<<i for i,row in enumerate(T)]
 for k in range(n):
  for i in range(n):
   if R[i]>>k&1:R[i]|=R[k]
 return all(row==(1<<n)-1 for row in R)

def main():
 checked=[]
 for h in range(3,128,2):
  m,q,lab=binary(h)
  for odd in (False,True):
   T,s=tournament(m,q,lab,odd)
   assert strong(T) and hp.dp(T,s)['fixed_hp']==h
   checked.append({'h':h,'n':len(T)})
 T,s=tournament(1,0,0,True)
 T[0]^=1<<1;T[1]^=1<<0
 assert strong(T) and hp.dp(T,s)['fixed_hp']==1
 checked.append({'h':1,'n':3})
 result={'strong_graph_checks':len(checked),'cases':checked,'h1_C3':{'outmasks':T,'sigma':s},'scope':'Constructed examples only; no exhaustive enumeration beyond five pairs.'}
 Path(__file__).with_name('strong_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Passed',len(checked),'literal strong-graph and fixed-path checks.')
if __name__=='__main__':main()
