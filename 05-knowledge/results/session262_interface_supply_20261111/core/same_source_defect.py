#!/usr/bin/env python3
"""Fixed labelled domination graph and base count, two exact defect branches."""
import json
from pathlib import Path
from response_kernel import Q,substitute,tt,born,ham

def join(A,B):
 n=len(A);m=len(B)
 return tuple(x|(((1<<m)-1)<<n) for x in A)+tuple(x<<n for x in B)
def tail(T,t):return born(join(T,tt(t)))
def values(n):
 x=2**n;z=3**n
 K=111*x*x//2-9*x+27
 ua=15*x*x+(682*x*z-44*z)//9+45*x//2-24
 ub=15*x*x+93*x*z+63*x-12*z-24
 return K,ua,ub,21*x,21*x//2

def main():
 out=[]
 for n in (3,4):
  A=tt(n);B=join((2,4,1),tt(n-3));TA,TB=substitute(A,B),substitute(B,A)
  K,ua,ub,ea,eb=values(n)
  if TA[0]!=TB[0] or TA[0].bit_count()!=n+3 or ham(TA)!=K or ham(TB)!=K:raise RuntimeError('source guard')
  for t in (0,1):
   ha,hb=ham(tail(TA,t)),ham(tail(TB,t))
   aa=K+2**t*ua-(ea if t==0 else 0);bb=K+2**t*ub-(eb if t==0 else 0)
   if (ha,hb)!=(aa,bb):raise RuntimeError('response')
   out.append({'n':n,'t':t,'base_H':K,'domination_count':1+2**(n+3+t),'H_AB':ha,'H_BA':hb,'gap':hb-ha})
 receipt={'status':'PASS','scope':'Four new two-branch literal controls; all-n proof is algebraic profile formula','controls':out}
 Path(__file__).with_name('same_source_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
