#!/usr/bin/env python3
"""Literal constructors for N(9,r,s), all r>=1 and s in {1,2,3}."""
from bilinear_kernel_certificate import substitute
from module_tools import *
import argparse,json
HERE=Path(__file__).resolve().parent

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def join(A,B):return tuple(row|(((1<<len(B))-1)<<len(A)) for row in A)+tuple(row<<len(A) for row in B)
def tail(B,O,t):
 n=len(B);out=list(B)+[O]+[0]*t
 for i in range(n):
  out[i]|=((1<<t)-1)<<(n+1)
  if not O>>i&1:out[i]|=1<<n
 for j in range(t):out[n+1+j]=(1<<n)|sum(1<<(n+1+k) for k in range(j+1,t))
 return tuple(out)
def realize(r,s):
 if r<1 or s not in (1,2,3):raise ValueError('require r>=1 and s=1,2,3')
 C=(2,4,1)
 if s==1:return join(C,tail(C,3,r))
 if s==2:return (32,9,3,69,15,30,55) if r==1 else substitute((0,),tt(r-1))
 return join(join(C,C),tail((0,),1,r+3))
def main():
 p=argparse.ArgumentParser();p.add_argument('--r',type=int);p.add_argument('--s',type=int);a=p.parse_args()
 if a.r is not None:
  if a.s is None:p.error('--s is required with --r')
  print(json.dumps({'r':a.r,'s':a.s,'rows':realize(a.r,a.s)}));return
 out=[]
 for r in range(1,5):
  for s in (1,2,3):
   T=realize(r,s);h=sum(m.endpoints(T)[-1]);target=1+2**s+9*2**(r+s)
   if h!=target:raise RuntimeError(('literal grid mismatch',r,s,h,target))
   out.append({'r':r,'s':s,'order':len(T),'H':h,'rows':T})
 result={'status':'PASS','scope':'All-r proofs in INCOMING_CONNECTIONS.md; twelve literal controls, no claim for s>=4','controls':out}
 HERE.joinpath('realize9_lower_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
