#!/usr/bin/env python3
"""Emit the literal three-parameter exchange-source family."""
import argparse,json
from response_kernel import substitute,tt,born,ham

def realize(a,b,t):
 if min(a,b)<1 or t<0:raise ValueError('a,b>=1 and t>=0 required')
 T=substitute(tt(a),tt(b));n=len(T)
 tail=tt(t)
 T=tuple(row|(((1<<t)-1)<<n) for row in T)+tuple(row<<n for row in tail)
 return born(T)

def formula(a,b,t):
 K=37*2**(a+b)-2**(a+1)-2**(b+1)+9
 U=62*2**a*3**b+10*2**(a+b)-4*2**a+23*2**b-4*3**b-8
 return K+2**t*U-(7*2**a if t==0 else 0)

def main():
 p=argparse.ArgumentParser();p.add_argument('--a',type=int,default=1);p.add_argument('--b',type=int,default=2);p.add_argument('--tail',type=int,default=0);p.add_argument('--check',action='store_true');args=p.parse_args()
 T=realize(args.a,args.b,args.tail);out={'rows':T,'order':len(T),'H_by_proved_formula':formula(args.a,args.b,args.tail)}
 if args.check:
  out['H_by_literal_DP']=ham(T)
  if out['H_by_literal_DP']!=out['H_by_proved_formula']:raise RuntimeError('mismatch')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
