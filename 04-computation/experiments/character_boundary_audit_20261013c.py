#!/usr/bin/env python3
"""Exact group-algebra audit of terminal character forcing and finite flux."""
from fractions import Fraction as F
import argparse,json

def ck(x,m):
 if not x:raise RuntimeError(m)
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args();N=2**16
 tau=[0]*(N+1);omega=[0]*(N+1)
 for n in range(2,N+1):
  z=n;t=o=0
  while z!=1:
   if z%2:z=(3*z+1)//2;o+=1
   else:z//=2
   t+=1;ck(t<1000000,'strict path cap')
  tau[n]=t;omega[n]=o
 controls=0
 for m in range(2,9):
  for a in range(m):
   for b in range(m):
    ex=[(a*t+b*o)%m for t,o in zip(tau,omega)]
    for n in range(2,N//2):
     target=n//2 if n%2==0 else (3*n+1)//2
     ck(ex[n]==(ex[target]+a+b*(n%2))%m,'nonterminal character law');controls+=1
    ck(ex[1]==0 and (ex[2]+a+b)%m==(2*a+b)%m,'terminal character defect');controls+=1
    if m>5:continue
    def add(x,y,scale=F(1),shift=0):
     z=x[:]
     for i,v in enumerate(y):z[(i+shift)%m]+=scale*v
     return z
    def prefix(end,s):
     v=[F(0)]*m
     for n in range(1,end+1):v[ex[n]]+=F(1,n**s)
     return v
    for cutoff in [2,3,4,7,8,15,16,31,32,63,64,127]:
     L=(2*cutoff-1)//3;J=L//2
     for s in [1,2,3]:
      Fn=prefix(cutoff,s);Fl=prefix(L,s);Fj=prefix(J,s)
      shift2=[F(0)]*m
      for n in range(2,cutoff+1,3):shift2[ex[n]]+=F(2,2*n-1)**s
      # Fl - a*2^-s Fj = b*(3/2)^s shifted_sum + (1-a*b).
      lhs=add(Fl,Fj,-F(1,2**s),a)
      rhs=add([F(0)]*m,shift2,F(3,2)**s,a+b)
      rhs[0]+=1;rhs[(2*a+b)%m]-=1
      ck(lhs==rhs,'finite group-algebra identity');controls+=1
      # General finite identity equals the root-shell numerator plus D(s)Fn.
      C=add(shift2,Fn,-F(1,3))
      K=add([F(0)]*m,C,F(3,2)**s,a+b);K[0]+=1;K[(2*a+b)%m]-=1
      delta=add([F(0)]*m,Fn);delta=add(delta,Fn,-F(1,2**s),a);delta=add(delta,Fn,-F(3,2)**s/F(3),a+b)
      sh=add(Fl,Fn,-1);sh=add(sh,add(Fn,Fj,-1),F(1,2**s),a)
      ck(add(K,delta,-1)==sh,'finite shell including D term');controls+=1
 out=dict(status='FINITE-EXACT',strict_source_universe=N,character_moduli=list(range(2,9)),fraction_identity_moduli=list(range(2,6)),controls=controls,terminal_defect='1-zeta^(2 alpha+beta)',hostile_omega_parity_at_1=dict(chi_1=1,branch_prediction=-1,defect=2))
 with open(args.output,'w')as f:json.dump(out,f,indent=2)
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
