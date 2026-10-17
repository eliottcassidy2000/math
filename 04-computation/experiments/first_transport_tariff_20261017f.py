#!/usr/bin/env python3
"""Exact 18-band tariff for the full first odd transport after braid folding."""
import argparse,json
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from third_fold_20261017f import ideal
from root_braid_fold_20261017f import coefficient,seed,need


def C(r,x):
 out=0;j=0;rho=r
 while 2*4**j*x/3<=4:
  if rho%3:
   a=2 if rho%3==1 else 1;s=((2**a*rho-1)//3)%3
   out+=(-1)**a*ideal(s,Q(2**a*4**j,3)*x)
  j+=1;rho=(4*rho+1)%9
 return out

def inverse(N,t):
 j=0;m=t;out=0
 while (2*m-1)//3<=4*N:
  if m%3:
   a=2 if m%3==1 else 1;s=(2**a*m-1)//3
   if s<=4*N:out+=(-1)**a*coefficient(N,s)
  j+=1;m=4*m+1
 return out

def run():
 rows=[[C(r,Q(9,2**(b+1))) for r in range(9)] for b in range(108)]
 need(all(rows[b]==rows[b%18] for b in range(108)),'six full periods of the exact band formula')
 need(all(sum(ideal(r,Q(a,b)) for r in range(3))==0 for a in range(1,9) for b in range(2,17)),'zero residue sum in every ideal column')
 for r in range(9):
  rr=r
  for j in range(9):rr=(4*rr+1)%9
  need(rr==r,'affine braid has order9 modulo9')
 tariff=sum(Q(sum(map(abs,rows[b])),8*2**b) for b in range(18))/(1-Q(1,2**18))
 need(tariff==Q(155896,87381),'full first-transport tariff')
 cases=0
 for N in range(1,257):
  forward=defaultdict(int)
  for s in range(1,4*N+1,2):
   if s%8==5:continue
   a=2 if s%8==1 else 1;forward[seed((3*s+1)//2**a)]+=(-1)**a*coefficient(N,s)
  for t in range(1,6*N+2,2):
   if t%8==5:continue
   need(inverse(N,t)==forward[t],'independent inverse-fibre and forward coefficient agreement');cases+=1
 native=json.loads(Path('04-computation/experiments/joint_logistic_current_20261017f/transport.json').read_text());errors=[]
 for row in native['runs']:
  stage=next(r for r in row['rows'] if r['depth']==1);mass=stage['residual_l1']+abs(stage['root_current'])
  errors.append(dict(N=row['N'],first_transport_l1=mass,error_numerator=87381*mass-155896*row['N']))
 return dict(status='PROVED scoped first-transport asymptotic; uniform iterated bound OPEN',band_coefficients=rows[:18],band_abs_sums=[sum(map(abs,r)) for r in rows[:18]],
   full_transport_coefficient=str(tariff),period=18,inverse_coefficient_checks=cases,denominator_factorization={'3':2,'7':1,'19':1,'73':1},finite_rows=errors)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();r=run();Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(r,sort_keys=True))
