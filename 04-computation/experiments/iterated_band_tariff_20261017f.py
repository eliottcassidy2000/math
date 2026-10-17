#!/usr/bin/env python3
"""Exact finite-depth signed transport bands and ternary fibre controls.

The period and fibre balance are proved inductively in the companion note.
This computes the rational leading l1 constants, not a uniform-depth bound.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from third_fold_20261017f import ideal
from root_braid_fold_20261017f import coefficient,seed,need


def tables(D):
 allrows=[[[ideal(r,Q(3,2**b)) for r in range(3)] for b in range(6)]]
 results=[]
 for d in range(D+1):
  mod=3**(d+1);period=2*mod
  if d:
   previous=allrows[-1];oldmod=mod//3;oldperiod=2*oldmod;rows=[]
   for b in range(period):
    row=[]
    for r in range(mod):
     out=0;j=0;rho=r
     while 2*j<=b:
      if rho%3:
       a=2 if rho%3==1 else 1;bb=b-(a-1)-2*j
       if bb>=0:
        s=((2**a*rho-1)//3)%oldmod
        out+=(-1)**a*previous[bb%oldperiod][s]
      rho=(4*rho+1)%mod;j+=1
     row.append(out)
    rows.append(row)
   allrows.append(rows)
  rows=allrows[d]
  need(all(sum(row[r+k*(mod//3)] for k in range(3))==0 for row in rows for r in range(mod//3)),'last ternary digit has exactly zero signed sum')
  # Abel summation on each three-lift fibre: polynomial divisible by 1-z^(mod/3).
  for row in rows:
   q=mod//3;quot=[row[r] for r in range(q)]+[row[r]+row[r+q] for r in range(q)]
   rec=quot+[0]*q
   for r,c in enumerate(quot):rec[r+q]-=c
   need(rec==row,'literal conductor-factor reconstruction')
  tariff=sum(Q(sum(map(abs,row)),2**b*4*2**d) for b,row in enumerate(rows))/(1-Q(1,2**period))
  need(tariff>0,'positive envelope')
  if d:need(tariff<=Q(results[-1]['tariff']),'one-step signed transport l1 monotonicity')
  results.append(dict(depth=d,residue_modulus=mod,binary_band_period=period,tariff=str(tariff),tariff_decimal=float(tariff),
    maximum_absolute_coefficient=max(abs(x) for row in rows for x in row),coarse_fourier_modes_annihilated=mod//3,signed_bands=rows))
 return allrows,results


def finite(N,D):
 weights=defaultdict(int)
 for s in range(1,4*N+1,2):
  if s%8!=5:weights[s]+=coefficient(N,s)
 masses=[sum(map(abs,weights.values()))]
 for d in range(D):
  nxt=defaultdict(int)
  for s,c in weights.items():
   a=2 if s%8==1 else 1;nxt[seed((3*s+1)//2**a)]+=(-1)**a*c
  weights=nxt;masses.append(sum(map(abs,weights.values())))
 return masses


def run(D):
 rows,result=tables(D)
 need(Q(result[0]['tariff'])==Q(16,7) and Q(result[1]['tariff'])==Q(155896,87381),'independent leading constants')
 # Evaluate a second full period directly from the unwrapped recurrence.
 for d in range(1,D+1):
  mod=3**(d+1);oldmod=mod//3;period=2*mod;oldperiod=2*oldmod
  for b in range(period,2*period):
   for r in range(mod):
    rho=r;out=0;j=0
    while 2*j<=b:
     if rho%3:
      a=2 if rho%3==1 else 1;bb=b-(a-1)-2*j
      if bb>=0:out+=(-1)**a*rows[d-1][bb%oldperiod][((2**a*rho-1)//3)%oldmod]
     j+=1;rho=(4*rho+1)%mod
    need(out==rows[d][b%period][r],'unwrapped second period')
 samples=[]
 for N in [1024,16384,262144]:
  mass=finite(N,D);samples.append(dict(N=N,l1_by_depth=mass,leading_error=[str(Q(v)-N*Q(row['tariff'])) for v,row in zip(mass,result)]))
 return dict(status='PROVED fixed-depth recurrence and tariffs; uniform-in-depth control OPEN',depths=result,samples=samples)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--depth',type=int,default=3);p.add_argument('--output',required=True);a=p.parse_args();r=run(a.depth)
 Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps(dict(depths=[{k:v for k,v in d.items() if k!='signed_bands'} for d in r['depths']],samples=r['samples']),sort_keys=True))
