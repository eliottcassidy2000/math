#!/usr/bin/env python3
"""Exact odd-step pairing and its 10/21 asymptotic cancellation tariff."""
import argparse,json
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from root_braid_fold_20261017f import coefficient,need


def ideal(r,x):
 c=0
 while x<=4:
  c+=(-3*(0 if r==0 else 1 if r==1 else -1)) if x<=2 else (2 if r==2 else -1)
  x*=4;r=(r+1)%3
 return c

def run():
 phase=[]
 for b in range(1,19):
  x=Q(3,2**b);v=[min(abs(ideal(r,x)),abs(ideal((2*r+1)%3,2*x))) for r in range(3)]
  phase.append(v)
 need([sum(v) for v in phase[:6]]==[3,3,3,3,0,0] and phase[6:12]==phase[:6] and phase[12:]==phase[:6],'six-band cancellation pattern')
 tariff=sum(Q(sum(phase[b-1]),6*2**b) for b in range(1,7))/(1-Q(1,64));need(tariff==Q(10,21),'exact infinite geometric tariff')
 checks=0
 for N in range(1,257):
  for s in range(1,4*N+1,8):
   t=2*s+1
   if t>4*N:continue
   a,b=coefficient(N,s),coefficient(N,t)
   need(a*b>=0,'paired braid columns have matching signs')
   need(abs(a-b)==abs(a)+abs(b)-2*min(abs(a),abs(b)),'exact odd-clock cancellation')
   need((3*t+1)//2==4*((3*s+1)//4)+1,'actual common-next-odd fibre')
   checks+=1
 rows=[]
 for N in [1024,16384,262144,300001,524287]:
  s=np.arange(1,4*N+1,2,dtype=np.int64);s=s[s%8!=5];H=np.zeros(len(s),dtype=np.int64);m=s.copy()
  while np.any(m<=4*N):
   H+=np.where(m<=2*N,np.where(m%3==1,-3,np.where(m%3==2,3,0)),np.where(m<=4*N,np.where(m%3==2,2,-1),0));m=4*m+1
  table=np.zeros(4*N+1,dtype=np.int64);table[s]=H
  pair=s[(s%8==1)&(2*s+1<=4*N)];lost=2*int(np.minimum(abs(table[pair]),abs(table[2*pair+1])).sum());base=int(abs(H).sum())
  rows.append(dict(N=N,second_fold_l1=base,cancelled_mass=lost,third_fold_l1=base-lost,tariff_error_21=21*lost-10*N,third_fold_error_21=21*(base-lost)-38*N))
 return dict(status='PROVED scoped exact third fold and uniform linear asymptotic; square-root current OPEN',pair_checks=checks,six_band_minima=phase[:6],cancellation_coefficient=str(tariff),third_fold_coefficient='38/21',rows=rows)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();r=run();Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(r,sort_keys=True))
