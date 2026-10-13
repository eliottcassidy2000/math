#!/usr/bin/env python3
"""Finite Dirichlet numerator and conditional tail audit; NUMERICAL, not intervals."""
import argparse,json,math
import numpy as np
from collatz_parity_modes_20261013b import resonance

def main():
 p=argparse.ArgumentParser();p.add_argument('--bits',required=True);p.add_argument('--exponent',type=int,default=26);p.add_argument('--output',required=True);args=p.parse_args()
 bits=np.fromfile(args.bits,dtype=np.uint8);N=2**args.exponent
 t=(1-2*np.unpackbits(bits,bitorder='little').astype(np.int16))[:N+1]
 reports=[]
 for f in [2.5,14.5,26.5,79.5,332.5]:
  m=resonance(f);s=1+complex(-m['decay'],2*np.pi*m['frequency'])/np.log(2)
  den=1+2**(-s)+(1.5**s)/3;dp=-np.log(2)*2**(-s)+np.log(1.5)*1.5**s/3
  total=0j;rows=[]
  # Direct original numerator E_N+R_N = sum_{n=2(3)}t(n)(n-1/2)^-s - F_N/3.
  # Independent shell formula uses only positive-index tail differences in F.
  cumulative={0:0j};endpoints={2**k for k in range(4,args.exponent+1)}
  for n in list(endpoints):endpoints.update([(2*n-1)//3,((2*n-1)//3)//2])
  for n in sorted(endpoints):
   prev=max(cumulative)
   val=cumulative[prev]
   for lo in range(prev+1,n+1,2**19):
    hi=min(n+1,lo+2**19);x=np.arange(lo,hi,dtype=np.float64);sg=t[lo:hi]
    z=np.exp(-s*np.log(x));val+=np.sum(sg*z)
   cumulative[n]=val
  direct=0j;prev=0
  for k in range(4,args.exponent+1):
   end=2**k
   for lo in range(prev+1,end+1,2**19):
    hi=min(end+1,lo+2**19);x=np.arange(lo,hi,dtype=np.float64);sg=t[lo:hi]
    z=-np.exp(-s*np.log(x))/3
    pick=(x.astype(np.int64)%3==2);z[pick]+=np.exp(-s*np.log(x[pick]-.5))
    direct+=np.sum(sg*z)
   L=(2*end-1)//3
   shell=(2/3)**s*((cumulative[end]-cumulative[L])+2**(-s)*(cumulative[end]-cumulative[L//2]))
   err=abs(shell-direct)
   if err>1e-8:raise RuntimeError(('shell identity',f,k,err))
   residue=-1.5**s*direct/dp
   amp=2*abs(residue)*2**((s.real-1)*29)
   # Conditional ONLY: |H(x)| <= C sqrt(x) for every real x >= end.
   # Each listed radius is per unit C, plus the unconditional R tail.
   beta=.5
   e_radius=(end**(beta-s.real)/3)*(1+abs(s)/(s.real-beta)) if s.real>beta else None
   r_radius=abs(s)/(2*s.real)*(1-1/(2*end))**(-s.real-1)*end**(-s.real)
   factor=2*abs(1.5**s/dp)*2**((s.real-1)*29)
   rows.append(dict(k=k,numerator=[direct.real,direct.imag],shell_error=err,amplitude29=amp,phase=float(np.angle(residue)),tail_amplitude_per_C=None if e_radius is None else factor*e_radius,R_tail_amplitude=factor*r_radius))
   prev=end
  reports.append(dict(base=f,s=[s.real,s.imag],root_residual=abs(den),rows=rows))
 with open(args.output,'w')as o:json.dump(dict(status='NUMERICAL; conditional tails are not established bounds for the infinite sequence',reports=reports),o,indent=2)
 for r in reports:
  print('base',r['base'])
  for row in r['rows']:
   if row['k'] in [10,14,18,22,args.exponent]:print(row)
if __name__=='__main__':main()
