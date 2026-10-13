#!/usr/bin/env python3
"""Arb interval certificate for ONE CONDITIONAL Collatz Dirichlet pole.
Hypothesis (UNPROVED): |H(x)| <= sqrt(x)*log2(x)/4 for all x >= 2^34.
The data are exact signed block moments of the finite ROOT parity census.
"""
import argparse,json,hashlib,time
import numpy as np
import mpmath as mp
from flint import arb,acb,ctx

def need(x,msg):
 if not x:raise RuntimeError(msg)
def main():
 p=argparse.ArgumentParser();p.add_argument('--raw-moments');p.add_argument('--moments',required=True);p.add_argument('--bits');p.add_argument('--census',required=True);p.add_argument('--output',required=True);args=p.parse_args()
 N=2**34;B=2**17;controls=0
 if args.raw_moments:
  raw=np.fromfile(args.raw_moments,dtype='<i8').reshape(-1,2);need(len(raw)==N//4096,'moment population')
  r=raw.reshape(-1,B//4096,2);v=np.empty((len(r),2),np.int64);v[:,0]=r[:,:,0].sum(axis=1);v[:,1]=(r[:,:,1]+np.arange(B//4096)[None,:]*4096*r[:,:,0]).sum(axis=1)
  if args.bits:
   t=1-2*np.unpackbits(np.fromfile(args.bits,dtype=np.uint8),bitorder='little').astype(np.int64);t[0]=0
   for j in range(min(256,len(t)//B)):
    a=t[j*B:(j+1)*B];need(v[j,0]==a.sum() and v[j,1]==np.dot(np.arange(B,dtype=np.int64),a),'independent raw bit moment');controls+=1
  np.savez_compressed(args.moments,moments=v,N=np.int64(N),B=np.int64(B))
 else:
  d=np.load(args.moments);v=d['moments'];need(int(d['N'])==N and int(d['B'])==B,'stored universe')
 census=json.load(open(args.census));need(census['N']==N,'census population');need(int(v[:,0].sum())+1==census['rows'][-1]['S'],'moment total');hN=abs(census['rows'][-1]['residues'][2][0])
 ctx.prec=160;mp.mp.dps=85
 theta=mp.log(3,2);defect=mp.mpf('26.5')*theta-mp.nint(mp.mpf('26.5')*theta);w=0j
 for j in range(1,121):w=mp.findroot(lambda z:mp.exp(-z)+mp.exp((theta-1)*z+2j*mp.pi*defect*j/120)-2,w)
 zz=1+(w+2j*mp.pi*mp.mpf('26.5'))/mp.log(2)
 c=acb(str(zz.real),str(zz.imag));radius=arb('1e-30');l2=arb(2).log();l15=(arb(3)/2).log()
 def D(s):return 1+acb(2)**(-s)+acb(arb(3)/2)**s/3
 def Dp(s):return -l2*acb(2)**(-s)+l15*acb(arb(3)/2)**s/3
 M=l2*l2*arb(2)**(-c.real+radius)+l15*l15*(arb(3)/2)**(c.real+radius)/3
 need(D(c).abs_upper()+M*radius*radius/2 < Dp(c).abs_lower()*radius,'Rouche root disk')
 need(Dp(c).abs_lower()>M*radius,'simple root')
 s=acb(arb(str(zz.real),'1e-30'),arb(str(zz.imag),'1e-30'));sigma=s.real
 L=(2*N-1)//3;J=L//2;start_long=(J+1+B-1)//B;start_short=(L+1+B-1)//B
 longsum=acb(0);shortsum=acb(0);started=time.monotonic()
 for j in range(start_long,len(v)):
  lo=j*B;center=arb(2*lo+B-1)/2;sm=int(v[j,0]);mu1=arb(2*int(v[j,1])-(B-1)*sm)/2
  power=acb(center)**(-s);val=power*(sm-s*mu1/center);longsum+=val
  if j>=start_short:shortsum+=val
  if j%32768==0:print('certified blocks',j,'seconds',round(time.monotonic()-started,1),flush=True)
 # n=N has even Terras stopping parity (N=2^34).
 last=acb(N)**(-s);longsum+=last;shortsum+=last
 def quadrature_error(cut,start):
  lo=start*B
  partial=(lo-cut-1)*arb(cut+1)**(-sigma)
  sum2=arb(B)*(B*B-1)/12
  # Tail of decreasing positive function bounds every second derivative.
  derivative_sum=arb(lo)**(-sigma-2)+arb(lo)**(-sigma-1)/(B*(sigma+1))
  return partial+(s*(s+1)).abs_upper()*sum2*derivative_sum/2
 errlong=quadrature_error(J,start_long);errshort=quadrature_error(L,start_short)
 numerator=acb(arb(2)/3)**s*(shortsum+acb(2)**(-s)*longsum)
 error=(arb(2)/3)**sigma*(errshort+arb(2)**(-sigma)*errlong)
 numerator+=acb(arb(0,error.upper()),arb(0,error.upper()))
 residue=-acb(arb(3)/2)**s*numerator/Dp(s)
 # H(N)=26869 from the exact all-prefix census; independent formula crosschecked there.
 eta=sigma-arb(1)/2
 etail=arb(hN)/3*arb(N)**(-sigma)+s.abs_upper()/(12*l2)*arb(N)**(-eta)*(arb(N).log()/eta+1/(eta*eta))
 rtail=s.abs_upper()/(2*sigma)*(1-arb(1)/(2*N))**(-sigma-1)*arb(N)**(-sigma)
 rho=(arb(3)/2)**sigma/Dp(s).abs_lower()*(etail+rtail)
 scale=2*arb(2**29)**(sigma-1)
 low=(scale*(residue.abs_lower()-rho)).lower();high=(scale*(residue.abs_upper()+rho)).upper()
 need(low>arb('0.10'),'conditional residue is nonzero with margin')
 etail_power=arb(hN)/3*arb(N)**(-sigma)+s.abs_upper()/2*arb(N)**(-eta)/eta
 rho_power=(arb(3)/2)**sigma/Dp(s).abs_lower()*(etail_power+rtail)
 power_interval=[str((scale*(residue.abs_lower()-rho_power)).lower()),str((scale*(residue.abs_upper()+rho_power)).upper())]
 energy_coeff=(2*(1-arb(2)**(-2*sigma-1))/(2*sigma+1)).sqrt()/(1-arb(2)**(-eta))
 etail_energy=arb(hN)/3*arb(N)**(-sigma)+s.abs_upper()/3*arb(N)**(-eta)*energy_coeff
 rho_energy=(arb(3)/2)**sigma/Dp(s).abs_lower()*(etail_energy+rtail)
 energy_interval=[str((scale*(residue.abs_lower()-rho_energy)).lower()),str((scale*(residue.abs_upper()+rho_energy)).upper())]
 need(scale*(residue.abs_lower()-rho_energy)>arb('0.30'),'energy-conditional residue margin')
 report=dict(energy_hypothesis='For every dyadic X >= 2^34, integral_X^(2X) H(floor(x))^2 dx <= 2 X^2',energy_hypothesis_amplitude_interval=energy_interval,status='CONDITIONAL: rigorously interval-checked implication, NOT a proof of the discrepancy hypothesis or a modal expansion',hypothesis='For every real x >= 2^34, |sum_{n<=x}(3*1_{n=2 mod3}-1)t(n)| <= sqrt(x)*log2(x)/4.',N=N,B=B,stronger_hypothesis='|H(x)| <= 1.5 sqrt(x) for every x >= 2^34',stronger_hypothesis_amplitude_interval=power_interval,independent_bit_moment_controls=controls,root_disk_center=str(c),root_disk_radius=str(radius),root_residual_upper=str(D(c).abs_upper()),simple_root_derivative_lower=str(Dp(s).abs_lower()),numerator=str(numerator),quadrature_disk_error=str(error),residue_before_future_tail=str(residue),future_residue_radius=str(rho),twice_residue_amplitude_at_log2n29_lower=str(low),twice_residue_amplitude_at_log2n29_upper=str(high),moments_sha256=hashlib.sha256(open(args.moments,'rb').read()).hexdigest())
 with open(args.output,'w')as f:json.dump(report,f,indent=2)
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
