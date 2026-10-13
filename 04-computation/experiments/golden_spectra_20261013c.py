#!/usr/bin/env python3
"""Numerical spectra on TWO declared golden rays; scalar model is a hypothesis."""
import argparse,json,numpy as np,mpmath as mp
from scipy.optimize import least_squares
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True);args=ap.parse_args();reports=[]
D=json.load(open(args.input));r=np.array(D['bins'],float);x=(np.log2(r[:,0])+np.log2(r[:,1]))/2
sel=(x>=14)&(x<22);x=x[sel];y=r[sel,3:5]/r[sel,2,None];lo=r[sel,0];hi=r[sel,1]
mp.mp.dps=30;theta=mp.log((1+mp.sqrt(5))/2,2)
for f in [1.5,11.5,24.5,37.5,47.5,60.5,83.5,96.5]:
 defect=f*theta-mp.nint(f*theta);w=0j
 for i in range(1,60):
  d=defect*i/59;w=mp.findroot(lambda w:.25*mp.exp(-w)+.75*mp.exp((theta-1)*w+2j*mp.pi*d)-1,w)
 pred=[float(f+w.imag/(2*mp.pi)),float(-w.real)]
 def residual(p,ray=0,ret=False):
  lam=-p[1]+2j*np.pi*p[0];s=1+lam/np.log(2);v=np.log(hi/lo)
  z=np.exp(lam*(np.log2(lo)-18))*np.expm1(s*v)/(s*np.expm1(v))
  B=np.c_[z.real,z.imag,np.ones(len(x))];c=np.linalg.lstsq(B,y[:,ray],rcond=None)[0];err=y[:,ray]-B@c
  return (float(np.hypot(*c[:2])),float(np.sqrt(np.mean(err**2)))) if ret else err
 fit=[]
 for ray in range(2):
  res=least_squares(lambda p:residual(p,ray),pred,bounds=([f-.2,-.2],[f+.2,.8]))
  fit.append([*res.x,*residual(res.x,ray,True)])
 reports.append(dict(base=f,scalar_model=pred,rays=[dict(frequency=float(v[0]),decay=float(v[1]),amplitude18=v[2],residual_rms=v[3])for v in fit]))
 print(f,'pred',pred,'fits',fit)

with open(args.output,'w')as fp:json.dump(dict(status='NUMERICAL; one-line fits do not authenticate scalar closure',interval=[14,22],rays=['n','n phi^2'],reports=reports),fp,indent=2)
