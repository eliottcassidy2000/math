#!/usr/bin/env python3
"""Describe four actual lifts, retaining both coefficients and old slack."""
import argparse,json
from math import gcd
from pathlib import Path
import numpy as np

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main(directory):
    profiles=[]
    for E in range(1,5):
        q,Q=5**E,5**(E+1)
        H=np.memmap(directory/f'lift-queue-{E}.bin',dtype='<i4',mode='r')
        J=np.memmap(directory/f'lift-queue-{E+1}.bin',dtype='<i4',mode='r')
        need(len(H)==q*(q-q//5) and len(J)==Q*(Q-Q//5),'complete fibres')
        index=np.arange(len(J),dtype=np.int64)
        row=index//Q;m=row+1+row//4;b=index%Q
        mm,bb=m%q,b%q;pi=(mm-mm//5-1)*q+bb
        D=8*5**(E-1);N=5**E-1
        old_scale=gcd(D,N);new_scale=gcd(5*D,5*N+4)
        K=new_scale*J.astype(np.int64)-5*old_scale*H[pi]
        values=np.unique(K)
        need(np.array_equal(values,np.arange(4-D,5*D,4)),'observed full interval of multiples of four')
        profile=dict(from_level=E,min_K=int(K.min()),max_K=int(K.max()),distinct_values=len(values),identity_K=int(K[0]))
        if E==4:
            examples=[]
            for a in (1,2):
                inv=pow(2**a,-1,Q);tm=3*m*inv%Q;tb=(3*b+1)*inv%Q
                target=(tm-tm//5-1)*Q+tb
                naive=K[target]-K-5*D*((tm==1)&(tb==0))+4*a
                bad=np.flatnonzero(naive<0);need(len(bad)>0,'discarded slack hostile')
                i=int(bad[0]);mi,bi=int(tm[i])%q,int(tb[i])%q
                j=(mi-mi//5-1)*q+bi
                R=sum((mi-1)%5**e==0 and bi%5**e==0 for e in range(1,E+1))
                sigma=old_scale*(int(H[j])-int(H[pi[i]]))-D*R+N*a
                need(int(naive[i])+5*sigma>=0,'retained slack repairs actual lift')
                examples.append(dict(a=a,failures_without_slack=len(bad),source=[int(m[i]),int(b[i])],
                                     target=[int(tm[i]),int(tb[i])],Ksource=int(K[i]),Ktarget=int(K[target[i]]),
                                     deficit_without_slack=int(-naive[i]),old_slack=sigma,repaired_slack=int(naive[i])+5*sigma))
            profile['omitted_slack_hostiles']=examples
        profiles.append(profile)
    return dict(profiles=profiles,scope='FINITE-EXACT lift profile through E=5; no induction theorem')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,default=Path('/tmp'));p.add_argument('--output',type=Path,default=Path(__file__).with_name('lift_profiles.json'));a=p.parse_args()
    report=json.dumps(main(a.directory),indent=2,sort_keys=True)+'\n';a.output.write_text(report);print(report,end='')
