import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, sys
T=int(float(sys.argv[1]))
P=np.zeros(T//3+10,dtype=np.int64)
for u in range(1,len(P)):
    P[u]= P[u//2] if u&1 else P[u//2]+u//2+((u//2)&1)
def rep(q,t,s,a,n,f,e):
    r=np.zeros(T+1,dtype=np.int32)
    X=1
    while True:
        base=t*X+a*P[X]+n+f*(X&1)
        if q*X+s+base>T: break
        ymax=(T-base)//(q*X+s)+1
        Y=np.arange(1,ymax+1,dtype=np.int64)
        F=(q*X+s)*Y+base+a*P[Y]+e*(Y&1)
        F=F[F<=T]
        np.add.at(r,F,1)
        X+=1
    return r
forms={'3-cycle (a=3)':(12,6,3,3,5,2,1),'pure hyperbola (a=0)':(12,6,3,0,5,2,1)}
for name,fm in forms.items():
    r=rep(*fm)
    print(name)
    for lo in [10**5,10**6,int(T/3)]:
        hi=min(2*lo,T)
        seg=r[lo|1:hi:2].astype(float)
        # local mean by residue mod 9 (odd numbers) to remove the mod-3 bias
        idx=np.arange(lo|1,hi,2)
        mean=np.zeros_like(seg); 
        for c in range(9):
            msk=(idx%9==c)
            if msk.any(): mean[msk]=seg[msk].mean()
        var=((seg-mean)**2).mean(); m=seg.mean()
        print(f'  N in [{lo},{hi}): mean r = {m:.3f}, var = {var:.3f}, var/mean = {var/m:.3f}, misses = {np.mean(seg==0):.4f}, (var/mean)/mean = {var/m/m:.3f}')
