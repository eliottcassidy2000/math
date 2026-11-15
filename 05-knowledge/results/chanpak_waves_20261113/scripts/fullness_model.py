"""Wave 25 (P73): Poisson occupancy model for the ridge band, calibrated on the exact order-10 census.
lambda_n = (classes with profile in band(0.5,0.5))/(admissible band points); occupancy ~ 1 - exp(-lambda/k)."""
import numpy as np, collections, math
CLASSES={10:9733056,11:903753248,12:154108311168}
def ridge(H,F,W):
    mh=np.average(H,weights=W); sh=np.sqrt(np.average((H-mh)**2,weights=W))
    mf=np.average(F,weights=W); beta=np.average((H-mh)*(F-mf),weights=W)/sh**2; alpha=mf-beta*mh
    s=np.sqrt(np.average((F-alpha-beta*H)**2,weights=W)); return mh,sh,alpha,beta,s
def band_points(n,mh,sh,alpha,beta,s,r=0.5):
    step=8 if n%2==0 else 4; eps=1 if n%2==0 else -1; tot=0
    for h in range(int(mh-r*sh)|1,int(mh+r*sh)+1,2):
        c=alpha+beta*h; lo=(int(c-r*s)//step+1)*step; hi=int(c+r*s)
        # count f in [lo,hi] step 'step' with f = eps*h mod 3
        for f0 in range(lo,lo+3*step,step):
            if (f0-eps*h)%3==0:
                tot+= max(0,(hi-f0)//(3*step)+1) if f0<=hi else 0; break
    return tot
def inband(H,F,mh,sh,alpha,beta,s,r=0.5): return (np.abs(H-mh)<=r*sh)&(np.abs(F-alpha-beta*H)<=r*s)
# exact order 10
d=np.loadtxt('profhf_10.txt',dtype=np.int64); H,F,C=d[:,0].astype(float),d[:,1].astype(float),d[:,2].astype(float)
p=ridge(H,F,C); m=inband(H,F,*p)
pts=band_points(10,*p); cls=C[m].sum(); occ=m.sum()/pts
lam=cls/pts; k=lam/(-math.log(1-occ))
print(f"order 10 exact: band points {pts}, classes in band {cls:.0f}, lambda {lam:.2f}, occupancy {occ:.3f} -> collision factor k = {k:.2f}")
for n,fn in ((10,'samp10.txt'),(11,'samp11.txt'),(12,'samp12.txt')):
    S=np.loadtxt(fn,dtype=np.int64); Hs,Fs=S[:,0].astype(float),S[:,1].astype(float)
    q=ridge(Hs,Fs,np.ones(len(Hs))); P=inband(Hs,Fs,*q).mean(); pts=band_points(n,*q)
    lam=CLASSES[n]*P/pts; pred=1-math.exp(-lam/k)
    print(f"order {n} (samples {len(Hs)}): ridge f = {q[2]:.0f} + {q[3]:.2f} h, residual sd {q[4]:.0f}; P(band) {P:.4f}; band points {pts}; lambda {lam:.1f}; predicted occupancy {pred:.4f}")
