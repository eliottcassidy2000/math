# Unit-scale local limit test for S = alpha H + beta P2 over uniform random tournaments (tiltsample q = 0.5).
import numpy as np, glob, sys
n=int(sys.argv[1]); al,be=(int(sys.argv[2]),int(sys.argv[3])) if len(sys.argv)>3 else (3,1)
def load(f):
    raw=open(f).read().split('\n'); raw=[l for l in raw if len(l.split())==3]
    return np.array([l.split()[:2] for l in raw],dtype=np.int64).reshape(-1,2)
A=np.concatenate([load(f) for f in sorted(glob.glob('u%d_*.txt'%n))])
S=al*A[:,0]+be*A[:,1]; H=A[:,0]
print('order %d, %d uniform samples; S = %dH + %dP2: mean %.0f, sd %.0f; log S mean %.3f sd %.3f, skew %.3f'
      %(n,len(S),al,be,S.mean(),S.std(),np.log(S).mean(),np.log(S).std(),((np.log(S)-np.log(S).mean())**3).mean()/np.log(S).std()**3))
lo,hi=np.percentile(S,[20,80]).astype(np.int64)
lat=6 if (n%2==1 and al%3==0) else 2
Sb=S[(S>=lo)&(S<=hi)]; c=np.bincount(Sb-lo,minlength=hi-lo+1)[:hi-lo+1]
N=np.arange(lo,hi+1); keep=(N%2==1) if lat==2 else (N%6==3)
N=N[keep]; c=c[keep]
print('bulk [%d, %d]: %d lattice points (step %d), mean count %.2f per point'%(lo,hi,len(N),lat,c.mean()))
# singular series: mean count by residue class, normalized
for m in (8,16,32,3,9,27,5,7,24):
    cls=N%m; vals=sorted(set(cls.tolist()))
    means=np.array([c[cls==v].mean() for v in vals]); means/=means.mean()
    print('  mod %2d: relative density by class min %.3f max %.3f   %s'%(m,means.min(),means.max(),' '.join('%d:%.2f'%(v,x) for v,x in zip(vals,means)) if m<=16 else ''))
# smoothness: divide out a smooth density (moving average over +-K points) times the mod-48 factor
m=48 if lat==2 else 144
cls=N%m; fac=np.zeros(m); 
for v in set(cls.tolist()): fac[v]=c[cls==v].mean()
fac/=fac[fac>0].mean()
base=c/np.where(fac[cls]>0,fac[cls],1)
K=60; ker=np.ones(2*K+1)/(2*K+1); smooth=np.convolve(base,ker,mode='same')
expect=smooth*fac[cls]
inner=slice(K,len(N)-K)
res=(c[inner]-expect[inner]); disp=(res**2).mean()/expect[inner].mean()
print('  after smooth density x mod-%d factor: dispersion index Var(c - expected)/mean = %.3f (Poisson = 1)'%(m,disp))
# lag correlations of residuals (fine structure)
r=res/np.sqrt(expect[inner])
for lag in (1,2,3,4,6,8,12):
    print('    residual autocorrelation at lag %2d points: %.3f'%(lag,np.corrcoef(r[:-lag],r[lag:])[0,1]))
