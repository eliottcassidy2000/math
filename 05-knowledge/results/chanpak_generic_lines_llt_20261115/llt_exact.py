# Exact class counts per line for S = 3H + P2 (orders 8-10): relative unit-scale fluctuation around a smooth density.
import numpy as np
for n in (8,9,10):
    S=np.loadtxt('S_%d.txt'%n,dtype=np.int64)
    lat=6 if n%2 else 2; off=3 if n%2 else 1
    lo,hi=np.percentile(S,[20,80]).astype(np.int64)
    lo+=(off-lo)%lat
    N=np.arange(lo,hi+1,lat)
    c=np.bincount(S[(S>=lo)&(S<=hi)]-lo,minlength=hi-lo+1)[::lat][:len(N)].astype(float)
    for K in (25,100):
        ker=np.ones(2*K+1)/(2*K+1); sm=np.convolve(c,ker,mode='same'); inner=slice(K,len(N)-K)
        rel=(c[inner]-sm[inner])/sm[inner]
        print('order %2d (lattice %d): %6d bulk lines, mean classes/line %8.1f; smoothing +-%3d lines: relative fluctuation sd %.3f, '
              'Poisson-equivalent dispersion %.2f, min/smooth %.3f'%(n,lat,len(N),c.mean(),K,rel.std(),((c[inner]-sm[inner])**2).mean()/sm[inner].mean(),(c[inner]/sm[inner]).min()))
    r=(c[inner]-sm[inner])/np.sqrt(sm[inner])
    print('           residual autocorrelation, lags 1..6 lines:',' '.join('%.2f'%np.corrcoef(r[:-L],r[L:])[0,1] for L in range(1,7)))
