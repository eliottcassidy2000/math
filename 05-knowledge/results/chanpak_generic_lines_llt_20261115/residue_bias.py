import numpy as np
def ptab(n):
    p=np.zeros(n,dtype=np.int64)
    for u in range(1,n):
        w=u>>1; p[u]=p[w] if u&1 else p[w]+w+(w&1)
    return p
P=ptab(1<<14)
def rev(x,a):
    r=0
    for i in range(a): r=(r<<1)|((x>>i)&1)
    return r
a=12; R=np.array([rev(x,a) for x in range(1<<a)])
xs=np.arange(1<<(a-1),1<<a); ys=np.arange(0,1<<(a-1))
X,Y=np.meshgrid(xs,ys,indexing='ij'); X=X.ravel(); Y=Y.ravel()
m=((X^Y)&1)==1; X=X[m]; Y=Y[m]
d=R[X]-R[Y]; delta=X-Y; dp=P[X]-P[Y]
rhs=2*delta+dp
sel=(d>8)&(rhs>0)&(rhs//d>=8)&(rhs//d<256)&(d!=delta)
frac=(rhs[sel]%d[sel])/d[sel]
print('generic cores in strip:',sel.sum())
h,_=np.histogram(frac,bins=16,range=(0,1)); print('distribution of (rhs mod d)/d in 16 bins:',h)
# residue of rhs mod d restricted to small set? check rhs mod 2 vs d parity, and gcd structure
print('fraction with d | rhs: %.5f   mean 1/d: %.5f'%(np.mean(rhs[sel]%d[sel]==0),np.mean(1.0/d[sel])))
print('rhs parity: %s ; d parity: %s'%(np.bincount(rhs[sel]%2),np.bincount(d[sel]%2)))
# d | rhs  vs  d | rhs/gcd: check whether gcd(d, rhs) tends to be large
g=np.gcd(d[sel],rhs[sel]); print('mean gcd(d,rhs): %.3f ; fraction gcd>1: %.3f'%(g.mean(),(g>1).mean()))
# exact identity check: is rhs = d * something + small? test rhs mod d against (delta - d) or p-structure
r2=(rhs[sel]-(delta[sel]))%d[sel]; print('fraction d | rhs - delta: %.5f'%np.mean(r2==0))
r3=(dp[sel])%d[sel]; print('fraction d | p(x)-p(y): %.5f'%np.mean(r3==0))
