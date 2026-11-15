# Corrected random-divisibility test for generic two-ended cores of G(u) = 2u + p(u) (c=2, k=1, phi=0).
# Collision (TWISTED_LATTICES eq. (1)):  d w + delta r + [w odd] - eps [r even] = ((c/k) 2^j + 2^j - 1) delta + p(x)-p(y)
# For each core and each (r, parity of w) the exact w is checked; prediction = sum over (r, e_w) of 1[w_real in range]/|d|.
import numpy as np
def ptab(n):
    p=np.zeros(n,dtype=np.int64)
    for u in range(1,n):
        w=u>>1; p[u]=p[w] if u&1 else p[w]+w+(w&1)
    return p
P=ptab(1<<15)
def rev(x,a):
    r=0
    for i in range(a): r=(r<<1)|((x>>i)&1)
    return r
c,k=2,1
for a,j in ((12,0),(13,0),(11,2),(12,2),(10,4)):
    R=np.array([rev(x,a) for x in range(1<<a)])
    xs=np.arange(1<<(a-1),1<<a); ys=np.arange(0,1<<(a-1))
    X,Y=np.meshgrid(xs,ys,indexing='ij'); X=X.ravel(); Y=Y.ravel()
    m=((X^Y)&1)==1; X=X[m]; Y=Y[m]
    d=R[X]-R[Y]; delta=X-Y; eps=(X&1)-(Y&1)
    gen=(d!=delta)&(np.abs(d)>8)
    X,Y,d,delta,eps=X[gen],Y[gen],d[gen],delta[gen],eps[gen]
    base=((c//k)*2**j+2**j-1)*delta+P[X]-P[Y]
    for Wb in (3,5,7):
        W=1<<Wb; act=0; pred=0.0
        for r in range(1<<j):
            for e in (0,1):
                num=base-delta*r-e+(eps*(1 if (r%2==0 and j>0) else 0))
                ok=(num%np.abs(d)==0)&(d!=0)
                w=np.where(ok,num//np.where(d==0,1,d),-1)
                act+=int((ok&(w>=W)&(w<2*W)&((w&1)==e)).sum())
                wr=num/d.astype(float)
                pred+=float((((wr>=W)&(wr<2*W))/np.abs(d)).sum())
        print('a=%2d j=%d W=2^%d: generic two-ended collisions actual %6d, random-divisibility prediction %9.1f, ratio %.3f'%(a,j,Wb,act,pred,act/max(pred,1e-9)))
