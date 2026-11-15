"""Wave 25 (HYP-9692 reduction): band geometry and target windows by order.
Band B_n(1): |h - m| <= sd_h, |f - (alpha + beta h)| <= s. A reader with (a, c), c/a in [1/3, 1), crosses the band
in an h-interval of length about 2 s/(beta - c/a) >= 2 s/beta; admissible points on its line recur every <= 24 a
in h. So a guaranteed point needs a <= a*(n) = s/(12 beta). Order n then covers targets N in
[L_n, a*(n) U_n], with L_n = alpha + (beta-1)(m - sd_h) and U_n = alpha + beta (m + sd_h), using readers of every
scale 1 <= a <= a*(n) whose log-spacing is below log(U_n/L_n)."""
import numpy as np
def geom(H,F,W):
    m=np.average(H,weights=W); sd=np.sqrt(np.average((H-m)**2,weights=W))
    mf=np.average(F,weights=W); beta=np.average((H-m)*(F-mf),weights=W)/sd**2; alpha=mf-beta*m
    s=np.sqrt(np.average((F-alpha-beta*H)**2,weights=W)); return m,sd,alpha,beta,s
rows=[]
d=np.loadtxt('profhf_10.txt',dtype=np.int64); rows.append((10,)+geom(d[:,0].astype(float),d[:,1].astype(float),d[:,2].astype(float)))
for n,fn in ((11,'samp11.txt'),(12,'samp12.txt')):
    S=np.loadtxt(fn,dtype=np.int64); rows.append((n,)+geom(S[:,0].astype(float),S[:,1].astype(float),np.ones(len(S))))
print(" n   mean h     sd h      ridge f = alpha + beta h      resid s   a*(n)   L_n          a* U_n       U_n/L_n")
prev=None
for n,m,sd,alpha,beta,s in rows:
    astar=s/(12*beta); L=alpha+(beta-1)*(m-sd); U=alpha+beta*(m+sd)
    print(f"{n:2d} {m:10.0f} {sd:9.0f}   {alpha:12.0f} + {beta:6.2f} h   {s:9.0f}  {astar:6.1f}  {L:12.4g}  {astar*U:12.4g}   {U/L:5.2f}")
