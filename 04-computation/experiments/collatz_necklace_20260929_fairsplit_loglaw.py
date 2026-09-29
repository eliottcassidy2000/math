"""j=3 fair cuts as a 2D cyclic lattice walk hitting the origin.
(a) numpy Monte Carlo of P(N_3>=1) and E[N_3 | N_3>=1] for m up to 1e5 (rho=1/2).
(b) exact second moment: E[N_3^2] = (m/C(3m,3x)) sum_{delta=0}^{m-1} sum_s C(delta,s)^3 C(m-delta,x-s)^3
    (two fair 3-cuts at 0 and delta force the six window counts s,x-s,s,x-s,s,x-s), and the
    second-moment lower bound P(N_3>=1) >= E[N]^2/E[N^2]."""
import numpy as np, sys
from math import lgamma, exp, log, comb, sqrt, pi
rng=np.random.default_rng(20260929)
def mc(m,x,T):
    K=3*m; X=3*x; hits=0; tot=0
    for _ in range(T):
        w=np.zeros(K,dtype=np.int8); w[rng.choice(K,X,replace=False)]=1
        F=np.concatenate(([0],np.cumsum(np.concatenate((w,w)))))
        r=np.arange(m)
        ok=(F[r+m]-F[r]==x)&(F[r+2*m]-F[r+m]==x)&(F[r+3*m]-F[r+2*m]==x)
        n=int(ok.sum()); hits+=(n>0); tot+=n
    return hits/T, tot/T, (tot/hits if hits else float('nan'))
print("(a) Monte Carlo, rho=1/2:")
for m,T in [(1000,3000),(3000,2000),(10000,1000),(30000,500),(100000,200)]:
    p,e,c=mc(m,m//2,T)
    print(f"  m={m}: P(N>=1)={p:.4f}  E[N]={e:.3f}  E[N|N>=1]={c:.2f}   P*log(m)={p*log(m):.3f}  P*sqrt(log m)={p*sqrt(log(m)):.3f}")
def lc(n,k):
    if k<0 or k>n: return None
    return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
print("\n(b) exact E[N_3], E[N_3^2], second-moment bound, rho=1/2:")
for m in [10,30,100,300,1000,3000,10000]:
    x=m//2; LC=lc(3*m,3*x)
    EN=m*exp(3*lc(m,x)-LC)
    # second moment
    tot=0.0
    for d in range(m):
        # sum_s C(d,s)^3 C(m-d,x-s)^3
        s_lo=max(0,x-(m-d)); s_hi=min(d,x)
        if s_lo>s_hi: continue
        s=np.arange(s_lo,s_hi+1)
        la=np.array([lc(d,int(t)) for t in s]); lb=np.array([lc(m-d,int(x-t)) for t in s])
        tot+=np.exp(3*la+3*lb-LC).sum()
    EN2=m*tot
    print(f"  m={m}: E[N]={EN:.4f}  E[N^2]={EN2:.4f}  E[N^2]/log m={EN2/log(m):.4f}  lower bound P>=E^2/E[N^2]={EN*EN/EN2:.4f}")
