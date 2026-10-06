"""Six-seven session (2026-10-06): companion of sixseven_20261006_lyapunov_adam.py. Run: python3 <this> n r starts seed (lowrank) or n starts seed (collect)."""
import sys, time, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sixseven_20261006_lyapunov_adam import make
n,r,starts,seed=map(int,sys.argv[1:5])
gg=make(n); rng=np.random.default_rng(seed); hits=0; best=-1; t0=time.time()
for s in range(starts):
    U=rng.standard_normal((n,r)); V=rng.standard_normal((n,r))
    A=U@V.T; c=np.sqrt(np.linalg.norm(A)); U/=c; V/=c
    mU=np.zeros_like(U); vU=np.zeros_like(U); mV=np.zeros_like(V); vV=np.zeros_like(V); b1,b2=0.9,0.999; bb=-1
    for t in range(1,1601):
        lr=0.02 if t<=800 else 0.006*(1-0.85*((t-800)/800)**2)
        A=U@V.T; g,G,ss,sk=gg(A); bb=max(bb,g)
        G=G-np.sum(G*A)*A
        GU=G@V; GV=G.T@U
        mU=b1*mU+(1-b1)*GU; vU=b2*vU+(1-b2)*GU*GU; mV=b1*mV+(1-b1)*GV; vV=b2*vV+(1-b2)*GV*GV
        U=U+lr*(mU/(1-b1**t))/(np.sqrt(vU/(1-b2**t))+1e-12); V=V+lr*(mV/(1-b1**t))/(np.sqrt(vV/(1-b2**t))+1e-12)
        A=U@V.T; c=np.sqrt(np.linalg.norm(A)); U/=c; V/=c
    if bb>1e-9: hits+=1
    best=max(best,bb)
print(f"n={n} rank<={r}: {starts} Adam runs, {hits} with gap>1e-9, best gap {best:.3e} ({time.time()-t0:.0f}s)",flush=True)
