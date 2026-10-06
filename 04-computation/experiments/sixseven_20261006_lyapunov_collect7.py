"""Six-seven session (2026-10-06): companion of sixseven_20261006_lyapunov_adam.py. Run: python3 <this> n r starts seed (lowrank) or n starts seed (collect)."""
import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sixseven_20261006_lyapunov_adam import make, adam_run
n,starts,seed=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
gg=make(n); rng=np.random.default_rng(seed); found=[]
for s in range(starts):
    b,A=adam_run(n,gg,rng)
    if b>1e-9:
        sv=np.linalg.svd(A,compute_uv=False); found.append((b,sv/sv[0]))
        print(f"gap {b:.3e}  normalised singular values {np.round(sv/sv[0],3)}",flush=True)
print(f"n={n} seed={seed}: {len(found)} counterexamples in {starts} runs")
if found:
    S=np.array([f[1] for f in found]); print("median profile:",np.round(np.median(S,axis=0),3)," max of the two smallest:",np.round(S[:,-2:].max(axis=0),3))
