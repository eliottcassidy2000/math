"""M(S) moment relaxation with MOSEK (SCS was non-monotone/inaccurate). Univariate: M~=max_{c in[-1,1]} min_v
(1-T_v(c))/2. Bisect s with MOSEK feasibility (crisp). Validate {1,2,3,4}->0.345, AP->sin^2(pi/14)=0.0495."""
import numpy as np, cvxpy as cp, warnings, time
from numpy.polynomial import chebyshev as C
warnings.filterwarnings("ignore")
def gcoef(v):
    g=-0.5*C.cheb2poly(np.eye(v+1)[v]); g[0]+=0.5; return g
def build(S,R):
    Vs=sorted(set(S)); y=cp.Variable(2*R+1); sp=cp.Parameter()
    cons=[y[0]==1, cp.bmat([[y[i+j] for j in range(R+1)] for i in range(R+1)])>>0,
          cp.bmat([[y[i+j]-y[i+j+2] for j in range(R)] for i in range(R)])>>0]
    for v in Vs:
        gc=gcoef(v); sz=max(1,R-(v+1)//2+1)
        cons.append(cp.bmat([[cp.sum([gc[k]*y[i+j+k] for k in range(len(gc))])-sp*y[i+j] for j in range(sz)] for i in range(sz)])>>0)
    return cp.Problem(cp.Minimize(0),cons), sp, y
def Mtilde(S,R,steps=30):
    prob,sp,y=build(S,R); lo,hi=0.0,0.9
    for _ in range(steps):
        mid=(lo+hi)/2; sp.value=mid
        try: prob.solve(solver=cp.MOSEK); feas=prob.status.startswith("optimal")
        except Exception: feas=False
        if feas: lo=mid
        else: hi=mid
    return lo, y.value
for S,tgt in [([1,2,3,4],1/5),([1,2,4],1/3),(list(range(1,14)),1/14)]:
    R=max(S); t0=time.time(); Mt,yv=Mtilde(S,R)
    true=np.sin(np.pi*tgt)**2
    # flat extension: rank of moment matrix
    Mmat=np.array([[yv[i+j] for j in range(R+1)] for i in range(R+1)]); ev=np.sort(np.abs(np.linalg.eigvalsh(Mmat)))[::-1]
    rank=int(np.sum(ev>1e-6*ev[0]))
    print(f"S(max={R}): M~={Mt:.6f} => M={np.arcsin(np.sqrt(max(Mt,0)))/np.pi:.6f} (true {tgt:.6f}=sin^2 {true:.6f}, gap {Mt-true:+.5f}); moment-matrix rank(atoms)={rank} [{time.time()-t0:.0f}s]")
