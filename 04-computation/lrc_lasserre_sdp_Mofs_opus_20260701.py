"""Tight SOS dual (single solve) for M~(S)=max_{c in[-1,1]} min_v g_v(c):
 min gamma s.t. gamma - sum_v lam_v(c) g_v(c) = sig0(c) + (1-c^2) sig1(c),  lam_v,sig0,sig1 SOS, sum_v lam_v(c)=1.
lam_v(c) c-dependent (degree 2e) => TIGHT. Returns M~ and the magic-function multipliers."""
import numpy as np, cvxpy as cp, warnings, time
from numpy.polynomial import chebyshev as C
warnings.filterwarnings("ignore")
def gcoef(v):
    g=-0.5*C.cheb2poly(np.eye(v+1)[v]); g[0]+=0.5; return g   # g_v monomial coeffs len v+1
def sos_dual(S,e):
    Vs=sorted(set(S)); dmax=max(Vs)
    Dprod=2*e+dmax                     # deg of lam_v*g_v
    m0=(Dprod+1)//2                    # sig0 half-deg
    D=2*m0
    gam=cp.Variable()
    Qv={v:cp.Variable((e+1,e+1),PSD=True) for v in Vs}
    Q0=cp.Variable((m0+1,m0+1),PSD=True); Q1=cp.Variable((m0,m0),PSD=True)
    def hankel(Q,k):  # coeff of c^k in z^T Q z, z=[1,c,..,c^{n-1}]
        n=Q.shape[0]; return cp.sum([Q[i,k-i] for i in range(max(0,k-(n-1)),min(k,n-1)+1)]) if 0<=k<=2*(n-1) else 0
    cons=[]
    # sum_v lam_v(c) = 1
    for k in range(2*e+1):
        cons.append(cp.sum([hankel(Qv[v],k) for v in Vs]) == (1.0 if k==0 else 0.0))
    # identity gamma[k==0] - sum_v (lam_v g_v)_k = (sig0)_k + (sig1(1-c^2))_k
    for k in range(D+1):
        lhs = (gam if k==0 else 0) - cp.sum([ cp.sum([gcoef(v)[b]*hankel(Qv[v],k-b) for b in range(len(gcoef(v))) if 0<=k-b]) for v in Vs])
        rhs = hankel(Q0,k) + (hankel(Q1,k) - hankel(Q1,k-2))
        cons.append(lhs==rhs)
    prob=cp.Problem(cp.Minimize(gam),cons); prob.solve(solver=cp.CLARABEL)
    return prob.value, prob.status
for S,tgt in [([1,2,3],1/4),([1,2,3,4],1/5),([1,2,4],1/3)]:
    true=np.sin(np.pi*tgt)**2
    for e in [1,2,3]:
        t0=time.time()
        try:
            val,st=sos_dual(S,e)
            print(f"S={S} e={e}: M~={val:.6f} ({st})  true sin^2={true:.6f} gap {val-true:+.6f}  M={np.arcsin(np.sqrt(max(val,0)))/np.pi:.5f} [{time.time()-t0:.1f}s]")
        except Exception as ex: print(f"S={S} e={e}: FAIL {ex}")
