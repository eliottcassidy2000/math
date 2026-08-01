"""FC(2) on the anti-symmetric eigenspace, with SCALE-CORRECTED residuals.

The moments int_0^{1/2} h^{2k} dt shrink like (max|h|)^{2k} whatever h is, so an
absolute residual measures DECAY, not cancellation.  The honest test uses the
relative moment
        rho_k = | int_0^{1/2} h^{2k} dt |  /  int_0^{1/2} |h|^{2k} dt   in [0,1],
which is 0 exactly when the k-th moment genuinely cancels and 1 when there is
no cancellation at all.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import quad

def make(D):
    odd=[j for j in range(1,D+1,2)]
    def h_of(c):
        def h(t):
            s=0j
            for i,j in enumerate(odd): s+=c[i]*t**j
            return s
        return h
    return odd, h_of

def rel_moments(c, D, kmax):
    odd,h_of=make(D); h=h_of(c)
    out=[]
    for k in range(1,kmax+1):
        num_r=quad(lambda t: (h(t)**(2*k)).real,0,0.5,limit=200)[0]
        num_i=quad(lambda t: (h(t)**(2*k)).imag,0,0.5,limit=200)[0]
        den  =quad(lambda t: abs(h(t))**(2*k),0,0.5,limit=200)[0]
        out.append((abs(complex(num_r,num_i))/den if den>0 else 0.0, complex(num_r,num_i), den))
    return out

for D,K in ((5,2),(5,3),(7,3),(9,4)):
    odd,_=make(D); N=len(odd)
    def resid(p):
        c=p[:N]+1j*p[N:]
        rm=rel_moments(c,D,K)
        o=[r[0] for r in rm]
        o.append(float(np.sum(np.abs(c)**2))-1.0)
        return o
    rng=np.random.default_rng(11); best=1e9; bx=None
    for _ in range(12):
        p0=rng.normal(size=2*N); p0/=np.linalg.norm(p0)
        r=least_squares(resid,p0,xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=800)
        v=float(np.max(np.abs(resid(r.x))))
        if v<best: best,bx=v,r.x
        if best<1e-9: break
    print(f"deg<={D}, k=1..{K}: best RELATIVE residual {best:.3e}  -> "
          f"{'SOLVED' if best<1e-8 else 'NOT FOUND'}")
    if bx is not None:
        c=bx[:N]+1j*bx[N:]
        rm=rel_moments(c,D,max(K,5))
        print("    relative moments rho_k, k=1..%d: %s" % (max(K,5),
              ", ".join(f"{r[0]:.3e}" for r in rm)))
