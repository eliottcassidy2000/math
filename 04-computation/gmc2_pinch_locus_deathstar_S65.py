# For the GMC radial D(r,t) = (1 - beta(r)t)^2 - 4 alpha(r) t^2, the fold/pinch events (D = d/dr D = 0)
# locate boxeph's arcs. Claim: pinches sit at r with alpha'(r)^2 = alpha(r) beta'(r)^2.
# For alpha = gamma r (const a,c), beta = b(r): alpha'=gamma, so pinch locus is  r*b'(r)^2 = gamma.
# This makes "stacked jumps" (two pinches sharing the same t-arc) an explicit coincidence -> high codim.
import numpy as np
def pinch_check(gamma, bcoef, label):
    # b(r)=sum bcoef[i] r^i ; b'(r); pinch locus g(r)=r*b'(r)^2 - gamma = 0
    b=np.polynomial.Polynomial(bcoef); bp=b.deriv()
    # g(r) = r*bp(r)^2 - gamma
    r=np.linspace(0.01,15,200000)
    g=r*bp(r)**2 - gamma
    # sign changes = pinch locations
    sc=np.where(np.diff(np.sign(g)))[0]
    roots=[r[i] for i in sc]
    # verify at a root: D and dDdr vanish simultaneously at the corresponding t
    verified=[]
    for rp in roots:
        bpv=bp(rp); bv=b(rp)
        if abs(bpv)<1e-9: continue
        tp=1.0/(bv-2*gamma/bpv) if abs(bv-2*gamma/bpv)>1e-9 else None
        if tp is None: continue
        D=(1-bv*tp)**2-4*gamma*rp*tp**2
        dD=2*(1-bv*tp)*(-bpv*tp)-4*gamma*tp**2   # d/dr D, with alpha'=gamma
        verified.append((round(rp,3),round(tp,4),abs(D)<1e-6,abs(dD)<1e-6))
    print(f"[{label}] gamma={gamma}, b={bcoef}: pinch roots of r*b'^2=gamma at r={[round(x,3) for x in roots]}")
    print(f"   (r_p, t_p, D==0, dD/dr==0): {verified}  -> folds sit exactly on the locus")
pinch_check(1.0,[0,0,1.0],"b=r^2")          # b'=2r, r*(2r)^2=4r^3=1 -> one pinch
pinch_check(1.0,[0,1.0,1.0],"b=r+r^2")      # b'=1+2r, r(1+2r)^2=1
pinch_check(-1.0,[0,1.0,-1.0],"b=r-r^2 (gamma<0, complex branch regime)")
print("\nSTACKED JUMPS = two distinct pinch roots r_p1 != r_p2 sharing the SAME t_p AND same rate.")
print("t_p = 1/(b(r_p) - 2*gamma/b'(r_p)). Equal-t_p for two roots is one equation on the coefficients")
print("(codim 1); equal RATE (e^{-r_p} arc modulus) is a second (codim 2) -> stacked jumps are a thin")
print("sub-stratum of the already-thin tie locus. This is a concrete handle for boxeph's S182 edge (3).")
