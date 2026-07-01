"""OBSERVER-LENS via the reproducing kernel at the marked origin z=1. For mu_{S,t}=uniform on {e^{2pi i v t}},
build orthonormal OPUC by Gram-Schmidt in L^2(mu), form CD kernel K_N(1,1)=sum|phi_j(1)|^2. Emptiness at z=1
(lonely) => LOW local density => LARGE K_N(1,1). Does max_t K track the lonely time / M(S)?"""
import numpy as np
def cd_kernel_at_1(S,t):
    S=np.array(S); N=len(S); z=np.exp(2j*np.pi*S*t); w=np.ones(N)/N  # atoms & weights
    # Gram-Schmidt on monomials 1,z,...,z^{N-1} in L^2(mu):  <f,g>=sum_k w_k f(z_k) conj(g(z_k))
    V=np.array([z**j for j in range(N)]).T   # Vandermonde N x N, row=atom, col=degree
    def ip(a,b): return np.sum(w*a*np.conj(b))
    phis=[]; K1=0.0
    for j in range(N):
        u=V[:,j].copy().astype(complex)
        for p in phis: u=u-ip(u,p)*p
        nrm=np.sqrt(ip(u,u).real)
        if nrm<1e-12: continue
        u=u/nrm; phis.append(u)
        val1=np.sum([ (1.0**0) ]) # phi_j(1): evaluate the poly with these coeffs at z=1 -- track via coeffs instead
    # simpler: evaluate phi_j(1) by re-doing GS symbolically on coefficient vectors
    def ipc(ca,cb):  # inner product of two polys given by coeff vectors, via atom values
        av=sum(ca[m]*z**m for m in range(len(ca))); bv=sum(cb[m]*z**m for m in range(len(cb)))
        return np.sum(w*av*np.conj(bv))
    coeffs=[]; K1=0.0
    for j in range(N):
        c=np.zeros(N,dtype=complex); c[j]=1
        for pc in coeffs: c=c-ipc(c,pc)*pc
        nrm=np.sqrt(ipc(c,c).real)
        if nrm<1e-12: continue
        c=c/nrm; coeffs.append(c)
        phij_1=sum(c[m]*(1.0**m) for m in range(N))  # phi_j(1)
        K1+=abs(phij_1)**2
    return K1
def M_grid(S,Q=20000):
    S=np.array(S); best=0
    for j in range(1,Q): 
        t=j/Q; d=np.min(np.abs(S*t-np.round(S*t)))
        if d>best: best=d
    return best
# For the AP, sweep t and see if K_N(1,1) is maximized at the lonely time 1/14
S=list(range(1,14)); ts=np.linspace(0.001,0.5,500)
Ks=[cd_kernel_at_1(S,t) for t in ts]
tmax=ts[int(np.argmax(Ks))]
print(f"AP {{1..13}}: t maximizing K_N(1,1) = {tmax:.4f}  (lonely time 1/14={1/14:.4f})   Kmax={max(Ks):.1f}")
# compare configs at their lonely times: K_N(1,1) vs M
print("\nconfig                         K_N(1,1)@lonely    M(S)")
for name,S,t in [("AP {1..13}",list(range(1,14)),1/14),
                 ("construction",list(range(1,13))+[182],14/183),
                 ("GW {1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24],1/14),
                 ("spread(random)",[1,17,38,59,88,123,167,201,244,290,331,377,399],0.234)]:
    K=cd_kernel_at_1(S,t); M=M_grid(S)
    print(f"  {name:28s} {K:10.2f}       {M:.4f}")
print("\n  => if max_t K_N(1,1) lands on the lonely time, the reproducing kernel at z=1 IS the observer-lens objective;")
print("     maximizing emptiness-at-basepoint = the lonely runner, in Szego/OPUC language.")
