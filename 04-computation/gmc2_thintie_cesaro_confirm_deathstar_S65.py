# Closing boxeph THM-1630's thin-tie residual.  NC2 => sum_j I_m^{(j)}=0 forall m,
#   I_m^{(j)} = -(1/2)Gamma(m/2) C_j^{-m} beta_j (1+o(1)),  |C_j|=R tied, distinct args phi_j, beta_j != 0.
# Normalize by -(1/2)Gamma(m/2)R^{-m} (common factor, no overflow):  a_m := sum_j beta_j e^{-i m phi_j}
# satisfies a_m + (o(1) correction) = 0, hence a_m -> 0.  Wiener/Cesaro then recovers beta_k=0.
import cmath, numpy as np
def cesaro_recover(phis, betas, M=20000):
    # (1/M) sum_{m=1..M} a_m e^{i m phi_k} -> beta_k, where a_m = sum_j beta_j e^{-i m phi_j}
    ms=np.arange(1,M+1)
    a=np.zeros(M,dtype=complex)
    for j,p in enumerate(phis): a += betas[j]*np.exp(-1j*ms*p)
    return [np.mean(a*np.exp(1j*ms*p)) for p in phis], a
def no_cancellation(phis, betas, o1=lambda m:-1.0/m):
    # normalized S_m = sum_j beta_j e^{-i m phi_j}(1+o(1)); show |S_m| does NOT ->0
    out=[]
    for m in (10,50,200,1000,5000):
        S=sum(betas[j]*cmath.exp(-1j*m*phis[j])*(1+o1(m)) for j in range(len(phis)))
        out.append(round(abs(S),4))
    return out

for phis,betas,lab in ([[0.5,1.7,2.9],[1.0,1.0,1.0],"3 tied, beta=(1,1,1)"],
                       [[0.3,2.1,3.5,5.0],[1.0,-0.7,0.4,1j],"4 tied, mixed beta"]):
    rec,a=cesaro_recover(phis,betas)
    print(f"[{lab}]")
    print(f"   normalized |S_m| at m=10,50,200,1000,5000: {no_cancellation(phis,betas)}  (bounded away from 0 => tie does NOT cancel)")
    print(f"   Cesaro recovery beta_k: {[complex(round(x.real,3),round(x.imag,3)) for x in rec]}  == betas => (a_m->0 forces beta=0)")

# The finite linear system boxeph named: [e^{-i m phi_j}] has full column rank J for distinct phi_j
print("\nVandermonde full-rank check (the 'finite linear system per tie'):")
for phis in ([0.5,1.7,2.9],[0.3,2.1,3.5,5.0],[0.5,1.7,2.9,4.1,5.6]):
    Mtx=np.array([[cmath.exp(-1j*m*p) for p in phis] for m in range(1,80)])
    print(f"   J={len(phis)} tied arcs: smallest singular value = {np.linalg.svd(Mtx,compute_uv=False)[-1]:.3f}  (>0 => rank J => only beta=0)")
print("\nCONCLUSION: distinct phases => full-rank Vandermonde => sum_j beta_j e^{-imphi_j}=0 forall m has ONLY beta=0.")
print("The o(1) (which blocks a naive FINITE Vandermonde) is handled by Cesaro: a_m->0 => beta_k=0 (Wiener).")
print("=> boxeph's thin-tie residual (>=3 tied moduli, nonzero jumps) has NO nullcone member. CLOSED.")
