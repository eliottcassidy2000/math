import numpy as np
print("=== LENS 2 (autocorrelation energy): merit factor / Littlewood / LRC disc_v are ONE object ===")
# aperiodic autocorrelation energy of a +-1 sequence (merit-factor denominator) = a Parseval energy
# of |S(theta)|^2 -- the SAME functional shape as LRC disc_v = sum_{m!=0}|chat_{mv}|^2.
def autocorr_energy(seq):
    N=len(seq); Cs=[sum(seq[n]*seq[n+k] for n in range(N-k)) for k in range(1,N)]
    E_spatial=sum(c*c for c in Cs)                       # sum_{k>=1} C_k^2
    # Parseval: (1/2pi)int |S|^4 = N^2 + 2 sum_{k>=1} C_k^2, S(t)=sum seq[n]e(nt)
    th=np.linspace(0,1,20000,endpoint=False)
    S=sum(seq[n]*np.exp(2j*np.pi*n*th) for n in range(N))
    E_parseval=(np.mean(np.abs(S)**4)-N**2)/2
    return E_spatial,E_parseval
seq=[1,1,1,-1,-1,1,-1]  # a short Barker-like sequence
Es,Ep=autocorr_energy(seq)
print(f"  +-1 seq {seq}: autocorr energy sum C_k^2 = {Es}, Parseval (int|S|^4 - N^2)/2 = {Ep:.3f}  (equal => same object)")
print("  So merit factor F = N^2/(2*energy), Littlewood L^1/L^4, and LRC disc_v are all L^p-energies of a")
print("  structured exponential sum -> one lens.")

print("\n=== LENS 4 (no-common-root / real-rootedness): GMC-Hermite = Lee-Yang/MSS interlacing ===")
def He_roots(m):
    # probabilists' Hermite via companion matrix of the 3-term recurrence He_{k+1}=x He_k - k He_{k-1}
    # build coefficients
    H=[np.array([1.0]),np.array([0.0,1.0])]
    for k in range(1,m): H.append(np.polynomial.polynomial.polysub(
        np.r_[0,H[k]], k*np.r_[H[k-1],0,0][:len(H[k])+1] if len(H[k-1])<=len(H[k])+1 else k*H[k-1]))
    # simpler: use numpy Hermite_e
    from numpy.polynomial import hermite_e as He
    c=[0]*m+[1]; return np.sort(He.hermeroots(c))
for m in (3,4):
    r1,r2=He_roots(m),He_roots(m+1)
    # interlacing: between consecutive roots of He_{m+1} lies exactly one root of He_m
    inter=all(r1[i]<r2[i+1]<r1[i+1] for i in range(len(r1)-1)) if len(r1)>=2 else True
    common=any(min(abs(r2-x))<1e-9 for x in r1)
    print(f"  He_{m} roots {np.round(r1,3)}, He_{m+1} roots {np.round(r2,3)}: interlace={inter}, share a root={common}")
print("  Consecutive Hermite are real-rooted and INTERLACE (no common root) -- the exact rigidity behind")
print("  Lee-Yang / de Bruijn-Newman (RH via real-rooted approximants) and MSS/Kadison-Singer (interlacing).")
