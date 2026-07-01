"""TASK 1: LRC atom pair-correlation. Atoms = units (Z/N)* (the AP lonely times). 1st moment (Fourier) =
Ramanujan sum c_N(k) (klein HYP-3793 = Lefschetz trace S23). 2nd moment = pair-correlation A(m)=#{unit pairs
diff m}; Wiener-Khinchin: A_hat(k) = |c_N(k)|^2 (the atoms' 2nd-moment power spectrum = SQUARED Ramanujan)."""
import numpy as np
from math import gcd
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def ramanujan(N,k): return sum(np.exp(2j*np.pi*k*a/N) for a in units(N))
def paircorr(N):
    U=units(N); A=np.zeros(N)
    for a in U:
        for b in U: A[(a-b)%N]+=1
    return A
for N in [14,10,6]:
    U=units(N); phi=len(U); A=paircorr(N)
    cN=np.array([ramanujan(N,k) for k in range(N)])
    Ahat=np.fft.fft(A)                     # discrete Fourier of the pair-correlation
    wk = np.allclose(Ahat, np.abs(cN)**2, atol=1e-6)
    print(f"N={N}: {phi} atoms (units {U})")
    print(f"  1st moment c_N(k) (Ramanujan, k=0..{N-1}): {[int(round(x.real)) for x in cN]}")
    print(f"  2nd moment / pair-corr A(m)=#unit-pairs diff m: {[int(x) for x in A]}")
    print(f"  Wiener-Khinchin  A_hat(k) == |c_N(k)|^2 (2nd moment = SQUARED Ramanujan)? {wk}")
    print(f"  A(0)=#atoms={int(A[0])}=phi; sum A = phi^2 = {int(A.sum())}; the m!=0 values = resonance-spacing multiset")
# The variance / 2nd central moment of the atom Fourier (spectral 2nd moment)
N=14; cN=np.array([ramanujan(N,k) for k in range(N)]); 
print(f"\nAtoms mod 14: sum_k |c_N(k)|^2 = {int(round((np.abs(cN)**2).sum()))} = N*phi = {14*6} (Parseval); mean |c|^2 = phi = 6.")
print("=> the atoms' pair-correlation is EXACTLY the squared-Ramanujan power spectrum: 1st moment=trace c_N(k), 2nd=|c_N(k)|^2.")
