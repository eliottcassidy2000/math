import numpy as np
from itertools import product
print("LRC AS A MOMENT NULLCONE. The danger comb 1_{D_i}(t) = sum_k hhat(k) e^{2pi i k v_i t},")
print("hhat(k) = sin(2 pi k lam)/(pi k), hhat(0)=2 lam. The UNCOVERED measure of the union of")
print("combs D_i is the CONSTANT Fourier term of prod_i (1 - 1_{D_i}):")
print("   uncovered = sum over (k_1..k_r) with sum k_i v_i = 0  of  prod_i (-c_{k_i}),")
print("where c_0 = -(1-hhat(0)) for the '1-' and c_k = hhat(k) for k!=0. This is a CHARGE-0")
print("(RESONANCE-LATTICE) sum of products -- structurally IDENTICAL to E[P^m] = sum over")
print("charge-0 reps of prod(coeffs). Frequency = charge; integral over circle kills nonzero")
print("frequency; resonance lattice {sum k_i v_i = 0} = the charge-0 lattice.")
print()
def hhat(k,lam):
    if k==0: return 2*lam
    return np.sin(2*np.pi*k*lam)/(np.pi*k)
def uncovered_direct(V,lam,Ngrid=200000):
    t=np.arange(Ngrid)/Ngrid
    safe=np.ones(Ngrid,bool)
    for v in V:
        d=np.abs(((v*t+0.5)%1)-0.5)  # ||v t||
        safe &= (d>=lam)
    return safe.mean()
def uncovered_resonance(V,lam,K=6):
    """sum over |k_i|<=K with sum k_i v_i = 0 of prod( (1-hhat(0)) if k=0 else -hhat(k) )...
       careful: prod(1 - 1_{D_i}) fourier: 1_{D_i} has coeff hhat(k) at freq k v_i (k incl 0).
       (1 - 1_{D_i}) has coeff [k=0](1) - hhat(k) => c^{(i)}_0 = 1 - hhat(0), c^{(i)}_k = -hhat(k)."""
    r=len(V)
    tot=0.0
    for ks in product(range(-K,K+1),repeat=r):
        if sum(k*v for k,v in zip(ks,V))!=0: continue
        term=1.0
        for k in ks:
            term *= (1-hhat(0,lam)) if k==0 else (-hhat(k,lam))
        tot+=term
    return tot
print("  VERIFY on small speed sets (lam=1/n):")
print("   V              lam     uncovered(direct)   uncovered(resonance sum, K=6)")
for V,lam in [([1,2],1/3),([1,2,3],1/4),([1,3],1/4),([1,2,3,4],1/5)]:
    d=uncovered_direct(V,lam); r=uncovered_resonance(V,lam,K=8)
    print(f"   {str(V):14s} {lam:.3f}   {d:.5f}            {r:.5f}   (match: {abs(d-r)<0.01})")
print()
print("  => the uncovered measure IS the resonance-lattice (charge-0) moment sum of prod hhat.")
print("  hhat(k)=0 <=> sin(2 pi k lam)=0 <=> 2 n | 2k for lam=1/n... i.e. the parity/resonance")
print("  annihilation (THM-1200). SIGNS of hhat = sin alternate => this is the BOTH-SIGNS")
print("  (charge-indefinite) hard case -- exactly why LRC(14) resists the positivity lever.")
