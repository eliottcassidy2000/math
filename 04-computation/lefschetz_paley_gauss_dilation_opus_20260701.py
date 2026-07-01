"""Lefschetz trick for the hard (free-Z2, p=3 mod4, Paley) side. Where Brouwer/symmetric-SOS fails, Lefschetz
COUNTS: (A) Paley skew spectrum = +-i*sqrt(p) = the GAUSS SUM (Frobenius/Weil eigenvalue = a Lefschetz trace);
(B) dilation phi_v: t->vt has L(phi_v)=1-v = -(resonance count) (the runner-return fixed points)."""
import numpy as np
def legendre(a,p):
    a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def paley_skew(p):  # S_ij = chi(j-i), circulant Paley matrix (skew for p=3 mod4)
    S=np.zeros((p,p))
    for i in range(p):
        for j in range(p): S[i,j]=legendre(j-i,p)
    return S
print("(A) PALEY skew spectrum = the GAUSS SUM +-i*sqrt(p) (Frobenius eigenvalue = Lefschetz trace):")
for p in [3,7,11,19,23]:
    S=paley_skew(p); ev=np.linalg.eigvals(S)
    ims=sorted(set(round(v.imag,4) for v in ev))
    skew_ok = np.allclose(S,-S.T)
    print(f"  p={p:>2} (=3 mod4): Paley matrix skew? {skew_ok}; eigenvalues (imag parts) = {ims}  vs +-sqrt(p)=+-{round(np.sqrt(p),4)}  => GAUSS SUM")
# (B) dilation Lefschetz on the circle
print("\n(B) dilation phi_v: t->vt on S^1: Lefschetz L=1-v = -(# fixed pts = resonances):")
for v in [2,3,7,13,14]:
    L=1-v; fixed=v-1
    print(f"  v={v:>2}: L(phi_v)=1-v={L}; #Fix (v t = t mod1 => t=k/(v-1)) = {fixed} resonances (signed by L)")
# The character-sum: singular-series moment = sum of Ramanujan/Gauss sums = a LEFSCHETZ trace (klein HYP-3793)
print("\n(C) Paley Cayley spectrum (Cayley of the sqrt(p) skew) -- CONCENTRATED (not spread like transitive roots):")
for p in [7,11]:
    S=paley_skew(p).astype(float); U=np.linalg.solve(np.eye(p)+S,np.eye(p)-S)
    ang=sorted(set(round(np.angle(z)/np.pi,3) for z in np.linalg.eigvals(U)))
    print(f"  p={p}: Cayley angles/pi = {ang}  (1 at 0 = fixed vertex; the rest at +-2 arctan(sqrt p)/pi -- the Gauss-sum pair)")
print("\n=> HARD-SIDE TRICK: the free Z2 kills the Brouwer/symmetric-SOS certificate, but the FROBENIUS/dilation")
print("   Lefschetz counts are EXACT arithmetic (Gauss sums sqrt(p), Ramanujan sums, resonance counts 1-v).")
print("   The Borsuk-Ulam obstruction (topological hardness) is DUAL to Weil/Lefschetz exactness (arithmetic).")
