"""mac-mini-S57 (LEM-007): the S-ARC OVERLAP FOURIER MASS LAW + the resonance mechanism.
Derived by expanding uncovered = prod_i(1-h(y-e_i x)), h=1_{[0,1/7)}:
  W-hat(m) = sum_{(n_i): sum n_i=0, sum n_i e_i=m} prod_i ghat_{n_i},  ghat_0=6/7, ghat_n=-(1-e(-n/7))/(2pi i n).
  => L_S-hat(m) = sum_{balanced on S} prod hhat (triple law: n_i d_ik + n_j d_jk = m). VERIFIED.
The (6/7)^{k-2} INACTIVE-ARC DAMPING (from the k-|S| zero entries, each ghat_0=6/7) is the
resonance: Var(W) leading = (6/7)^{2(k-2)} * pair-resonance-kernel ~ (6/7)^{2(k-2)}*R2*V1phi,
cutting the naive R2*V1phi by 15x (the bulk of the 96% IE cancellation). c_k=(6/7)^{2(k-2)}*V1phi."""
import numpy as np
from collections import Counter
L=1/7
def hhat(n): return L if n==0 else (1-np.exp(-2j*np.pi*n*L))/(2j*np.pi*n)
def Lhat_formula(E3,m,NN=60):
    ei,ej,ek=E3;tot=0j
    for ni in range(-NN,NN+1):
        for nj in range(-NN,NN+1):
            nk=-ni-nj
            if ni*ei+nj*ej+nk*ek==m: tot+=hhat(ni)*hhat(nj)*hhat(nk)
    return tot
def triple_ov(E3,GX=1500,GY=3000):
    xs=(np.arange(GX)+0.5)/GX;ys=(np.arange(GY)+0.5)/GY;out=np.zeros(GX)
    for xi,xx in enumerate(xs):
        cov=np.ones(GY,bool)
        for e in E3: cov&=(((ys-(e*xx)%1.0)%1.0)<L)
        out[xi]=cov.mean()
    return xs,out
print("Triple-overlap law L_ijk-hat(m) = sum_{balanced} hhat*hhat*hhat:")
for E3 in [(0,1,3),(1,4,9)]:
    xs,ov=triple_ov(E3)
    ok=all(abs(np.mean(ov*np.exp(-2j*np.pi*m*xs))-Lhat_formula(E3,m))<2e-3 for m in range(4))
    print(f"  E={E3}: E[L]={ov.mean():.5f} (L^3={L**3:.5f}); formula matches m=0..3: {ok}")
V1phi=2/(3*343)-1/49**2
print(f"\nresonance damping: c_k=(6/7)^(2(k-2))*V1phi (k=11) = {(6/7)**18*V1phi:.3e}; empirical ~5.6e-5")
print(f"naive R2*V1phi (block)= {770*V1phi:.4f}; damped (6/7)^18* = {(6/7)**18*770*V1phi:.4f}; true Var(W)=0.047")
