"""
opus-2026-06-07-S706 : cross-correlation is the ADJOINT of convolution, and this
duality unifies the repo's clock/shell/converse/additive-energy objects.

User seed: for complex-valued functions the cross-correlation operator is the adjoint
of the convolution operator. Integrate both, their relation, and how they extend the repo.

CORE IDENTITIES (verified on Z/N with Hermitian inner product <a,b>=sum conj(a)*b):
  convolution        (h*g)(x)   = sum_y h(y) g(x-y)
  cross-correlation  (h@g)(x)   = sum_y conj(h(y)) g(x+y)
  ADJOINT:  <h*a, b> = <a, h@b>   (so C_h^* = R_h)
  REFLECTION LINK:  h@g = (sigma h-bar) * g ,  (sigma f)(x)=f(-x)   [the S702 antipodal involution]

REPO MAP (the synthesis):
  * convolution  1_S * 1_S = SUMSET multiplicities (sums v_i+v_j)  = LRC SHELL face (mod 2n-1, S700)
  * correlation  1_S @ 1_S = DIFFERENCE multiplicities (v_i-v_j)   = LRC CLOCK face (mod n,   S701)
    and 1_S @ 1_S = 1_{sigma S} * 1_S, so clock=shell composed with sigma (S702 antipodal) = the ADJOINT.
  * tournament CONVERSE T->T* = adjoint of the circulant adjacency-convolution; self-converse
    (worry-set, THM-402) = self-adjoint-up-to-multiplier; skew-adjacency S satisfies S*=-S.
  * unit distance U(P) = autocorrelation 1_P@1_P summed on the unit sphere; additive energy
    E(S)=||1_S@1_S||^2 = sum_xi |hat 1_S(xi)|^4  (Wiener-Khinchin / S599 Cayley bridge).
  * all diagonalised by Fourier; ADJOINT = CONJUGATE the symbol (S599g spectral unification).
"""
import numpy as np
from itertools import combinations
from collections import Counter
rng=np.random.default_rng(0)

def conv(h,g):            # circular convolution on Z/N
    return np.fft.ifft(np.fft.fft(h)*np.fft.fft(g))
def corr(h,g):            # circular cross-correlation (h@g)(x)=sum_y conj(h(y)) g(x+y)
    return np.fft.ifft(np.conj(np.fft.fft(h))*np.fft.fft(g))
def ip(a,b):              # Hermitian inner product
    return np.vdot(a,b)   # = sum conj(a)*b

print("="*76)
print("(1) ADJOINT IDENTITY  <h*a,b> = <a,h@b>  on Z/N  (complex, random)")
print("="*76)
for N in [7,12,27]:
    h=rng.standard_normal(N)+1j*rng.standard_normal(N)
    a=rng.standard_normal(N)+1j*rng.standard_normal(N)
    b=rng.standard_normal(N)+1j*rng.standard_normal(N)
    lhs=ip(conv(h,a),b); rhs=ip(a,corr(h,b))
    # reflection link: h@g = (sigma h-bar)*g
    sig_hbar=np.conj(np.roll(h[::-1],1))   # (sigma f)(x)=f(-x): reverse then shift; conj
    link=np.max(np.abs(corr(h,b)-conv(sig_hbar,b)))
    print(f"  N={N:2d}: |<h*a,b>-<a,h@b>| = {abs(lhs-rhs):.2e}   "
          f"|h@b-(sigma h-bar)*b| = {link:.2e}")
print("  => cross-correlation IS the adjoint of convolution; the link is the antipodal sigma (S702).")

print("\n"+"="*76)
print("(2) LRC speed set: convolution=SUMSET (shell, 2n-1) vs correlation=DIFFs (clock, n)")
print("="*76)
def sets_for(S,N):
    ind=np.zeros(N);
    for v in S: ind[v%N]=1
    return ind
# an LRC config (movers), n=6 so clock mod n=6, shell mod 2n-1=11
S=[1,2,3,5,12]; n=6
print(f"  S={S}, n={n}: clock modulus n={n}, shell modulus 2n-1={2*n-1}")
# exact sumset / difference multiplicities (integers)
sums=Counter((a+b) for a in S for b in S)
diffs=Counter((a-b) for a in S for b in S)
print(f"  SUMSET (convolution 1_S*1_S) support sizes: {len(sums)} distinct sums; "
      f"shell-partners a+b==0 mod {2*n-1}: "
      f"{sorted((a,b) for a in S for b in S if a<b and (a+b)%(2*n-1)==0)}")
print(f"  DIFFSET (correlation 1_S@1_S) is symmetric (d and -d): "
      f"{ {d:c for d,c in sorted(diffs.items()) if d>0} }")
# verify convolution gives sumset multiplicities, correlation gives difference mult (mod big N)
N=64; ind=sets_for(S,N)
cv=np.real(conv(ind,ind)); cr=np.real(corr(ind,ind))
chk_sum=all(abs(cv[x%N]-sums.get(x,0))<1e-9 for x in range(-20,40))
chk_dif=all(abs(cr[x%N]-diffs.get(x,0))<1e-9 for x in range(-20,40))
print(f"  conv==sumset-mult: {chk_sum};  corr==diff-mult: {chk_dif}")

print("\n"+"="*76)
print("(3) ADDITIVE ENERGY two ways: ||1_S@1_S||^2 (space) == sum |hat 1_S|^4 (Fourier)")
print("="*76)
for S in [[1,2,3,5,12],[0,1,2,3,4],[1,3,4,9,11]]:
    N=64; ind=sets_for(S,N)
    auto=np.real(corr(ind,ind))                  # autocorrelation = difference mult
    E_space=int(round(np.sum(auto**2)))          # = sum_x r(x)^2, r=#{(a,b):a-b=x}
    F=np.fft.fft(ind); E_four=np.real(np.sum(np.abs(F)**4))/N
    E_quad=sum(1 for a in S for b in S for c in S for d in S if a+b==c+d)  # direct E(S)
    print(f"  S={S}: E(space ||auto||^2)={E_space}  E(Fourier sum|F|^4/N)={E_four:.3f}  "
          f"E(direct #a+b=c+d)={E_quad}")
print("  => additive energy = ||autocorrelation||^2 = 4th Fourier moment (Wiener-Khinchin).")
print("     This is exactly the repo's unit-distance/Cayley additive-energy object (S599).")

print("\n"+"="*76)
print("(4) TOURNAMENT: converse = adjoint of circulant adjacency; self-converse=self-adjoint")
print("="*76)
def round_tour(m):           # Cay(Z/m, H), H={1..(m-1)/2}; A[i,j]=1 iff (j-i) in H
    half=(m-1)//2; H=set(range(1,half+1))
    A=np.zeros((m,m))
    for i in range(m):
        for d in H: A[i,(i+d)%m]=1
    return A,H
for m in [5,7]:
    A,H=round_tour(m)
    AT=A.T                                   # converse = transpose = adjoint (real)
    # converse connection set = -H mod m
    negH=set((-h)%m for h in H)
    # circulant symbol = DFT of first row; eigenvalues of A
    # ADJOINT = CONJUGATE THE SYMBOL: symbol of A is the DFT of the connection indicator 1_H;
    # symbol of A^T (=converse, conn set -H) must equal its conjugate, per frequency.
    w=np.exp(-2j*np.pi/m)
    symA =np.array([sum(w**((h*k)) for h in H) for k in range(m)])
    symAT=np.array([sum(w**(((-h)%m)*k) for h in H) for k in range(m)])
    sym_conj=np.max(np.abs(symAT-np.conj(symA)))
    # skew-adjacency S = A - A^T : check S^* = -S
    Sk=A-A.T; skew=np.max(np.abs(Sk.T+Sk))
    # self-converse via a multiplier lambda: lambda*H == -H (mod m)
    selfconv=any(sorted((l*h)%m for h in H)==sorted(negH) for l in range(1,m) if np.gcd(l,m)==1)
    print(f"  C_{m}: converse conn set -H={sorted(negH)} (vs H={sorted(H)}); "
          f"symbol(A^T)=conj(symbol(A)): {sym_conj:.1e}; "
          f"skew S*=-S: {skew:.1e}; self-converse(multiplier): {selfconv}")
print("  => converse T* = A^T = adjoint of the adjacency-convolution; self-converse worry-set")
print("     (THM-402) = self-adjoint-up-to-multiplier; the {+-1} skew-adjacency is skew-Hermitian.")
print("\nDONE.")
