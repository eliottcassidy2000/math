#!/usr/bin/env python3
"""
klein-2026-07-08-S194: the EXACT Fourier transform of the uncovered-measure
function 𝒲 -- the shared 𝒲̂-decay constant (mac-mini THM-664 handoff (a)),
made EXACT and a-priori.

𝒲(phi_1,...,phi_{k-1}) = uncovered measure of the phase set {0, phi_1,...,phi_{k-1}}
                       = int_0^1 prod_{i=0}^{k-1} (1 - 1[phi_i in (t-1/7, t)]) dt
(phi_0 = 0 pinned).  W = Sigma_gaps (gap - 1/7)_+ = this integral (a point t is
'uncovered' iff the 1/7-arc just before it holds no phase).

DERIVED CLOSED FORM (to verify):
  𝒲̂(n) = (-1)^r (6/7)^{(k-1)-r} [prod_{i: n_i != 0} b0(n_i)] * (1[sigma=0] - c(sigma))
  r = #{i: n_i != 0},  sigma = sum_i n_i,
  b0(m) = (e(m/7)-1)/(2 pi i m)  (b0(0)=1/7),  |b0(m)| = |sin(pi m/7)|/(pi|m|),
  c(s)  = (1 - e(-s/7))/(2 pi i s)  (c(0)=1/7),  |c(s)| = |sin(pi s/7)|/(pi|s|).
Sanity: n=0 => (6/7)^{k-1} * (1 - 1/7) = (6/7)^k  (mac-mini's 𝒲̂(0)).  Decay:
prod 1/(pi|n_i|), zero if any 7|n_i, extra 1/|sigma| unless sigma=0.

This script: (1) builds 𝒲 on a grid for k-1=2,3 phases, FFTs it, and compares
𝒲̂(n) to the closed form; (2) checks Parseval sum|𝒲̂|^2 = int 𝒲^2 = E[W^2].
"""
import numpy as np
import cmath

INV7 = 1.0/7.0

def W_uncovered(phases):
    """uncovered measure of a phase set on the circle (Sigma (gap-1/7)_+)."""
    p = np.sort(np.mod(np.asarray(phases,float),1.0))
    g = np.diff(p); g = np.append(g, 1.0 - p[-1] + p[0])
    return np.maximum(g - INV7, 0.0).sum()

def b0(m):
    if m == 0: return 1.0/7.0
    return (cmath.exp(2j*np.pi*m/7.0) - 1.0)/(2j*np.pi*m)

def c_of(s):
    if s == 0: return 1.0/7.0
    return (1.0 - cmath.exp(-2j*np.pi*s/7.0))/(2j*np.pi*s)

def What_formula(n, k):
    n = np.asarray(n)
    nz = [int(v) for v in n if v != 0]
    r = len(nz)
    sigma = int(n.sum())
    prod = 1.0+0j
    for v in nz: prod *= b0(v)
    bal = (1.0 if sigma == 0 else 0.0) - c_of(sigma)
    return ((-1)**r) * (6.0/7.0)**((k-1)-r) * prod * bal

def build_W_grid(kphases, N):
    """𝒲 on N^{kphases} grid; kphases = k-1."""
    axes = [ (np.arange(N))/N for _ in range(kphases) ]
    grids = np.meshgrid(*axes, indexing='ij')
    shape = grids[0].shape
    Wg = np.zeros(shape)
    it = np.nditer(grids[0], flags=['multi_index'])
    for _ in it:
        idx = it.multi_index
        ph = [0.0] + [grids[d][idx] for d in range(kphases)]
        Wg[idx] = W_uncovered(ph)
    return Wg

for k in (3, 4):
    kp = k-1
    N = 84 if kp==2 else 42   # multiple of 7 helps; keep k=4 grid modest
    print(f"\n===== k={k}  ({kp} phase vars), grid N={N} =====")
    Wg = build_W_grid(kp, N)
    Fh = np.fft.fftn(Wg)/(N**kp)     # Fh[n] = (1/N^kp) sum W e(-n.x)  ~ 𝒲̂(n)
    # compare on a set of small n
    def freq(idx): return tuple(int(v) if v<=N//2 else int(v)-N for v in idx)
    maxerr=0.0; checks=0; samples=[]
    rng = range(-3,4)
    import itertools
    for nn in itertools.product(rng, repeat=kp):
        # locate FFT index
        fidx = tuple(v % N for v in nn)
        num = Fh[fidx]
        ana = What_formula(nn, k)
        err = abs(num-ana)
        maxerr=max(maxerr,err); checks+=1
        if abs(nn[0])<=2 and (kp<3 or abs(nn[1])<=1):
            if len([x for x in nn if x!=0])<=2:
                samples.append((nn, num, ana, err))
    print(f"  checked {checks} freqs |n_i|<=3;  MAX |FFT - formula| = {maxerr:.2e}")
    print(f"  W_hat(0) FFT={Fh[tuple([0]*kp)].real:.6f}  formula={What_formula([0]*kp,k).real:.6f}  (6/7)^k={(6/7)**k:.6f}")
    # Parseval
    lhs = (np.abs(Fh)**2).sum()
    rhs = (Wg**2).mean()
    print(f"  Parseval: sum|W_hat|^2 = {lhs:.6f}   int W^2 = {rhs:.6f}   (E[W^2] moment)")
    for (nn,num,ana,err) in samples[:5]:
        print(f"    n={nn}: FFT={num:+.5f}  formula={ana:+.5f}")
print("\nIf MAX err ~ grid-discretization (small) and W_hat(0)=(6/7)^k, the closed form is CONFIRMED.")
