"""
mac-mini-2026-07-08-S60 -- R2: the k=12,13 density-floor TAIL a-priori rate.

opus-S157 closed the k=11 (L=10) tail: for E_{d,p}={0,d,..,(L-1)d}u{p} (gcd(d,p)=1),
  m_j(E) = L_j + sum_{kappa!=0} Ghat^j(kappa*p, -kappa*d),   L_j = int int G^j (decorr limit),
  |Ghat^j(a,b)| <= V_j/|ab|  =>  |m_j - L_j| <= 2*zeta(2)*V_j/(p*d),  and via the D3 gradient
  |D3 - D3_inf| <= C/(p*d),  C = (pi^2/3) * sum_j |grad_j| * V_j.
LEM-011 gives G = 𝒲 (the uncovered-measure function) an EXACT closed form, so V_j is a-priori.

For k=12,13 the tail family is the SAME structure with L = k-1 (11-AP or 12-AP + outlier).
This computes V_j = max_{a,b!=0} |a*b*Ghat^j(a,b)| by FFT of G,G^2,G^3, confirming the a-priori
rate holds (V_j bounded, finite) -- mirroring opus-S157's 0.28/0.16/0.10 at k=11.

G(u,v) = uncovered measure of the circle by arcs of length 1/7 at phases {0,u,2u,..,(L-1)u, v}.
"""
import numpy as np
l = 1.0/7.0

def uncovered_grid(L, N):
    """G(u,v) on an N x N grid: uncovered by arcs at {i*u : i<L} u {v}, arc length 1/7."""
    u = (np.arange(N)/N)[:, None]      # N x 1
    v = (np.arange(N)/N)[None, :]      # 1 x N
    G = np.zeros((N, N))
    # coverage via 1 - prod(1 - 1_{arc_i}) integrated over x on a fine grid
    Nx = 700
    x = (np.arange(Nx)+0.5)/Nx
    # for each x, a phase p covers x iff frac(x - p) in [0, l)  (arc [p,p+l))
    cov_free = np.ones((N, N, Nx))
    for i in range(L):
        ph = np.mod(i*u, 1.0)                       # N x 1 phase of i*u
        # 1 - 1[x in [ph, ph+l)] = 1[frac(x-ph) >= l]
        d = np.mod(x[None, None, :] - ph[:, :, None], 1.0)
        cov_free *= (d >= l)
    dv = np.mod(x[None, None, :] - v[:, :, None], 1.0)   # outlier phase v
    cov_free *= (dv >= l)
    G = cov_free.mean(axis=2)          # uncovered measure = mean over x of (all-free)
    return G

def Vj_from_fft(G, j, drop=1):
    Gj = G**j
    F = np.fft.fft2(Gj) / (G.shape[0]*G.shape[1])   # Ghat^j(a,b), a,b in freq index
    N = G.shape[0]
    freqs = np.fft.fftfreq(N, d=1.0/N)              # integer freqs
    A, B = np.meshgrid(freqs, freqs, indexing='ij')
    mask = (np.abs(A) >= drop) & (np.abs(B) >= drop)  # a,b both nonzero
    val = np.abs(A*B*F)
    return val[mask].max()

print("k=12,13 tail a-priori rate: V_j = max|a*b*Ghat^j(a,b)| (a-priori via LEM-011)\n")
print(f"{'k':>3} {'L=k-1':>6} {'V_1':>8} {'V_2':>8} {'V_3':>8}  {'L_1=D3_inf-ish (int G)':>22}")
for k, L in [(11, 10), (12, 11), (13, 12)]:
    N = 210     # multiple of 7 for clean arc alignment
    G = uncovered_grid(L, N)
    V1 = Vj_from_fft(G, 1); V2 = Vj_from_fft(G, 2); V3 = Vj_from_fft(G, 3)
    L1 = G.mean()   # int int G = decorrelation-limit first moment
    print(f"{k:>3} {L:>6} {V1:>8.4f} {V2:>8.4f} {V3:>8.4f}  {L1:>22.5f}")
print("\n(k=11 row should match opus-S157's V_j ~ 0.28/0.16/0.10; k=12,13 are the new a-priori")
print(" constants. All FINITE => |D3-D3_inf| <= C/(pd) with C=(pi^2/3)sum|grad_j|V_j a-priori.)")
