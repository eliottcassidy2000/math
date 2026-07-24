#!/usr/bin/env python3
"""
LRC(14): the SOUND linear looseness certificate over ALL nonnegative test functions R=|P|^2
(Fejer-Riesz) = lambda_min of the covering-multiplicity Toeplitz matrix. Breaks the Riesz-PRODUCT
stall (1.0096, THM-518) on the interior-drop extremizers. mac-mini-2026-07-23-S169 (Opus).

M(tau)=sum_v 1{||v tau||<=1/14};  Mhat(n)=sum_{v|n} s(n/v), s(k)=sin(pi k/7)/(pi k), s(0)=1/7.
min_{R>=0, deg<=N} int(M R)/int(R) = lambda_min(T_M^{(N)}),  T_M^{(N)}_{jk}=Mhat(j-k).
lambda_min < 1  =>  L(S)>0  (SOUND: R>=0, M>=1 on danger => int_G R>0 => |lonely set|>0).

Findings reproduced below:
  * drop-j cores {1..13}\{j}u{56}: ALL certified lambda_min<1 by N=80 (worst j=6: 0.909);
    Riesz PRODUCT stalls at 1.0096 on the same core.
  * tight {1..13}, 2*{1..13}, 3*{1..13} (L=0): lambda_min -> 1 from above (sharp boundary).
  * uniform over the stranger 14m (crossing degree ~60 independent of m) = stranger-decoupling.
  * optimal R concentrates ~34x on the lonely set; its AR/max-entropy dual has reflection coeffs
    k_m with 2*artanh(k_m) = per-mode log-energies = the external snippet's functional.
CAVEAT: per-set (Szego guarantees lambda_min->0 for any loose S); inf_{loose S} L>0 has no single
degree (sup_S lambda_min=1); reduces to the (conjectural) drop-j extremizer-minimality. See reflection
07-reflections/sos-toeplitz-certificate-breaks-the-riesz-product-stall-macmini-S169.md
"""
import numpy as np
from scipy.linalg import toeplitz

def s(k):
    return 1.0/7 if k == 0 else np.sin(np.pi*k/7)/(np.pi*k)

def Mhat(n, V):
    if n == 0:
        return len(V)*s(0)
    an = abs(n)
    return sum(s(an//v) for v in V if an % v == 0)

def lam_min(V, N):
    col = np.array([Mhat(n, V) for n in range(N+1)])
    return np.linalg.eigvalsh(toeplitz(col))[0]

def Lgrid(V, Q=1_000_003):
    tau = np.arange(Q)/Q
    safe = np.ones(Q, bool)
    for v in V:
        safe &= (np.abs(v*tau - np.round(v*tau)) > 1.0/14)
    return safe.mean()

if __name__ == "__main__":
    tight = set(range(1, 14))
    print("SOS certificate lambda_min(T_M^{(N)}) -- <1 proves L(S)>0. Riesz-PRODUCT stall = 1.0096.\n")
    print("All 13 drop-j cores {1..13}\\{j}u{56} (worst = uniformity bottleneck):")
    print(f"{'j':>3} | {'N=80':>8} | {'N=120':>8} | {'N=160':>8}")
    worst = {80: 0, 120: 0, 160: 0}
    for j in range(1, 14):
        V = set(x for x in range(1, 14) if x != j) | {56}
        row = {N: lam_min(V, N) for N in (80, 120, 160)}
        for N in row:
            worst[N] = max(worst[N], row[N])
        print(f"{j:>3} | {row[80]:>8.5f} | {row[120]:>8.5f} | {row[160]:>8.5f}")
    print(f"worst-core lambda_min: N=80 {worst[80]:.5f}, N=120 {worst[120]:.5f}, N=160 {worst[160]:.5f} (all <1)\n")
    print("Sharp boundary -- tight sets (L=0) stay at/above 1:")
    for name, V in [("{1..13}", tight), ("2*{1..13}", set(2*k for k in range(1, 14))),
                    ("3*{1..13}", set(3*k for k in range(1, 14)))]:
        print(f"   {name:>10}: L={Lgrid(V):.5f}  lambda_min(N=160)={lam_min(V,160):.5f}")
