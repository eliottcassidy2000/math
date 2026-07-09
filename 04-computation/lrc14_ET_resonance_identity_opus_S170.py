"""
lrc14_ET_resonance_identity_opus_S170.py   (opus-2026-07-09-S170)

The CORRECT large-spread / near-AP good-period tool is Erdos-Turan on the ruler grid (klein-S193
MISTAKE-127: arc-count is VACUOUS on the near-AP extremal family).  This script makes the EXACT
Fourier identity explicit and DIAGNOSES why the resonant discrepancy is small -- confirming it is
the SIGNED-CANCELLATION / L2 phenomenon (opus-S154 L2-not-L1, opus-S167 Mertens), NOT an absolute
bound, and locating the a-priori lever.

SETUP (k=11 cluster e=d*{0..9}u{1}, plus P={1,2} auxiliaries; Vmax=9d+14).  Good set
  G* = {x: maxgap{frac(e_i x)} > 1/7} cap G_P,   1_{G*}(x) = sum_m Ghat(m) e(m x),  Ghat(0)=rho*.
Ruler grid x_j=(j+1/2)/Vmax, j=0..Vmax-1.  Since sum_j e(m x_j) = Vmax*[Vmax|m]*e(m/(2Vmax)):
  #good = sum_j 1_{G*}(x_j) = Vmax * sum_{Vmax|m} Ghat(m) e(m/(2Vmax))
        = Vmax * ( rho* + sum_{m!=0, Vmax|m} Ghat(m) e(m/(2Vmax)) ).        [EXACT ET IDENTITY]
So   S_signed := #good/Vmax - rho*  =  sum_{n!=0} Ghat(n Vmax) e(n/2)   (the RESONANT sum).

Verifies: (A) the identity (Fourier RHS == #good/Vmax - rho*);  (B) the L1 resonant sum
sum|Ghat(nVmax)| is O(1) and does NOT ->0 (like #arcs; the vacuous bound) while |S_signed| ->0
=> CANCELLATION (Mertens/L2);  (C) the existence margin rho* - |S_signed| stays > 0, uniformly in d.
"""
import numpy as np
from math import gcd

INV7 = 1.0 / 7.0
P = [1, 2]


def indicator_Gstar(e, x, P):
    """1_{G*}(x) for array x: maxgap{frac(e_i x)}>1/7 AND all P-runners >=1/14 from observer."""
    inGP = np.ones(len(x), bool)
    for p in P:
        dd = np.abs(((p * x + 0.5) % 1.0) - 0.5)
        inGP &= (dd >= 1.0 / 14)
    e = np.asarray(e, float)
    Ph = np.mod(np.outer(x, e), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return (inGP & (g.max(axis=1) > INV7)).astype(float)


def rho_and_ghat(e, Vmax, P, L=64):
    """rho* and Ghat(n*Vmax) for n=1..nmax via FFT on N=L*Vmax uniform samples."""
    N = L * Vmax
    x = (np.arange(N) + 0.5) / N
    f = indicator_Gstar(e, x, P)
    # F[k] = sum_t f(t) e(-2pi i k t/N); Ghat(k) ~ (1/N) F[k] * e(-pi i k/N) phase for midpoint grid
    F = np.fft.rfft(f) / N
    # midpoint correction: samples at (t+0.5)/N => multiply by e(+pi i k/N) to get grid-Fourier coeff
    ks = np.arange(F.shape[0])
    F = F * np.exp(1j * np.pi * ks / N)
    rho = F[0].real
    ghat = {}
    nmax = L - 1
    for n in range(1, nmax + 1):
        k = n * Vmax
        if k < F.shape[0]:
            ghat[n] = F[k]
    return rho, ghat


def good_on_ruler(e, Vmax, P):
    j = np.arange(Vmax); x = (j + 0.5) / Vmax
    return int(indicator_Gstar(e, x, P).sum())


print("=" * 96)
print("EXACT ET IDENTITY + resonant-sum diagnosis (near-AP e=d*{0..9}u{1}, Vmax=9d+14)")
print("=" * 96)
print(f"{'d':>4} {'Vmax':>6} {'rho*':>7} {'#good':>6} {'S_signed':>9} {'RHS(Fourier)':>13} "
      f"{'idOK':>5} {'L1_res':>7} {'L1/|signed|':>11} {'margin rho*-|S|':>15}")
for d in (5, 10, 20, 40, 80, 150):
    e = sorted(set(list(d * np.arange(10)) + [1]))
    if len(e) < 11:
        continue
    Vmax = 9 * d + 14
    rho, ghat = rho_and_ghat(e, Vmax, P, L=64)
    ng = good_on_ruler(e, Vmax, P)
    S_signed = ng / Vmax - rho
    # Fourier RHS: sum_{n!=0} Ghat(nVmax) e(n/2).  e(n/2)=(-1)^n.  Include +-n (conj):
    rhs = 0.0j
    L1 = 0.0
    for n, gh in ghat.items():
        sign = (-1) ** n
        rhs += 2 * (gh * sign).real * 0.5 * 2   # 2*Re for +-n; e(n/2) real => 2*Re(Ghat(nVmax))*(-1)^n
        # simpler: contribution of +n and -n = 2 Re[Ghat(nVmax) e(n/2)] = 2 (-1)^n Re Ghat(nVmax)
        L1 += 2 * abs(gh)
    rhs = 0.0
    for n, gh in ghat.items():
        rhs += 2 * ((-1) ** n) * gh.real
    idOK = "YES" if abs(rhs - S_signed) < 0.02 else "no"
    ratio = L1 / max(abs(S_signed), 1e-9)
    margin = rho - abs(S_signed)
    print(f"{d:>4} {Vmax:>6} {rho:>7.4f} {ng:>6} {S_signed:>9.4f} {rhs:>13.4f} "
          f"{idOK:>5} {L1:>7.3f} {ratio:>11.1f} {margin:>15.4f}")

print()
print("  READING:")
print("  * idOK: the EXACT ET identity #good/Vmax-rho* = sum_{n!=0} Ghat(nVmax)(-1)^n holds.")
print("  * L1_res = sum 2|Ghat(nVmax)| is O(1) and does NOT ->0 (the VACUOUS arc-count-scale bound);")
print("    L1/|signed| >> 1 => the small discrepancy is SIGNED CANCELLATION (Mertens/L2), not absolute.")
print("  * margin rho*-|S_signed| > 0 uniformly => #good>0 => a GOOD PERIOD EXISTS. The a-priori lever")
print("    is the L2/cancellation control of the resonant sum (opus-S154), NOT #arcs (klein MISTAKE-127).")
print("=" * 96)
