"""
lrc14_ocf_shell_transport_opus_S172.py   (opus-2026-07-09-S172)

TRANSPORT the tournament Walsh-OCF clean truncation to the covering master law, to make the S171
absolute residual bound R_abs < (6/7)^k A-PRIORI.

LEM-011 (klein-S194) exact closed form (k coords, phi_0=0 pinned, n in Z^{k-1}):
  What(n) = (-1)^r (6/7)^{k-1-r} [prod_{n_i!=0} b0(n_i)] (1[sigma=0] - c(sigma)),
  b0(m)=(e(m/7)-1)/(2pi i m),  c(sigma)=(1-e(-sigma/7))/(2pi i sigma), c(0)=1/7, sigma=sum n_i,
  r=#nonzeros,  |What(n)| <= (6/7)^{k-1-r} prod 1/(pi|n_i|) min(6/7, 1/(pi|sigma|)),  What=0 if 7|n_i.
Geometric per-coord factor (7/6)/pi = 0.371 < 1  => the r-shells decay geometrically (the OCF truncation).

E_grid[W]-(6/7)^k = R = sum_{n!=0, Vmax|n.e} What(n)   (THM-664).  This decomposes R_abs by r-shell:
  R_abs = sum_r R_abs^(r),  R_abs^(r) = sum_{r nonzeros, Vmax|n.e} |What(n)|.
Goal: show R_abs = R_abs^(1)+R_abs^(2)+R_abs^(3)+tail < (6/7)^k, dominated by low r, tail geometric.

This computes, for dissociated / 7-structured / tight-AP k=13 clusters at hard-window Vmax: the r-shell
R_abs^(r) (r=1,2,3) via the exact closed form over |n_i|<=H, the shell ratios (geometric?), the total
vs (6/7)^13, and the MULTI-DIM R_abs vs the 1-D FFT R (S171: does the |.|-per-n bound already beat main,
or is intra-frequency cancellation needed?).
"""
import numpy as np, itertools, cmath
from math import gcd

K = 13
MAIN = (6.0 / 7.0) ** K
TWOPI = 2 * np.pi


def b0(m):
    if m == 0:
        return 1.0 / 7
    return (cmath.exp(1j * TWOPI * m / 7) - 1) / (TWOPI * 1j * m)


def cfun(s):
    if s == 0:
        return 1.0 / 7
    return (1 - cmath.exp(-1j * TWOPI * s / 7)) / (TWOPI * 1j * s)


def What(nz):
    """nz = list of (coord_index_ignored, value) nonzero entries; returns What(n)."""
    r = len(nz)
    sigma = sum(v for _, v in nz)
    prod = 1.0 + 0j
    for _, v in nz:
        prod *= b0(v)
    ind = (1.0 if sigma == 0 else 0.0) - cfun(sigma)
    return ((-1) ** r) * (6.0 / 7) ** (K - 1 - r) * prod * ind


def shell_Rabs(E, Vmax, rmax=3, H=16):
    """R_abs^(r) for r=1..rmax over |n_i|<=H, restricted to Vmax|n.e, 7 not| n_i.  E has e_0=0."""
    e = E[1:]                       # e_1..e_{k-1} (e_0=0 pinned)
    dim = len(e)                    # = k-1 = 12
    vals = [v for v in range(-H, H + 1) if v != 0 and v % 7 != 0]
    shells_abs = {}
    shells_sig = {}
    for r in range(1, rmax + 1):
        tot_abs = 0.0
        tot_sig = 0.0 + 0j
        for coords in itertools.combinations(range(dim), r):
            for combo in itertools.product(vals, repeat=r):
                ne = sum(c * e[i] for c, i in zip(combo, coords))
                if ne % Vmax != 0:
                    continue
                w = What(list(zip(coords, combo)))
                tot_abs += abs(w)
                tot_sig += w
        shells_abs[r] = tot_abs
        shells_sig[r] = tot_sig
    return shells_abs, shells_sig


def W_series(E, x):
    Earr = np.asarray(E, float)
    Ph = np.mod(np.outer(x, Earr), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - 1.0 / 7, 0, None).sum(axis=1)


def fft_R(E, Vmax, L=64):
    N = L * Vmax
    x = np.arange(N) / N
    F = np.fft.rfft(W_series(E, x)) / N
    Rabs = 0.0; Rsig = 0.0
    for n in range(1, L):
        m = n * Vmax
        if m < F.shape[0]:
            Rabs += 2 * abs(F[m]); Rsig += 2 * F[m].real
    return Rabs, Rsig


print("=" * 100)
print(f"OCF SHELL TRANSPORT: R_abs by r-shell (exact LEM-011 closed form).  main=(6/7)^{K}={MAIN:.5f}")
print("=" * 100)
families = {
    "dissociated Sidon   ": [0, 1, 3, 7, 12, 20, 30, 44, 60, 78, 95, 115, 140],
    "dissociated random  ": sorted(set([0] + list(np.random.RandomState(5).choice(range(1, 130), 12, replace=False)))),
    "7-structured =1 mod7": [7 * i + (0 if i == 0 else 1) for i in range(13)],  # 0 then 1 mod7... fix below
    "tight AP {0..12}    ": list(range(13)),
}
# fix 7-structured: all diffs 0 mod 7 => e_i = 7*i (co-offsets), shift so e_0=0
families["7-structured =0 mod7"] = [7 * i for i in range(13)]
del families["7-structured =1 mod7"]

for name, E in families.items():
    E = sorted(set(E))
    if len(E) != 13:
        print(f"  {name}: bad size {len(E)}"); continue
    spread = max(E)
    Vmax = spread + 1
    sh_abs, sh_sig = shell_Rabs(E, Vmax, rmax=3, H=16)
    tot_abs = sum(sh_abs.values())
    tot_sig = sum(sh_sig.values()).real
    fRabs, fRsig = fft_R(E, Vmax)
    print(f"\n  {name}  (spread={spread}, Vmax={Vmax})")
    print(f"    R_abs^(1)={sh_abs[1]:.5f}  R_abs^(2)={sh_abs[2]:.5f}  R_abs^(3)={sh_abs[3]:.5f}  "
          f"ratio 2/1={sh_abs[2]/max(sh_abs[1],1e-12):.3f} 3/2={sh_abs[3]/max(sh_abs[2],1e-12):.3f}")
    print(f"    multi-dim  R_abs(r<=3) = {tot_abs:.5f}  ({tot_abs/MAIN:.3f} main) | R_signed(r<=3) = {tot_sig:.5f}")
    print(f"    1-D FFT    R_abs       = {fRabs:.5f}  ({fRabs/MAIN:.3f} main) | R_signed        = {fRsig:.5f}")
    print(f"    multi-dim R_abs < main? {'YES' if tot_abs < MAIN else 'NO'}   (tight-AP expected NO = the boundary)")

print()
print("  READING: shell ratios <1 => geometric OCF truncation (transport works). multi-dim R_abs<main")
print("  => the a-priori |.|-per-n bound ALONE closes it (no intra-frequency cancellation needed).")
print("  tight AP is the boundary (R_abs ~ main). Dissociated: R_abs << main a-priori.")
print("=" * 100)
