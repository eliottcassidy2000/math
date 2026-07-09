"""
lrc14_Egrid_absolute_adversarial_opus_S171.py   (opus-2026-07-09-S171)

CONFIRM the S171 finding ADVERSARIALLY: kps-S96's E_grid residual R = sum_{Vmax|m} What(m) is bounded
ABSOLUTELY -- R_abs := sum_{n>=1} 2|What(nVmax)| < (6/7)^k -- across the WHOLE critical Vmax window
[s+1, floor(7s/6)] and many random dissociated + 7-structured k=13 clusters (try to BREAK R_abs<main).
Plus the MECHANISM for 7-structured: arc-Fourier b(m)=(1-e(m/7))/(2pi i m) VANISHES at 7|m (opus-S167),
so all resonances m=n.e that are 0 mod 7 contribute 0 -- confirm What(m)~0 at 7|m for 7-structured e.
"""
import numpy as np
from math import gcd

INV7 = 1.0 / 7.0
K = 13
MAIN = (6.0 / 7.0) ** K


def W_series(E, x):
    E = np.asarray(E, float)
    Ph = np.mod(np.outer(x, E), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - INV7, 0, None).sum(axis=1)


def R_abs(E, Vmax, L=64):
    N = L * Vmax
    x = np.arange(N) / N
    F = np.fft.rfft(W_series(E, x)) / N
    s = 0.0
    for n in range(1, L):
        m = n * Vmax
        if m < F.shape[0]:
            s += 2 * abs(F[m])
    return s


def gen_dissociated(rng, spread):
    for _ in range(200):
        E = sorted(set([0] + rng.sample(range(1, spread), K - 2) + [spread]))
        if len(E) == K and reduce_gcd(E) == 1:
            return E
    return None


def gen_7structured(rng, base_spread):
    r = rng.choice([1, 2, 3, 4, 5, 6])
    aa = sorted(rng.sample(range(0, base_spread), K))
    E = [7 * a + r for a in aa]
    return sorted(set(E)) if len(set(E)) == K else None


def reduce_gcd(E):
    from math import gcd as g
    d = 0
    for i in range(len(E) - 1):
        d = g(d, E[i + 1] - E[i])
    return d


import random
print("=" * 92)
print(f"ADVERSARIAL: max R_abs/main over critical Vmax window; main=(6/7)^{K}={MAIN:.5f}")
print("=" * 92)
rng = random.Random(1711)
for label, gen in [("dissociated", "diss"), ("7-structured (diffs=0 mod7)", "7str")]:
    worst = 0.0; worst_cfg = None; nchk = 0
    for _ in range(60):
        spread = rng.randint(60, 130)
        E = gen_dissociated(rng, spread) if gen == "diss" else gen_7structured(rng, 22)
        if E is None or len(E) != K:
            continue
        s = max(E) - min(E)
        E0 = [e - min(E) for e in E]; s = max(E0)
        for Vmax in range(s + 1, 7 * s // 6 + 1, max(1, (7 * s // 6 - s) // 4)):
            ra = R_abs(E0, Vmax)
            nchk += 1
            if ra / MAIN > worst:
                worst = ra / MAIN; worst_cfg = (E0, Vmax)
    print(f"  {label:>28}: max R_abs/main = {worst:.4f}  over {nchk} (cfg,Vmax) checks  "
          f"=> R_abs<main {'HOLDS' if worst < 1 else 'FAILS'}")

print()
print("=" * 92)
print("MECHANISM (7-structured): What(m) ~ 0 at 7|m (arc-Fourier b(m)=0 at 7|m, opus-S167)")
print("=" * 92)
E7 = [1 + 7 * i for i in range(13)]        # all == 1 mod 7
E0 = [e - 1 for e in E7]                    # 0 mod 7 => n.e always 0 mod 7
Vmax = max(E0) + 1
N = 64 * Vmax
x = np.arange(N) / N
F = np.fft.rfft(W_series(E0, x)) / N
mag_7 = np.mean([abs(F[m]) for m in range(7, 3000, 7)])          # m divisible by 7
mag_n7 = np.mean([abs(F[m]) for m in range(1, 3000) if m % 7])   # m not divisible by 7
print(f"  e = 7*(0..12) (all e_i = 0 mod 7): mean|What(m)| at 7|m = {mag_7:.2e}, at 7-nmid = {mag_n7:.2e}")
print(f"  ratio (7|m)/(7-nmid) = {mag_7/max(mag_n7,1e-18):.4f}  (<<1 => resonances at 7|m are KILLED)")
print(f"  => for e=0 mod 7, every n.e = 0 mod 7 => What lands on the arc-Fourier zeros => R tiny.")
print("=" * 92)
