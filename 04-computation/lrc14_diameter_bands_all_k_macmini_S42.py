"""
mac-mini-2026-07-07-S42 (HYP-4838 pending) -- THE DIAMETER BANDS FOR EVERY LEG k=8..13:
extending monad-S2's tail-diameter theorem (HYP-4817, k=13) to the k=8..12 G_P-coupled
legs via THM-530's union bound.

ARCHITECTURE. The k<13 legs consume G2(P,E) = meas(G_P ∩ Good_E) >= m_P (hlarge).
Union bound: G2 >= meas(G_P) + mu_{1/7}(E) - 1 >= mu_{1/7}(E) - thr_k, where
thr_k = 1 - min_{|P|=13-k} meas(G_P) (THM-530 exact minima, re-verified by monad-S1).
So    mu_{1/7}(E) >= thr_k + m_P   ==>   G2 >= m_P   (the hlarge bar).
Tail-diameter (monad-S2 / kps-S59): mu_{1/7}(E) >= mu_{1/7}(AP_{D+1}), D = diam(E).
==> DIAMETER BAND for leg k:  D0(k) := max{D : mu(AP_{D+1}) >= thr_k + m_P};
    every k-point cluster with diam <= D0(k) discharges its hlarge case.
    (k-point cluster has diam >= k-1, so the band is [k-1, D0(k)] -- may be empty.)

This script: (1) exact bars thr_k + m_P as Fractions; (2) exact mu(AP_m) crossings via
monad's Farey-cell engine (imported); (3) the exact band table; (4) numeric residual
probes: adversarial min mu over k-point families with diam > D0(k) -- the per-k residual
slack; (5) the k=8 razor check exact.
"""
import sys, random
sys.path.insert(0, '04-computation')
from fractions import Fraction as F
import numpy as np

from lrc14_tail_diameter_theorem_monad_S2 import mu_ap_exact, M_P  # exact engine (credit monad-S2)

# THM-530 / monad-S1 exact per-|P| minima of meas(G_P), |P| = 1..5 (k = 12..8)
GP_MIN = {1: F(6,7), 2: F(66,91), 3: F(55,91), 4: F(1979,4004), 5: F(2243,5880)}
THR = {13: F(0)}
for psz, gp in GP_MIN.items():
    THR[13-psz] = 1 - gp

BAR = {k: THR[k] + M_P for k in THR}   # mu needed: union bound delivers G2 >= mu - thr_k >= m_P

print("=== exact bars: mu_{1/7}(E) >= thr_k + m_P discharges leg k via union bound ===")
for k in sorted(BAR):
    print(f"  k={k:2d}: thr={THR[k]} = {float(THR[k]):.5f}   bar=thr+m_P = {float(BAR[k]):.5f}")

print("\n=== exact mu(AP_m) ladder (monad-S2 engine) and the per-k crossings ===")
mu_cache = {}
def mu(m):
    if m not in mu_cache:
        mu_cache[m] = mu_ap_exact(m)
    return mu_cache[m]

# compute the ladder over the relevant range
for m in range(8, 36):
    v = mu(m)
    tags = [f"k={k}" for k in sorted(BAR) if k <= 12 and v >= BAR[k]]
    print(f"  mu(AP_{m:2d}) = {float(v):.5f}   clears bars: {','.join(tags) if tags else '-'}")

print("\n=== THE DIAMETER BAND TABLE ===")
print(f"{'k':>3s} {'bar':>8s} {'D0(k)':>6s} {'band [k-1, D0]':>15s} {'band size':>9s}")
bands = {}
for k in range(8, 13):
    D0 = None
    for m in range(8, 60):
        if mu(m) >= BAR[k]:
            D0 = m - 1
        else:
            if D0 is not None and m > k:   # ladder is decreasing past small m
                break
    bands[k] = D0
    lo = k - 1
    if D0 is not None and D0 >= lo:
        print(f"{k:3d} {float(BAR[k]):8.5f} {D0:6d}   [{lo},{D0}]        {D0-lo+1:9d}")
    else:
        print(f"{k:3d} {float(BAR[k]):8.5f} {'--':>6s}   EMPTY BAND {'':9s}")
# k=13 for reference (monad-S2): D0 = 75
print(f" 13 {float(BAR[13]):8.5f} {75:6d}   [12,75]      {75-12+1:9d}   (monad-S2 exact)")

print("\n=== k=8 razor check (exact fractions) ===")
for m in (8,9,10,11):
    v = mu(m)
    print(f"  mu(AP_{m}) = {v} = {float(v):.6f}  {'>=':>2s} bar_8={float(BAR[8]):.6f}? {v >= BAR[8]}")

# ---------------- residual probes ----------------
print("\n=== residual probes: adversarial min mu over k-point families with diam > D0(k) ===")
GRID = 60000
xs = (np.arange(GRID)+0.5)/GRID
def mu_np(E):
    ph = np.mod(np.outer(xs, np.array(E,float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return float((g.max(axis=1) > 1/7).mean())

random.seed(42)
for k in range(8, 13):
    D0 = bands[k]
    if D0 is None or D0 < k-1:
        # band empty: residual = ALL diameters; probe from diam = k-1 up
        Dmin = k - 1
    else:
        Dmin = D0 + 1
    best, bestE = 9.9, None
    for trial in range(8):
        # start: stretched k-point sets with diam in [Dmin, 3*Dmin]
        D = random.randint(Dmin, max(Dmin+10, 3*Dmin))
        interior = sorted(random.sample(range(1, D), k-2))
        E = [0] + interior + [D]
        cur = mu_np(E)
        for step in range(300):
            i = random.randrange(1, k-1)   # keep endpoints (diam pinned >= Dmin)
            cand = E.copy(); cand[i] = random.randint(1, E[-1]-1)
            cand = sorted(set(cand))
            if len(cand) != k or (cand[-1]-cand[0]) < Dmin: continue
            cv = mu_np(cand)
            if cv < cur: E, cur = cand, cv
        if cur < best: best, bestE = cur, E
    slack = best/float(BAR[k])
    print(f"  k={k:2d} residual diam>={Dmin:3d}: adversarial min mu ~ {best:.4f} "
          f"(bar {float(BAR[k]):.4f}, {slack:.2f}x) at {bestE}")
print("\nNOTE: residual mins are NUMERIC lower-envelope probes (grid 60k, 8x300 descent);")
print("bands are EXACT (Farey-cell engine).")
