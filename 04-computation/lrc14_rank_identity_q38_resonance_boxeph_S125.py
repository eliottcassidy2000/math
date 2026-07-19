#!/usr/bin/env python3
"""
boxeph-2026-07-19-S125 — pushing codex-S16's rank identity onto q=38, via the resonance debt.

Owner: push the rank identity onto q=38, think resonance hypothesis.

codex-S16 rank identity (exact inclusion-exclusion): for danger combs D_v(lam)={t:||vt||<=lam},
active set S(t)={v: t in D_v}, r(t)=max(|S(t)|-1,0),
    mu(uncovered) = mu([0,1)) - sum_v mu(D_v) + integral r(t) dt.
mu(D_v(lam)) = 2*lam for every v. At lam=3/38 with 12 speeds: sum_v mu(D_v)=12*6/38=36/19.
So COVERING at 3/38 (mu(uncovered)=0, the hole is the single point t*) forces the COVER-DEBT
    integral r(t) dt = 36/19 - 1 = 17/19   (the overlap that must be supplied).

RESONANCE DEBT (opus-S531): mu(LONELY at lam) = sum_{k: sum k_i v_i = 0} prod ghat(k_i), with
ghat(0)=1-2lam and ghat(k) = sin(2 pi k lam)/(pi k). At lam=3/38: ghat(k)=sin(3 pi k/19)/(pi k),
which VANISHES at k==0 mod 19 -- so the 3/38 comb spectrum lives in the 1/19 alphabet (apex-19 of
38=2*19; parallels LRC(14)'s 1/7). The debt is written in multiples of 1/19; mod-19-aligned
resonances contribute NOTHING.

RESONANCE HYPOTHESIS (3/38): the cover-debt 17/19 must be supplied by pairwise (and higher)
resonances. This script:
 (A) verifies the rank-identity accounting at 3/38 (cover-debt 17/19) on a fine grid;
 (B) confirms the apex-19 spectral vanishing;
 (C) localizes the cover-debt: |S(t)| peaks at small-q resonances (t=p/q, q small) -- the debt
     concentrates at the MEDIUM/small moduli, NOT at q=38, which is exactly where the competing
     deeper holes (S124's needles) form;
 (D) the resonance-debt ratio |DEBT|/CREDIT for {1..12} at its critical 1/13 vs band-filled
     families -- the AP is the unique clean debt=credit (mu=0) config, and it sits at 1/13.
"""
from fractions import Fraction as F
from math import sin, pi, gcd


def dist_int(x):
    return abs(x - round(x))


def active(sp, t, lam):
    return [v for v in sp if dist_int(v * t) <= lam + 1e-12]


LAM = 3.0 / 38.0
base = list(range(1, 13))
band = [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35]   # covering, band-filled mod 38, pair (3,35)

print("=" * 74)
print("(A) rank identity at lam=3/38: cover-debt = sum mu(D_v) - 1 (when covered)")
print("=" * 74)
N = 200000
for name, sp in [("{1..12}", base), ("band-filled", band)]:
    # fine-grid estimate of mu(union) and integral r
    grid_pts = N
    cov = 0
    sumS = 0
    for i in range(grid_pts):
        t = (i + 0.5) / grid_pts
        k = len(active(sp, t, LAM))
        if k >= 1:
            cov += 1
        sumS += k
    mu_union = cov / grid_pts
    int_r = sumS / grid_pts - mu_union          # integral(|S|-1 over covered) = integral|S| - mu(union)
    sum_mu = 12 * 2 * LAM
    uncovered = 1 - mu_union
    print(f"  {name:>12}: sum mu(D_v)={sum_mu:.4f} (=36/19={36/19:.4f})  mu(union)={mu_union:.4f}  "
          f"uncovered={uncovered:.4f}  int r={int_r:.4f}  (17/19={17/19:.4f})")

print()
print("=" * 74)
print("(B) apex-19: the 3/38-comb spectrum ghat(k)=sin(3*pi*k/19)/(pi*k) vanishes at k==0 mod 19")
print("=" * 74)
for k in [1, 2, 6, 12, 18, 19, 20, 38, 19 * 2]:
    val = sin(3 * pi * k / 19) / (pi * k)
    print(f"    ghat({k:>2}) = {val:+.6f}   {'<-- ZERO (k==0 mod 19)' if k % 19 == 0 else ''}")

print()
print("=" * 74)
print("(C) cover-debt localization: |S(t)| at resonances t=p/q (small q) -- where the debt lives")
print("=" * 74)
print("    (a high |S| at t=p/q means many combs coincide there = a resonance carrying overlap)")
for name, sp in [("{1..12}", base), ("band-filled", band)]:
    peaks = []
    for q in range(2, 40):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            t = p / q
            k = len(active(sp, t, LAM))
            if k >= 3:
                peaks.append((k, q, p))
    peaks.sort(reverse=True)
    top = peaks[:8]
    print(f"  {name:>12}: top overlap resonances (|S|, q, p): "
          f"{[(k, q, p) for k, q, p in top]}")
    # where is the DEEPEST hole (min over t of |S|=0 region)? report smallest q giving an empty point
    print(f"                max overlap |S|={peaks[0][0] if peaks else 0} at q={peaks[0][1] if peaks else '-'}")

print()
print("=" * 74)
print("(D) resonance-debt ratio |DEBT|/CREDIT (opus-S531) at radius lam, Fourier truncation")
print("=" * 74)


def ghat(k, lam):
    if k == 0:
        return 1 - 2 * lam
    return sin(2 * pi * k * lam) / (pi * k)


def lonely_measure(sp, lam, K=60):
    # mu(LONELY) = sum_{k: sum k_i v_i = 0} prod ghat(k_i); truncate |k_i|<=K, up to pairwise+triple
    # r=0 (credit):
    n = len(sp)
    credit = ghat(0, lam) ** n
    # r=2 pairwise debt: for each pair i<j and integer a!=0 with a*v_i - a'*v_j... simplest resonance
    # k_i v_i + k_j v_j = 0 => k_i = m*v_j/g, k_j=-m*v_i/g, g=gcd(v_i,v_j)
    debt = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            g = gcd(sp[i], sp[j])
            ki0, kj0 = sp[j] // g, -sp[i] // g
            for m in range(1, K):
                ki, kj = m * ki0, m * kj0
                if abs(ki) > K or abs(kj) > K:
                    break
                term = ghat(0, lam) ** (n - 2) * ghat(ki, lam) * ghat(kj, lam)
                debt += 2 * term   # +m and -m
    return credit, debt


for name, sp in [("{1..12}", base), ("band-filled", band), ("{1..11,24}", list(range(1, 12)) + [24])]:
    for lam, lname in [(1 / 13, "1/13"), (3 / 38, "3/38"), (2 / 25, "2/25")]:
        credit, debt = lonely_measure(sp, lam)
        ratio = abs(debt) / credit if credit else float('nan')
        print(f"  {name:>12} at lam={lname:>5}: credit={credit:.5f} pairwise-debt={debt:+.5f} "
              f"|debt|/credit={ratio:.4f}")
    print()
print("DONE.")
