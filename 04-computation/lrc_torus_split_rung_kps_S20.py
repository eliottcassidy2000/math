#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20 (HYP-4247): THE COUPLED-TORUS SPLIT RUNG -- verification.

Coupled 2-torus system (the J-K lift-limit object, HYP-4262 dictionary):
  base speeds w_i (teeth ||w_i t||), lifted slopes r_i > 0 with couplings a_i
  (teeth ||r_i theta + a_i t||).  Value M2 = sup_{(t,theta)} min over all teeth.

THE RUNG: if the system is rho-covered at every (t, theta) with rho <= 1/12,
then 2*rho*(#lifted) >= 1.
  Proof shape: citation on the base gives t0 with base margins >= 1/12 > rho;
  at t0 the lifted combs alone must rho-cover the theta-circle; on a D-grid
  each lifted comb hosts <= (2 rho) D + 3 r_i points (the rho-parametric sharp
  visit count); D <= sum => D (1 - 2 rho l) <= 3 sum(r) => contradiction for
  D large unless 2 rho l >= 1.

COROLLARY (the (G)2 l<=6 kill): at rho = 2/25, every coupled 2-torus system
with l <= 6 lifted has a 2/25-clear point (t, theta) => M2 >= 2/25 => the
open gap (1/13, 2/25) is EMPTY of <= 6-lifted coupled values.

This script verifies:
  (1) the rho-parametric count (2 rho) D + 3 w on random inhomogeneous combs;
  (2) the corollary: random coupled systems with l <= 6 always have a
      2/25-clear point (searched constructively);
  (3) the wall arithmetic at 2/25 and 3/38 (l >= 7 both);
  (4) sharpness data at l = 7: shifted distinct-frequency comb systems'
      minimum theta-uncovered fraction (can 7 lifted combs theta-cover?).
"""
import random
from fractions import Fraction

random.seed(20260706)

def dist_to_int(x):
    f = x - int(x)
    if f < 0: f += 1.0
    return min(f, 1.0 - f)

# ---------- (1) the rho-parametric count ----------
def check_count(trials=20000):
    worst = -1e9; bad = 0
    for _ in range(trials):
        w = random.randint(1, 40)
        D = random.randint(1, 300)
        rho = random.uniform(0.005, 0.5)
        c = random.uniform(-3, 3)
        cnt = 0
        for j in range(D):
            x = c + w * j / D
            if dist_to_int(x) < rho:
                cnt += 1
        bound = 2 * rho * D + 3 * w
        slack = bound - cnt
        if slack < worst or worst < -1e8: worst = slack if worst < -1e8 else min(worst, slack)
        if cnt > bound + 1e-9: bad += 1
    return bad, worst

bad, worst = check_count()
print(f"(1) rho-parametric count (2rho)D+3w: violations = {bad}/20000 (min slack {worst:.3f})")
assert bad == 0

# ---------- (2) the l<=6 clear-point corollary at rho = 2/25 ----------
RHO = 2.0 / 25.0

def clear_point_search(base, lifted, rho=RHO, tgrid=4000, thgrid=4000):
    """Search (t, theta) with all base ||w t|| >= rho and lifted ||r th + a t|| >= rho.
    Constructive echo of the proof: t from base-citation candidates, theta swept."""
    # candidate t: maximize base margin over a fine grid (citation guarantees >= 1/12)
    best_t, best_m = None, -1.0
    for jt in range(tgrid):
        t = (jt + 0.5) / tgrid
        m = min(dist_to_int(w * t) for w in base) if base else 1.0
        if m > best_m:
            best_m, best_t = m, t
    if best_m < rho:
        return None  # base itself rho-covered everywhere (shouldn't happen, citation)
    t = best_t
    for jth in range(thgrid):
        th = (jth + 0.5) / thgrid
        ok = all(dist_to_int(r * th + a * t) >= rho for (r, a) in lifted)
        if ok:
            return (t, th, best_m)
    return None

fails = 0; trials2 = 400
for _ in range(trials2):
    nb = random.randint(0, 6)
    nl = random.randint(1, 6)
    if nb + nl > 12: nb = 12 - nl
    base = random.sample(range(1, 60), nb) if nb else []
    lifted = [(random.randint(1, 30), random.randint(-30, 30)) for _ in range(nl)]
    res = clear_point_search(base, lifted)
    if res is None:
        fails += 1
        print("   FAIL:", base, lifted)
print(f"(2) l<=6 clear-point (2/25) on random coupled systems: {trials2 - fails}/{trials2} found, {fails} failures")
assert fails == 0

# ---------- (3) wall arithmetic ----------
for c, q in [(2, 25), (3, 38)]:
    rho = Fraction(c, q)
    lmin = 1
    while 2 * rho * lmin < 1:
        lmin += 1
    print(f"(3) wall at rho = {c}/{q}: 2*rho*l >= 1 forces l >= {lmin}  (2*{c}/{q}*6 = {float(2*rho*6):.4f} < 1: {2*rho*6 < 1})")
    assert lmin == 7

# ---------- (4) sharpness at l = 7: can 7 lifted combs theta-cover? ----------
def uncovered_fraction(lifted, t, rho=RHO, grid=20000):
    unc = 0
    for j in range(grid):
        th = (j + 0.5) / grid
        if all(dist_to_int(r * th + a * t) >= rho for (r, a) in lifted):
            unc += 1
    return unc / grid

# random 7-lifted systems with random couplings: min over sampled t of covered?
best_cover = 1.0; cover_hits = 0
for trial in range(200):
    lifted = [(random.randint(1, 25), random.randint(-25, 25)) for _ in range(7)]
    # a system "covers at t0" if uncovered fraction ~ 0
    t = random.random()
    uf = uncovered_fraction(lifted, t, grid=4000)
    if uf < best_cover: best_cover = uf
    if uf == 0.0: cover_hits += 1
print(f"(4) l=7 random systems: min theta-uncovered fraction over 200 samples = {best_cover:.5f}; full covers seen = {cover_hits}")

# an engineered cover attempt: frequencies 1..7 with tuned shifts (greedy)
freqs = [1, 2, 3, 4, 5, 6, 7]
def try_engineer(seed_shifts=None, iters=3000):
    shifts = seed_shifts or [random.random() for _ in freqs]
    def unc(shs):
        grid = 3000; u = 0
        for j in range(grid):
            th = (j + 0.5) / grid
            if all(dist_to_int(f * th + s) >= RHO for f, s in zip(freqs, shs)):
                u += 1
        return u / grid
    cur = unc(shifts)
    for _ in range(iters):
        k = random.randrange(7)
        old = shifts[k]
        shifts[k] = random.random()
        new = unc(shifts)
        if new <= cur: cur = new
        else: shifts[k] = old
    return cur

eng = min(try_engineer() for _ in range(3))
print(f"(4b) engineered cover attempt (freqs 1..7, tuned shifts): min uncovered = {eng:.5f}")
print("     (>0 means even l=7 struggles at these frequencies; the rung only needs l>=7 NECESSARY)")
print("ALL CHECKS PASSED")
