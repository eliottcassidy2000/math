#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S9 -- HYP-4332: locate the Newman content in CircleClearFloor.

CircleClearFloor (the crux's last covering obligation): l distinct positive
integer frequencies w_1<...<w_l, combs {tau : ||w_i tau - phi_i|| < rho},
rho = 2/25.  Do they leave an uncovered gap (a clear point at level >= rho)?

CLAIM (this session): SPLIT by shifts.
 (A) UNSHIFTED (phi_i = 0) = LRC(l): there is tau with min_i ||w_i tau|| >=
     1/(l+1).  For l <= 11: 1/(l+1) >= 1/12 > 2/25 => clearance >= 1/12-2/25
     = 1/300.  So the unshifted floor is a ONE-LINE CITATION, not Newman.
 (B) SHIFTED (phi_i free): sLRC is FALSE (from n=5), so shifts may let distinct
     frequencies cover.  This is the genuine Newman/DMNR content -- IF the
     shifts are free.  In the lift application the shifts are CONSTRAINED
     (phi_i = a_i t_0 from one base coupling), not free.

TESTS:
 T1  confirm 1/(l+1) > 2/25 for l <= 11 (the unshifted citation bound).
 T2  UNSHIFTED covering floor: for random distinct-frequency sets, the best
     clear point (max_tau min_i ||w_i tau||) -- confirm >= 1/(l+1) (LRC), and
     that combs at 2/25 never cover for l <= 11.
 T3  SHIFTED-ADVERSARIAL: can FREE shifts phi_i make distinct frequencies cover
     the circle at 2/25 for l <= 11?  (hill-climb the shifts to minimize the
     uncovered fraction.)  If YES -> the lemma needs constrained shifts; the
     open content is real.  If NO even with free shifts -> the whole
     CircleClearFloor is citation-adjacent.
 T4  CONSTRAINED shifts (phi_i = a_i * t_0, one parameter t_0 + integer a_i):
     the actual lift structure -- does it ever cover?  (this is what the
     application needs; my S6b says the 2-torus is subsumed anyway.)
"""
from fractions import Fraction as F
from math import gcd
import random, time
random.seed(9)

T0 = time.time()
def log(m=""): print(m, flush=True)
RHO = 2/25

def dist(x):
    x = x - round(x)
    return abs(x)

def clear_level(freqs, phases, N=4000):
    """max_tau min_i ||w_i tau - phi_i|| on a grid (a lower bound on the true max)."""
    best = 0.0
    for j in range(N):
        tau = j / N
        m = 1.0
        for w, ph in zip(freqs, phases):
            d = dist(w * tau - ph)
            if d < m:
                m = d
                if m <= best: break
        if m > best: best = m
    return best

def uncovered_frac(freqs, phases, N=6000):
    """fraction of the circle uncovered by radius-RHO combs."""
    unc = 0
    for j in range(N):
        tau = j / N
        if all(dist(w * tau - ph) >= RHO for w, ph in zip(freqs, phases)):
            unc += 1
    return unc / N

log("T1: 1/(l+1) vs 2/25 = 0.08 (unshifted LRC citation bound):")
for l in range(7, 14):
    log(f"   l={l}: 1/(l+1) = {1/(l+1):.5f}  {'> 2/25 (clears, CITATION)' if 1/(l+1) > RHO else '<= 2/25 (LRC insufficient)'}")

log("\nT2: UNSHIFTED clear level (best tau, phi=0) for random distinct freq sets:")
for l in (7, 9, 11, 12):
    mn = 1.0
    for _ in range(200):
        freqs = sorted(random.sample(range(1, 60), l))
        if gcd(*freqs) != 1: continue
        cl = clear_level(freqs, [0]*l)
        mn = min(mn, cl)
    log(f"   l={l}: min clear level over samples = {mn:.5f}  (LRC bound 1/(l+1) = {1/(l+1):.5f}; "
        f"combs cover at 2/25? {'NO' if mn > RHO else 'maybe'})")

log("\nT3: SHIFTED-ADVERSARIAL -- can FREE shifts make distinct freqs COVER at 2/25?")
log("    (hill-climb phases to MINIMIZE uncovered fraction; 0 => covered)")
for l in (7, 9, 11, 12, 13):
    best_cover = 1.0
    for restart in range(40):
        freqs = sorted(random.sample(range(1, 40), l))
        if gcd(*freqs) != 1: continue
        ph = [random.random() for _ in range(l)]
        cur = uncovered_frac(freqs, ph, 2000)
        for step in range(60):
            i = random.randrange(l)
            old = ph[i]; ph[i] = random.random()
            nu = uncovered_frac(freqs, ph, 2000)
            if nu <= cur: cur = nu
            else: ph[i] = old
        best_cover = min(best_cover, cur)
        if cur == 0: break
    log(f"   l={l}: min uncovered fraction over freq-sets & shifts = {best_cover:.4f}  "
        f"{'<< COVERS with free shifts!' if best_cover < 0.005 else 'floor remains (Newman)'}")

log(f"\nREADING:")
log(f" - l <= 11 unshifted: CITATION (1/(l+1) >= 1/12 > 2/25), clearance >= 1/300. Not Newman.")
log(f" - the OPEN content = shifted distinct-frequency covering (T3); locate whether free")
log(f"   shifts cover (=> the lemma needs the constrained/lift shift structure).")
log(f"[t = {time.time()-T0:.0f}s]")
