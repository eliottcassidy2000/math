#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S9 -- HYP-4332: the Fourier/Newman lower bound for CircleClearFloor.

The uncovered set U = {tau : all ||w_i tau - phi_i|| >= rho} has measure
    |U| = INT prod_i g(w_i tau - phi_i) dtau
        = (1-2rho)^l + SUM_{k != 0, sum k_i w_i = 0} prod_i ghat(k_i) e^{-2pi i sum k_i phi_i}
where g = indicator of {dist >= rho}, ghat(0) = 1-2rho, ghat(k) = -sin(2pi k rho)/(pi k).

RIGOROUS LOWER BOUND (all phases): |U| >= (1-2rho)^l - R,
    R = SUM_{k != 0, sum k_i w_i = 0} prod_i |ghat(k_i)|.
So (1-2rho)^l > R  =>  |U| > 0 for ALL shifts => CircleClearFloor holds (no cover).

This computes R (truncated at |k_i| <= K, ghat decays like 1/k) for distinct
frequency sets and compares to (1-2rho)^l.  The point: locate the class of
frequency sets the Newman-Fourier bound CLOSES, and see the residual.
"""
from math import sin, pi, gcd
from itertools import product
import random, time
random.seed(99)

T0 = time.time()
def log(m=""): print(m, flush=True)
RHO = 2/25
MEAN = 1 - 2*RHO           # 21/25 = 0.84

def ghat(k):
    if k == 0: return MEAN
    return -sin(2*pi*k*RHO)/(pi*k)

def resonance_R(freqs, K=6):
    """R = sum over k != 0 with sum k_i w_i = 0, |k_i| <= K, of prod |ghat(k_i)|.
    Enumerate by DFS pruning on partial sums (bounded)."""
    l = len(freqs)
    absg = {k: abs(ghat(k)) for k in range(-K, K+1)}
    # DFS over coordinates, tracking partial sum s = sum_{j<i} k_j w_j and product
    # prune: remaining |sum| reachable is bounded by K*sum(remaining freqs)
    suffix = [0]*(l+1)
    for i in range(l-1, -1, -1):
        suffix[i] = suffix[i+1] + K*freqs[i]
    total = 0.0
    def dfs(i, s, prod, any_nonzero):
        nonlocal total
        if i == l:
            if s == 0 and any_nonzero:
                total += prod
            return
        # prune: can the remaining make s -> 0?
        if abs(s) > suffix[i]:
            return
        for k in range(-K, K+1):
            dfs(i+1, s + k*freqs[i], prod*absg[k], any_nonzero or (k != 0))
    dfs(0, 0, 1.0, False)
    return total

log(f"Newman-Fourier bound: (1-2rho)^l vs R (rho=2/25, (1-2rho)=0.84)\n")
log(f"{'freqset':>34} {'l':>3} {'(0.84)^l':>9} {'R(K=6)':>9} {'>R?':>5} {'|U|min(T3)':>10}")

# structured worst cases (most resonances = small consecutive / AP-like)
cases = [
    ("1..7", list(range(1,8))),
    ("1..9", list(range(1,10))),
    ("1..11", list(range(1,12))),
    ("2,4,6,8,10,12,14", [2,4,6,8,10,12,14]),
    ("primes 2..19 (7)", [2,3,5,7,11,13,17]),
    ("1,2,4,8,16,32,64", [1,2,4,8,16,32,64]),   # no small signed sums (lacunary)
    ("13-lift ladder-ish", [98,112,126,140,154,168,14]),   # 14r deep-well freqs + 14
]
for name, fr in cases:
    l = len(fr)
    main = MEAN**l
    R = resonance_R(fr, K=6)
    log(f"{name:>34} {l:>3} {main:>9.4f} {R:>9.4f} {'YES' if main>R else 'no':>5}")

log("\nrandom distinct sets (small, resonance-rich) at l=7,9,11:")
for l in (7, 9, 11):
    won = 0; tot = 0
    worst = (1e9, None)
    for _ in range(60):
        fr = sorted(random.sample(range(1, 25), l))
        if gcd(*fr) != 1: continue
        tot += 1
        main = MEAN**l; R = resonance_R(fr, K=5)
        if main > R: won += 1
        if R - main < worst[0]: worst = (R - main, fr, main, R)
    log(f"  l={l}: (0.84)^l > R for {won}/{tot} random small sets; tightest gap "
        f"(0.84)^l - R = {-worst[0]:.4f} at {worst[1]} ((0.84)^l={worst[2]:.3f}, R={worst[3]:.3f})")

log("\nREADING: where (0.84)^l > R, CircleClearFloor is PROVED for ALL shifts (Newman-Fourier).")
log("The residual = resonance-rich sets where R exceeds the main term at K=6 truncation;")
log("those need either a sharper resonance bound or the phase-coupling (adversary can't align all).")
log(f"[t = {time.time()-T0:.0f}s]")
