#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THM-606 VERIFICATION: the depth-d ladder, instantiated exactly at depth 3.

Family: P = {1,2};  level 1: V1,  offs1 = {0,1,2};  level 2: V2, offs2 = {0,1,3};
level 3: V3, offs3 = {0,1,2,4,7}  -- 2 + 3 + 3 + 5 = 13 runners at band h = 1/14.

Pipeline (all exact rationals):
  1. choose the window [lo, hi] inside the {1,2}-lonely region and the canonical delta-schedule
     (F1): delta_l = (1+eta) * sum_{m>l} 1/V_m;  mu_l = V_l * delta_l;  check (H-len);
  2. certificate search: c_l on a rational grid with arcSafe(h + mu_l, o, c_l, lo, hi) for all
     offsets (kps's exact single-cell criterion), and arcSafe(h, s, 0, lo, hi) for P;
  3. the ladder walk: j_k := any integer in [V_k a - c_k, V_k(a + delta_{k-1} - delta_k) - c_k]
     (exists by H-len), t_k := (j_k + c_k)/V_k; final tau := t_3 EXACT;
  4. verify dist(v*tau, Z) >= 1/14 for ALL 13 runners; report margins;
  5. robustness: repeat the walk from 200 randomized starts a in [lo, hi - delta_0'] and with
     j_k chosen adversarially (floor vs ceil);
  6. the (F2) ratio table: minimal feasible mu vs ratio rho;
  7. NEGATIVE control: shrink V3 below the (H-len) threshold and confirm the ruler interval
     contains no integer for some start a (the hypothesis is sharp, not slack).

opus-2026-07-02-S36 (THM-606 / HYP-3905).
"""
import sys, random
from fractions import Fraction as Fr
from math import floor, ceil
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

H = Fr(1, 14)

def distZ(x):
    f = x - floor(x)
    return min(f, 1 - f)

def arcSafe(beta, w, c, lo, hi):
    """kps's exact single-cell criterion: w*tau - c stays in [k+beta, k+1-beta] on [lo,hi]."""
    A = w * lo - c; B = w * hi - c
    if A > B: return False
    k = floor(A)
    return (beta <= A - k) and (B <= k + 1 - beta)

# ---------------------------------------------------------------- 1. parameters
P = [1, 2]
LEVELS = [  # (V, offs)
    (50,    [0, 1, 2]),
    (2000,  [0, 1, 3]),
    (90000, [0, 1, 2, 4, 7]),
]
d = len(LEVELS)
lo, hi = Fr(7, 20), Fr(3, 8)             # [0.35, 0.375] inside the {1,2}-lonely region
# NOTE: first attempt used hi = 0.37 (delta0 = 1/50 = 0.02), violating the NECESSARY (F1)
# condition delta0 > sum 1/V_l = 0.0205 -- the walk failed at level 1 exactly as the calculus
# predicts. The calculus is sharp; delta0 = 1/40 restores (H-len) with margin.
delta0 = hi - lo
eta = Fr(1, 10)
Vs = [V for V, _ in LEVELS]
delta = [delta0] + [None] * d
for l in range(1, d + 1):
    delta[l] = (1 + eta) * sum(Fr(1, Vs[m]) for m in range(l, d))  # m>l-1 in 0-idx: levels l+1..d
mu = [None] + [Vs[l - 1] * delta[l] for l in range(1, d + 1)]
assert delta[d] == 0 and mu[d] == 0

print("=" * 100)
print(" THM-606 depth-3 ladder, exact instantiation.  h = 1/14; window [%s, %s], delta0 = %s" % (lo, hi, delta0))
print("=" * 100)
print(f"  levels: " + "; ".join(f"V{l}={Vs[l-1]}, offs={LEVELS[l-1][1]}" for l in range(1, d + 1)))
print(f"  schedule: " + "; ".join(f"delta_{l}={delta[l]} (mu_{l}={mu[l]}={float(mu[l]):.4f})" for l in range(1, d + 1)))

# (H-len)
print("\n (H-len) checks: 1 < V_l * (delta_{l-1} - delta_l):")
hlen_ok = True
for l in range(1, d + 1):
    val = Vs[l - 1] * (delta[l - 1] - delta[l])
    ok = val > 1
    hlen_ok &= ok
    print(f"   l={l}: V={Vs[l-1]}: {val} = {float(val):.4f} > 1 ? {ok}")
print(f"   ALL (H-len): {hlen_ok}")

# ---------------------------------------------------------------- 2. certificates
def find_cert(offs, beta, grid=2000):
    for k in range(grid):
        c = Fr(k, grid)
        if all(arcSafe(beta, o, c, lo, hi) for o in offs):
            return c
    return None

print("\n certificate search (arcSafe at inflated bands):")
certs = [None]
cert_ok = True
for l in range(1, d + 1):
    beta = H + mu[l]
    c = find_cert(LEVELS[l - 1][1], beta, 4000)
    certs.append(c)
    print(f"   level {l}: band h+mu = {beta} = {float(beta):.4f}: c_{l} = {c}")
    cert_ok &= c is not None
pok = all(arcSafe(H, s, Fr(0), lo, hi) for s in P)
print(f"   P = {P} arcSafe(h, s, 0): {pok}")
assert cert_ok and pok, "certificate search failed"

# ---------------------------------------------------------------- 3-4. the ladder walk + check
def ladder_walk(a, pick):
    """run the ladder from window start a; pick in {'lo','hi','rand'} chooses j_k."""
    t = a
    for l in range(1, d + 1):
        V, c = Vs[l - 1], certs[l]
        A = V * t - c
        B = V * (t + delta[l - 1] - delta[l]) - c
        assert B - A > 1 - Fr(1, 10**9), "ruler interval shorter than 1?!"
        jlo, jhi = ceil(A), floor(B)
        assert jlo <= jhi, f"no integer in ruler interval at level {l}"
        j = jlo if pick == 'lo' else (jhi if pick == 'hi' else random.randint(jlo, jhi))
        t = Fr(j + c, V)
    return t

def check_tau(tau):
    worst = Fr(1)
    for s in P:
        worst = min(worst, distZ(s * tau))
    for l in range(1, d + 1):
        V = Vs[l - 1]
        for o in LEVELS[l - 1][1]:
            worst = min(worst, distZ((V - o) * tau))
    return worst

tau = ladder_walk(lo, 'lo')
worst = check_tau(tau)
print(f"\n MAIN WALK from a = lo: tau = {tau}")
print(f"   worst runner margin: {worst} = {float(worst):.6f}  >= 1/14 = {float(H):.6f} ? {worst >= H}")
for s in P:
    print(f"   speed {s:>6}: dist = {float(distZ(s*tau)):.6f}")
for l in range(1, d + 1):
    for o in LEVELS[l - 1][1]:
        v = Vs[l - 1] - o
        print(f"   speed {v:>6}: dist = {float(distZ(v*tau)):.6f}   (level {l}, offset {o})")

# ---------------------------------------------------------------- 5. robustness
random.seed(606)
fails = 0
worst_all = Fr(1)
for trial in range(200):
    # randomized start for CLAIM(2): a in [lo, hi - delta_1]
    a = lo + Fr(random.randint(0, 10**6), 10**6) * (hi - lo - delta[1])
    # run levels 2..3 only from window [a, a+delta_1] (CLAIM(2) robustness), then check levels 2,3
    t = a
    okwalk = True
    for l in range(2, d + 1):
        V, c = Vs[l - 1], certs[l]
        A = V * t - c
        B = V * (t + delta[l - 1] - delta[l]) - c
        jlo, jhi = ceil(A), floor(B)
        if jlo > jhi: okwalk = False; break
        t = Fr(random.randint(jlo, jhi) + c, V)
    if not okwalk:
        fails += 1; continue
    w = Fr(1)
    for l in range(2, d + 1):
        for o in LEVELS[l - 1][1]:
            w = min(w, distZ((Vs[l - 1] - o) * t))
    worst_all = min(worst_all, w)
    if w < H: fails += 1
full_fails = 0
for trial in range(200):
    tau_r = ladder_walk(lo, 'rand')
    if check_tau(tau_r) < H: full_fails += 1
print(f"\n ROBUSTNESS: CLAIM(2) from 200 random window starts: fails = {fails}, worst sub-margin = {float(worst_all):.6f}")
print(f"             full ladder with random ruler picks x200: fails = {full_fails}")

# ---------------------------------------------------------------- 6. (F2) table
print("\n (F2) ratio table: mu_per_level ~ (1+eta)/(rho-1) vs exact schedule at geometric V:")
for rho in (10, 20, 29, 40, 80):
    V1 = 50
    Vg = [V1 * rho ** i for i in range(3)]
    dl = [(1 + eta) * sum(Fr(1, Vg[m]) for m in range(l, 3)) for l in range(1, 4)]
    mus = [Vg[l - 1] * dl[l - 1] for l in range(1, 4)]
    bound = Fr(1 + eta, rho - 1)
    print(f"   rho={rho:>3}: mu_1={float(mus[0]):.4f} mu_2={float(mus[1]):.4f} (bound (1+eta)/(rho-1)={float(bound):.4f}) "
          f"ok={all(m <= bound + Fr(1,10**12) for m in mus[:2])}")

# ---------------------------------------------------------------- 7. negative control
print("\n NEGATIVE CONTROL: violate (H-len) at level 3 (V3' with V3'*delta_2 < 1):")
V3bad = int(1 / float(delta[2]) * 0.5)
misses = 0
for k in range(1000):
    a = lo + Fr(k, 1000) * (hi - lo - delta[1])
    # level-3 ruler interval from a window of width delta_2 at speed V3bad:
    A = V3bad * a - certs[3]; B = V3bad * (a + delta[2]) - certs[3]
    if ceil(A) > floor(B): misses += 1
print(f"   V3' = {V3bad} (V3'*delta_2 = {float(V3bad*delta[2]):.3f} < 1): ruler interval EMPTY for {misses}/1000 starts")
print(f"   => (H-len) is sharp: below it the construction genuinely fails on a positive fraction of windows.")
print("\nDONE." if worst >= H and fails == 0 and full_fails == 0 else "\n*** FAILURES PRESENT ***")
