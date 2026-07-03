# -*- coding: utf-8 -*-
"""
kps-2026-07-03-S25. Ground-truth probe for the regime-C drifting floor.

Question: for a 7-far block {w1..w7} cited against k=6 bounded runners (B=22),
window L = 2*delta, delta = (13-k)/(14*(k+1)*B), does the Hunter one-pair ledger
    safe >= L - singles_sum + firstpair_credit
have a chance of being positive over ADVERSARIAL window positions?

We measure, over many window offsets tau:
  - singles_sum(tau) = sum_w |window ∩ danger(w)|   (danger(w) = {t: dist(w t,Z)<1/14})
  - for EACH consecutive/any pair, pairmass(tau) = |window ∩ danger(wi) ∩ danger(wj)|
  - ledger_onepair(tau) = L - singles_sum + max_pair pairmass   (best single pair)
  - ledger_allpairs(tau) = L - singles_sum + sum_{path} pairmass (full Hunter, path = sorted)
Report the MIN over tau (adversarial). Positive min => closable at that config.
Sweep w1 and span to map the closable region.
"""
import numpy as np

h = 1.0/14.0           # danger radius
k = 6; B = 22
delta = (13-k)/(14*(k+1)*B)
L = 2*delta
print(f"window L = {L:.6f}  (delta={delta:.6f})")

def dist_to_Z(x):
    x = x - np.floor(x)
    return np.minimum(x, 1.0-x)

def danger_mask(w, ts):
    return dist_to_Z(w*ts) < h

def measure_config(ws, n_tau=400, n_fine=20000):
    """min over tau of the one-pair and all-pair Hunter ledgers."""
    ws = sorted(ws)
    # fine grid inside a window, reused per tau
    s = np.linspace(0, L, n_fine, endpoint=False)
    dt = L/n_fine
    min_one = 1e9; min_all = 1e9
    # adversarial tau: sweep window start over a full period-ish range
    taus = np.linspace(0, 1.0, n_tau, endpoint=False)
    for tau in taus:
        ts = tau + s
        masks = [danger_mask(w, ts) for w in ws]
        singles = sum(m.sum() for m in masks) * dt
        # pair masses (all pairs)
        pair_best = 0.0
        path_sum = 0.0
        for i in range(len(ws)-1):
            pm = (masks[i] & masks[i+1]).sum()*dt
            path_sum += pm
        # best single pair over ALL pairs (we can reorder, so max gap pair too)
        for i in range(len(ws)):
            for j in range(i+1, len(ws)):
                pm = (masks[i] & masks[j]).sum()*dt
                if pm > pair_best: pair_best = pm
        one = L - singles + pair_best
        allp = L - singles + path_sum
        min_one = min(min_one, one); min_all = min(min_all, allp)
    return min_one, min_all

configs = {
    "consec 23..29 (tight,small)":      [23,24,25,26,27,28,29],
    "consec 100..106":                  list(range(100,107)),
    "span60 @100 (100,110,..,160)":     [100,110,120,130,140,150,160],
    "span300 @100":                     [100,150,200,250,300,350,400],
    "consec 1000..1006":                list(range(1000,1007)),
    "span300 @1000":                    [1000,1050,1100,1150,1200,1250,1300],
    "consec 4000..4006":                list(range(4000,4007)),
    "consec 7400..7406 (regime A edge)":list(range(7400,7407)),
    "span300 @7400":                    [7400,7450,7500,7550,7600,7650,7700],
    "span1500 @7400 (spread,large)":    [7400,7650,7900,8150,8400,8650,8900],
    "span1500 @100 (spread,small)":     [100,350,600,850,1100,1350,1600],
}
print(f"\n{'config':40s} {'min_1pair':>12s} {'min_allpair':>12s}  {'w1':>6s}")
for name, ws in configs.items():
    m1, ma = measure_config(ws, n_tau=300, n_fine=8000)
    flag = "  CLOSABLE" if ma > 0 else ""
    flag1 = " [1p+]" if m1 > 0 else ""
    print(f"{name:40s} {m1:12.6f} {ma:12.6f}  {ws[0]:6d}{flag}{flag1}")

# fee vs floor scaling summary
print(f"\nL/49 (floor density scale) = {L/49:.6f}")
for w in [23,100,1000,3769,7400,22556]:
    fee = 3.0/w   # ~ 7 runners * 3/(7w)
    print(f"  w1={w:6d}: fee~3/w1={fee:.6f}   6L/49(klein budget)={6*L/49:.6f}   fee<budget? {fee < 6*L/49}")
