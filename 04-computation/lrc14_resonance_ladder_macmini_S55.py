"""
mac-mini-2026-07-07-S55 (HYP-5177) -- THE RESONANCE-LADDER LAW (the coprime lens).

CLAIMS (verify each):
 L1 (config-level dilation invariance, 2-line proof): the substitution y = cx is
    measure-preserving and maps the config of cE at x to the config of E at y; hence
    EVERY config functional -- including per-collapse-level window mass -- is dilation
    invariant.  Verification: the S54 window formula evaluated for cE at the x-rational
    p/(den) with den = c*q/gcd-structure reproduces E's window values scaled correctly.
 L2 (the ladder): a config-level-q collapse of cE lives at x-rationals with denominator
    Q = c*q/gcd(c, q*)-structure; concretely for x = p/Q primitive, the config of cE
    collapses mod q = Q/gcd(c, Q) * (gcd part of numerator...) -- operationally:
    cx mod 1 = (c p mod Q)/Q reduced; config-denominator q = Q/gcd(cp, Q).
 L3 (accounting incompleteness): the q <= 6 x-denominator window sum UNDERCOUNTS mu
    for c-adic families: the missing mass sits at x-denominators in (6, 6c] -- locate
    it for the 3-adic cascade E[U]-minimizer and the parity record (2-adic).
 L4 (the retro-explanation): the S41 Farey-shell 'generic' mass of the 3-adic cascade
    = ladder mass at denominators 9, 12, 15, 18 -- quantify the recovery.
"""
import numpy as np
from math import gcd
import random as rnd
rnd.seed(55)

TH = 1/7

def Wq_formula_at(E, p, Q):
    """leading-order window width of family E at the x-rational p/Q (any Q)."""
    # config at x = p/Q + t: frac(e (p/Q)) = (e p mod Q)/Q; class labels c_i = e*p mod Q
    classes = {}
    for e in E:
        c = (e * p) % Q
        classes.setdefault(c, []).append(e)
    occ = sorted(classes)
    r = len(occ)
    tot = 0.0
    for sign in (+1, -1):
        best = 0.0
        for s in range(r):
            c1 = occ[s]; c2 = occ[(s+1) % r]
            G = (c2 - c1) % Q
            if r == 1: G = Q
            gapval = G/Q - TH
            if gapval <= 0: continue
            left, right = classes[c1], classes[c2 if r > 1 else c1]
            D = (max(left) - min(right)) if sign > 0 else (max(right) - min(left))
            if D <= 0: continue
            best = max(best, gapval / D)
        tot += best
    return tot

def total_window_sum(E, Qmax):
    """sum of formula windows over ALL primitive x-rationals p/Q, Q <= Qmax."""
    tot = 0.0; per_Q = {}
    for Q in range(1, Qmax+1):
        s = 0.0
        for p in range(0, Q):
            if Q > 1 and gcd(p, Q) != 1: continue
            if Q == 1 and p != 0: continue
            s += Wq_formula_at(E, p, Q)
        per_Q[Q] = s; tot += s
    return tot, per_Q

def mu_measured(E, grid=400000):
    xs = (np.arange(grid)+0.5)/grid
    ph = np.mod(np.outer(xs, np.array(E, float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:,0]+1-ph[:,-1])[:,None]], axis=1)
    return float((g.max(axis=1) > TH).mean())

print("=== L1/L2: dilation ladder -- windows of cE vs windows of E ===")
AP = list(range(1, 14))
for c in (2, 3):
    cE = [c*e for e in AP]
    tot_E, perQ_E = total_window_sum(AP, 6)
    tot_cE, perQ_cE = total_window_sum(cE, 6*c)
    print(f"  c={c}: sum windows E (Q<=6) = {tot_E:.4f}; sum windows {c}E (Q<={6*c}) = {tot_cE:.4f}"
          f"  (mu = {mu_measured(AP):.4f} both -- invariance {'OK' if abs(tot_E-tot_cE)<1e-6 else 'CHECK'})")
    print(f"    {c}E per-Q mass: " + " ".join(f"Q{Q}:{perQ_cE[Q]:.4f}" for Q in sorted(perQ_cE) if perQ_cE[Q] > 1e-9))
    # predicted ladder: window of E at p/q reappears for cE at denominators c*q (coprime part)
    pred = {}
    for q in range(1, 7):
        if perQ_E[q] > 1e-9:
            pred[c*q] = perQ_E[q]  # mass moves q -> c*q wholesale when gcd(c,q)=1... report
    print(f"    ladder prediction (mass q -> {c}q): " + " ".join(f"Q{Q}:{v:.4f}" for Q,v in sorted(pred.items())))

print("\n=== L3/L4: where the cascades hide their window mass ===")
FAMS = {
 '3-adic EU-min': [0,30,36,45,50,54,60,63,70,72,81,90,108],
 'parity record': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'AP13': AP,
}
for name, E in FAMS.items():
    mu = mu_measured(E)
    t6, p6 = total_window_sum(E, 6)
    t18, p18 = total_window_sum(E, 18)
    hidden = {Q: v for Q, v in p18.items() if Q > 6 and v > 1e-4}
    print(f"  {name:16s}: mu = {mu:.4f}; window sum Q<=6 = {t6:.4f} ({t6/mu*100:.0f}% of mu); "
          f"Q<=18 = {t18:.4f} ({t18/mu*100:.0f}%)")
    if hidden:
        print(f"      hidden ladder mass: " + " ".join(f"Q{Q}:{v:.4f}" for Q, v in sorted(hidden.items())))
