"""
mac-mini-2026-07-07-S54 (HYP-5147) -- the EXACT per-q WINDOW FORMULA (kps-S72 handoff
(i): the per-q lemma needs the exact window width from residue data + intra-class drift).

SETUP. Near x = p/q + t (q <= 6, gcd(p,q)=1, small t): frac(e_i x) = frac(c_i/q + e_i t),
c_i = e_i p mod q. Points with the same c cluster; cluster c occupies (to first order)
[c/q + t min_{class c} e, c/q + t max_{class c} e] (t > 0; mirrored for t < 0).
Let the OCCUPIED residues be c_(1) < ... < c_(r) with cyclic residue gaps
G_s = c_(s+1) - c_(s) (integers summing to q).  The gap between cluster s and s+1:
    gap_s(t) = G_s/q - t * (maxE(c_(s)) - minE(c_(s+1)))   for t > 0
             (edge drift rate D_s^+ = maxE of the left cluster minus minE of the right;
              for t < 0 the roles mirror: D_s^- = maxE(right) - minE(left) ... care).
maxgap(x) = max_s gap_s(t) to first order.  The q-window (one p):
    w(p,q) = |{t : max_s (G_s/q - t D_s) > 1/7}| summed over both signs
           = sum over sign of max_s (G_s/q - 1/7)_+ / D_s^{sign}   ... with the max
    taken correctly: the window endpoint is where the LAST gap above 1/7 drops below:
    t*_{sign} = max_s (G_s/q - 1/7)_+ / D_s^{sign}, so w(p,q) = t*_+ + t*_-.
FORMULA:  W_q(E) ~ sum_{p in (Z/q)^*} [max_s (G_s/q - 1/7)_+/D_s^+ + max_s (...)/D_s^-].

THIS SCRIPT: (1) implement the formula; (2) compare against DIRECT measurement
(attribute good-x mass to its governing q via nearest small-denominator rational --
kps's attribution) at k = 13 and k = 8 for the AP and a bank; targets: kps-S72's
measured W_q(AP_13) = 0.065/0.086/0.054/0.078/0.016 (q = 2..6); (3) the q=2 case
structure: W_2 formula explicitly + the extremal comparison over the bank (is the AP
minimal? what carries the inequality: residue-miss (G >= 2) vs drift (D)?).
"""
import numpy as np
from math import gcd
from itertools import combinations
import random as rnd
rnd.seed(54)

TH = 1/7

def Wq_formula(E, q):
    """leading-order window mass at resonance q."""
    tot = 0.0
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        classes = {}
        for e in E:
            c = (e * p) % q
            classes.setdefault(c, []).append(e)
        occ = sorted(classes)
        r = len(occ)
        # cyclic gaps and edge drifts, both time directions
        for sign in (+1, -1):
            best = 0.0
            for s in range(r):
                c1, c2 = occ[s], occ[(s+1) % r]
                G = (c2 - c1) % q
                if G == 0: G = q  # r == 1: single cluster, gap q... handle r=1: gap = q
                if r == 1: G = q
                gapval = G/q - TH
                if gapval <= 0: continue
                left, right = classes[c1], classes[c2 if r > 1 else c1]
                if sign > 0:
                    # t>0: left cluster's right edge moves at max(left); right cluster's
                    # left edge at min(right): closing rate = max(left) - min(right)
                    D = max(left) - min(right)
                else:
                    D = max(right) - min(left)
                if D <= 0:
                    # gap OPENS in this direction: window extends to the cell boundary
                    # -- cap at the Farey cell scale (conservative: skip; flag)
                    continue
                best = max(best, gapval / D)
            tot += best
    return tot

def Wq_measured(E, qmax=6, grid=400000):
    """attribute each good x to the smallest q <= 6 with a rational p/q within its
    'collapse range' -- practical attribution: nearest rational with q <= 6."""
    xs = (np.arange(grid)+0.5)/grid
    ph = np.mod(np.outer(xs, np.array(E, float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:,0]+1-ph[:,-1])[:,None]], axis=1)
    good = g.max(axis=1) > TH
    # nearest small-denominator rational attribution
    best_q = np.full(grid, 99, dtype=int)
    best_d = np.full(grid, 9.9)
    for q in range(1, qmax+1):
        for p in range(0, q+1):
            if q > 1 and gcd(p, q) != 1: continue
            d = np.abs(xs - p/q)
            upd = d < best_d
            best_d[upd] = d[upd]; best_q[upd] = q
    W = {}
    for q in range(1, qmax+1):
        W[q] = float((good & (best_q == q)).mean())
    return W, float(good.mean())

print("=== AP_13: formula vs measured vs kps-S72 ===")
AP13 = list(range(1, 14))
Wm, mu = Wq_measured(AP13)
print(f"  measured (nearest-rational attribution): q=1..6 -> "
      f"{[round(Wm[q],4) for q in range(1,7)]}; sum = {sum(Wm.values()):.4f} = mu = {mu:.4f}")
print(f"  kps-S72 reported W_q(AP_13) q=2..6: 0.065/0.086/0.054/0.078/0.016")
for q in range(2, 7):
    f = Wq_formula(AP13, q)
    print(f"  q={q}: formula = {f:.4f}   measured = {Wm[q]:.4f}")

print("\n=== bank comparison at q=2 (the extremal question) ===")
BANK = {
 'AP13': AP13,
 'GW': list(range(1,12))+[13,24],
 'parity record': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'miss-2 (all odd)': [1,3,5,7,9,11,13,15,17,19,21,23,25],
 'random': sorted(rnd.sample(range(1,200),13)),
}
for name,E in BANK.items():
    f2 = Wq_formula(E, 2)
    Wm2, mu2 = Wq_measured(E)
    print(f"  {name:16s}: W_2 formula = {f2:.4f}, measured = {Wm2[2]:.4f}  (mu = {mu2:.4f})")
print("\nq=2 structure note: residues mod 2 split E into odd/even classes at x ~ p/2;")
print("miss (all one parity) => G = 2, gap 1 - 1/7 huge => wide window; full => G = 1,")
print("gap 1/2 - 1/7, drift = cross-parity edge rates. The extremal comparison is")
print("(3/14)/D_AP vs family values -- printed above for the dichotomy check.")
