#!/usr/bin/env python3
"""
kps-2026-07-07-S57 -- a cleaner reduction of Route 1's density floor: the mean MAX-GAP.

opus-S131's Paley-Zygmund reduction: mu_1/7(E) = P_x(maxgap{frac(e_i x)} > 1/7) >= E[U],
E[U] = E_x[ sum_gaps (gap-1/7)_+ ], open lemma inf_E E[U] > 0 (obstructed by triple+ overlaps).

CREATIVE SIMPLIFICATION (kps-S57).  Two elementary reductions to a SINGLE order statistic:

 (1) POINTWISE:  U(x) = sum_gaps (gap-1/7)_+ >= (maxgap-1/7)_+ >= maxgap - 1/7.
     So  E[U] >= E[maxgap] - 1/7.   (connects to opus's PZ)

 (2) REVERSE-MARKOV (bypasses E[U] entirely):  maxgap in [0,1], so for a=1/7,
        mu_1/7 = P(maxgap > 1/7) >= (E[maxgap] - 1/7)/(1 - 1/7) = (7/6)(E[maxgap] - 1/7).

Either way the density floor  mu_1/7(E) > 0  reduces to the MEAN MAX-GAP bound

        inf_E  E[maxgap]  >  1/7,

which is EMPIRICALLY COMFORTABLE: inf E[maxgap] ~ 0.20 (adversarial descent), margin ~0.06
above 1/7 = 0.1429 -- NOT razor-thin (unlike the tight raw-loneliness M=1/14).  This is one
order statistic with a comfortable margin, vs E[U]'s full inclusion-exclusion.

Rigorous partials toward inf E[maxgap] > 1/7 (both fall a little short -- the remaining content
is 'max beats the length-biased/origin gap', a three-distance/discrepancy statement):
  * E[maxgap] >= E[sum g_i^2] >= 1/k = 1/13   (Cauchy-Schwarz; length-biased mean).
  * E[maxgap] >= E[gap_0] = E[min frac] + 1 - E[max frac] ~ 2/(k+1) = 1/7   (origin gap).
"""
import numpy as np, random

def stats(E, res=8000):
    xs = (np.arange(res) + 0.5) / res
    ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0); ph.sort(axis=1)
    gaps = np.empty_like(ph); gaps[:, :-1] = np.diff(ph, axis=1); gaps[:, -1] = (ph[:, 0] + 1.0) - ph[:, -1]
    mg = gaps.max(axis=1)
    U = np.sum(np.maximum(gaps - 1.0/7, 0.0), axis=1)
    return (mg > 1.0/7).mean(), mg.mean(), U.mean(), (gaps**2).sum(axis=1).mean(), (ph[:,0] + 1.0 - ph[:,-1]).mean()

print("Route-1 density floor via the mean max-gap (1/7 = %.5f)" % (1.0/7))
print(" family                  mu_1/7  E[maxgp] E[U]   (7/6)(Emg-1/7) E[sg^2] E[gap0]")
for nm, E in [("AP {1..13}", list(range(1,14))),
              ("2*AP (M=1/14 tight)", [2*i for i in range(1,14)]),
              ("prim-sat 2*{1..12}+13", [2,4,6,8,10,12,13,14,16,18,20,22,24]),
              ("min-E[maxgap] adversarial", [2,6,8,10,11,12,14,16,18,20,22,26,42]),
              ("random", [3,17,29,44,58,71,89,103,120,137,155,170,190])]:
    mu, mg, eu, sg, g0 = stats(E)
    print("  %-24s %.4f  %.4f  %.4f  %.4f       %.4f  %.4f" %
          (nm, mu, mg, eu, (7.0/6)*(mg-1.0/7), sg, g0))

# adversarial min of E[maxgap]
random.seed(7); best = 1.0; bestE = None
for trial in range(40):
    H = random.choice([16,20,26,40,70]); E = sorted(random.sample(range(1,H+1),13)); cur = stats(E, 2500)[1]
    for _ in range(70):
        i = random.randrange(13); new = random.randint(1, random.choice([20,40,80]))
        if new in E: continue
        cand = sorted(set(E[:i]+E[i+1:]+[new]))
        if len(cand) < 13: continue
        c = stats(cand, 2500)[1]
        if c < cur - 1e-5: E = cand; cur = c
    if cur < best: best = cur; bestE = E
print()
print("adversarial inf E[maxgap] = %.5f  (> 1/7 = %.5f, margin %+.5f)" % (best, 1.0/7, best - 1.0/7))
print("extremal:", bestE)
print()
print("=> DENSITY FLOOR mu_1/7 > 0  <=  inf_E E[maxgap] > 1/7 (mean max-gap, margin ~0.06).")
print("   Elementary reverse-Markov replaces the E[U] inclusion-exclusion; cleaner target.")
