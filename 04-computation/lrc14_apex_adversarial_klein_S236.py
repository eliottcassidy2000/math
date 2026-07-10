#!/usr/bin/env python3
"""
Adversarial mu_inf minimization at k=7 + the consecutive apex at k=8..13
(klein-2026-07-10-S236, THM-687(D) data). Imports the exact evaluator from
lrc14_two_scale_cluster_limit_klein_S236. Output:
05-knowledge/results/lrc14_apex_adversarial_klein_S236.out
"""
import sys, random
sys.path.insert(0, "04-computation")
from lrc14_two_scale_cluster_limit_klein_S236 import mu_inf
from fractions import Fraction as F
random.seed(236)

P = [1, 2, 3, 4, 5, 6]
cands = []
for d in range(1, 6):
    cands.append(("AP d=%d" % d, [j * d for j in range(7)]))
for d in range(1, 5):
    cands.append(("AP d=%d detuned" % d, [j * d for j in range(6)] + [6 * d + 1]))
cands.append(("geometric", [0, 1, 2, 4, 8, 16, 32]))
cands.append(("pairs", [0, 1, 7, 8, 14, 15, 21]))
for i in range(30):
    E = sorted(random.sample(range(0, 31), 7))
    cands.append(("rand%d" % i, [e - E[0] for e in E]))
best = None
for name, E in cands:
    mi, mP = mu_inf(P, E)
    if best is None or mi < best[0]:
        best = (mi, name, E)
    if float(mi) < 0.08:
        print(f"  low: {name:16s} E={E}  mu_inf = {float(mi):.6f}")
print(f"\nminimum over {len(cands)} candidate shapes: mu_inf = {best[0]} = "
      f"{float(best[0]):.6f}  ({best[1]}, E={best[2]})")
for k in range(8, 14):
    P2 = list(range(1, 14 - k)) or [1]
    mi, mP = mu_inf(P2, list(range(k)))
    print(f"k={k}: P={P2} E=0..{k-1}: mu_inf = {float(mi):.6f} "
          f"{'>0' if mi > 0 else 'ZERO'}")
