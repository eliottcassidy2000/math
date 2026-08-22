"""THE MONEY EXPERIMENT.  gamma* = 0.5979874356654402 is the proved floor rate and
capacity permits every gamma > gamma* at EVERY R (gamma_c(R) < gamma* for all R).
Does the CONSTRUCTION close at a rate just above gamma*, with BOUNDED D0?

If gamma = 0.598 closes at R = 8,16,32,64 with small D0, then C <= 1.598 < 8/5 = 1.6
for those epochs -- an improvement on the repo's standing upper bound, and evidence
that the archimedean floor log_5(5 phi^2) is TIGHT.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # repo 04-computation (was a deleted session worktree)
import amm12592_gamma35_beam_deathstar as beam
from fractions import Fraction as F
import time

RATES = [("3/5   (control)", 3, 5),
         ("0.599", 599, 1000),
         ("0.598", 299, 500),
         ("0.5980", 1495, 2500),
         ("100/167 simplest in gap", 100, 167),
         ("0.59799 (just above gamma*)", 59799, 100000)]

def eff(R, g1, g2, D0):
    return max(F((g1 * (R + i)) // g2 + D0, R + i) for i in range(R))

print("rate                          R    result                         D0  eff.rate")
for (name, g1, g2) in RATES:
    for R in (8, 16, 32, 64):
        hit = None
        t0 = time.time()
        for D0 in range(0, 4):
            for (bw, ctrl, span) in [(250,2,2),(800,2,2),(800,2,3),(2000,2,2)]:
                sol, msg = beam.solve(R, g1=g1, g2=g2, D0=D0, beam=bw, ctrl=ctrl, span=span)
                if sol and beam.verify(R, sol, g1=g1, g2=g2, D0=D0):
                    hit = (D0, bw, ctrl, span); break
                if time.time() - t0 > 240: break
            if hit or time.time() - t0 > 240: break
        if hit:
            print(f"{name:28s} {R:4d}  CLOSES + verified            D0={hit[0]}  "
                  f"{float(eff(R,g1,g2,hit[0])):.6f}", flush=True)
        else:
            print(f"{name:28s} {R:4d}  did not close                 -", flush=True)
