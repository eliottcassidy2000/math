"""THE D0 TRAP.  d_i = floor(gamma(R+i)) + D0, so the EFFECTIVE rate at epoch R is
        gamma_eff = max_i d_i/(R+i),
which exceeds gamma by up to D0/R.  At R=32 with D0=4 that is +0.125 -- enormous.
So "gamma = 0.48 closes at R=32 with D0=4" does NOT witness C <= 1.48; the honest
rate it witnesses is gamma_eff.  Only as R -> infinity does a bounded D0 become free.

Correct statement of the constant: C = 1 + gamma requires ONE gamma closing at EVERY
R with BOUNDED D0.  Since gamma_c(R) increases to gamma*, the binding regime is large R.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # repo 04-computation (was a deleted session worktree)
import amm12592_gamma35_beam_deathstar as beam
from gammac import gamma_c, GSTAR
from fractions import Fraction as F

def eff_rate(R, g1, g2, D0):
    return max(F((g1 * (R + i)) // g2 + D0, R + i) for i in range(R))

def closes(R, g1, g2, D0, bw=800):
    for (b, c, s) in [(250,2,2),(bw,2,2),(bw,2,3)]:
        sol, _ = beam.solve(R, g1=g1, g2=g2, D0=D0, beam=b, ctrl=c, span=s)
        if sol and beam.verify(R, sol, g1=g1, g2=g2, D0=D0): return True
    return False

print("Demonstrating the trap on the cases just found:")
for (R, g1, g2, D0) in [(32,1,2,3), (32,49,100,4), (32,48,100,4), (64,3,5,0), (32,3,5,0)]:
    er = eff_rate(R, g1, g2, D0)
    print(f"  R={R:3d} gamma={g1}/{g2}={g1/g2:.4f} D0={D0}:  EFFECTIVE rate = {float(er):.6f}"
          f"   (gamma_c({R}) = {gamma_c(R, lo=0.02, hi=0.75):.6f})")

print("\nHONEST construction threshold: min over (gamma, D0) that CLOSE of the effective rate.")
print("   R    gamma_c(R)    best effective rate achieved    witness (g1/g2, D0)")
DEN = 240
for R in (8, 16, 32, 64):
    gc = gamma_c(R, lo=0.02, hi=0.75)
    best = None
    for D0 in range(0, 5):
        for k in range(int(0.30 * DEN), int(0.62 * DEN) + 1):
            er = eff_rate(R, k, DEN, D0)
            if best is not None and er >= best[0]: continue
            if closes(R, k, DEN, D0):
                best = (er, f"{k}/{DEN}", D0)
    if best:
        print(f"{R:5d}   {gc:.6f}      {float(best[0]):.6f}"
              f"                     {best[1]}, D0={best[2]}", flush=True)
    else:
        print(f"{R:5d}   {gc:.6f}      none closed", flush=True)
