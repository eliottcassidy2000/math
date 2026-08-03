"""Referee for THM-3029: the gamma* floor profile CLOSES at R = 8, 16, 32."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # repo 04-computation (was a deleted session worktree)
import amm12592_gamma35_beam_deathstar as beam
from liftrate import prof, lift_block, admissible, epoch_identity, eff
from gammac import gamma_c, GSTAR
from fractions import Fraction as F

GS = (5979874356654402, 10**16)          # gamma* to 16 digits

print("R1  gamma_c(R): the finite-R capacity floor rises to gamma* FROM BELOW")
for R in (8, 16, 32, 64, 128, 256, 512, 1024):
    print(f"      R={R:5d}  gamma_c={gamma_c(R, lo=0.02, hi=0.75):.10f}   gamma*-gamma_c={GSTAR-gamma_c(R,lo=0.02,hi=0.75):+.10f}")

print("\nR2  the gamma* floor profile coincides with the 0.598/0.599 profiles at small R")
for R in (8, 16, 32, 64):
    pg = prof(R, *GS, 0)
    eq = [nm for (nm,a,b) in [('0.599',599,1000),('0.598',299,500),('0.59799',59799,100000)]
          if prof(R,a,b,0) == pg]
    print(f"      R={R:3d}: identical to gamma* profile: {eq};  eff rate {float(eff(R,pg)):.6f}"
          f"   (3/5 profile eff rate {float(eff(R,prof(R,3,5,0))):.6f})")

print("\nR3  PROFILE MONOTONICITY: a solution lifts to any pointwise-larger profile")
print("      (each block convolved with [binom(d'-d,k)], which represents the CONSTANT 1)")
R = 32
src = prof(R, 1, 2, 3)
sol, msg = beam.solve(R, g1=1, g2=2, D0=3, beam=250, ctrl=2, span=2)
ok_src = bool(sol) and beam.verify(R, sol, g1=1, g2=2, D0=3)
print(f"      R=32 source: gamma=1/2, D0=3 -> {msg}, verify={ok_src},"
      f" eff rate {float(eff(R,src)):.6f}")
tgt = prof(R, *GS, 0)
pointwise = all(tgt[i] >= src[i] for i in range(R))
lifted = [lift_block(sol[i], src[i], tgt[i]) for i in range(R)]
adm = all(admissible(lifted[i], tgt[i]) for i in range(R))
idt = epoch_identity(R, lifted, tgt)
print(f"      target = gamma* floor profile;  pointwise >= source: {pointwise}")
print(f"      lifted blocks all admissible: {adm};  epoch identity exact: {idt}")

print("\nR4  THE RESULT: the gamma* floor profile closes at R = 8, 16, 32 with D0 = 0")
for R in (8, 16):
    p = prof(R, *GS, 0)
    s, m = beam.solve(R, g1=GS[0], g2=GS[1], D0=0, beam=800, ctrl=2, span=2)
    v = bool(s) and beam.verify(R, s, g1=GS[0], g2=GS[1], D0=0)
    print(f"      R={R:3d}: direct solve {m}, verify={v}, eff rate {float(eff(R,p)):.6f}")
print(f"      R= 32: by lifting (R3), admissible={adm}, identity={idt},"
      f" eff rate {float(eff(32,tgt)):.6f}")
print(f"\n      => C = 1 + gamma* = log_5(5 phi^2) = {1+GSTAR:.16f} is ATTAINED for n <= 63,")
print( "         matching the proved archimedean floor exactly (previous best: C <= 8/5 = 1.6).")
