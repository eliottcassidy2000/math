#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S78 -- THE OPEN CORE IS FINITE: the MSS velocity bound + the modern LRC frontier.

INSPIRATION (web, July 2026) the repo was unaware of (grep: zero mentions of Malikiosis/Santos/Schymura/
Rosenfeld/product bound):

  Malikiosis-Santos-Schymura 2024 (arXiv:2411.06903, "Linearly-exponential checking is enough"):
    UNCONDITIONALLY, to verify LRC for n+1 runners it suffices to check velocity tuples with all
    velocities <= C(n+1,2)^(n-1) <= n^(2n).  (n = # nonzero velocities.)  [improves Tao 2018 n^{O(n^2)}]

  Frontier (proven, elementary prime-filtering + this bound, NO analysis):
    Barajas-Serra 2008: 7 runners.  Rosenfeld 2025: 8 runners.  2025-26: 9 and 10 runners.

  Rosenfeld's method (8 runners, arXiv:2509.14111): Malikiosis product bound (upper bound on a
    counterexample's velocity product) + PRIME FILTERING (Lemma 6/7: for a prime p meeting a
    structural covering condition mod (k+1)p, p divides every counterexample product) -> the forced
    prime product exceeds the MSS bound -> contradiction.  27 primes {31..163}, backtracking search.

CONSEQUENCE FOR LRC14 (14 runners = 13 nonzero speeds):
  a counterexample has ALL 13 speeds <= C(14,2)^12 = 91^12 ~ 3.2e23.  FINITE.
  => the repo's 15-session premise "THM-523 does NOT bound speeds -> UNBOUNDED far configs are the real
     difficulty" is OBSOLETE. There are no unbounded speeds in a counterexample; the huge-speed analytic
     work (equidistribution HYP-3786/3788, signed correction HYP-3787, >=7-huge cross-harmonic) is really
     EFFECTIVIZING the finite ceiling 3e23 down to the searchable range, NOT handling infinity.

  The repo's covering reduction (THM-523) = the MSS/Rosenfeld reduction to the gcd=1 hard core.
  The repo's band-prime reduction (HYP-3750, primes {17,19,23}) = Rosenfeld's PRIME FILTERING (Lemma 6/7).
  => the repo already has the SOTA tool; the frontier template says LEAN INTO prime-filtering (elementary),
     complementing (or instead of) the analytic route.
"""
from math import comb, log10, lcm

n_runners = 14
n_speeds = n_runners - 1                       # 13
# MSS: verify LRC for (n+1) runners <=> check velocities <= C(n+1,2)^(n-1), n = # speeds
B = comb(n_runners, 2) ** (n_speeds - 1)       # C(14,2)^12 = 91^12
print("="*84)
print("MSS UNCONDITIONAL VELOCITY BOUND for LRC14 (14 runners, 13 speeds)")
print("="*84)
print(f"  suffices to check velocities <= C(14,2)^12 = 91^12 = {B}")
print(f"  ~ 10^{log10(B):.2f}   (looser n^2n = 13^26 ~ 10^{26*log10(13):.1f})")
print(f"  => a counterexample (if any) has all 13 speeds in [1, {B:.3g}]. FINITE.")

print()
print("="*84)
print("THE FINITE OPEN CORE vs the repo's thresholds (all << the MSS ceiling)")
print("="*84)
rows = [
    ("construction patch n(n-1)", n_runners*(n_runners-1), "forced single-patch multiple (S73)"),
    ("lcm(13,14)", lcm(13, 14), "= n(n-1); the one-slot cover of q=13 AND q=14"),
    ("lazy-cut bounded ceiling", n_runners*(n_runners-1), "rigorous at n=12 (HYP-3782)"),
    ("signed-correction thresh (max)", 259, "per-core N/(3(1-2r)|L_C|) (HYP-3787)"),
    ("MSS unconditional ceiling", B, "ALL counterexample speeds below this (finite)"),
]
for name, val, note in rows:
    print(f"  {name:32s} {val:<26} {note}")
print()
print("  => the interval (182, 3.2e23] is a FINITE (astronomically large) window the analytic tools")
print("     must clear -- effectivize the ceiling, do not prove asymptotic equidistribution.")

print()
print("="*84)
print("REFRAMED PROOF MAP for LRC14 (following the Rosenfeld/MSS template)")
print("="*84)
print("  [MSS]     speeds bounded: all 13 speeds <= 91^12 (finite).")
print("  [THM-523] covering reduction: non-covering killed by q-witness; reduce to covering 13-sets.")
print("  [HYP-3750]band-prime filtering {17,19,23}: the repo's analog of Rosenfeld Lemma 6/7 --")
print("            LEAN IN: check the Lemma-6 covering condition per band prime (elementary, SOTA-proven).")
print("  [S77]     safe-band residue frame + deep-well isolation: bulk M>=0.108>>1/14 (slack to 1/14).")
print("  [lazy-cut]search the bounded residual (target M<1/14 for speed, HYP-3792 slack lever).")
print("  [S73-75]  effectivize the (182, ceiling] window (equidistribution/signed correction, now with a")
print("            FINITE upper limit 91^12, so single/few large speeds are cleanly bounded).")
print()
print("HONEST: MSS makes the open core FINITE, not FEASIBLE (91^12 is unsearchable). It corrects the")
print("'unbounded is the difficulty' framing, connects the repo to the frontier (Rosenfeld prime-filtering =")
print("HYP-3750), and redirects the huge-speed effort to EFFECTIVE bounds. It does not by itself close LRC14.")
print("DONE.")
