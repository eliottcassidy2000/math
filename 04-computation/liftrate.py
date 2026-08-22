"""liftrate: profile / lift / admissibility / epoch-identity toolkit (THM-3029).

RECONSTRUCTED 2026-08-03 (boxeph). The original module lived only in the
deleted session worktree /tmp/math-wt-coinC2/04-computation and was never
committed, breaking the reproduction chain of THM-3029's referee
amm12592_floor_rate_attained_thm3029.py. The definitions below are taken
verbatim from the committed amm12592_profile_monotonicity_thm3029.py, whose
module body (kept here, at module level, as in the original) prints exactly
the first 7 lines of the committed expected output
05-knowledge/results/amm12592_floor_rate_attained_thm3029.out -- that
byte-level match is the validation of this reconstruction.

Semantics (THM-2966 / THM-3002 / THM-3026):
  A block at degree d is a coefficient vector (delta_k)_{k=0..d} in the basis
      B_{d,k}(x) = x^{d-k} (1-x)^k,
  ADMISSIBLE iff  |delta_k| <= binom(d,k)  and  delta_k == binom(d,k) (mod 2)
  (Lucas-box capacity + parity).
  prof(R, g1, g2, D0)[i] = floor((g1/g2)*(R+i)) + D0,  i = 0..R-1.
  lift_block(delta, d, dp): convolution with [binom(dp-d,k)]_k, the admissible
  block representing the CONSTANT 1 (x + (1-x) = 1); re-expresses the same
  polynomial at degree dp >= d, still admissible (THM-3026 (L)+(M)).
  epoch_identity(R, sol, d): the exact epoch-closure identity of THM-3002 (*):
      sum_{i=0}^{R-1} x^i Delta_i(x) == (1-x)^{R-1}   in Z[x],
      Delta_i(x) = sum_k delta_{i,k} B_{d_i,k}(x).
  eff(R, d) = max_i d_i/(R+i) as an exact Fraction (the effective rate;
  guards against the D0 slack trap of THM-3029 sec. 2).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amm12592_gamma35_beam_deathstar as beam
from math import comb
from fractions import Fraction as F

def prof(R, g1, g2, D0): return [(g1 * (R + i)) // g2 + D0 for i in range(R)]

def conv(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b): r[i + j] += u * v
    return r

def lift_block(delta, d, dp):
    """lift an admissible block at d to degree dp >= d"""
    if dp == d: return list(delta)
    return conv(delta, [comb(dp - d, k) for k in range(dp - d + 1)])

def admissible(delta, d):
    if len(delta) - 1 > d: return False
    return all(abs(v) <= comb(d, k) and (v - comb(d, k)) % 2 == 0 for k, v in enumerate(delta))

def epoch_identity(R, sol, d):
    acc = [0] * (6 * R + 16)
    for i, de in enumerate(sol):
        for k, v in enumerate(de):
            if v:
                for t, c in enumerate(beam.basis_poly(d[i], k)): acc[i + t] += v * c
    return beam.trim(acc) == beam.trim(beam.qpow(R - 1))

def eff(R, d): return max(F(d[i], R + i) for i in range(R))

# ---------------------------------------------------------------------------
# Module-level demonstration (as in the original lost module): lines 1-7 of
# the committed expected output of THM-3029's referee. Runs at import time.
# ---------------------------------------------------------------------------
print("R=32: solve at gamma=1/2, D0=3, then LIFT onto sub-3/5 profiles.")
R = 32
src = prof(R, 1, 2, 3)
sol, msg = beam.solve(R, g1=1, g2=2, D0=3, beam=250, ctrl=2, span=2)
print(f"   source solve: {msg}, verify={beam.verify(R, sol, g1=1, g2=2, D0=3) if sol else False}")
print(f"   source profile eff rate = {float(eff(R, src)):.6f}")

for (nm, g1, g2, D0) in [("0.599", 599, 1000, 0), ("0.598", 299, 500, 0),
                         ("0.59799", 59799, 100000, 0), ("gamma* exactly-ish", 5979875, 10000000, 0)]:
    tgt = prof(R, g1, g2, D0)
    if any(tgt[i] < src[i] for i in range(R)):
        bad = [i for i in range(R) if tgt[i] < src[i]][:5]
        print(f"   {nm:20s}: NOT pointwise larger (first bad rows {bad}) -- lift unavailable")
        continue
    lifted = [lift_block(sol[i], src[i], tgt[i]) for i in range(R)]
    ok_adm = all(admissible(lifted[i], tgt[i]) for i in range(R))
    ok_id = epoch_identity(R, lifted, tgt)
    print(f"   {nm:20s}: lifted -> all blocks admissible: {ok_adm};  epoch identity: {ok_id};"
          f"  eff rate = {float(eff(R, tgt)):.6f}")
