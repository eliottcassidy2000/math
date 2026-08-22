"""PROFILE MONOTONICITY (THM-3026 lemma L) turns a search result into a THEOREM.

If an epoch closes with degree profile (d_i), it closes with ANY pointwise-larger
profile (d'_i >= d_i): each block lifts by convolving its deltas with
[binom(d'-d,k)]_k, which is the admissible block representing the CONSTANT 1
(since x + (1-x) = 1).  The epoch identity sum_i x^i Delta_i = q^{R-1} is unchanged.

So the beam FAILING at a larger profile while SUCCEEDING at a smaller one is a pure
search artefact.  Exploit it: lift the gamma=1/2 solutions onto sub-3/5 profiles.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # repo 04-computation (was a deleted session worktree)
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
