#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S23: THE Z_2 INVOLUTION LENS on the LRC tight locus --
connecting to the tournament project's self-complementary (SC) structure.

DICHOTOMY SYNTHESIS (user steer: add/mult, odd/even, pos/neg, rat/irrat):
The tournament half of THIS project is governed by the Z_2 COMPLEMENT involution
(T <-> T^op, self-complementary tournaments, merged metagraph G_n/Z_2, THM-024).
The LRC tight locus is protected by -1 in (Z/13)* : the involution sigma(i)=-i=13-i,
under which the AP {1..12} is SET-INVARIANT.  These are the SAME Z_2.  At t=1/13 the
witness points {i/13} are reflection-symmetric through 1/2 (i/13 <-> 1 - i/13), a
'self-complementary' point configuration.  The reflection sign-BALANCES the
resonance/safe sum (pos/neg dichotomy).

We test, sharply:
 (1) are ALL tight (dilated-AP) families sigma-invariant (as residue sets)?
 (2) SIGN-BALANCE: for sigma-invariant families the safe-measure resonance sum
     pairs k <-> -k; does the witness point set's reflection symmetry force a
     balanced (real, sign-paired) safe sum?
 (3) does sigma-invariance SHARPEN tight-vs-loose beyond residue-pinning?  I.e.
     among residue-pinned families, do sigma-asymmetric ones have systematically
     higher M (looser)?  (RP+SV was necessary-not-sufficient, mac-mini S13.)
 (4) the tournament bridge: is the tight AP the 'self-complementary' fixed point
     of the LRC Z_2, mirroring SC tournaments as fixed points of the complement?
"""
from fractions import Fraction
import numpy as np, random

P = 13
LO, HI = Fraction(1, 13), Fraction(2, 25)

def M_exact(v, Qc=None):
    S = int(sum(abs(x) for x in v)); Q = Qc or min(4*S, 3*max(abs(x) for x in v)+40)
    va = np.array(v, dtype=np.int64); bn, bd = 0, 1
    for q in range(2, Q+1):
        a = np.arange(1, q, dtype=np.int64); r = np.outer(va, a) % q
        d = np.minimum(r, q-r); bq = int(d.min(axis=0).max())
        if bq*bd > bn*q: bn, bd = bq, q
    return Fraction(bn, bd)

def residues(v): return sorted(x % P for x in v)
def sigma_invariant_residues(v):
    """Is the residue SET invariant under i -> -i mod 13?"""
    R = set(x % P for x in v)
    return R == set((-x) % P for x in R)

AP = list(range(1, 13))

print("=== (1) are all dilated-AP (tight) families sigma-invariant? ===", flush=True)
for d in range(1, 13):
    v = [d*k for k in range(1, 13)]
    inv = sigma_invariant_residues(v)
    print(f"  d={d:2d}: residues sigma-invariant = {inv}  (M={float(M_exact(v)):.5f})", flush=True)

print(flush=True)
print("=== (2) sign-balance: witness reflection symmetry of the AP at t=1/13 ===", flush=True)
# points {i/13}: reflection through 1/2 sends i/13 -> (13-i)/13.  AP set-invariant.
pts = sorted((i % P) for i in AP)
refl = sorted((P - i) % P for i in AP)
print(f"  AP residues        : {pts}", flush=True)
print(f"  reflected (13 - i) : {refl}", flush=True)
print(f"  reflection-symmetric point set: {pts == refl}", flush=True)
print("  => resonance sum pairs k <-> -k with conjugate Fourier terms => REAL, sign-balanced.", flush=True)

print(flush=True)
print("=== (3) does sigma-invariance sharpen tight-vs-loose among residue-pinned families? ===", flush=True)
# residue-pinned = residues are a full system {1..12} (permutation).  Build such families
# with lifts; split by whether the SPEED set (not just residues) is sigma-related, and by
# whether the WITNESS-scale structure is reflection-symmetric.  Compare M distributions.
random.seed(20260706)
# We test: among residue-pinned lifts, is there ANY sigma-'balanced' loose family that is
# near-tight?  And is the AP the unique sigma-invariant + minimal-height family?
def is_dilated_AP(v):
    sv = sorted(v); d = sv[0]
    return d > 0 and all(sv[k] == d*(k+1) for k in range(12))
near_tight_sigma = []; near_tight_asym = []
count = 0
for _ in range(40000):
    # residue-pinned lift: v_i = i + 13 c_i, c_i in {0,1,2}  (keeps residue i)
    c = [random.randint(0, 2) for _ in range(12)]
    v = [AP[i] + 13*c[i] for i in range(12)]
    M = M_exact(v); count += 1
    if M < Fraction(1, 12):   # near-tight (below 1/12, in the danger zone toward the gap)
        # is the LIFT sigma-symmetric?  c_i vs c_{residue -i}: pair i with (13-i)-residue runner
        # runner with residue i is index i-1; its sigma-partner has residue 13-i = index 12-i
        sym = all(c[i] == c[11-i] for i in range(12))   # reflection of the lift vector
        (near_tight_sigma if sym else near_tight_asym).append((float(M), v, c))
print(f"  scanned {count} residue-pinned lifts; near-tight (M<1/12): "
      f"{len(near_tight_sigma)} sigma-symmetric-lift, {len(near_tight_asym)} asymmetric-lift", flush=True)
if near_tight_sigma:
    ms = [x[0] for x in near_tight_sigma]
    print(f"    sigma-symmetric-lift near-tight M range: [{min(ms):.5f}, {max(ms):.5f}]", flush=True)
if near_tight_asym:
    ma = [x[0] for x in near_tight_asym]
    print(f"    asymmetric-lift near-tight M range:      [{min(ma):.5f}, {max(ma):.5f}]", flush=True)
# the sharp question: are the ONLY M=1/13 the sigma-symmetric ones (dilations)?
tights = [(float(M_exact([AP[i]+13*c[i] for i in range(12)])), c)
          for c in [tuple(random.randint(0,2) for _ in range(12)) for _ in range(3000)]]
n13 = [c for (m,c) in tights if abs(m - 1/13) < 1e-9]
print(f"  among 3000 lifts, {len(n13)} at M=1/13; all reflection-symmetric lifts? "
      f"{all(all(c[i]==c[11-i] for i in range(12)) for c in n13) if n13 else 'n/a'}", flush=True)

print(flush=True)
print("=== (4) the tournament bridge: AP = self-complementary fixed point of LRC Z_2 ===", flush=True)
print("  SC tournaments = fixed points of complement (T=T^op up to iso), have -1 anti-aut (THM-024).", flush=True)
print("  AP {1..12} = fixed set of sigma:i->-i in (Z/13)*, tight at the roots of unity.", flush=True)
print("  Both: the order-2 element of the governing group pins the extremal/fixed configuration.", flush=True)
print(flush=True)
print("VERDICT: if all tight = sigma-invariant AND asymmetric lifts are systematically looser", flush=True)
print("=> sigma-invariance (LRC self-complementarity) is a genuine tight-locus condition,", flush=True)
print("the same Z_2 that governs the tournament SC structure -- one involution, two halves.", flush=True)
