#!/usr/bin/env python3
"""residue6_closure_assembly_macmini_S118.py -- mac-mini-2026-07-16-S118.
CLOSING NEGATIVE RESIDUE 6: assemble the universal bound q(a,b,c) <= 47/100 (THM-904 (3))
from the ANOVA decomposition of codex's 84-weight beta + explicit tails.

Structure: q(a,b,c) = beta_0 + sum_{pairs} D_{ij}(ratio) + T(a,b,c) where
  - singleton channels integrate to ZERO for every speed (sectors exactly uniform);
  - D(p:q) = integral of the pair channel at primitive ratio p:q (EXACT breakpoint sweep);
  - T = the zero-marginal triple channel; |T| <= tail bound via relation-lattice Fourier.
Empirics turn into the theorem: if beta_0 + (sum of the three largest COMPATIBLE pair
values) + |T|_max < 47/100 outside an explicit finite box covered by codex's scan, (3)
closes and with it -F_6(E) <= 47/490 < 0.097 -- the sole limiting sign of THM-891.
"""
import sys, os
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "04-computation")
sys.stdout.reconfigure(line_buffering=True)
import importlib.util as ilu
spec = ilu.spec_from_file_location("cert", "04-computation/lrc14_residue6_triple_certificate_codex_S18.py")
cert = ilu.module_from_spec(spec)
try:
    spec.loader.exec_module(cert)
except SystemExit:
    pass

beta = cert.beta
SECTIONS = range(7)

# ANOVA channels
mean = sum(beta(a, b, c) for a in SECTIONS for b in SECTIONS for c in SECTIONS) / Fr(343)
sing = [sum(beta(s, b, c) for b in SECTIONS for c in SECTIONS) / Fr(49) - mean for s in SECTIONS]
pairch = {}
for s1 in SECTIONS:
    for s2 in SECTIONS:
        pairch[(s1, s2)] = (sum(beta(s1, s2, c) for c in SECTIONS) / Fr(7)
                            - mean - sing[s1] - sing[s2])
tripch = {}
for s1 in SECTIONS:
    for s2 in SECTIONS:
        for s3 in SECTIONS:
            tripch[(s1, s2, s3)] = (beta(s1, s2, s3) - mean - sing[s1] - sing[s2] - sing[s3]
                                    - pairch[(s1, s2)] - pairch[(s1, s3)] - pairch[(s2, s3)])
print(f"beta_0 (Haar/generic value) = {mean} = {float(mean):.6f}")
print(f"singleton channel max |.|  = {float(max(abs(v) for v in sing)):.6f} (integrates to 0 for every speed)")
print(f"margin to 47/100 at generic: {float(Fr(47,100) - mean):.6f}")

# exact pair-channel integral at primitive ratio (p, q): joint sectors of ({px},{qx})
def D_pair(p, q):
    # breakpoints of sec(px), sec(qx): x = k/(7p), k/(7q)
    pts = sorted(set([Fr(k, 7 * p) for k in range(7 * p + 1)] + [Fr(k, 7 * q) for k in range(7 * q + 1)]))
    tot = Fr(0)
    for i in range(len(pts) - 1):
        x = (pts[i] + pts[i + 1]) / 2
        s1 = int((x * p % 1) * 7); s2 = int((x * q % 1) * 7)
        tot += pairch[(s1, s2)] * (pts[i + 1] - pts[i])
    return tot

print("\npair-channel integrals D(p:q) (exact), primitive p<q:")
vals = []
for q in range(2, 41):
    for p in range(1, q):
        if gcd(p, q) == 1:
            vals.append((abs(D_pair(p, q)), p, q))
vals.sort(reverse=True)
for v, p, q in vals[:8]:
    print(f"   D({p}:{q}) = {float(D_pair(p,q)):+.6f}  |.| = {float(v):.6f}")
big = [(float(D_pair(p, q)), p, q) for v, p, q in vals[:20]]
maxD = float(vals[0][0])
# decay check
import math
far = max(float(v) for v, p, q in vals if p + q > 20)
print(f"   max |D| = {maxD:.6f}; max |D| with p+q > 20: {far:.6f} (decay visible)")

# triple channel: worst-case bound = max over states (pointwise) and relation-tail estimate
maxT_point = max(abs(v) for v in tripch.values())
l1_slice = max(sum(abs(tripch[(s1, s2, s3)]) for s3 in SECTIONS) / 7 for s1 in SECTIONS for s2 in SECTIONS)
print(f"\ntriple channel: pointwise max = {float(maxT_point):.6f}")

# exact q check against codex's scan max
def q_exact(a, b, c):
    pts = sorted(set([Fr(k, 7 * a) for k in range(7 * a + 1)] +
                     [Fr(k, 7 * b) for k in range(7 * b + 1)] +
                     [Fr(k, 7 * c) for k in range(7 * c + 1)]))
    tot = Fr(0)
    for i in range(len(pts) - 1):
        x = (pts[i] + pts[i + 1]) / 2
        tot += beta(int((x * a % 1) * 7), int((x * b % 1) * 7), int((x * c % 1) * 7)) * (pts[i + 1] - pts[i])
    return tot
q147 = q_exact(1, 4, 7)
print(f"\nreproduction: q(1,4,7) = {q147} = {float(q147):.6f} (codex: 81/175 = {float(Fr(81,175)):.6f})")

# assembly probe: decompose q(1,4,7)
d14, d17, d47 = D_pair(1, 4), D_pair(1, 7), D_pair(4, 7)
trip147 = q147 - mean - d14 - d17 - d47
print(f"decomposition: beta0 {float(mean):+.5f} + D(1:4) {float(d14):+.5f} + D(1:7) {float(d17):+.5f}"
      f" + D(4:7) {float(d47):+.5f} + T {float(trip147):+.5f}")

# the assembled sufficient condition: worst three compatible pairs + worst triple remainder
print("\nASSEMBLY: q <= beta0 + D1 + D2 + D3 + T. Need < 0.47.")
print(f"   crude worst case: {float(mean) + 3 * maxD + float(maxT_point):.4f} (pointwise T -- too crude)")
top3 = sum(v for v, p, q in vals[:3])
print(f"   top-3 pair values (incompatible speeds in general): beta0 + top3 + T_point = "
      f"{float(mean) + float(top3) + float(maxT_point):.4f}")
print(f"   REALISTIC route: T's integral obeys the relation tail (BV-Fourier ~ C/|k1 k2 k3|);")
print(f"   empirical T at (1,4,7): {float(trip147):+.5f}; compute T over the scan to calibrate:")
worstT = 0.0; worst_at = None
for (a, b, c) in [(1,4,7),(1,2,3),(1,2,4),(2,3,5),(1,3,4),(1,6,7),(3,4,7),(1,4,5),(2,5,7),(1,2,7),(1,3,7),(2,4,7)]:
    t = float(q_exact(a, b, c) - mean - D_pair(a//gcd(a,b), b//gcd(a,b)) * 0 - D_pair(*sorted((a//gcd(a,b), b//gcd(a,b)))) - D_pair(*sorted((a//gcd(a,c), c//gcd(a,c)))) - D_pair(*sorted((b//gcd(b,c), c//gcd(b,c)))))
    if abs(t) > abs(worstT): worstT, worst_at = t, (a, b, c)
    print(f"      T{(a,b,c)} = {t:+.6f}")
print(f"   worst sampled T = {worstT:+.6f} at {worst_at}")
print("\nVERDICT INPUTS: beta0, max pair D, pair-decay, worst T -- see analysis in canon file.")
print("DONE")
