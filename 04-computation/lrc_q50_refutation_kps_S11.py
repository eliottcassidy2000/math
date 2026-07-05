#!/usr/bin/env python3
"""
lrc_q50_refutation_kps_S11.py  --  kind-pasteur-2026-07-05-S11, HYP-4137.

REFUTES the Q50 conjecture / TemplateDichotomy at bound 50 (HYP-4119 mac-mini,
HYP-4127 kps-S10) by CONSTRUCTING an explicit Fin-13 family that satisfies EVERY
hypothesis of `TemplateDichotomy` (LRCTemplateSurface.lean) yet has NO rational
2/25-witness at any denominator s <= 50 -- and is NOT tight-shaped.

Mechanism (the MISTAKE-096 phenomenon at the profile level): the profile pins
residues only mod q <= 25; a witness at a modulus whose prime-power factors are
all <= 25 ("pinned-only") is HEIGHT-INDEPENDENT, but a witness at a "free"
modulus (27,29,31,32,37,41,43,47,49) depends on residues the profile does not
control.  A NEEDFREE shape (no pinned-only witness in [26,50]) can be lifted to
high height and, by CRT on the free residues, have EVERY free modulus in [26,50]
PINNED simultaneously -- killing all witnesses <= 50.  The true witness modulus
grows with height (~ c*ln(height)); NO fixed bound closes the loose branch.

This does NOT disprove LRC(14): the family has a real witness at t=a/53 (margin
>= 2/25), so `TightLooseDichotomy` (the REAL-valued loose branch, the actual
reduction lrc14_of_dichotomy_and_corner) is UNAFFECTED.  Only the bound-50
"template" refinement is refuted.

All verification is INDEPENDENT (fresh implementations of every predicate).
"""
from math import gcd
from functools import reduce
import sys

def lcm(a, b): return a * b // gcd(a, b)
L = reduce(lcm, range(2, 26))           # lcm(2..25) = 26,771,144,400

def mu(s): return (2 * s + 24) // 25     # ceil(2s/25)
def distZ(x, q):                         # distance of x/q to nearest integer, times q
    r = x % q
    return min(r, q - r)

# ---- the exact TemplateWitness / loose-branch predicate (fresh implementation) ----
def has_template_witness(base, s):
    """Is there k with mu(s) <= (v*k) mod s <= s-mu(s) for all base v?  (== 2/25
    margin at t=k/s).  Ranges k over 1..s-1 (k=0 is degenerate)."""
    m = mu(s)
    for k in range(1, s):
        if all(m <= (v * k) % s <= s - m for v in base):
            return k
    return None

def min_witness_modulus(base, hi):
    for s in range(2, hi + 1):
        if has_template_witness(base, s) is not None:
            return s
    return None

# ---- exact profile predicates (fresh; match the Lean defs) ----
def covering_family(v):                  # CoveringFamily: q|some v_i for q in 2..14
    return all(any(x % q == 0 for x in v) for q in range(2, 15))
def compressed(v):                       # forall i exists j!=i : |v_i| <= 13|v_j|
    return all(any(j != i and abs(v[i]) <= 13 * abs(v[j]) for j in range(len(v)))
               for i in range(len(v)))
def primitive(v):
    return reduce(gcd, [abs(x) for x in v]) == 1
def pinning_le25(base):                  # F5: pinned (no +-1 witness) at every q<=25
    for q in range(2, 26):
        if any(x % q == 0 for x in base):
            continue
        for b in range(1, q // 2 + 1):
            if gcd(b, q) != 1:
                continue
            if not any(x % q in (b, q - b) for x in base):
                return False
    return True
def is_tight_shape(base):                # exists c>=2 : every |v|=c*j, j in 1..12
    a = [abs(x) for x in base]
    g = reduce(gcd, a)
    for c in range(2, g + 1):
        if g % c: continue
        if all((x % c == 0 and 1 <= x // c <= 12) for x in a):
            return c
    return None

# =====================================================================
#  BUILD the counterexample from the NEEDFREE seed F1
# =====================================================================
seed = [1, 3, 4, 7, 10, 11, 13, 16, 17, 19, 23, 36]     # F1 base, NEEDFREE
assert min_witness_modulus(seed, 25) is None, "seed must have no witness <=25"

free_pp = [32, 27, 49, 29, 31, 37, 41, 43, 47]          # free/lifted prime powers in [26,50]
def allowed(fv, q):
    g = gcd(L, q); base = fv % g
    return [r for r in range(q) if r % g == base]

import random
random.seed(20260705)
def find_pin(q):
    al = [allowed(fv, q) for fv in seed]
    for _ in range(3_000_000):
        res = [random.choice(a) for a in al]
        # PIN q  <=>  no 2/25 witness at q for this residue multiset
        m = mu(q); good = False
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            if all(distZ(r * a, q) >= m for r in res):
                good = True; break
        if not good:
            return res
    return None

pin = {}
for q in free_pp:
    r = find_pin(q)
    assert r is not None, f"failed to pin {q}"
    pin[q] = r

# CRT-combine per runner over pairwise-coprime prime powers
pplist = [32, 27, 25, 49, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
M_full = reduce(lambda a, b: a * b, pplist)
def crt(pairs):
    R, M = 0, 1
    for r, m in pairs:
        g = gcd(M, m); assert (r - R) % g == 0
        t = ((r - R) // g * pow(M // g, -1, m // g)) % (m // g)
        R += M * t; M = M // g * m
    return R % M
base_res = []
for i in range(12):
    fv = seed[i]
    pairs = [(pin[pp][i] % pp, pp) if pp in pin else (fv % pp, pp) for pp in pplist]
    base_res.append(crt(pairs))

# --- SINGLE-SCALE realization: all 12 base runners in [M_full, 12 M_full] ---
# multipliers 1..12 (distinct) give ratio 12 -- as single-scale as the AP {1..12}
mult = list(range(1, 13))
base = sorted(base_res[i] + mult[i] * M_full for i in range(12))
assert len(set(base)) == 12
# --- istar (13th runner): covers 14, is the argmax, within 13x of base max ---
bmax = max(base)
istar = bmax
while not (istar > bmax and istar % 14 == 0 and istar <= 13 * bmax):
    istar += 1
V = base + [istar]                       # Fin 13 family
ISTAR_IDX = 12

# =====================================================================
#  INDEPENDENT VERIFICATION
# =====================================================================
print("=" * 70)
print("Q50 / TemplateDichotomy(bound 50) REFUTATION -- explicit counterexample")
print("=" * 70)
print(f"L = lcm(2..25) = {L}")
print(f"M_full = product of prime powers (2^5 3^3 5^2 7^2 11..47) : {len(str(M_full))} digits")
print(f"family height ~ 10^{len(str(max(V)))-1}")
print(f"base multipliers (single-scale, ratio {max(mult)}/{min(mult)}={max(mult)//min(mult)}): {mult}")

checks = {}
checks["all nonzero"] = all(x != 0 for x in V)
checks["CoveringFamily (q|some v, q=2..14)"] = covering_family(V)
checks["compressed (all i ex j: |v_i|<=13|v_j|)"] = compressed(V)
checks["primitive (gcd=1)"] = primitive(V)
checks["istar is argmax"] = all(abs(V[i]) <= abs(V[ISTAR_IDX]) for i in range(13))
checks["distinct 12 base + istar"] = len(set(V)) == 13
base12 = [V[i] for i in range(13) if i != ISTAR_IDX]
checks["F5 pinning mod <=25 (base)"] = pinning_le25(base12)
checks["NOT tight-shaped (base)"] = is_tight_shape(base12) is None

# the CRUX: no template witness at ANY s <= 50 for the 12 base runners
wit50 = [s for s in range(2, 51) if has_template_witness(base12, s) is not None]
checks["NO 2/25-witness at any s<=50 (base)"] = (wit50 == [])

print("\n--- profile hypotheses of TemplateDichotomy ---")
for k, ok in checks.items():
    print(f"  [{'PASS' if ok else 'FAIL'}] {k}")

print(f"\n  witness moduli s<=50 for base: {wit50}  (must be empty)")
mw = min_witness_modulus(base12, 400)
print(f"  TRUE minimum witness modulus (scan 2..400): {mw}")
# confirm it IS a real loose family (TightLooseDichotomy loose branch holds)
kk = has_template_witness(base12, mw) if mw else None
print(f"  real witness exists at t={kk}/{mw}  => TightLooseDichotomy loose branch HOLDS")

allok = all(checks.values())
print("\n" + ("*** ALL HYPOTHESES PASS + NO WITNESS <=50 : TemplateDichotomy(50) REFUTED ***"
              if allok else "!!! some check failed -- NOT a valid counterexample !!!"))
print(f"\nbase runners (12): {base12}")
print(f"istar: {istar}")
sys.exit(0 if allok else 1)
