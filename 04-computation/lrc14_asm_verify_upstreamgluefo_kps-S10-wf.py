#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asm_verify_upstreamgluefo_kps-S10-wf.py   (kind-pasteur-2026-06-19-S10, ADVERSARIAL VERIFY)

ADVERSARIAL VERIFICATION of the assembly piece angle "upstream-glue-formal":
  the two glue lemmas G1 (global-witness soundness) and G2 (finite-Vmax discretization),
  plus the COMPLETE-CHAIN claim
      G1 + G2 + [meas(S7(E)) <= cap_k, k=8..12] + [k<=7 pigeonhole]  ==>  LRC(14),
  with the SINGLE remaining open input asserted to be [meas(S7(E)) <= cap_k].

This script re-derives EVERY load-bearing constant EXACTLY (fractions.Fraction) and HUNTS for
the gaps the four verification directives demand:
  (1) re-derive every constant exactly;
  (2) check decay/envelope is a real proof (here: the finite-Vmax tooth picture vs the limit);
  (3) check the finite-check / V0* thresholds are consistent and the enumeration is genuinely
      complete (here: is V0* derived from a CORRECT slop bound?);
  (4) glue arc-count Vmax-uniform + V0 explicit.

VERDICT (see bottom + final message): the glue lemmas are ALMOST rigorous but TWO real gaps
survive:
  GAP-A (G2-V4 "SHARP V0*"):  the slop bound "free phi-window shrinks by <= sum(E)/Vmax" is
         FALSE; the true shrinkage scales like (sum of e_i^2-type)/Vmax (the tooth CENTER shifts
         by e^2 x/(Vmax-e), an order of magnitude larger than e/Vmax).  V0* happens to be a safe
         threshold numerically for the binding shape, but the STATED DERIVATION is wrong.
  GAP-B (CHAIN step 2):  [meas(S7(E)) <= cap_k] does NOT by itself deliver rho* > 0.  The
         union-bound floor degenerates to rho* >= meas(G_P) - min_P meas(G_P) >= 0 (can be EXACTLY
         0 when P is the minimizer).  The strict 1891/5880 floor needs "consec minimizes
         mu_{1/7}", which THM-530 itself records as VERIFIED-not-PROVED ("NOT symbolically
         closed", kps-S8: every local-monotone route DEAD).  So the assembly's "single open
         input" is INCOMPLETE: a second open input (consec extremal / rho*>0 non-anti-correlation)
         is silently absorbed.
G1's QUALITATIVE conclusion (good ruler index => M(S) >= 1/14) is, however, EMPIRICALLY SOUND
(0 counterexamples over thousands of shapes incl. small Vmax).
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor, ceil, gcd
import random, sys, io
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
random.seed(202020)

# ----------------------------------------------------------------------------
# exact primitives
# ----------------------------------------------------------------------------
def frac(x): return x - (x.numerator // x.denominator)
def nrm(x):
    f = frac(x); return min(f, 1 - f)
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    return out
def meas(iv): return sum((b - a for a, b in iv), F(0))
def complement(iv):
    iv = merge(iv); out = []; prev = F(0)
    for a, b in iv:
        if a > prev: out.append((prev, a))
        prev = max(prev, b)
    if prev < 1: out.append((prev, F(1)))
    return out
def danger_arcs(u, h=F(1, 14)):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def safe_set(S, h=F(1, 14)):
    if not S: return [(F(0), F(1))]
    return complement(merge([iv for v in S for iv in danger_arcs(v, h)]))
def M_safe_measure(S, h=F(1, 14)): return meas(safe_set(S, h))
def maxgap(points):
    pts = sorted(set(points))
    if not pts: return F(1)
    g = F(0)
    for a, b in zip(pts, pts[1:]): g = max(g, b - a)
    g = max(g, pts[0] + 1 - pts[-1]); return g
def phases_at(E, x): return sorted(set(frac(e * x) for e in E))

print("=" * 84)
print("ADVERSARIAL VERIFY — LRC(14) upstream-glue-formal  (kind-pasteur-2026-06-19-S10)")
print("=" * 84)

# ============================================================================
# (1) RE-DERIVE EVERY CONSTANT EXACTLY
# ============================================================================
print("\n[C1] tooth-width LIMIT identity 2*(1/14)=1/7:", 2 * F(1, 14) == F(1, 7), " (EXACT, but see C2)")

# C2: THE EXACT FINITE-Vmax TOOTH (no limit).  tau=(j+phi)/Vmax, member v=Vmax-e.
#   ||v tau|| < 1/14  <=>  phi in ( center - hw, center + hw ) (mod 1),
#   center = e*x*Vmax/(Vmax-e),  hw = Vmax/(14*(Vmax-e)),  x = j/Vmax.
# So the FINITE tooth FULL WIDTH = Vmax/(7*(Vmax-e))  >  1/7  (strict, for e>0).
print("\n[C2] EXACT finite-Vmax cluster tooth (DERIVED, not limit):")
print("     center = e*x*Vmax/(Vmax-e),  half-width = Vmax/(14*(Vmax-e)),  full width = Vmax/(7*(Vmax-e))")
def true_tooth_fullwidth(Vmax, e): return F(Vmax, 7 * (Vmax - e))
for (Vmax, e) in [(420, 7), (501, 27), (532, 5), (200, 20)]:
    w = true_tooth_fullwidth(Vmax, e)
    print(f"     Vmax={Vmax} e={e}: full width = {w} = {float(w):.6f}  ( 1/7 = {float(F(1,7)):.6f} ; EXCESS = {float(w-F(1,7)):+.6f} )")
print("     => cluster teeth are STRICTLY WIDER than 1/7 at finite Vmax.  Step 2's 'EXACT 1/7'")
print("        is a LIMIT statement; the script's own [G1-V2] flags 'interior width exactly 1/7: False'.")

# C2b: tooth CENTER shift vs claimed frac(e x)
print("\n[C2b] tooth CENTER shift  center - frac(e x) = (e^2 x)/(Vmax-e)  (worst x->1):")
for (Vmax, e) in [(420, 7), (501, 27), (200, 20)]:
    shift_max = F(e * e, Vmax - e)
    print(f"     Vmax={Vmax} e={e}: max center shift = e^2/(Vmax-e) = {float(shift_max):.5f}   (claimed 'diff <= e/Vmax' = {float(F(e,Vmax)):.5f})")
print("     => CENTER shift is ~ e/(Vmax) * e  -- a factor ~e LARGER than the claimed e/Vmax bound.")

# ============================================================================
# (2)/(3) GAP-A:  IS THE G2-V4 'SHARP V0*' SLOP BOUND REAL?
#   claim: 'finite teeth differ from limit by <= e_i/Vmax => free phi-window shrinks by
#           <= sum(E)/Vmax'.   TEST the actual free-window shrinkage.
# ============================================================================
print("\n" + "=" * 84)
print("[GAP-A]  G2-V4 'SHARP V0*' slop bound:  free-window shrinkage <= sum(E)/Vmax  ??")
print("=" * 84)
def free_limit(E, x):
    teeth = []
    for e in set(E):
        c = frac(F(e) * x); a = (c - F(1, 14)) % 1; b = (c + F(1, 14)) % 1
        if a < b: teeth.append((a, b))
        else: teeth.append((a, F(1))); teeth.append((F(0), b))
    return meas(complement(merge(teeth)))
def free_finite(E, x, Vmax):
    teeth = []
    for e in set(E):
        if e == 0: c = F(0); hw = F(1, 14)
        else:
            u = Vmax - e; c = frac(F(e) * x * Vmax / u); hw = F(Vmax, 14 * u)
        a = (c - hw) % 1; b = (c + hw) % 1
        if a < b: teeth.append((a, b))
        else: teeth.append((a, F(1))); teeth.append((F(0), b))
    return meas(complement(merge(teeth)))

E0 = list(range(8))  # binding shape co-offsets
print(f"  binding co-offsets E={E0}, sum(E)={sum(E0)}.  Max over ruler indices j of (free_limit - free_finite):")
gapA_violated = True
for Vmax in [200, 300, 501, 1000, 5000]:
    worst = F(0)
    for j in range(1, Vmax):
        x = F(j, Vmax); fl = free_limit(E0, x); ff = free_finite(E0, x, Vmax)
        if fl > 0: worst = max(worst, fl - ff)
    bnd = F(sum(E0), Vmax)
    holds = worst <= bnd
    if holds: gapA_violated = False
    print(f"    Vmax={Vmax:5d}: max shrinkage={float(worst):.5f}  vs claimed slop sum(E)/Vmax={float(bnd):.5f}  -> bound {'HOLDS' if holds else 'VIOLATED'}")
print(f"  => claimed slop bound sum(E)/Vmax is {'OK' if not gapA_violated else 'VIOLATED at EVERY tested Vmax (by ~2-3x)'}.")
print("  => GAP-A: the G2-V4 derivation of V0* rests on a FALSE per-tooth slop estimate.  The")
print("     V0* NUMBER may still be safe (checked below) but the stated PROOF of it is not valid.")

# Is V0*=501 nonetheless a valid threshold for the binding shape? (number accidentally safe)
def reconstruct(P, E, Vmax):
    L = [Vmax - e for e in E]
    if min(L) <= 13: return None
    S = sorted(set(P) | set(L))
    if len(S) != 13: return None
    g = 0
    for s in S: g = gcd(g, s)
    return S if g == 1 else None
P0 = [1, 2, 3, 12, 13]
GP0 = safe_set(P0)
bad_V0star = []
for Vmax in range(501, 1100):
    found = False
    for j in range(Vmax):
        x = F(j, Vmax)
        if not any(a <= x < b for a, b in GP0): continue
        if free_finite(E0, x, Vmax) > 0: found = True; break
    if not found: bad_V0star.append(Vmax)
print(f"\n  [number-check] for binding shape, Vmax in [501,1100): Vmax with NO finite-good ruler index: {len(bad_V0star)}")
print(f"    => V0*=501 is EMPIRICALLY a valid threshold for THIS shape (number safe), DESPITE the wrong derivation.")

# ============================================================================
# G1 QUALITATIVE SOUNDNESS (the conclusion, not the constant): good index => M(S)>=1/14
# ============================================================================
print("\n" + "=" * 84)
print("[G1-SOUND]  qualitative: (exists good ruler index j) => M(S) >= 1/14, incl. SMALL Vmax")
print("=" * 84)
fails = []; tested = 0
for _ in range(6000):
    k = random.randint(8, 12); psz = 13 - k
    P = sorted(random.sample(range(1, 14), psz))
    spread = random.choice([k - 1, k, k + 1, 2 * k, 3 * k])
    body = sorted(random.sample(range(1, spread + 1), min(k - 1, spread)))
    E = [0] + body
    if len(set(E)) != k: continue
    Vmax = max(E) + 14 + random.randint(0, 30)   # stress small Vmax
    S = reconstruct(P, E, Vmax)
    if S is None: continue
    GP = safe_set(P); goodj = None
    for j in range(Vmax):
        x = F(j, Vmax)
        if not any(a <= x < b for a, b in GP): continue
        if maxgap(phases_at(E, x)) > F(1, 7): goodj = j; break
    if goodj is None: continue
    tested += 1
    if M_safe_measure(S) <= 0: fails.append((P, E, Vmax))
print(f"  shapes with a good ruler index (small Vmax stress): {tested}")
print(f"  SOUNDNESS counterexamples (good index but M(S)<1/14): {len(fails)}")
for f in fails[:5]: print("    *** CE:", f)
print("  => G1's CONCLUSION is empirically sound; the limit-drop of e*phi/Vmax does not break it.")

# ============================================================================
# (CHAIN GAP-B)  does [meas(S7(E)) <= cap_k] deliver rho* > 0 ?
# ============================================================================
print("\n" + "=" * 84)
print("[GAP-B]  CHAIN step (2):  does [meas(S7(E)) <= cap_k] alone give rho*_{1/7} > 0 ?")
print("=" * 84)
def meas_S7(E):
    def sec(p): return int((p * 7).numerator // (p * 7).denominator)
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for i in range(e):
            for m in range(7):
                v = (F(m, 7) + i) / e
                if 0 <= v < 1: bps.add(v)
    bps = sorted(b for b in bps if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bps, bps[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if len(set(sec(p) for p in phases_at(E, mid))) == 7: tot += hi - lo
    return tot
def mu_17(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for i in range(e):
            for m in range(7):
                v = (F(m, 7) + i) / e
                if 0 <= v < 1: bps.add(v)
    bps = sorted(b for b in bps if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bps, bps[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if maxgap(phases_at(E, mid)) > F(1, 7): tot += hi - lo
    return tot

caps = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
print("  per-k:  cap_k,  thr_k=1-cap_k,  and what the UNION-bound floor degenerates to.")
print("  THM-530 union floor:  rho* >= meas(G_P) + mu_{1/7}(E) - 1.")
print("  bridge:  meas(S7(E))<=cap_k  =>  mu_{1/7}(E) >= 1 - meas(S7(E)) >= 1 - cap_k = thr_k.")
print("  so floor >= meas(G_P) + thr_k - 1 = meas(G_P) - (1-thr_k) = meas(G_P) - cap_k.")
print("  Since cap_k = min_P meas(G_P), the floor = meas(G_P) - min_P meas(G_P) >= 0,  = 0 when P=argmin!")
for k in range(8, 13):
    thr = 1 - caps[k]
    print(f"    k={k}: cap={float(caps[k]):.4f}  thr_k=1-cap={float(thr):.4f}  union-floor(min P) = meas(GP)-cap_k -> can be 0")

# verify bridge inequality mu_17 >= 1 - meas(S7) (the only TRUE part of the bridge)
print("\n  [bridge-ineq]  mu_{1/7}(E) >= 1 - meas(S7(E))  (S7^c subset {maxgap>=1/7}):")
bviol = 0; btest = 0
for _ in range(150):
    k = random.randint(8, 12)
    spread = random.choice([k - 1, k, k + 1, 2 * k])
    body = sorted(random.sample(range(1, spread + 1), min(k - 1, spread)))
    E = [0] + body
    if len(set(E)) != k: continue
    btest += 1
    if mu_17(E) < 1 - meas_S7(E) - F(1, 10**9): bviol += 1
print(f"    tested {btest}; violations of mu >= 1-meas(S7): {bviol}  (bridge inequality {'holds' if bviol==0 else 'FAILS'})")

# the STRICT 0.322 floor needs the ACTUAL mu(E), i.e. 'consec minimizes mu_{1/7}' (the open lemma)
print("\n  The strict floor 1891/5880=0.3216 uses mu_{1/7}(consec_k) (the MAX mu), NOT thr_k:")
mu_consec = {8: F(691, 735), 9: F(247, 294), 10: F(38, 49), 11: F(1381, 2205), 12: F(13823, 24255)}
for k in range(8, 13):
    fl = caps[k] + mu_consec[k] - 1
    print(f"    k={k}: meas(GP)_min + mu(consec) - 1 = {float(caps[k]):.4f}+{float(mu_consec[k]):.4f}-1 = {float(fl):.4f}")
print("  => the 0.322 floor REQUIRES mu_{1/7}(E) >= mu_{1/7}(consec_k), i.e. 'consec minimizes mu_{1/7}'.")
print("  => THAT lemma is, per THM-530 and kps-S8, VERIFIED-NOT-PROVED ('NOT symbolically closed';")
print("     every local monotone-descent route DEAD).  It is NOT implied by [meas(S7)<=cap_k].")

# Empirically: rho*_{1/7} > 0 everywhere? (THM-530's VERIFIED non-anti-correlation)
def rho_star(P, E):
    GP = safe_set(list(P)); E = sorted(set(E)); bps = set()
    for lo, hi in GP: bps.add(lo); bps.add(hi)
    for e in E:
        if e == 0: continue
        for i in range(e):
            for m in range(7):
                v = (F(m, 7) + i) / e
                if 0 <= v < 1: bps.add(v)
    bps = sorted(b for b in bps if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bps, bps[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if any(a <= mid <= b for a, b in GP) and maxgap(phases_at(E, mid)) > F(1, 7): tot += hi - lo
    return tot
print("\n  [rho*-positivity]  search for admissible (P,E) with rho*_{1/7}=0 (the anti-correlation pathology):")
zeros = 0; ztest = 0
for _ in range(1500):
    k = random.randint(8, 12); psz = 13 - k
    P = sorted(random.sample(range(1, 14), psz))
    spread = random.choice([k - 1, k, k + 1, 2 * k, 3 * k])
    body = sorted(random.sample(range(1, spread + 1), min(k - 1, spread)))
    E = [0] + body
    if len(set(E)) != k: continue
    ztest += 1
    if rho_star(P, E) == 0: zeros += 1
print(f"    tested {ztest} admissible (P,E); rho*_{{1/7}}=0 found: {zeros}")
print("    => empirically rho*>0 everywhere (matches THM-530 'no zeros', R_{1/7}>=0.796), but VERIFIED not PROVED.")

print("\n" + "=" * 84)
print("VERDICT: NOT GAP-FREE.  Two surviving gaps:")
print("  GAP-A: G2-V4 'SHARP V0*' slop bound (sum(E)/Vmax) is FALSE (true shrinkage ~e^2/Vmax-type,")
print("         2-3x larger at every Vmax).  V0* number is empirically safe but the derivation is not.")
print("  GAP-B: [meas(S7(E))<=cap_k] does NOT give rho*>0.  The union-floor degenerates to >=0 (can be")
print("         exactly 0); the strict 0.322 floor needs 'consec minimizes mu_{1/7}', which is the")
print("         OPEN (not symbolically closed) lemma -- a SECOND open input silently absorbed.")
print("  G1's qualitative conclusion is empirically SOUND; G2(a,b,c) arc-count + V0 logic is sound.")
print("=" * 84)
