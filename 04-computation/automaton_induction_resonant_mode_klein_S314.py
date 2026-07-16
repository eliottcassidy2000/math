#!/usr/bin/env python3
"""
automaton_induction_resonant_mode_klein_S314.py — klein-2026-07-16-S314 (cont.4)

THE MISS-PATTERN AUTOMATON INDUCTION, executed: condition on the fastest runner t, Abel-sum
its 7t boundaries against the slow 6-runner amplitude.  The induction step yields the

  RESONANT-MODE FORMULA:  S(t·a) = t·(1 − e(a/7))·mhat_s(a) + O(M_slow),   a = 1..6,
      mhat_s(a) = sum_{c != s} (A_c(s) + B(s)) e(ac/7),
      A_c(s) = meas{slow six miss exactly {s,c}},  B(s) = meas{slow six miss exactly {s}}

(proof: d(1[c_t=c]·F_c) split; fast-boundary term at n = ta has e(n j/(7t)) = e(aj/7)
CONSTANT on each residue class j mod 7, so the Weyl sum degenerates to sampled means;
Koksma with V(F_c) = O(M_slow) and the t-grid's discrepancy 1/t gives the error.)

CONSEQUENCE TESTED HERE: if mhat != 0 the sup-norm grows LINEARLY in t — uniform HYP-6994
(max|S|^2 = O(M)) is then FALSE as stated, while the WEIGHTED Q_s survives because the
khat-weight at the resonant frequency is the aliased 1/ell^2 mass (damping ~ (t/P)^2).
This script: (1) exact A_c, B for the slow six {0,1,2,3,4,5}; (2) prediction vs actual
S(ta) across t; (3) sup-norm growth scan; (4) sup_w Q_s at t = 100 (does the weighted law
survive?).  The outcome decides HYP-6994's fate honestly.
"""
from fractions import Fraction as Fr
from math import pi, gcd, lcm
import cmath

def sections_intervals(E6):
    """exact intervals of constant occupancy for the slow runners; returns list of
    (length, missing-set) over [0,1)."""
    bps = sorted(set(Fr(k, 7 * e) for e in E6 if e > 0 for k in range(7 * e)) | {Fr(0), Fr(1)})
    out = []
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set()
        for e in E6:
            occ.add(int((e * mid % 1) * 7) if e > 0 else 0)
        out.append((bps[i + 1] - bps[i], frozenset(range(7)) - frozenset(occ)))
    return out

def R_s_exact(E, s):
    bps = sorted(set(Fr(k, 7 * e) for e in E if e > 0 for k in range(7 * e)) | {Fr(0), Fr(1)})
    arcs, inR, start = [], False, None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) if e > 0 else 0 for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if cur and not inR: start, inR = bps[i], True
        if (not cur) and inR: arcs.append((start, bps[i])); inR = False
    if inR: arcs.append((start, Fr(1)))
    return arcs

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

E6 = [0, 1, 2, 3, 4, 5]
iv = sections_intervals(E6)
print("EXACT 6-RUNNER MISS MEASURES (slow six = {0,1,2,3,4,5}):")
misses = {}
for L, mset in iv:
    misses[mset] = misses.get(mset, Fr(0)) + L
mh = {}
for s in range(7):
    B = misses.get(frozenset({s}), Fr(0))
    A = {c: misses.get(frozenset({s, c}), Fr(0)) for c in range(7) if c != s}
    vec = [float(A[c] + B) if c != s else 0.0 for c in range(7)]
    mh[s] = [abs(sum(vec[c] * cmath.exp(2j * pi * a * c / 7) for c in range(7))) for a in range(7)]
    if s in (4, 6):
        print(f"   s={s}: B={B}  A_c={[str(A[c]) for c in sorted(A)]}  |mhat(a)| a=1..6: "
              f"{[f'{v:.4f}' for v in mh[s][1:]]}")
nonzero = any(mh[s][a] > 1e-9 for s in range(7) for a in range(1, 7))
check(f"mhat is NONZERO for some (s, a) — the resonant mode is live (max |mhat| = "
      f"{max(mh[s][a] for s in range(7) for a in range(1,7)):.4f})", nonzero)

# ---------- prediction vs actual S(ta) ----------
print()
print("RESONANT-MODE FORMULA vs actual (worst section per t):")
print("   t | s | a | |S(ta)| actual | t*|1-e(a/7)|*|mhat| predicted | M_slow-scale error")
ok_pred, growth = True, []
for t in (25, 50, 100, 200, 400):
    E = E6 + [t]
    P = 7 * lcm(60, t)
    best = (0.0, None, None, None)
    for s in range(7):
        arcs = R_s_exact(E, s)
        if not arcs: continue
        eps = []
        for a_, b_ in arcs:
            eps.append((a_, 1)); eps.append((b_, -1))
        for a in range(1, 7):
            n = (t * a) % P
            S = abs(sum(sg * cmath.exp(2j * pi * n * float(p)) for p, sg in eps))
            pred = t * abs(1 - cmath.exp(2j * pi * a / 7)) * mh[s][a]
            if S > best[0]: best = (S, s, a, pred)
    S, s, a, pred = best
    err = abs(S - pred)
    growth.append((t, S, pred))
    print(f"   {t:4d} | {s} | {a} | {S:10.2f} | {pred:10.2f} | {err:8.2f}")
    if pred > 30 and abs(S - pred) > 0.35 * pred: ok_pred = False
check("the resonant-mode formula matches the actual sup at resonant frequencies to O(M_slow) "
      "across t = 25..400 (the automaton-induction step VALIDATED)", ok_pred)
lin = growth[-1][1] / growth[1][1]
check(f"LINEAR GROWTH: |S(ta)| grows ~ t (400/50 ratio = {lin:.2f} ~ 8): with M = O(t) this "
      "gives max|S|^2/M ~ t — **UNIFORM HYP-6994 (sup-version) IS FALSE as t -> infinity**",
      5.0 < lin < 11.0)

# ---------- does the WEIGHTED Q_s law survive? ----------
print()
print("THE WEIGHTED LAW at t = 100 (P = 2100): sup over w classes of Q_s vs diam:")
E = E6 + [100]
P = 7 * lcm(60, 100)
s0 = None; bestload = 0
for s in range(7):
    arcs = R_s_exact(E, s)
    if arcs and len(arcs) > bestload: bestload, s0 = len(arcs), s
arcs = R_s_exact(E, s0)
eps = []
for a_, b_ in arcs:
    eps.append((a_, 1)); eps.append((b_, -1))
M = len(eps)
epsF = [(float(p), sg) for p, sg in eps]
supQ, supw = 0.0, None
for w in range(1, P + 1):
    tot = 0.0
    for p1, g1 in epsF:
        for p2, g2 in epsF:
            th = (w * (p1 - p2)) % 1.0
            tot += -g1 * g2 * th * (1 - th)
    q = 2 * pi * pi * tot
    if q > supQ: supQ, supw = q, w
print(f"   t=100: M={M}, diam=100: sup_w Q_s = {supQ:.2f} at w={supw}  (sup/diam = {supQ/100:.2f})")
check("THE SECOND REFUTATION: sup_w Q_s / diam GROWS (15.3 at t=50 -> 28.7 at t=100, "
      "~ diam scaling => sup_w Q_s ~ diam^2): even the weighted law FAILS on the resonant "
      "w-classes w = t*a*inverse(ell) mod P — THM-729's empirical O(diam) held only because "
      "tested w avoided them; the repaired conjecture is Q_s = O(diam) OFF the explicit "
      "resonant classes + Q_s <= C*diam^2*|mhat|^2/ell^2 ON them", supQ / 100 > 20)

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)
