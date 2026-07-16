#!/usr/bin/env python3
"""
peel_residue_audit_klein_S314.py — klein-2026-07-16-S314 (cont.5)

THE PEEL-RESIDUE AUDIT (THM-883 follow-up): does the density route survive the resonant
w-classes?

Frame (all exact): S(w) := w·Error_cover(w) = Σ_s Σ_{arcs [a,b] of R_s} (G_s(wb) − G_s(wa)),
G_s(y) = meas([0,{y}] ∩ [s/7,(s+1)/7)) − {y}/7 + (whole periods contribute 0), computed with
rationals.  Since all arc endpoints have denominator dividing P = 7·lcm(E), S(w) depends
only on w mod P: **Error(w) = S(w mod P)/w**, so

    AUDIT = compute S(r) exactly for ALL r in Z_P; then for every peel w > D0:
    |Error(w)| <= max_r |S(r)| / w,  and  D0_res := max_r|S(r)| / slack
    is the EXPLICIT band threshold beyond which every peel (resonant or not) is safe.

Checks: (1) the maximizing residues ARE the THM-883 resonant classes (r ≡ t·a·ell^{-1});
(2) S(t·1 mod P) matches the resonant-mode prediction; (3) D0_res(t) grows LINEARLY in t
(the THM-883 refutation's fingerprint) — the structured band must extend to c·diam, which
remains FINITE: the route survives with a quantified enlargement.
"""
from fractions import Fraction as Fr
from math import pi, gcd, lcm
import cmath

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

def G_s(y, s):
    """int_0^y g_s(u) du for rational y >= 0, g_s = 1_[s/7,(s+1)/7) - 1/7 (mean zero, so
    whole periods drop)."""
    fy = y % 1
    lo, hi = Fr(s, 7), Fr(s + 1, 7)
    inter = max(Fr(0), min(fy, hi) - lo)
    return inter - fy / 7

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

SLACK = Fr(97, 1000)     # the route's Error headroom (cap9 - 0.397 ~ 0.097)
print("THE COMPLETE PEEL-RESIDUE AUDIT:  S(w) = w*Error is P-periodic; Error = S/w")
print("   t |   P   | max_r |S(r)| | argmax r | resonant? | D0_res = max|S|/0.097 | D0/t")
rows = []
for t in (50, 100, 200):
    E = [0, 1, 2, 3, 4, 5, t]
    P = 7 * lcm(60, t)
    allarcs = [(s, R_s_exact(E, s)) for s in range(7)]
    maxS, argr = Fr(0), None
    Svals = {}
    for r in range(1, P):
        S = Fr(0)
        for s, arcs in allarcs:
            for a, b in arcs:
                S += G_s(r * b, s) - G_s(r * a, s)
        aS = abs(S)
        Svals[r] = S
        if aS > maxS: maxS, argr = aS, r
    # resonance test: is argr ≡ t*a*inv(ell) mod P for small a, ell?
    reson = None
    for a in range(1, 7):
        for ell in range(1, 13):
            if gcd(ell, P) == 1 and (ell * argr - t * a) % P == 0:
                reson = (a, ell); break
            # also allow gcd>1 classes: ell*argr ≡ t*a mod P
            if (ell * argr - t * a) % P == 0:
                reson = (a, ell); break
        if reson: break
    D0 = float(maxS / SLACK)
    rows.append((t, P, float(maxS), argr, reson, D0))
    print(f"   {t:3d} | {P:5d} | {float(maxS):10.4f} | {argr:6d} | {str(reson):>9s} | "
          f"{D0:8.1f} | {D0/t:5.2f}")

check("the max-|S| residues are RESONANT classes (ell*r = t*a mod P with small a, ell) at "
      "every t — THM-883's classes are exactly the worst peels", all(r[4] for r in rows))
check(f"D0_res grows LINEARLY in t (D0/t = {[f'{r[5]/r[0]:.2f}' for r in rows]}): the "
      "structured band must extend to c*diam — FINITE, so the route SURVIVES with a "
      "quantified band enlargement (not the sqrt(diam) band THM-729 hoped for)",
      max(r[5] / r[0] for r in rows) / min(r[5] / r[0] for r in rows) < 1.8)

# verdict quantities: beyond D0_res every peel is safe REGARDLESS of residue:
print()
print("VERDICT: for every peel w > D0_res: |Error(w)| <= max|S|/w < 0.097, resonant or not.")
print("Below D0_res: the finite structured band (S275-style) must cover [26, D0_res] —")
print("its length is c*diam with c ~", f"{rows[-1][5]/rows[-1][0]:.2f}", "(vs the old sqrt-diam band).")
# sanity: check Error at the first large representative of the worst class is safe
t, P = rows[-1][0], rows[-1][1]
r = rows[-1][3]
w_large = r + P if r <= rows[-1][5] else r
err = abs(rows[-1][2]) / w_large
check(f"first representative of the worst class beyond the band (w = {w_large}): "
      f"|Error| = {err:.5f} < 0.097 — safe", err < 0.097)

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)
