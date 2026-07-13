#!/usr/bin/env python3
"""
lrc14_multirunner_argA_exact_opus_S270.py
=========================================
opus-2026-07-13-S270.  Answer to kps-S127 cont.70's ASK (HYP-6248): the crack body
{1..11,13,84} falls between the runner-1 lemma's Argument A (measure: |S_rest cap D_1| <= 1/14
via the smallest even speed alone -- FAILS since |S_rest| = 0.0666 < 1/14 = 0.0714) and
Argument B (equidistribution: spread rest -- STRESSED by the 9-consecutive run 2..11).
kps's proposed fix: MULTI-RUNNER Argument A -- carve D_1 by SEVERAL small runners:

    |S_rest cap D_1|  <=  |D_1 cap S_{w1} cap ... cap S_{wk}|   for any {w1..wk} subset rest,

and the right-hand side is EXACTLY COMPUTABLE in rational arithmetic (finitely many rational
intervals).  Multi-runner Arg A closes the body as soon as some prefix P satisfies

    |D_1 cap S_P|  <  |S_rest|          (then S_rest cannot be a subset of D_1  =>  LRC holds).

This script computes, in EXACT Fraction arithmetic (no grids, no floats in the logic):
  1. |S_rest| exactly (rest = {2..11,13,84});
  2. the carving curve |D_1 cap S_2 cap ... cap S_j| as small runners are added (ascending
     order, plus a greedy order), each value an exact rational;
  3. the minimal prefix that crosses below |S_rest| -- the certificate kps asked for;
  4. kps's specific suggestion |S_2 cap S_4 cap D_1|;
  5. the surplus L = |S_rest| - |S_rest cap D_1| (must equal the body's true good-set measure;
     cross-checks klein's grid value L ~ 0.0054);
  6. the same carving on the DEEP WELL {1..12,182} core peel for contrast, and the closed-form
     pattern of |D_1 cap S_{2..j}| (the multi-runner threshold curve replacing (s-1)/(7s)).

PERSPECTIVE NOTE (owner's frame, this session): peeling runner 1 and carving D_1 by the rest is
the CORE-endpoint of the same peel family whose FAR-endpoint is klein's THM-731/732 (peel the
far element, certify on its fine grid).  One identity family, 13 peels = 13 perspectives; the
runner-1 lemma and the far-element certificate are its two ends.  See the S270 reflection.
"""
from fractions import Fraction as F

THR = F(1, 14)

# ---------------- exact interval engine on [0,1) ----------------

def normalize(ivs):
    """sort + merge overlapping/touching intervals"""
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return [(a, b) for a, b in out]

def intersect(A, B):
    """intersection of two normalized interval lists"""
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo:
            out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def measure(A):
    return sum(b - a for a, b in A)

def safe_set(w):
    """S_w = {t in [0,1): ||w t|| >= 1/14} as exact intervals"""
    return [((k + THR) / w, (k + 1 - THR) / w) for k in range(w)]

def danger_set(w):
    """D_w = {t: ||w t|| < 1/14}"""
    out = []
    for k in range(w + 1):
        a, b = F(k, w) - THR / w, F(k, w) + THR / w
        out.append((max(a, F(0)), min(b, F(1))))
    return normalize(out)

def inter_all(sets):
    cur = [(F(0), F(1))]
    for s in sets:
        cur = intersect(cur, s)
    return cur

# ---------------- the crack body ----------------

BODY = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 84]
REST = [w for w in BODY if w != 1]

print("=" * 96)
print("MULTI-RUNNER ARGUMENT A -- exact rational carving of D_1 on the crack body {1..11,13,84}")
print("(kps-S127 cont.70 ask, HYP-6248 -> HYP-6505).  All values exact Fractions; floats for display.")
print("=" * 96)

S_rest_ivs = inter_all(safe_set(w) for w in REST)
S_rest = measure(S_rest_ivs)
print("\n|S_rest| (rest = {2..11,13,84})  =  %s  =  %.6f    (kps grid 2e6 gave 0.0666)"
      % (S_rest, float(S_rest)))
print("   S_rest has %d intervals" % len(S_rest_ivs))

D1 = danger_set(1)
S_rest_cap_D1 = intersect(S_rest_ivs, D1)
m_full = measure(S_rest_cap_D1)
L = S_rest - m_full
print("|S_rest cap D_1|                 =  %s  =  %.6f" % (m_full, float(m_full)))
print("L = |S_rest| - |S_rest cap D_1|  =  %s  =  %.6f    (true good-set measure; klein ~0.0054)"
      % (L, float(L)))
print("Arg A (single, s=2) threshold 1/14 = %.6f  vs |S_rest| = %.6f  -> Arg A FAILS (confirmed)"
      % (float(F(1, 14)), float(S_rest)))

# ---------------- carving curve, ascending order ----------------

print("\n" + "-" * 96)
print("CARVING CURVE (ascending order): m_j = |D_1 cap S_2 cap ... cap S_j|  vs target |S_rest|")
print("-" * 96)
print("%8s %28s %12s %10s  %s" % ("add w", "m_j (exact)", "m_j", "target", "status"))
cur = D1
added = []
crossed_at = None
for w in REST:
    cur = intersect(cur, safe_set(w))
    added.append(w)
    m = measure(cur)
    status = "BELOW target -- multi-runner Arg A CLOSES" if m < S_rest else ""
    if m < S_rest and crossed_at is None:
        crossed_at = list(added)
    print("%8d %28s %12.6f %10.6f  %s" % (w, str(m), float(m), float(S_rest), status))

print("\nMinimal ascending prefix that closes: P = {%s}" %
      (", ".join(map(str, crossed_at)) if crossed_at else "NONE -- full rest needed"))

# ---------------- kps's specific suggestion ----------------

m24 = measure(inter_all([D1, safe_set(2), safe_set(4)]))
print("\nkps's suggested pair {2,4}:  |D_1 cap S_2 cap S_4| = %s = %.6f  (target %.6f)  -> %s"
      % (m24, float(m24), float(S_rest), "closes" if m24 < S_rest else "NOT enough"))

# ---------------- greedy order ----------------

print("\n" + "-" * 96)
print("GREEDY CARVING (each step: add the rest-speed that minimizes the measure)")
print("-" * 96)
cur = D1
pool = list(REST)
greedy_seq = []
while pool:
    best = min(pool, key=lambda w: measure(intersect(cur, safe_set(w))))
    cur = intersect(cur, safe_set(best))
    pool.remove(best)
    greedy_seq.append(best)
    m = measure(cur)
    mark = "  <-- BELOW target" if m < S_rest else ""
    print("   add %3d  ->  m = %-24s = %.6f%s" % (best, str(m), float(m), mark))
    if m < S_rest:
        break
print("Greedy closing set: {%s}" % ", ".join(map(str, greedy_seq)))

# ---------------- the closed-form threshold curve ----------------

print("\n" + "-" * 96)
print("THE MULTI-RUNNER THRESHOLD CURVE  T_j = |D_1 cap S_2 cap ... cap S_j|  (consecutive 2..j),")
print("the replacement for Arg A's single-runner (s-1)/(7s):")
print("-" * 96)
cur = D1
for j in range(2, 15):
    cur = intersect(cur, safe_set(j))
    m = measure(cur)
    print("   j = %2d :  T_j = %-22s = %.6f   (1/14 = %.6f)" % (j, str(m), float(m), float(F(1,14))))

# ---------------- contrast: deep well core peel ----------------

print("\n" + "-" * 96)
print("CONTRAST -- deep well {1..12,182}, same core peel (rest = {2..12,182}):")
REST_DW = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182]
S_rest_dw_ivs = inter_all(safe_set(w) for w in REST_DW)
S_rest_dw = measure(S_rest_dw_ivs)
m_dw = measure(intersect(S_rest_dw_ivs, D1))
print("   |S_rest| = %s = %.6f ;  |S_rest cap D_1| = %.6f ;  L = %.6f"
      % (S_rest_dw, float(S_rest_dw), float(m_dw), float(S_rest_dw - m_dw)))
cur = D1
for w in REST_DW:
    cur = intersect(cur, safe_set(w))
    m = measure(cur)
    if m < S_rest_dw:
        print("   ascending carving crosses below target after adding w = %d  (m = %.6f < %.6f)"
              % (w, float(m), float(S_rest_dw)))
        break
else:
    print("   ascending carving NEVER crosses (needs full rest)")

print("\ndone.")
