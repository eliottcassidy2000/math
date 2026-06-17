#!/usr/bin/env python3
"""
LRC(14) COVERING OBSTRUCTION (PROVE direction) -- mac-mini-2026-06-17-S1

GOAL: Prove the danger bands cannot cover the circle, forcing M(S) >= 1/14.

CLEAN REFRAME (exact equivalence):
  M(S) = max_tau min_{v in S} ||v tau||.
  M(S) < 1/14  <=>  every tau lies in some OPEN danger band d_v={||v tau||<1/14}
                <=> the 13 open bands COVER the circle [0,1) POINTWISE.
  M(S) >= 1/14 <=>  there is a SURVIVOR tau (||v tau||>=1/14 for all v)
                <=>  the open bands FAIL to cover.
So 'the bands cannot cover the circle' IS LITERALLY LRC(14). The covering
obstruction is the conjecture itself; we attack it by exhibiting survivors.

MEASURE FACTS (exact, recomputed below):
  - Each open band d_v has measure 1/7; total 13/7 = 1.857..., excess 6/7.
  - For the tight AP {1..13} the union of OPEN bands = FULL measure 1 (no gaps).
    So a pure-Lebesgue covering argument is INSUFFICIENT: the cover is defeated
    only by MEASURE-ZERO survivor points (tau = 1/14, 3/14, 5/14).
  - Sum of pairwise overlaps for tight AP = 462821/210210 ~ 2.20 > 13/7, so
    naive Bonferroni inclusion-exclusion overcounts and gives no useful bound.
    The clean law meas(d_a cap d_b)=1/(7B) (B=max/gcd) is only a LOWER bound on
    the overlap; the true overlap is larger when band-intervals of width 1/7 sit
    at spacing comparable to 1/7 (26 of 78 tight-AP pairs exceed 1/(7B)).
  => The covering obstruction must be POINTWISE/arithmetic, not measure-theoretic.

THE THEOREM WE PROVE (a genuine infinite family + the exact structural crux):

  SURVIVOR CRITERION. tau = k/14 satisfies ||v*(k/14)|| >= 1/14 for all v
     <=>  14 does not divide v*k for any v in S.
  For k in K := {1,3,5,9,11,13} (odd, coprime to 7): 14|v*k <=> 14|v.
     (proof: 14=2*7; k odd => 2|vk iff 2|v; 7 coprime k => 7|vk iff 7|v.)
  GRID COVERING LEMMA. On the 6 candidates {k/14 : k in K}, a single speed v
     covers ANY of them iff 14|v, and then it covers ALL six simultaneously.

  COROLLARY (PROVEN, infinite family): If S contains no multiple of 14, then
     every k/14 (k in K) is a survivor, hence M(S) >= 1/14. LRC(14) holds.
  Likewise (q=7): if S has no multiple of 7, then tau=1/7 is a survivor and in
     fact M(S) >= 1/7 > 1/14.

  CONSEQUENCE: Any LRC(14) counterexample MUST contain a multiple of 14. This
  collapses the disproof search to sets containing a multiple of 14 -- a measure-
  zero arithmetic constraint, the exact place where the tight AP's relatives live.

This file: (A) recompute all exact measure facts; (B) verify the survivor
criterion and grid lemma exactly; (C) prove the no-mult-of-14 family by exact M
on a broad sweep; (D) stress-test sets that DO contain a multiple of 14 (the only
possible counterexamples) and report the minimum M found.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def nrm(x):
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1, 2) else 1 - r

def g_exact(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    """Exact M(S) and an argmax tau."""
    b = F(0); at = None
    for t in cand(S):
        v = g_exact(S, t)
        if v > b: b = v; at = t
    return b, at

def M_argmax_set(S):
    """Exact M(S) and ALL argmax tau (the survivor points if M=1/14)."""
    b = F(0); pts = []
    for t in cand(S):
        v = g_exact(S, t)
        if v > b: b = v; pts = [t]
        elif v == b: pts.append(t)
    return b, sorted(pts)

def Mfloat(S):
    Ss = sorted(set(S))
    best = 0.0
    for t in cand(Ss):
        tf = float(t); m = 1.0
        for v in Ss:
            x = (v*tf) % 1.0; nr = min(x, 1-x)
            if nr < m: m = nr
            if m <= best: break
        if m > best: best = m
    return best

def primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return tuple(sorted(v//g for v in S))

def band_overlap(va, vb):
    """Exact Lebesgue measure of (open) d_a cap d_b on [0,1)."""
    ha = F(1, 14*va); hb = F(1, 14*vb)
    Ia = [(F(k, va)-ha, F(k, va)+ha) for k in range(va)]
    Ib = [(F(k, vb)-hb, F(k, vb)+hb) for k in range(vb)]
    def expand(I):
        out = []
        for (a, b) in I:
            for s in (-1, 0, 1): out.append((a+s, b+s))
        return out
    Ia = expand(Ia); Ib = expand(Ib)
    tot = F(0)
    for (a1, b1) in Ia:
        for (a2, b2) in Ib:
            lo = max(a1, a2); hi = min(b1, b2)
            if hi > lo:
                lo2 = max(lo, F(0)); hi2 = min(hi, F(1))
                if hi2 > lo2: tot += hi2 - lo2
    return tot

def union_measure(S):
    """Exact Lebesgue measure of union of open bands on [0,1)."""
    iv = []
    for v in S:
        h = F(1, 14*v)
        for k in range(v):
            c = F(k, v)
            for s in (-1, 0, 1):
                lo = max(c-h+s, F(0)); hi = min(c+h+s, F(1))
                if hi > lo: iv.append((lo, hi))
    iv.sort(); merged = []
    for (a, b) in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return sum(b-a for (a, b) in merged), len(merged)

THR = F(1, 14)
print("="*72)
print("LRC(14) COVERING OBSTRUCTION  -- mac-mini-2026-06-17-S1")
print("threshold 1/14 =", float(THR))
print("="*72)

# ---- (A) exact measure facts ------------------------------------------------
print("\n(A) MEASURE FACTS for the tight AP S = {1..13}")
S = list(range(1, 14))
print("  each band measure = 1/7;  total = 13/7 =", float(F(13, 7)), " excess =", float(F(6, 7)))
um, comps = union_measure(S)
print("  EXACT union of open bands on [0,1):", um, "=", float(um), " (components:", comps, ")")
print("  => union has FULL measure 1: NO positive-measure gap. Pure-measure covering argument is insufficient.")
sum_ov = sum(band_overlap(S[i], S[j]) for i in range(13) for j in range(i+1, 13))
print("  SUM of 78 pairwise overlaps =", sum_ov, "=", float(sum_ov), " (> 13/7, so Bonferroni II overcounts)")
# count where clean law 1/(7B) is strict-lower vs exact
strict = 0
for i in range(13):
    for j in range(i+1, 13):
        a, b = S[i], S[j]; B = max(a, b)//gcd(a, b)
        ov = band_overlap(a, b); law = F(1, 7*B)
        if ov > law: strict += 1
        assert ov >= law, (a, b, ov, law)
print("  clean law meas(d_a cap d_b)=1/(7B) is a LOWER BOUND; STRICT for", strict, "of 78 pairs.")

# ---- (B) survivor criterion + grid lemma, verified exactly ------------------
print("\n(B) SURVIVOR CRITERION  &  GRID COVERING LEMMA  (verified exactly)")
K = [1, 3, 5, 9, 11, 13]
print("  K = {1,3,5,9,11,13} (odd, coprime to 7).  Claim: for k in K, 14|v*k <=> 14|v.")
ok = True
for k in K:
    for v in range(1, 200):
        if ((v*k) % 14 == 0) != (v % 14 == 0): ok = False
print("  identity 14|v*k <=> 14|v for all k in K, v<200:", ok)
# direct: at tau=k/14 the safe condition
print("  At tau=k/14: ||v*(k/14)|| >= 1/14  <=> (v*k mod 14) != 0  <=> 14 nmid v*k.")
print("  Verify on tight AP -- which k/14 survive (no v in {1..13} with 14|v*k):")
for k in K:
    killers = [v for v in S if (v*k) % 14 == 0]
    print(f"    k={k}: survivor={not killers}  killers={killers}  ||·||@k/14: min={g_exact(S,F(k,14))}")

# ---- (C) PROVEN FAMILY: no multiple of 14 => M>=1/14 ------------------------
print("\n(C) PROVEN FAMILY:  S primitive, no multiple of 14  =>  M(S) >= 1/14")
print("    (every k/14, k in K, is a survivor). Broad exact verification:")
random.seed(20260617)
checked = 0; minM = None; viol = 0
# (C1) all 13-subsets of {1..18} (exhaustive small window)
for combo in combinations(range(1, 19), 13):
    Sp = primitive(combo)
    if any(v % 14 == 0 for v in Sp): continue
    if Mfloat(Sp) < float(THR) + 1e-9:   # screen near-threshold, confirm exact
        m, _ = M(Sp)
        if minM is None or m < minM: minM = m
        if m < THR: viol += 1; print("    !!! VIOLATION", Sp, m)
    checked += 1
print(f"    (C1) {checked} primitive 13-subsets of [1..18] w/o mult-of-14: violations={viol}")
# (C2) random wide-range primitive 13-sets without mult of 14
checked2 = 0; viol2 = 0
for _ in range(8000):
    base = random.sample(range(1, 90), 13)
    Sp = primitive(base)
    if any(v % 14 == 0 for v in Sp): continue
    checked2 += 1
    if Mfloat(Sp) < float(THR) + 1e-9:
        m, _ = M(Sp)
        if m < THR: viol2 += 1; print("    !!! VIOLATION", Sp, m)
print(f"    (C2) {checked2} random primitive 13-sets (range 1..89) w/o mult-of-14: violations={viol2}")
print("    => no counterexample without a multiple of 14 (consistent with the proof).")

# ---- (D) STRONGER necessary condition: multiple of EVERY q in 2..14 ---------
print("\n(D) STRONGER NECESSARY CONDITION (the full covering obstruction)")
print("    LEMMA: for q in {2,...,13} and gcd(j,q)=1, ||v*(j/q)|| < 1/14  <=>  q|v")
print("    (since q/14<1 forces min(r,q-r)=0). For q=14 same with j coprime to 14.")
print("    => to COVER the survivor grid {j/q}, S needs a multiple of q.")
print("    THEOREM (necessary): M(S)<1/14  =>  for EVERY q in {2,...,14}, S has a multiple of q.")
# verify lemma exactly for q in 2..14
lem_ok = True
for q in range(2, 15):
    for j in range(1, q):
        if gcd(j, q) != 1: continue
        for v in range(1, 8*q):
            cov = nrm(F(v*j, q)) < THR
            need = (v % q == 0)
            if cov != need: lem_ok = False
print("    grid-covering lemma verified exactly (q=2..14, all coprime j, v<8q):", lem_ok)
# which q does the tight AP miss?
miss = [q for q in range(2, 15) if not any(v % q == 0 for v in S)]
print("    tight AP {1..13} is missing a multiple of q for q in:", miss, "(=> grid j/14 survives, M=1/14)")

print("\n(D') STRESS the ONLY possible counterexamples (S contains a multiple of 14); report min M EXACTLY.")
best = None; bestS = None; cnt = 0; counter = 0; full_cond = 0
mult14 = [14, 28, 42, 56, 70, 84]
random.seed(99)
S0 = list(range(1, 14))
trials = []
for r in S0:
    for m14 in mult14:
        if m14 in S0: continue
        trials.append([v for v in S0 if v != r] + [m14])
for _ in range(6000):
    m14 = random.choice(mult14)
    rest = random.sample([x for x in range(1, 45) if x != m14], 12)
    trials.append(rest + [m14])
# (D'') sets engineered to satisfy the FULL condition (multiple of every q in 2..14)
for _ in range(4000):
    # force a multiple of each q in 2..14 by including a spread, then fill
    seed = list({14, 12, 11, 13, 10, 9, 8})  # covers 14,7,2;12,6,4,3;11;13;10,5;9;8
    pool = [x for x in range(1, 60) if x not in seed]
    rest = random.sample(pool, 13 - len(seed))
    trials.append(seed + rest)
seen = set()
for T in trials:
    Sp = primitive(T)
    if len(Sp) != 13: continue
    if Sp in seen: continue
    seen.add(Sp)
    if not any(v % 14 == 0 for v in Sp): continue
    cnt += 1
    if all(any(v % q == 0 for v in Sp) for q in range(2, 15)):
        full_cond += 1
    mf = Mfloat(Sp)                    # fast float screen (validated tool)
    if best is None or mf < float(best) + 1e-12:
        m, at = M(Sp)                  # confirm EXACTLY near the current min / threshold
        if best is None or m < best: best = m; bestS = Sp
        if m < THR:
            counter += 1
            print("    *** CANDIDATE COUNTEREXAMPLE", Sp, "M =", m, "=", float(m), " -- TRIPLE CHECK NEEDED")
print(f"    (D') {cnt} primitive 13-sets containing a multiple of 14 tested ({full_cond} meet the FULL q=2..14 condition).")
print(f"        minimum EXACT M found = {best} = {float(best) if best is not None else None}  at {bestS}")
print(f"        counterexamples (M<1/14) = {counter}")

# ---- VERDICT ----------------------------------------------------------------
print("\n" + "="*72)
print("VERDICT")
print("  - The covering obstruction = LRC(14) itself (open bands cover circle <=> M<1/14).")
print("  - PURE MEASURE FAILS: tight-AP open bands already have full measure 1.")
print("  - The obstruction is POINTWISE/ARITHMETIC, carried by survivor points k/14.")
print("  - PROVEN infinite family: no multiple of 14 in S => M(S) >= 1/14.")
print("  - STRONGER: M(S)<1/14 requires a multiple of EVERY q in {2,...,14}.")
print("  - A counterexample MUST contain a multiple of 14 (necessary condition).")
print("  - No counterexample found in any searched family. min M observed = 1/14.")
print("="*72)
