#!/usr/bin/env python3
"""
lrc14_residual_bounded-speed-reduction_kps-S3-wf

ANGLE: "bounded-speed-reduction" for the residual case S3 of LRC(14).
(All set-minus written 'S minus v' in prose to avoid escape-sequence warnings.)

LRC(14) <=> M(S) >= 1/14 for every primitive covering 13-set S.
Case split (exhaustive, all PROVED elsewhere except S3):
  S1 (single-large, k<=1)            : PROVED
  S2 (clustered, Vmax < 13*Vmin)     : PROVED by LEMMA 1 (global witness arc J)
  S3 (Vmax >= 13*Vmin)               : THE ONLY OPEN CASE.

C(S) := [exists v in S: W(S-minus-v) > 1/(7v)]  =>  M(S) >= 1/14  (PROVED implication, THM-526).
Goal of this angle: turn S3 into a FINITE check by an EXPLICIT bounded-speed threshold V0*.

This script establishes / refutes the building blocks of that reduction, with EXACT
Fractions throughout.  Every numeric claim below is reproduced by running this file.

RESULTS (each tagged PROVED / VERIFIED / REFUTED):

 [R1] PROVED  - SCALE-INVARIANCE.  For any positive integer g, scaling S -> g*S maps
       G_{gS}(tau) = G_S(g*tau) (a g-to-1 cover), hence W(gA) = W(A)/g and Vmax(gS)=g Vmax(S),
       so the via-Vmax ratio  W(S-minus-Vmax)*7*Vmax  is INVARIANT under common scaling.  =>  WLOG gcd(S)=1.

 [R2] REFUTED - "W(S-minus-Vmax) has a universal positive lower bound" is FALSE.
       Counterexample family A={t,2t,...,12t}: W(A) = W({1..12})/t -> 0.  So a bounded-speed
       reduction CANNOT be obtained by lower-bounding W(S-minus-v) by a constant.  The via-v
       criterion is only scale-invariant as a RATIO W(S-minus-v)*7v, not as a width.

 [R3] REFUTED - The naive UNION-BOUND on any single window proving a whole-S witness FAILS.
       In the first gap of Vmin (width 6/(7 Vmin)) even "slow" H-runners (u < 13 Vmin) spin up
       to ~11 times, so the tooth count is large and sum of danger measures exceeds the window.
       Hence no closed-form measure proof of a one-window whole-S witness exists.

 [V1] VERIFIED - FIRST-GAP-OF-Vmin WHOLE-SET WITNESS.  In EVERY tested S3 set, the open arc
       I = (1/(14 Vmin), 13/(14 Vmin)) contains a point tau* SAFE FOR ALL 13 runners, giving
       M(S) >= 1/14 directly.  0 failures.  (This is LEMMA 1's arc, but we only require its
       intersection with G_{S-minus-Vmin} to be nonempty, not the full band-fit.)  NOT PROVED.

 [V2] VERIFIED - CRITERION C(S) (via best v) holds with ratio max_v W(S-minus-v)*7v >= 2.8 in the
       worst observed S3 set.  0 failures.  NOT PROVED.

 [PARTIAL/PROVED] - BOUNDED-PATTERN REDUCTION.  Combining [R1] + the dominant-large theorem
       (proved, mac-mini-S6) we PROVE: the S3 residual reduces to primitive covering 13-sets in
       which NOT all 12 non-max runners are <= 13 (i.e. some non-max runner is itself large).
       We quantify the explicit threshold and explain WHY the reduction does not close to a
       FINITE check: the residual family is infinite even mod scaling (e.g. {1,...,12,V} as V->inf
       is primitive, S3, and the 12 non-max runners ARE bounded -> handled; but {t,...,12t,V}
       with t,V both growing is primitive, S3, non-max runners UNbounded -> NOT a finite family).

CONCLUSION (honest): the bounded-speed-reduction CLOSES the sub-regime where the 12 non-max
runners are bounded (dominant-large, already proved).  For the genuine S3 residual it provides
NO finite reduction, because (a) widths are not bounded below [R2] and (b) no single-window
measure bound proves the whole-set witness [R3].  The first-gap-of-Vmin whole-set witness [V1]
is the cleanest TRUE statement and the right target for a future proof, but it is an irreducibly
Diophantine (three-distance / simultaneous) statement, not reachable by the bounded-speed angle.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

C = F(1, 14)

# ------------------------------------------------------------------ exact tools
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_for(tau, S):
    return all(min((u * tau) % 1, 1 - (u * tau) % 1) >= C for u in S)

def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def safe_arcs(A):
    """maximal safe arcs (lo,hi) on [0,1) for runner set A at level 1/14."""
    iv = []
    for v in set(A):
        hw = F(C, v); w = 2 * hw
        for j in range(v):
            c = F(j, v); a = (c - hw) % 1; b = a + w
            if b <= 1: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b - 1))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]: merged[-1] = (merged[-1][0], b)
        else: merged.append((a, b))
    arcs = []; prev = F(0)
    for a, b in merged:
        if a > prev: arcs.append((prev, a))
        prev = b
    if prev < 1: arcs.append((prev, F(1)))
    return arcs

def Wwidth(A):
    a = safe_arcs(A)
    if not a: return F(0)
    ws = [hi - lo for lo, hi in a]
    if a[0][0] == 0 and a[-1][1] == 1 and len(a) > 1:
        ws.append(a[0][1] + (1 - a[-1][0]))      # circular wrap
    return max(ws)

# exact M via THM-524 candidate set
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): Cc.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): Cc.add(F(k, d)); k += 1
    Cc.add(F(1, 2)); return Cc

def Mval(S):
    return max(min(nrm(x * t) for x in S) for t in cand(S))

def rand_cov_S3(rng, mx=200):
    S = set()
    order = [11, 13, 9, 8, 7, 5, 4, 3, 2, 12, 10, 6, 14]; rng.shuffle(order)
    for q in order:
        if any(x % q == 0 for x in S): continue
        S.add(q * rng.randint(1, mx // q))
    S = sorted(S)
    while len(S) < 13:
        S.append(rng.randint(2, mx)); S = sorted(set(S))
    return sorted(set(S))[:13]

# ====================================================================
def main():
    print("=" * 76)
    print("[R1] SCALE-INVARIANCE of the via-Vmax ratio  W(S\\{Vmax})*7*Vmax")
    print("=" * 76)
    rng = random.Random(2)
    allok = True
    for _ in range(8):
        base = sorted(random.sample(range(1, 40), 13))
        g = rng.choice([1, 2, 3, 5, 7])
        S = [g * x for x in base]
        A = [u for u in S if u != max(S)]; Ab = [u for u in base if u != max(base)]
        r1 = Wwidth(A) * 7 * max(S); r2 = Wwidth(Ab) * 7 * max(base)
        allok &= (r1 == r2)
    print(f"  via-Vmax ratio invariant under common scaling on 8 random pairs: {allok}")
    print("  PROVED: G_{gS}(tau)=G_S(g tau) => W(gA)=W(A)/g, Vmax scales by g => ratio fixed.")
    print("  CONSEQUENCE: WLOG gcd(S)=1 (primitive). [But primitivity does NOT bound speeds.]")

    print("\n" + "=" * 76)
    print("[R2] REFUTED: W(S\\{Vmax}) is NOT bounded below (kills the naive reduction).")
    print("=" * 76)
    for t in [1, 5, 20, 50]:
        A = [k * t for k in range(1, 13)]
        W = Wwidth(A)
        print(f"  A=t*{{1..12}}, t={t:3d}: W(A)={float(W):.6f}  W(A)*t={float(W*t):.6f} (constant)")
    print("  => W(A)=W({1..12})/t -> 0.  No constant lower bound on W exists.")

    print("\n" + "=" * 76)
    print("[R3] REFUTED: one-window union bound for a whole-S witness.")
    print("=" * 76)
    print("  Window I = first gap of Vmin = (1/(14 Vmin), 13/(14 Vmin)), |I| = 6/(7 Vmin).")
    print("  A 'slow' H-runner u<13 Vmin spins u*|I| = 6u/(7Vmin) < 78/7 ~ 11.1 times in I.")
    print("  So even slow runners contribute ~11 wide teeth; sum of danger measures > |I|.")
    print("  Hence NO closed-form union/measure bound proves a one-window whole-S witness.")
    # demonstrate the failure numerically on the worst sampled set
    Sdemo = [25, 90, 99, 143, 172, 192, 195, 242, 266, 280, 294, 310, 348]
    Vmin = min(Sdemo); H = [u for u in Sdemo if u != Vmin]
    I = F(6, 7 * Vmin)
    ub = sum((F(6 * u, 7 * Vmin) + 1) * F(1, 7 * u) for u in H)  # union bound on danger meas in I
    print(f"  demo S={Sdemo}: |I|={float(I):.5f}, union-bound danger meas in I = {float(ub):.3f}"
          f"  ({'exceeds' if ub > I else 'below'} |I|).")

    print("\n" + "=" * 76)
    print("[V1] VERIFIED: first-gap-of-Vmin contains a WHOLE-S safe point (M(S)>=1/14 directly).")
    print("=" * 76)
    rng = random.Random(41); tested = 0; fail = 0; worst = F(1); worstS = None
    for _ in range(3000):
        S = rand_cov_S3(rng)
        if len(S) != 13 or not covering(S): continue
        if max(S) < 13 * min(S): continue           # skip S2
        tested += 1
        Vmin = min(S); H = [u for u in S if u != Vmin]
        lo = F(1, 14 * Vmin); hi = F(13, 14 * Vmin)
        tot = F(0); pt = None
        for a, b in safe_arcs(H):
            l = max(a, lo); r = min(b, hi)
            if r > l:
                tot += r - l
                if pt is None: pt = (l + r) / 2
        if tot == 0:
            fail += 1; continue
        if not safe_for(pt, S):
            print("  BUG: witness not whole-S-safe", S, pt)
        if tot < worst: worst = tot; worstS = S
    print(f"  S3 sets tested: {tested};  first-gap-of-Vmin witness FAILURES: {fail}")
    print(f"  worst overlap meas( I cap G_H ) = {float(worst):.6f} (>0); worst S = {worstS}")
    print("  STATUS: VERIFIED (0 failures). NOT PROVED (irreducibly Diophantine).")

    print("\n" + "=" * 76)
    print("[V2] VERIFIED: criterion C(S) via best v, worst-case ratio over S3.")
    print("=" * 76)
    rng = random.Random(61); tested = 0; worst = F(99); worstS = None; cfail = 0
    for _ in range(500):
        S = rand_cov_S3(rng)
        if len(S) != 13 or not covering(S): continue
        if max(S) < 13 * min(S): continue
        tested += 1
        best = F(0)
        for v in S:
            A = [u for u in S if u != v]
            r = Wwidth(A) * 7 * v
            if r > best: best = r
        if best <= 1: cfail += 1
        if best < worst: worst = best; worstS = S
    print(f"  S3 tested: {tested};  C(S) (best v) FAILURES (ratio<=1): {cfail}")
    print(f"  worst best-v ratio max_v W(S\\v)*7v = {float(worst):.4f} (C holds iff >1)")
    print(f"  at S = {worstS}")

    print("\n" + "=" * 76)
    print("[PARTIAL] BOUNDED-PATTERN reduction: what the angle DOES and does NOT close.")
    print("=" * 76)
    # dominant-large worst-core width over bounded non-max pool B (reproduce the proved constant)
    for B in [13, 16]:
        mn = F(1); arg = None
        for A in combinations(range(1, B + 1), 12):
            W = Wwidth(list(A))
            if 0 < W < mn: mn = W; arg = A
        print(f"  B={B}: w_min(B)={mn}={float(mn):.6f} at {arg};  via-Vmax discharges all V>{F(1,7)/mn}.")
    print("  PROVED sub-regime: non-max runners all <=13 and V>=53 => C(S) via v=Vmax. (mac-mini-S6)")
    print("  RESIDUAL after that: some non-max runner is itself large. By [R2] its width is")
    print("  unbounded-small, and {t,2t,...,12t,V} (t,V both free) is primitive+S3 with UNbounded")
    print("  non-max runners => the residual is INFINITE even mod scaling. NO finite check exists")
    print("  along this angle.")

    print("\n" + "=" * 76)
    print("CROSS-CHECK: exact M(S) >= 1/14 on the sampled S3 sets (sanity).")
    print("=" * 76)
    rng = random.Random(77); tested = 0; mfail = 0; minM = F(1)
    for _ in range(400):
        S = rand_cov_S3(rng, mx=120)
        if len(S) != 13 or not covering(S): continue
        if max(S) < 13 * min(S): continue
        tested += 1
        m = Mval(S)
        if m < C: mfail += 1
        if m < minM: minM = m
    print(f"  exact-M S3 sets tested: {tested};  M(S) < 1/14 count: {mfail};  min M = "
          f"{minM} = {float(minM):.5f}  (threshold 1/14 = {float(C):.5f})")


if __name__ == "__main__":
    main()
