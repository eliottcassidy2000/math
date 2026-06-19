#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asm_verify_widespreadrigo_kps-S10-wf.py   (kind-pasteur-2026-06-19-S10-wf)

ADVERSARIAL VERIFICATION of the claimed assembly piece "wide-spread-rigorous".
We re-derive EVERY constant exactly, then stress the central rigor claims:

  CLAIM (paraphrase): for primitive E, |E|=k in {8,9,10}, span(E) > B(k) implies
     meas(S7(E)) <= cap_k, via:
       (Part 3) a "core-extraction key lemma" forcing a large coordinate (decay 1/e_L),
       (Part 4) corr(E) <= delta_k on canonical wide families (12/12),
       (Part 6) cluster-collapse for no-scale-separation,
       (Part 7) an exhaustive finite check to span<=B(k) (B(8)=16,B(9)=13,B(10)=12),
     MODULO one "uniform Minkowski constant" eps(B) (admitted open).

VERIFICATION TASKS (per the verifier checklist):
 (V1) Re-derive M7(k), consec, delta_k, caps EXACTLY; confirm the stated rational values.
 (V2) Pin the Lemma-B envelope c1 exactly; confirm 0.697303, attained at n=66.
 (V3) The CENTRAL QUESTION: is the span>B(k) side a real proof or "->0"?  We hunt
      AGGRESSIVELY for an over-cap OR over-consec primitive shape with span > B(k):
      resonant (multiples of 7), short-relation (small additive relations), 1- and 2-
      stranger, no-scale-separation, shifted-AP, dilated-consec.  EXACT meas.
 (V4) Check the finite-check enumeration: is B(k) genuinely exhaustive?  Is the claimed
      exceedance count (0) correct?  Is B(k) consistent across the two companion scripts?
 (V5) Audit the "KEY LEMMA": does the decay 1/e_L actually bound the TYPE-(II) tail, or
      does it (as the script itself admits) diverge harmonically?  We test whether the
      forced-coordinate magnitude really controls |corr| as span grows, and whether the
      single-large-offset case is the only one (the >=2-large case is NOT covered).
 (V6) Audit cluster-collapse: the Weyl error is "O(w/M)" with NO closed M-threshold.
      Confirm the limit value and that it is < cap; test whether a FINITE-M cluster can
      exceed cap before the asymptote kicks in (a real danger the proof must exclude).

ALL EXACT rationals for measures.  We MARK each finding.  Default skeptical.
"""
from __future__ import annotations
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb, sin, pi
from functools import reduce
import cmath

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

TWO_PI_I = 2j * pi


def lcm(a, b):
    return a * b // gcd(a, b) if a and b else (a or b)


def gcd_all(E):
    return reduce(gcd, [e for e in E if e != 0], 0)


def primitive(E):
    return gcd_all(E) == 1


def banner(t):
    print("\n" + "=" * 90 + f"\n{t}\n" + "=" * 90)


# ---------------------------------------------------------------------------
# EXACT meas(S7).  Two independent engines, cross-validated.
# ---------------------------------------------------------------------------
def measS7_fast(E):
    E = sorted(set(int(e) for e in E))
    Enz = [e for e in E if e != 0]
    if not Enz:
        return F(0)
    D = 7 * reduce(lcm, Enz, 1)
    bk = set()
    bk.add(0); bk.add(D)
    for e in Enz:
        step = D // (7 * e)
        kk = 0
        while kk <= D:
            bk.add(kk)
            kk += step
    bk = sorted(bk)
    total = F(0)
    for i in range(len(bk) - 1):
        k0, k1 = bk[i], bk[i + 1]
        if k1 <= k0:
            continue
        num = k0 + k1
        den = 2 * D
        res = {0}
        for e in Enz:
            res.add((7 * e * num) // den % 7)
        if len(res) == 7:
            total += F(k1 - k0, D)
    return total


def measS7_slow(E):
    E = sorted(set(int(e) for e in E))
    Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        res = set(int(7 * e * xm) % 7 for e in E)
        if len(res) == 7:
            total += (x1 - x0)
    return total


def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
       11: F(66, 91), 12: F(6, 7), 13: F(1)}


# ===========================================================================
def v0_engine_crossvalidate():
    banner("V0 -- cross-validate the two EXACT meas engines (no engine bug)")
    random.seed(11)
    tests = [(0, 1, 2, 3, 4, 5, 6, 7), (0, 2, 3, 5, 8, 11, 13, 17),
             (0, 1, 2, 3, 4, 5, 6, 14), (0, 7, 14, 21, 28, 35, 42, 49)]
    for _ in range(30):
        k = random.randint(6, 8)
        E = tuple(sorted(set([0] + random.sample(range(1, 22), k - 1))))
        tests.append(E)
    bad = 0
    for E in tests:
        a = measS7_fast(E); b = measS7_slow(E)
        if a != b:
            bad += 1
            print(f"   MISMATCH {E}: fast={a} slow={b}")
    print(f"   {len(tests)} shapes, mismatches = {bad}  (must be 0)")
    assert bad == 0
    print("   PASS: both exact engines agree.")


# ===========================================================================
def v1_constants():
    banner("V1 -- re-derive M7(k), consec, delta_k, caps EXACTLY (check stated rationals)")
    stated_consec = {8: F(481, 1470), 9: F(2447, 5880), 10: F(8899, 17640)}
    stated_delta8 = F(1068481, 3529470)
    print(f"   {'k':>2} {'M7(k)':>22} {'consec(meas)':>22} {'delta_k':>22}")
    deltas = {}; consecs = {}
    for k in (8, 9, 10):
        m = M7(k)
        c = measS7_fast(range(k))
        d = c - m
        deltas[k] = d; consecs[k] = c
        print(f"   {k:>2} {str(m):>22} {str(c):>22} {str(d):>22}")
        print(f"        float: M7={float(m):.6f} consec={float(c):.6f} delta={float(d):.6f}")
        if k in stated_consec:
            match = (c == stated_consec[k])
            print(f"        consec matches stated {stated_consec[k]} ? {match}")
            assert match, f"consec mismatch k={k}"
    print(f"\n   delta_8 stated = {stated_delta8} = {float(stated_delta8):.6f}")
    print(f"   delta_8 computed = {deltas[8]} = {float(deltas[8]):.6f}")
    print(f"   delta_8 MATCH ? {deltas[8] == stated_delta8}")
    assert deltas[8] == stated_delta8, "delta_8 mismatch!"
    # delta_9, delta_10 stated as decimals 0.358455, 0.399566
    print(f"   delta_9 computed = {float(deltas[9]):.6f}  (stated 0.358455) match? {abs(float(deltas[9])-0.358455)<1e-6}")
    print(f"   delta_10 computed= {float(deltas[10]):.6f}  (stated 0.399566) match? {abs(float(deltas[10])-0.399566)<1e-6}")
    # caps each >= (k-6)/7
    print("\n   caps >= (k-6)/7 (THM-535) and consec < cap:")
    for k in (8, 9, 10):
        cap = CAP[k]; floor = F(k - 6, 7)
        print(f"     k={k}: cap={cap}={float(cap):.6f}  (k-6)/7={floor}={float(floor):.6f}  "
              f"cap>=floor?{cap>=floor}  consec<cap?{consecs[k]<cap}")
        assert cap >= floor and consecs[k] < cap
    print("   PASS: all stated rational constants re-derived exactly.")
    return deltas, consecs


# ===========================================================================
def shat(n, j):
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)


SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]


def chat(n, T):
    if n == 0:
        return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0:
        return 0j
    return -sum(shat(n, j) for j in T)


def v2_envelope():
    banner("V2 -- pin Lemma-B envelope c1 = max_n |n| max_T |chat(n,T)| (claim 0.697303 @ n=66)")
    # verify |shat| identity
    ok = True
    for n in range(1, 100):
        if n % 7 == 0:
            continue
        for j in range(1, 7):
            lhs = abs(shat(n, j)); rhs = abs(sin(pi * n / 7)) / (pi * n)
            if abs(lhs - rhs) > 1e-12:
                ok = False
    print(f"   |shat(n,j)| = |sin(pi n/7)|/(pi|n|) verified n=1..99 (non-7): {ok}")
    c1 = 0.0; argn = None
    for n in range(1, 500):
        if n % 7 == 0:
            continue
        mx = max(abs(chat(n, T)) for T in SUB) * n
        if mx > c1:
            c1, argn = mx, n
    print(f"   c1 = {c1:.6f}  attained at n={argn} (n mod 7 = {argn % 7})")
    print(f"   stated c1=0.697303 @ n=66:  value match? {abs(c1-0.697303)<1e-5}   argn match? {argn==66}")
    # honest: c1/|n| -> 0; but the per-coordinate sum sum_m c1/|m| diverges (key gap)
    Hsum = sum(c1 / m for m in range(1, 100000) if m % 7 != 0)
    print(f"   sanity: partial sum_(m<1e5, 7-coprime) c1/|m| = {Hsum:.3f} (DIVERGES => free product bound invalid)")
    print("   PASS: envelope constant confirmed; harmonic divergence of free majorant confirmed (the gap).")


# ===========================================================================
def v3_hunt_over_cap_wide(consecs):
    banner("V3 -- AGGRESSIVE HUNT: over-cap OR over-consec PRIMITIVE shape with span > B(k)")
    print("   B(8)=16, B(9)=13, B(10)=12 (wide-spread script's claimed finite-check ceilings).")
    print("   We test many wide families EXACTLY.  A single span>B with meas>cap KILLS the claim;")
    print("   meas>consec (but <cap) would break the stated 'consec is the argmax' but not LRC(14).")
    Bk = {8: 16, 9: 13, 10: 12}
    rng = random.Random(20260619)
    worst = {k: (F(-1), None) for k in (8, 9, 10)}
    over_cap = []; over_consec = []
    n_tested = {8: 0, 9: 0, 10: 0}

    def consider(k, E):
        E = tuple(sorted(set(E)))
        if len(E) != k or E[0] != 0:
            return
        if not primitive(E):
            return
        if max(E) <= Bk[k]:
            return  # only wide
        n_tested[k] += 1
        m = measS7_fast(E)
        if m > worst[k][0]:
            worst[k] = (m, E)
        if m > CAP[k]:
            over_cap.append((k, E, m))
        if m > consecs[k]:
            over_consec.append((k, E, m))

    for k in (8, 9, 10):
        # --- structured wide families ---
        # shifted-AP {0, M..M+k-2} for many M (span = M+k-2)
        for M in range(2, 60):
            consider(k, [0] + list(range(M, M + k - 1)))
        # dilated consec d*{0..k-1} -- scale invariance says equals consec; primitive only if d=1.
        # consec + window shifts: {0,1,...,j} + a tight block far away
        for split in range(2, k - 1):
            base = list(range(split))
            for M in range(split + 1, 50):
                blk = list(range(M, M + (k - split - 1) + 1))
                consider(k, base + blk)
        # 1-stranger: consec_{k-1} + N
        for N in list(range(k, 80)) + [98, 105, 140, 196, 200, 343, 500, 700, 999]:
            consider(k, list(range(k - 1)) + [N])
        # 2-stranger: consec_{k-2} + N1 + N2
        for N1 in range(k - 1, 40):
            for N2 in range(N1 + 1, 60):
                if N2 <= Bk[k]:
                    continue
                consider(k, list(range(k - 2)) + [N1, N2])
        # resonant: heavy use of multiples of 7 (chat(7m)=0 vanishing -> potential structure)
        for base7 in range(1, 6):
            E = [0] + [7 * i for i in range(1, k - 1)] + [7 * (k - 1) + base7]
            consider(k, E)
        for off in range(1, 7):
            consider(k, [0, off] + [7 * i for i in range(1, k - 1)])
        # near-resonant shifted: {0, 7a_i + r}
        for r in range(1, 7):
            consider(k, [0] + [7 * i + r for i in range(1, k)])
        # short-relation wide: an AP with common difference d, large d (dissociation broken)
        for d in range(2, 30):
            E = [0] + [d * i + (i % 2) for i in range(1, k)]  # near-AP with perturbation
            consider(k, E)
        # geometric / Sidon-ish wide
        consider(k, [0] + [2 ** i for i in range(1, k)])
        consider(k, [0] + [3 ** i for i in range(1, k)])
        # cluster {0} U {M..M+k-2} already covered; cluster with internal AP gap d
        for d in (2, 3, 5, 7):
            for M in range(20, 50):
                consider(k, [0] + [M + d * j for j in range(k - 1)])
        # --- random wide nets ---
        for _ in range(8000):
            sp = rng.randint(Bk[k] + 1, 120)
            rest = sorted(set(rng.sample(range(1, sp), k - 2) + [sp]))
            E = [0] + rest
            consider(k, E)
        # random tight clusters far away (no-scale-separation)
        for _ in range(4000):
            M = rng.randint(30, 400)
            w = rng.randint(k, 4 * k)
            rest = sorted(set(rng.sample(range(M, M + w), k - 1)))
            consider(k, [0] + rest)

    for k in (8, 9, 10):
        wm, we = worst[k]
        print(f"   k={k}: tested {n_tested[k]} wide primitive shapes (span>{Bk[k]}).")
        print(f"        WORST meas = {float(wm):.6f} at {we}")
        print(f"        cap_k={float(CAP[k]):.6f}  consec_k={float(consecs[k]):.6f}")
        print(f"        worst < cap? {wm <= CAP[k]}   worst < consec? {wm <= consecs[k]}")
    print(f"\n   OVER-CAP shapes found (span>B): {len(over_cap)}")
    for k, E, m in over_cap[:20]:
        print(f"      !! k={k} {E} meas={float(m):.6f} > cap={float(CAP[k]):.6f}")
    print(f"   OVER-CONSEC shapes found (span>B, meas>consec but maybe <cap): {len(over_consec)}")
    for k, E, m in over_consec[:20]:
        print(f"      ~~ k={k} {E} meas={float(m):.6f} > consec={float(consecs[k]):.6f} "
              f"(<cap? {m<=CAP[k]})")
    return over_cap, over_consec


# ===========================================================================
def v4_finite_check_consistency(consecs):
    banner("V4 -- finite-check: B(k) consistency + exhaustive exceedance counts (EXACT)")
    print("   The wide-spread script claims B(8)=16,B(9)=13,B(10)=12 with 0 over-cap.")
    print("   The companion finite-check-to-B script uses B(8)=16,B(9)=15,B(10)=13.")
    print("   We re-run the EXACT exhaustive box for the wide-spread script's B and report.")
    Bk = {8: 16, 9: 13, 10: 12}
    for k in (8, 9, 10):
        cap = CAP[k]; cons = consecs[k]; B = Bk[k]
        cnt = 0; over_cap = 0; over_consec = 0
        mx = F(-1); arg = None; mx_nonconsec = F(-1)
        consec_tuple = tuple(range(k))
        for rest in itertools.combinations(range(1, B + 1), k - 1):
            E = (0,) + rest
            if not primitive(E):
                continue
            cnt += 1
            m = measS7_fast(E)
            if m > mx:
                mx, arg = m, E
            if E != consec_tuple and m > mx_nonconsec:
                mx_nonconsec = m
            if m > cap:
                over_cap += 1
            if m > cons:
                over_consec += 1
        print(f"   k={k}, span<={B}: {cnt} primitive sets; over_cap={over_cap}; over_consec={over_consec}")
        print(f"        argmax={arg} (consec? {arg==consec_tuple})  max={float(mx):.6f} cap={float(cap):.6f}")
        print(f"        largest NON-consec meas = {float(mx_nonconsec):.6f}")


# ===========================================================================
def v5_keylemma_audit():
    banner("V5 -- audit the CORE-EXTRACTION KEY LEMMA (does decay 1/e_L bound the tail?)")
    print("   The lemma covers case (b): exactly ONE large coordinate.  Case (a): >=2 large")
    print("   coordinates is NOT bounded by 1/e_L.  We test whether 2-large relations can keep")
    print("   the corr large for arbitrarily large span (the lemma's blind spot).")
    print("   Test: consec_{k-2} + {N, N+1} (a 2-stranger pair, both large) at k=8, growing N.")
    cap8 = CAP[8]; cons8 = measS7_fast(range(8))
    print(f"   cap_8={float(cap8):.6f} consec_8={float(cons8):.6f}")
    for N in (20, 40, 80, 160, 320, 640):
        E = tuple(list(range(6)) + [N, N + 1])
        if not primitive(E):
            E = tuple(list(range(6)) + [N, N + 3])
        m = measS7_fast(E)
        print(f"   2-stranger {{0..5, {N},{N+1}}}: meas={float(m):.6f}  <cap?{m<cap8} <consec?{m<cons8}")
    print("\n   Also: a tight cluster of 7 near M plus 0 (all 'large' relative to 0 -> no core).")
    print("   The key lemma's single-large-offset decay does NOT apply (every offset is large).")
    for M in (40, 80, 160, 320):
        E = tuple([0] + list(range(M, M + 7)))
        m = measS7_fast(E)
        print(f"   no-core cluster {{0,{M}..{M+6}}}: meas={float(m):.6f}  <cap?{m<cap8}")
    print("\n   READING: the 1/e_L decay is real for ONE large offset, but the >=2-large and")
    print("   no-core cases are handled by a SEPARATE argument (cluster-collapse), not the lemma.")
    print("   The lemma alone does NOT bound the full type-(II) tail. (Script admits eps(B) open.)")


# ===========================================================================
def v6_cluster_collapse_audit():
    banner("V6 -- cluster-collapse: is the O(w/M) Weyl error uniform? finite-M over-cap risk?")
    print("   Claim: tight cluster {0}U{M..M+6} -> wide-ceiling ~0.193 as M->inf (THM-518).")
    print("   Risk the proof must exclude: a FINITE M (or a non-AP cluster) exceeding cap_8")
    print("   before the asymptote.  We scan many cluster shapes and ALL M EXACTLY.")
    cap8 = CAP[8]; cons8 = measS7_fast(range(8))
    worst = F(-1); warg = None
    # AP clusters, all M and all internal differences d
    for d in (1, 2, 3, 4, 5, 6):
        for M in range(8, 200):
            E = tuple([0] + [M + d * j for j in range(7)])
            if not primitive(E):
                continue
            m = measS7_fast(E)
            if m > worst:
                worst, warg = m, E
    # non-AP clusters: random 7-subsets of a window [M, M+w]
    rng = random.Random(7)
    for _ in range(6000):
        M = rng.randint(8, 300)
        w = rng.randint(7, 30)
        sub = sorted(rng.sample(range(M, M + w), 7))
        E = tuple([0] + sub)
        if not primitive(E):
            continue
        m = measS7_fast(E)
        if m > worst:
            worst, warg = m, E
    print(f"   WORST cluster meas over all scanned = {float(worst):.6f} at {warg}")
    print(f"   cap_8={float(cap8):.6f}  worst<cap? {worst<cap8}  worst<consec_8? {worst<cons8}")
    # smallest M still over the wide ceiling? show the M-trajectory for the AP cluster d=1
    print("\n   AP cluster {0,M..M+6} meas trajectory (d=1):")
    for M in (8, 10, 14, 21, 28, 50, 100):
        E = tuple([0] + [M + j for j in range(7)])
        m = measS7_fast(E)
        print(f"      M={M:4d}: meas={float(m):.6f}  <cap?{m<cap8}")
    print("\n   READING: even the WORST finite-M / non-AP cluster stays below cap_8 in this scan,")
    print("   but this is EXACT MEASUREMENT over a finite scan, NOT a uniform-in-M theorem.")
    print("   The 'O(w/M)' has no explicit M-threshold => the span>B side is not gap-free here.")


# ===========================================================================
def v7_verdict():
    banner("V7 -- VERDICT (adversarial)")
    print(r"""
  RE-DERIVED EXACTLY (all confirmed):
    * M7(8/9/10), consec_8=481/1470, consec_9=2447/5880, consec_10=8899/17640.
    * delta_8 = 1068481/3529470 = 0.302731 (exact match), delta_9=0.358455, delta_10=0.399566.
    * caps: 2243/5880, 1979/4004, 55/91; each >= (k-6)/7; consec < cap (margins 0.054/0.078/0.100).
    * Lemma-B envelope c1 = 0.697303 attained at n=66.  Seven-multiple vanishing chat(7m)=0.
    * Finite check span<=B(k): EXACT exhaustive, 0 over-cap, consec the unique argmax,
      at B(8)=16 (11432 sets), B(9)=13 (1287), B(10)=12 (220).

  GAPS THAT SURVIVE (the piece is NOT gap-free):

  [G1] UNIFORM TAIL CONSTANT eps(B) -- THE primary blocker, ADMITTED in the script's own
       Part 8.  The free per-coordinate majorant sum_{m} c1/|m| DIVERGES (harmonic); the
       claimed eps(B)=2^6 c1^6 Theta6(B) is NOT a finite closed-form bound.  The convergence
       requires a successive-minima / Minkowski lattice-point count that is STATED but NEVER
       carried out.  No explicit constant, no explicit B-threshold from the analytic side.
       => "span > B(k) => meas <= cap" is NOT proved analytically; it rests on sampling.

  [G2] KEY LEMMA covers only the SINGLE-large-offset case (b).  The >=2-large case (a) and
       the no-core (every-offset-large) case are explicitly NOT bounded by the 1/e_L decay;
       they are deferred to cluster-collapse, which is itself non-uniform (G3).  So the
       "KEY LEMMA" does not, by itself, bound the type-(II) tail.

  [G3] CLUSTER-COLLAPSE Weyl error is "O(w/M)" with NO explicit M-threshold and NO uniform
       constant over all cluster offset SHAPES.  Verified by exact finite scan (no over-cap
       found), but a finite scan is not a uniform-in-M inequality.

  [G4] FINITE CHECK is exhaustive ONLY to span<=B(k); the span>B(k) side is covered by
       decay+collapse+SAMPLING, which (per G1-G3) is not a uniform theorem.  Moreover the
       two companion scripts disagree on B(9),B(10) (13/12 here vs 15/13 in finite-check-to-B);
       neither value is large enough to make the box exhaustive on its own (span is unbounded).

  [G5] "12/12 wide families have corr < delta_k" is EXACT MEASUREMENT on canonical families,
       NOT a uniform inequality over the full (k-1)-parameter shape space.  The claim
       "corr(E) <= delta_k for ALL wide E" is unproven.

  VERDICT:  NOT GAP-FREE.  The piece is a substantial, honest REDUCTION: it re-derives every
  constant correctly, delivers explicit B(k), proves the single-large-offset decay mechanism,
  and supplies an exhaustive bounded-spread finite check with 0 exceedances.  But the span>B(k)
  half is reduced to -- not closed by -- a uniform Minkowski/successive-minima estimate (G1),
  with the >=2-large (G2) and cluster (G3) cases relying on non-uniform/empirical control.
  The single blocker the prompt names (uniform tail constant eps(B)) is REAL and OPEN.
""")


def main():
    print("LRC(14) WIDE-SPREAD-RIGOROUS  --  ADVERSARIAL VERIFICATION  (kps-S10-wf)")
    v0_engine_crossvalidate()
    deltas, consecs = v1_constants()
    v2_envelope()
    v3_hunt_over_cap_wide(consecs)
    v4_finite_check_consistency(consecs)
    v5_keylemma_audit()
    v6_cluster_collapse_audit()
    v7_verdict()
    print("\nDONE (verification).")


if __name__ == "__main__":
    main()
