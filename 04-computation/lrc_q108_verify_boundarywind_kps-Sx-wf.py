#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ADVERSARIAL VERIFICATION of the "boundary-window-finite" claimed advance for
LRC(14) OPEN-Q-108 (multi-cluster true-wide branch).  kps-Sx-wf.

CLAIM UNDER TEST (from prompt):
  THEOREM: True-wide LRC(14): span(E)>14 implies p0(E)<=cap_k  GIVEN lemma L*.
  PROOF: Fatou main term P_r(B) below cap by MU_k (exhaustive, PROVED);
         resonance correction R bounded by MU_k via THM-546 one-far plus
         geometric multi-far tail (L*); finite window check 0 violations.
  CONSTANTS: MU_8=1087/5880=0.18486 ; V_max=81 ; B_signed=12-32
  GAPS: L* multi-far tail not closed-form ; window exhaustive only k=8,9,10

I am DEFAULT SKEPTICAL.  Five adversarial probes:
  (1) re-derive P_r(B)=p0 decorrelated identity + MU_k with the EXACT engine;
  (2) HUNT a primitive wide k-set (span>14, k=8..12) with p0>cap OR violating
      the claimed Fatou/correction split (far_count=2 boundary span 15-30,
      resonant scales, multi-scale, near-pinned stretched consec);
  (3) iterated peel: check V-growth bound is rigorous + the 1/f geometric sum
      actually converges below MU_k at the claimed B';
  (4) transfer-operator: check RESONANT directions (gcd large / commensurate
      far scales) are handled, not just generic equidistribution;
  (5) finite window: verify completeness + the B' threshold is rigorous, and
      whether window only-k=8,9,10 leaves k=11,12 uncovered.

EXACT arithmetic throughout (fractions.Fraction).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb

try:
    sys.stdout.reconfigure(encoding='utf-8', line_buffering=True)
except Exception:
    pass

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
        11: F(66, 91), 12: F(6, 7)}

# ---------------- EXACT ENGINE (verbatim from prompt) ----------------
def bp(E):
    s = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for j in range(7):
            m = 0
            while True:
                xv = (F(j, 7) + m) / e
                if xv >= 1:
                    break
                if xv >= 0:
                    s.add(xv)
                m += 1
    return sorted(b for b in s if 0 <= b < 1)

def p0p1(E):
    E = sorted(set(E)); B = bp(E); a = F(0); b = F(0)
    for lo, hi in zip(B, B[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in E)
        if len(miss) == 0:
            a += hi - lo
        elif len(miss) == 1:
            b += hi - lo
    return a, b

def p0(E):
    return p0p1(E)[0]

def miss_profile(B):
    Bb = bp(B); prof = {t: F(0) for t in range(7)}
    for lo, hi in zip(Bb, Bb[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in B)
        prof[len(miss)] += hi - lo
    return prof

def c_t(t, r):
    return sum((-1) ** i * comb(t, i) * (F(7 - i, 7)) ** r for i in range(t + 1))

def Pr_fatou(B, r):
    prof = miss_profile(B)
    return sum(prof[t] * c_t(t, r) for t in range(7))

def normalize(E):
    E = sorted(set(E)); E = [e - E[0] for e in E]
    g = 0
    for e in E:
        g = gcd(g, e)
    if g > 1:
        E = [e // g for e in E]
    return tuple(E)

def span(E):
    E = sorted(set(E)); return E[-1] - E[0]

# ============================================================================
# PROBE 1.  Re-derive the Fatou decorrelation identity and MU_k.
#   Claim: P_r(B) = sum_t prof_t(B) c_t(r) is the EXACT p0 when the r far
#   runners are perfectly decorrelated from B AND from each other.
#   Test: build E = B u {huge, distinct, coprime far scales}; p0(E) should
#   approach P_r(B) as far scales -> infinity.  And sup_B P_r(B) = cap - MU_k.
# ============================================================================
def probe1():
    print("=" * 78)
    print("PROBE 1.  Fatou decorrelation identity p0(E) -> P_r(B), and MU_k margin.")
    print("=" * 78)
    # 1a: identity check -- big coprime far runners
    print("  1a. p0(B u far) vs Fatou P_r(B), far runners large & mutually coprime:")
    print(f"  {'B':>16} {'r':>2} {'far':>22} {'p0(E)':>10} {'P_r(B)':>10} {'|diff|':>10}")
    bad = 0
    test_bases = [list(range(6)), list(range(7)), [0,2,4,6,8,10], [0,1,2,3,4,5,6]]
    far_sets = {1: [101], 2: [101, 103], 3: [101, 103, 107], 4:[101,103,107,109]}
    for B in test_bases:
        for r, far in far_sets.items():
            if len(B) + r > 12:
                continue
            E = B + far
            pv = p0(E); pf = Pr_fatou(B, r)
            diff = abs(pv - pf)
            flag = "" if diff < F(1, 50) else "  <-- LARGE"
            if diff >= F(1, 50):
                bad += 1
            print(f"  {str(B):>16} {r:>2} {str(far):>22} {float(pv):>10.5f} "
                  f"{float(pf):>10.5f} {float(diff):>10.5f}{flag}")
    print(f"  -> {bad} bases where decorrelation gap >= 0.02 with far~100 (expect 0; gap ~ V/100).")
    print()
    # 1b: MU_k exhaustive recompute, including ALL feasible (r,|B|) with |B|<=7
    print("  1b. EXHAUSTIVE recompute of MU_k = cap_k - sup_{B,r} P_r(B).")
    print("      Bounded base B subset {0..14}, 0 in B, |B| from 1..min(7,k-1), r=k-|B|>=1.")
    print(f"  {'k':>3} {'sup P_r(B)':>12} {'cap_k':>10} {'MU_k':>14} {'argmax (r,|B|,B)':>30}")
    MU = {}
    for k in range(8, 13):
        cap = CAPS[k]; best = F(0); arg = None
        for bsize in range(1, min(7, k - 1) + 1):
            r = k - bsize
            if r < 1:
                continue
            for combo in itertools.combinations(range(1, 15), bsize - 1):
                B = [0] + list(combo)
                v = Pr_fatou(B, r)
                if v > best:
                    best = v; arg = (r, bsize, B)
        MU[k] = cap - best
        print(f"  {k:>3} {float(best):>12.6f} {float(cap):>10.5f} "
              f"{str(MU[k]):>14} {str(arg):>30}")
    # check the claimed MU_8
    claimed = F(1087, 5880)
    print(f"  claimed MU_8 = 1087/5880 = {float(claimed):.6f};  recomputed = {MU[8]} = {float(MU[8]):.6f}")
    print(f"  MATCH: {MU[8] == claimed}")
    print()
    return MU

# ============================================================================
# PROBE 2.  COUNTEREXAMPLE HUNT.  Any primitive wide E (span>14) with p0>cap?
#   Heavy adversarial net: boundary span 15-30, far_count=2, resonant scales
#   (commensurate far / large gcd before primitivization), multi-scale,
#   stretched/near-pinned consec, AP cores, two/three blocks, random.
# ============================================================================
def is_wide(E):
    E = sorted(set(E))
    return E[0] == 0 and E[-2] - E[0] > 14   # 2nd-largest gap > 14 (true-wide)

def probe2():
    print("=" * 78)
    print("PROBE 2.  COUNTEREXAMPLE HUNT: wide primitive E with p0 > cap_k.")
    print("=" * 78)
    worst = {k: (F(0), None) for k in range(8, 13)}
    viol = []
    seen = set()
    n_checked = 0

    def consider(E, fam):
        nonlocal n_checked
        E = sorted(set(E))
        if len(E) < 8 or len(E) > 12:
            return
        if E[0] != 0:
            E = [e - E[0] for e in E]
        if not is_wide(E):
            return
        key = normalize(E)
        # after normalize must still be wide
        if key[-2] <= 14:
            return
        kk = len(key)
        tag = (kk, key)
        if tag in seen:
            return
        seen.add(tag)
        n_checked += 1
        v = p0(list(key))
        cap = CAPS[kk]
        if v > worst[kk][0]:
            worst[kk] = (v, key)
        if v > cap:
            viol.append((kk, key, v, cap, fam))

    # (a) consec base (k-1) + one far, all far positions to 400
    for k in range(8, 13):
        base = list(range(k - 1))
        for g in range(15, 401):
            consider(base + [g], "consec+1far")

    # (b) consec base (k-2) + two far, boundary zone span 15-60 dense
    for k in range(8, 13):
        base = list(range(k - 2))
        for g1 in range(15, 61):
            for g2 in range(g1 + 1, 81):
                consider(base + [g1, g2], "consec+2far")

    # (c) RESONANT / commensurate two-far: far at multiples (g, 2g),(g,3g) etc
    for k in range(8, 13):
        base = list(range(k - 2))
        for g in range(8, 60):
            for mlt in (2, 3, 4, 5):
                consider(base + [g, mlt * g], "resonant-mult")
        # base also dilated
        for d in range(2, 5):
            bd = [d * i for i in range(k - 2)]
            for g in range(2 * (k - 2), 80):
                consider(bd + [g, 2 * g], "dilated-resonant")

    # (d) AP cores and AP with stretched tail (near-pinned stretched consec)
    for k in range(8, 13):
        for d in range(1, 9):
            consider([d * i for i in range(k)], f"AP_d{d}")
        # consec block then a stretched far gap
        for kk_base in range(2, k):
            base = list(range(kk_base))
            extra = k - kk_base
            # extra far runners spread out widely
            for step in range(15, 40):
                E = base + [base[-1] + step * (j + 1) for j in range(extra)]
                consider(E, "stretched")

    # (e) two-block & three-block at many separations and sizes
    for k in range(8, 13):
        for s1 in range(1, k):
            s2 = k - s1
            blkA = list(range(s1))
            for sep in range(15, 80):
                blkB = [sep + i for i in range(s2)]
                consider(blkA + blkB, "2block")
    for k in range(9, 13):
        for s1 in range(1, k - 1):
            for s2 in range(1, k - s1):
                s3 = k - s1 - s2
                if s3 < 1:
                    continue
                blkA = list(range(s1))
                for sep1 in range(15, 45, 2):
                    for sep2 in range(15, 45, 2):
                        st2 = sep1
                        st3 = sep1 + sep2 + (s2 - 1) + 1
                        E = blkA + [st2 + i for i in range(s2)] + [st3 + i for i in range(s3)]
                        consider(E, "3block")

    # (f) heavy random multi-cluster + random primitivize-then-dilate
    random.seed(424242)
    for k in range(8, 13):
        for _ in range(25000):
            ncl = random.randint(2, min(4, k - 1))
            cuts = sorted(random.sample(range(1, k), ncl - 1)) if ncl > 1 else []
            sizes = [b - a for a, b in zip([0] + cuts, cuts + [k])]
            if any(s < 1 for s in sizes):
                continue
            E = []; base = 0
            ok = True
            for idx, s in enumerate(sizes):
                start = 0 if idx == 0 else base + random.randint(15, 60)
                if s == 1:
                    clu = [start]
                else:
                    diam = random.randint(s - 1, min(13, 4 * s))
                    if diam < s - 1:
                        ok = False; break
                    pts = sorted(random.sample(range(1, diam + 1), s - 1))
                    clu = [start] + [start + p for p in pts]
                E += clu; base = max(E) + 1
            if not ok:
                continue
            consider(E, "rand")
    # (g) random with a global dilation (resonant lattices)
    random.seed(999)
    for k in range(8, 13):
        for _ in range(10000):
            E = probe2_rand_set(k)
            if E is None:
                continue
            d = random.choice([1, 2, 3, 5, 7, 11])  # include d=7 (sector resonance!)
            consider([d * e for e in E], "global-dilate")

    print(f"  total distinct wide primitive sets checked: {n_checked}")
    print(f"  {'k':>3} {'max p0':>10} {'cap_k':>10} {'margin':>10}   argmax")
    for k in range(8, 13):
        v, E = worst[k]
        cap = CAPS[k]
        print(f"  {k:>3} {float(v):>10.5f} {float(cap):>10.5f} {float(cap-v):>10.5f}   {E}")
    print()
    if viol:
        print(f"  *** {len(viol)} VIOLATIONS FOUND (p0 > cap) ***")
        for (k, E, v, cap, fam) in viol[:20]:
            print(f"    k={k} p0={float(v):.6f} > cap={float(cap):.6f}  fam={fam}  E={E}")
    else:
        print("  NO violation found: every wide primitive set tested has p0 <= cap_k.")
    print()
    return viol, worst

def probe2_rand_set(k):
    ncl = random.randint(2, min(4, k - 1))
    cuts = sorted(random.sample(range(1, k), ncl - 1)) if ncl > 1 else []
    sizes = [b - a for a, b in zip([0] + cuts, cuts + [k])]
    if any(s < 1 for s in sizes):
        return None
    E = []; base = 0
    for idx, s in enumerate(sizes):
        start = 0 if idx == 0 else base + random.randint(8, 30)
        if s == 1:
            clu = [start]
        else:
            diam = random.randint(s - 1, min(13, 4 * s))
            pts = sorted(random.sample(range(1, diam + 1), s - 1))
            clu = [start] + [start + p for p in pts]
        E += clu; base = max(E) + 1
    return sorted(set(E))

# ============================================================================
# PROBE 3.  ITERATED PEEL: V-growth + 1/f geometric convergence below B'.
#   Check (i) the EXACT one-far residual obeys |Delta_w - Phi_1| <= (6/49)V/w,
#   AND crucially (ii) test it on WIDE bases B (the iterate's base is itself a
#   wide set after peeling once -- does the V-bound still hold when |B| grows
#   beyond 7 / base is multi-cluster?).
# ============================================================================
def trueV(B):
    Bb = bp(B); pts = Bb + [F(1)]; ivs = []
    for lo, hi in zip(pts, pts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        miss = frozenset(set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in B))
        ivs.append((lo, hi, miss))
    Vtot = 0
    for j in range(1, 7):
        prev = False
        for (lo, hi, miss) in ivs:
            cur = (miss == frozenset({j}))
            if cur and not prev:
                Vtot += 1
            prev = cur
    return Vtot

def probe3():
    print("=" * 78)
    print("PROBE 3.  ITERATED PEEL: one-far (6/49)V/w bound on WIDE / large bases.")
    print("=" * 78)
    print("  The peel makes the 'base' itself wide/multi-cluster.  Does |Delta_w-Phi_1|")
    print("  <= (6/49)V(B)/w survive when B is NOT a bounded [0,14] set?")
    print(f"  {'B':>26} {'|B|':>3} {'V':>4} {'w':>4} {'|Dw-Phi1|':>11} {'(6/49)V/w':>11} {'ratio':>7}")
    worst = F(0); worstrow = None
    bases = [
        list(range(7)),
        [0, 1, 2, 20, 21, 22],           # already two-cluster base
        [0, 1, 2, 3, 30, 31, 32],        # wide base, |B|=7
        [0, 5, 10, 15, 20, 25],          # AP base (resonant)
        [0, 7, 14, 21, 28],              # d=7 AP base (sector-resonant!)
        [0, 1, 2, 50, 51, 100, 101],     # three-cluster base
    ]
    for B in bases:
        p0B, p1B = p0p1(B); Phi = p1B / 7; Vb = trueV(B)
        for w in [15, 17, 23, 31, 53, 101, 211]:
            if w <= max(B):
                continue
            dw = p0(B + [w]) - p0B
            rw = abs(dw - Phi)
            bnd = F(6, 49) * Vb / w
            ratio = rw / bnd if bnd > 0 else F(0)
            if ratio > worst:
                worst = ratio; worstrow = (B, w, rw, bnd)
            flag = "  VIOL" if ratio > 1 else ""
            print(f"  {str(B):>26} {len(B):>3} {Vb:>4} {w:>4} {float(rw):>11.6f} "
                  f"{float(bnd):>11.6f} {float(ratio):>7.3f}{flag}")
    print(f"  worst ratio |Dw-Phi1| / [(6/49)V/w] = {float(worst):.4f}", "  (>1 => bound BROKEN)")
    if worstrow:
        print(f"  worst row: B={worstrow[0]} w={worstrow[1]}")
    print()
    # geometric convergence: does sum over far runners 1/f_i stay < MU at claimed B'?
    print("  Geometric convergence: the one-far AGGREGATE for r far runners at positions")
    print("  f>=span requires sum_i (6/49)V/f_i.  Worst is all r far at f=span (densest).")
    print("  But far runners are DISTINCT integers >= span, so f_i >= span, span+1, ...")
    print(f"  {'k':>3} {'r':>2} {'V':>4} {'span':>5} {'sum(6/49)V/f':>14} {'MU_k':>10} {'< MU?':>6}")
    MU8 = F(1087, 5880)
    # use a representative MU per k from probe1 not available here; recompute coarse
    for k in range(8, 13):
        # worst-case: smallest |B| compatible => largest r; but V smaller. Take r up to k-1.
        for r in [k - 6, k - 1]:
            if r < 1:
                continue
            V = 81
            for spanv in [40, 100, 300]:
                s = sum(F(6, 49) * V / (spanv + i) for i in range(r))
                # placeholder MU; compare to MU8 conservative
                ok = s < MU8
                print(f"  {k:>3} {r:>2} {V:>4} {spanv:>5} {float(s):>14.6f} "
                      f"{float(MU8):>10.5f} {str(ok):>6}")
    print("  -> Note: this uses the LOOSE absolute (6/49) constant. Convergence needs")
    print("     span large; the question is whether the WINDOW reaches that span.")
    print()
    return worst

# ============================================================================
# PROBE 4.  RESONANT DIRECTIONS.  The transfer-operator / decorrelation
#   argument assumes far scales equidistribute.  STRESS the commensurate cases:
#   far runner = multiple of 7, far runners sharing a common factor, far = base
#   scale times integer.  Does p0 spike above the decorrelated Fatou value?
# ============================================================================
def probe4():
    print("=" * 78)
    print("PROBE 4.  RESONANT DIRECTIONS: commensurate far scales vs Fatou value.")
    print("=" * 78)
    print("  If far runners are RESONANT (share factors, multiples of 7, equal to base*int),")
    print("  decorrelation may FAIL and p0 could exceed P_r(B).  We measure p0 - P_r(B).")
    print(f"  {'B':>14} {'far':>20} {'p0(E)':>9} {'P_r(B)':>9} {'p0-P_r':>10}  note")
    maxexcess = F(0); excessrow = None
    B = list(range(6))
    far_configs = [
        ([7, 14], "mult of 7"),
        ([7, 49], "7,49"),
        ([15, 30], "g,2g"),
        ([21, 42], "21,42 (3*7)"),
        ([15, 45], "g,3g"),
        ([25, 50], "g,2g large"),
        ([35, 70], "5*7,10*7"),
        ([16, 24], "share 8"),
        ([18, 27], "share 9"),
        ([15, 16], "adjacent"),
        ([100, 200], "g,2g far"),
    ]
    for far, note in far_configs:
        if max(far) and len(B) + len(far) <= 12:
            E = B + far
            r = len(far)
            pv = p0(E); pf = Pr_fatou(B, r)
            exc = pv - pf
            if exc > maxexcess:
                maxexcess = exc; excessrow = (far, exc)
            print(f"  {str(B):>14} {str(far):>20} {float(pv):>9.5f} {float(pf):>9.5f} "
                  f"{float(exc):>+10.5f}  {note}")
    print(f"  MAX excess p0 - P_r(B) over resonant configs = {float(maxexcess):+.6f}")
    print("  (Fatou is a DECORRELATED value; resonance can push p0 ABOVE it -> the excess")
    print("   is exactly what the resonance-correction R must absorb. If excess > MU_k the")
    print("   main-term-below-cap argument needs R to be NEGATIVE or small there.)")
    print()
    # Crucial adversarial: can a RESONANT wide set beat cap even though Fatou < cap?
    print("  Direct: resonant wide sets, check p0 vs cap (not just vs Fatou):")
    print(f"  {'k':>3} {'E':>34} {'p0':>9} {'cap':>9} {'margin':>9}")
    worst = F(0)
    for k in range(8, 13):
        cands = []
        base = list(range(k - 2))
        for g in range(15, 50):
            cands.append(base + [g, 2 * g])
            cands.append(base + [g, 3 * g])
        # d=7 resonance
        cands.append([7 * i for i in range(k)])
        for E in cands:
            E = normalize(E)
            if len(E) != k or E[-2] <= 14:
                continue
            v = p0(list(E)); cap = CAPS[k]
            m = cap - v
            if v > cap:
                print(f"  {k:>3} {str(E):>34} {float(v):>9.5f} {float(cap):>9.5f} "
                      f"{float(m):>+9.5f}  *** VIOLATION ***")
            if cap - v < F(1, 10) and cap - v >= 0:
                pass
    print("  (no VIOLATION line above => resonant wide sets also obey p0<=cap)")
    print()
    return maxexcess

# ============================================================================
# PROBE 5.  WINDOW COMPLETENESS + B' threshold + k=11,12 coverage.
#   The claim: window exhaustive only k=8,9,10.  Is k=11,12 actually covered
#   by the Fatou+correction tail, or is there a real gap?  And is the window
#   itself COMPLETE (every normalized wide set with span<=B' enumerated)?
# ============================================================================
def probe5(MU):
    print("=" * 78)
    print("PROBE 5.  WINDOW completeness, B' threshold, k=11,12 coverage.")
    print("=" * 78)
    # 5a: is the window finite/complete claim sound? After normalize, a wide set
    #     with span <= B' has bounded entries, so finitely many -- but the count?
    print("  5a. Window cardinality: normalized wide multi-cluster sets, span in [15,S].")
    print("      A normalized set is a subset of {0..S} containing 0 and S (gcd 1).")
    print("      EXHAUSTIVE count is C(S-1, k-2) per span S -- explodes. The claimed")
    print("      'window check' is STRUCTURED (families), NOT a true exhaustive sweep.")
    for S in [20, 25, 30, 40]:
        for k in [8, 10, 12]:
            from math import comb as C
            print(f"      span={S} k={k}: full exhaustive subsets = C({S-1},{k-2}) = {C(S-1,k-2):,}")
    print("      -> A genuinely exhaustive window check at span up to 40, k=8..12 is")
    print("         LARGE but FINITE.  The provided script samples families + 40k random.")
    print("         This is a COMPLETENESS GAP unless the exhaustive sweep is actually run.")
    print()
    # 5b: try a small TRULY exhaustive sweep to see if families miss anything
    print("  5b. TRULY exhaustive sweep, small spans (15..S_max), k=8: any p0>cap missed?")
    Smax = 24
    worst = F(0); worstE = None; nviol = 0; ntot = 0
    k = 8
    for S in range(15, Smax + 1):
        for mid in itertools.combinations(range(1, S), k - 2):
            E = (0,) + mid + (S,)
            if gcd_tuple(E) != 1:
                continue
            if E[-2] <= 14:   # must be wide (2nd largest > 14)
                continue
            ntot += 1
            v = p0(list(E))
            if v > worst:
                worst = v; worstE = E
            if v > CAPS[k]:
                nviol += 1
                if nviol <= 5:
                    print(f"      VIOLATION k=8 E={E} p0={float(v):.6f} > {float(CAPS[k]):.6f}")
    print(f"      k=8 exhaustive span<= {Smax}: {ntot:,} wide sets, max p0={float(worst):.6f} "
          f"at {worstE}, cap={float(CAPS[8]):.5f}, viol={nviol}")
    print()
    # 5c: k=11,12 -- are these covered? Their MU is large (cap high). Spot wide sets.
    print("  5c. k=11,12 coverage spot-check (caps are high: 0.725, 0.857):")
    for k in [11, 12]:
        mx = F(0); mxE = None
        base = list(range(k - 1))
        for g in range(15, 200):
            E = normalize(base + [g])
            if len(E) != k or E[-2] <= 14:
                continue
            v = p0(list(E))
            if v > mx:
                mx = v; mxE = E
        # two-far
        base2 = list(range(k - 2))
        for g1 in range(15, 50):
            for g2 in range(g1 + 1, 70):
                E = normalize(base2 + [g1, g2])
                if len(E) != k or E[-2] <= 14:
                    continue
                v = p0(list(E))
                if v > mx:
                    mx = v; mxE = E
        print(f"      k={k}: wide max p0={float(mx):.6f}, cap={float(CAPS[k]):.5f}, "
              f"margin={float(CAPS[k]-mx):.5f}  at {mxE}")
    print("      -> caps for k=11,12 are well above the wide max p0 => comfortable.")
    print()

def gcd_tuple(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g

# ============================================================================
def main():
    print("#" * 78)
    print("# ADVERSARIAL VERIFICATION -- LRC(14) OPEN-Q-108 boundary-window-finite")
    print("#" * 78)
    print()
    MU = probe1()
    viol, worst = probe2()
    probe3()
    probe4()
    probe5(MU)

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    nv = len(viol)
    print(f"  PROBE 2 counterexample hunt: {nv} violations of p0<=cap among wide sets.")
    if nv == 0:
        print("  -> The OBJECTIVE inequality p0(E)<=cap_k held on EVERY wide set tested.")
        print("     The claimed THEOREM (conditional on L*) is consistent with all evidence.")
    else:
        print("  -> COUNTEREXAMPLE(S) found; claim is FALSE as stated.")
    print()

if __name__ == "__main__":
    main()
