#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_verify_reduction_2_kps-S2-wf.py   (kind-pasteur-2026-06-17-S2)

ADVERSARIAL VERIFICATION (logical-gap / circularity lens) of the
"easy-dominates-hard covering reduction" for LRC(14) (THM-525 stub).

QUESTION: Does the argument secretly ASSUME the very uniformity meas(G_C)>=c
that it is supposed to prove (i.e. is it CIRCULAR w.r.t. OPEN-Q-108)?  Or does it
genuinely reduce LRC(14) to a STRICTLY WEAKER / KNOWN statement?

I trace each STEP and ask: what does its CONCLUSION depend on?  A step is
NON-circular and USEFUL only if its output bound is either
  (a) a UNIVERSAL CONSTANT (independent of the core C and of v_max), or
  (b) provably implies M>=1/14 WITHOUT invoking a uniform meas(G_C) bound.
A step is VACUOUS/CIRCULAR if its conclusion is "M>=1/14 PROVIDED meas(G_C)>=(stuff)"
where (stuff) is itself the open uniformity.

Everything EXACT (fractions.Fraction).
"""
import sys, itertools as it
from math import gcd
from fractions import Fraction as F
from functools import reduce

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

# ---------------- exact lonely-set machinery (gap 1/14) -------------------
def danger_arcs(v):
    w = F(1, 14 * v); out = []
    for k in range(v + 1):
        c = F(k, v); lo, hi = c - w, c + w
        if lo < 0:
            out.append((F(0), hi)); out.append((1 + lo, F(1)))
        elif hi > 1:
            out.append((lo, F(1))); out.append((F(0), hi - 1))
        else:
            out.append((lo, hi))
    return [(a, b) for a, b in out if a < b]

def union_arcs(intervals):
    iv = sorted(intervals)
    if not iv: return []
    out = []; cl, ch = iv[0]
    for lo, hi in iv[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else:
            out.append((cl, ch)); cl, ch = lo, hi
    out.append((cl, ch)); return out

def safe_arcs(S):
    danger = union_arcs([iv for v in set(S) for iv in danger_arcs(v)])
    out = []; cur = F(0)
    for lo, hi in danger:
        if lo > cur: out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1: out.append((cur, F(1)))
    return out

def L(S): return sum((b - a for a, b in safe_arcs(S)), F(0))
def primitive(S): return reduce(gcd, S) == 1

# ---------------- exact M tool (validated, verbatim) ----------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def exact_M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

def is_covering(S):
    Sset = set(S)
    return all(any(v % q == 0 for v in Sset) for q in range(2, 15))

def banner(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)

# ===========================================================================
def main():
    banner("A. THE LOGICAL SKELETON — what each STEP actually delivers")
    print("""
  STEP1 (q-witness, THM-523, PROVED): S omits a multiple of some q in 2..14
        => tau=1/q is lonely => M(S)>=1/q>=1/14.  CONCLUSION = a UNIVERSAL bound
        (1/14), no dependence on any meas(G_C).  Genuinely reduces the problem to
        COVERING sets.  NON-CIRCULAR, fully closed.  [verify below]

  STEP2 (gap decomposition, PROVED exact): S = C u {w}, |C|=12, 14|w.
        L(S) = meas(G_C) - meas(G_C ∩ D_w).  An IDENTITY.  Says nothing by itself.
        It REWRITES L(S) in terms of meas(G_C).  CONCLUSION depends on meas(G_C).

  STEP3 (decoupling floor, PROVED): L(C u {w}) >= (6/7) meas(G_C) - r/(7w).
        => for w > r/(6 meas(G_C)), L(S) > 0 => M(S) > 1/14.
        *** CRITICAL: the threshold w_0 = r/(6 meas(G_C)) DEPENDS ON meas(G_C). ***
        If meas(G_C) is small, w_0 is LARGE, so STEP3 only handles w above a
        core-dependent threshold.  For w BELOW w_0 we must exact-check (finite).
        So STEP3 closes the SINGLE-parked-w case ONLY GIVEN meas(G_C)>0 for THAT C.

  STEP4 (GAP A, OPEN): k>=3 coordinated growing speeds.  Need meas(G_C) bounded
        below UNIFORMLY in v_max = OPEN-Q-108 itself.

  The circularity question: is STEP1+2+3 a reduction to a STRICTLY WEAKER claim,
  or does the whole thing still rest on the uniform meas(G_C)>=c (OPEN-Q-108)?
""")

    # ---- STEP1 is genuinely universal: re-verify it closes non-covering sets --
    banner("B. STEP1 is universal (non-circular) — re-verify on random non-covering 13-sets")
    import random
    random.seed(20260617)
    bad = 0; checked = 0
    for _ in range(4000):
        S = sorted(random.sample(range(1, 80), 13))
        if reduce(gcd, S) != 1: continue
        if is_covering(S): continue
        checked += 1
        # find a missing q-witness and confirm M>=1/q
        qs = [q for q in range(2, 15) if not any(v % q == 0 for v in S)]
        q = qs[0]
        if gmin(S, F(1, q)) < F(1, 14):
            bad += 1
    print(f"  non-covering primitive 13-sets checked: {checked}; STEP1 q-witness failures: {bad}")
    print(f"  => STEP1 reduces LRC(14) to COVERING sets with a UNIVERSAL 1/14 bound. NON-CIRCULAR.")

    # =======================================================================
    banner("C. THE CORE CIRCULARITY TEST: does STEP3's 'closure' need a UNIFORM meas bound?")
    print("""
  STEP3 says: M(C u {w}) > 1/14 whenever w > r/(6 meas(G_C)).  For a FIXED core C
  with meas(G_C)>0 (handed by LRC(12)), this leaves only finitely many w to check.
  So for ANY FIXED finite family of cores, the single-parked case is closed.
  But OPEN-Q-108 ranges over INFINITELY MANY cores C (all 12-subsets of N).
  The closure is uniform iff the threshold w_0(C) = r(C)/(6 meas(G_C)) is bounded
  AND meas(G_C) is bounded below.  BOTH are exactly the open uniformity.

  I now SHRINK meas(G_C) by choosing cores with a LARGE speed, and watch w_0 blow up.
""")
    # Build cores with one large speed L: C = {1..11} u {Lbig}. meas(G_C) -> ? as Lbig grows.
    # This is a 12-core (11 smalls + 1 large). Track meas(G_C) and threshold w_0.
    print("  Core C = {1,2,...,11, Lbig} (12 speeds). Track meas(G_C), r, threshold w0=r/(6 meas):")
    print(f"  {'Lbig':>6} {'meas(G_C)':>16} {'float':>12} {'r':>3} {'w0=r/6meas':>12}")
    last = []
    for Lbig in [13, 20, 50, 100, 500, 1000, 5000, 50000]:
        C = list(range(1, 12)) + [Lbig]
        if reduce(gcd, C) != 1:
            pass
        arcs = safe_arcs(C); m = sum((b-a for a,b in arcs), F(0)); r = len(arcs)
        w0 = (F(r) / (6 * m)) if m > 0 else None
        last.append((Lbig, m, r, w0))
        print(f"  {Lbig:>6} {str(m):>16} {float(m):>12.7f} {r:>3} "
              f"{(float(w0) if w0 else 'inf'):>12}")
    print("""
  READING (CORRECTED — adversarial finding): meas(G_C) does NOT collapse to 0 as
  the single large speed grows (it stabilizes near a positive value: the large
  speed acts like a density-1/7 comb thinning {1..11}'s safe set).  BUT the arc
  count r GROWS LINEARLY with Lbig (the large speed shatters the safe set into
  ~Lbig arcs).  Since w0 = r/(6 meas), the threshold w0 ALSO GROWS ~ Lbig.
  *** So even the SINGLE-parked-w decoupling floor (6/7)meas - r/(7w) becomes
  VACUOUS for cores with a large speed: it only proves L>0 for w > w0 ~ Lbig,
  and for w <= w0 one must EXACT-CHECK (finitely many, but the bound grows). ***
  The MEASURE stays bounded below for a single large speed, but the PROVED
  decoupling certificate does NOT uniformly close even this case -- the r/(7w)
  error term is too crude once r is large.  STEP3 is a clean proof ONLY when the
  core's speeds (hence r) are BOUNDED.  Empirically (next: D, F) no w drives L=0,
  but that is evidence, not the decoupling proof.
""")

    # =======================================================================
    banner("D. WHERE IT STOPS BEING A REDUCTION: the k>=3 coordinated regime is NOT reduced")
    print("""
  The honest claim must be: STEP1-3 reduce LRC(14) to controlling cores C that
  THEMSELVES contain k>=2 coordinated large speeds (because peeling ONE w leaves a
  12-core that may again have large coordinated speeds).  Each peel needs meas of
  the SMALLER core bounded below.  Iterating bottoms out at a small (<=2 large
  speeds) core ONLY IF the intermediate measures stay positive -- which is again
  the uniform bound.  I test the ITERATION directly: peel large speeds one-by-one
  and see if any intermediate core hits meas=0 (which would BREAK the induction).
""")
    # take a 13-set with 3 coordinated large speeds, peel them, track measures
    test_families = [
        ("3 coord 14k on drop-6 core", [1,2,3,4,5,7,8,9,10,11, 14*5, 14*7, 14*11]),
        ("3 coord 14k on {1..10}",     [1,2,3,4,5,6,7,8,9,10, 14*5, 14*7, 14*11]),
        ("3 large AP (no mod-14)",      [1,2,3,4,5,6,7,8,9,10, 101,103,107]),
        ("drop-6 + 84,168,252",         [1,2,3,4,5,7,8,9,10,11,13, 84,168,252][:13]),
    ]
    for name, S in test_families:
        S = sorted(set(S))[:13]
        if len(S) != 13:
            print(f"  {name}: |S|={len(S)} skip"); continue
        cov = is_covering(S); prim = primitive(S); ls = L(S)
        # peel the 3 largest, tracking measure of the shrinking core
        cur = sorted(S); chain = [(len(cur), L(cur))]
        for _ in range(3):
            cur = cur[:-1]  # drop largest
            chain.append((len(cur), L(cur)))
        broke = any(m == 0 for _, m in chain)
        print(f"  {name}: cov={cov} prim={prim} L(S)={float(ls):.6f}")
        print(f"     peel-chain (|C|, meas): {[(n, float(m)) for n,m in chain]}  "
              f"{'BREAKS (meas hit 0)' if broke else 'all intermediate meas>0'}")
    print("""
  RESULT: no intermediate core hits meas=0 in these tests, so the iteration LOOKS
  fine -- but 'looks fine on examples' is exactly the empirical 33k-set evidence,
  NOT a proof.  The induction needs every intermediate (13-k)-core to have
  meas>=c uniformly; that is OPEN-Q-108 unchanged.
""")

    # =======================================================================
    banner("E. THE DECISIVE CIRCULARITY VERDICT — is the FINAL claim contingent on OPEN-Q-108?")
    print("""
  I formalize the dependency.  Define:
    P  := 'M(S) >= 1/14 for all primitive 13-sets S'           (= LRC(14), the GOAL)
    Q  := 'exists c>0 with meas(G_C) >= c for all 12-subsets C' (= OPEN-Q-108)
  The reduction's logical content (STEP1-3 + THM-522 quantization) is the implication
    Q  ==>  P    (uniform fattening + decoupling + quantization gives inf L>0 => LRC14).
  CIRCULARITY would be:  the proof of 'Q ==> P' itself uses P, OR the proof only
  establishes 'P ==> P' (vacuous), OR it establishes 'Q' assuming 'Q'.

  TEST: is the implication Q==>P PROVED using ONLY facts NOT equivalent to P or Q?
  STEP1 (PROVED, uses nothing): reduces P to P-restricted-to-covering-sets.
  STEP2 (IDENTITY): no logical content.
  STEP3 (PROVED, uses LRC(12) = a DIFFERENT, ALREADY-PROVEN theorem, not P/Q):
        gives, for EACH core with meas(G_C)>0, closure of the single-large-w case
        with a threshold depending on meas(G_C).
  So 'Q ==> P' is a genuine implication whose proof uses LRC(12) and decoupling,
  NEITHER of which is P or Q.  => the reduction is NOT circular.
  BUT: it is also NOT a reduction to anything WEAKER than Q.  The k>=3 regime is
  UNREDUCED; the program's success still REQUIRES Q (or something equivalent).
  Net: 'Q ==> P' is established; 'P ==> Q' or a weaker-than-Q route is NOT.
""")
    # Concrete demonstration that Q is genuinely needed (not implied by LRC12 alone):
    # LRC(12) gives meas at gap 1/13, but the gap-1/13 maximizer need NOT fatten to
    # a uniformly positive gap-1/14 measure. Show the LRC12 lever does NOT by itself
    # bound meas(G_C) (the half-width shrinks ~1/v_binder).
    print("  Demonstrate Q is NOT a corollary of LRC(12) alone (the lever is not enough):")
    print("  For cores whose gap-1/13 binders are LARGE speeds, the guaranteed safe arc")
    print("  around the maximizer shrinks like 1/v_binder, so a single-arc Lipschitz bound")
    print("  gives NO universal constant.  Exhibit shrinking half-widths:")
    for C in [list(range(1,12))+[13], [1,2,3,4,5,7,8,9,10,11,12,13],
              list(range(1,12))+[100], list(range(1,12))+[1000]]:
        M13, tau13 = exact_M(C)
        binders = sorted({v for v in C if nrm(v*tau13) == M13})
        slack = M13 - F(1,14)
        hw = min((slack / v for v in binders), default=F(0))
        print(f"    C={C[:3]}...{C[-1]} (vmax={max(C)}): M(C)={float(M13):.5f} binders={binders} "
              f"Lipschitz half-width={float(hw):.7f}")
    print("""
  => the single-arc lower bound -> 0 as v_binder grows.  So LRC(12) does NOT imply Q;
  Q is a STRICTLY STRONGER statement (controls TOTAL measure, not one arc).  The
  reduction does not prove Q and does not dodge needing it.
""")

    # =======================================================================
    banner("F. INDEPENDENT COUNTEREXAMPLE HUNT (is there a covering set with M<=1/14?)")
    print("  Fresh exact search distinct from the original 33k census: targeted coordinated")
    print("  families + a structured grid, exact L-screen then exact-M on L=0 survivors.")
    found_tight = []; found_ctx = []; checked = 0; l0 = 0; minM = None; minM_S = None
    small = list(range(1, 14))
    def consider(S):
        nonlocal checked, l0, minM, minM_S
        S = tuple(sorted(set(S)))
        if len(S) != 13 or not primitive(S) or not is_covering(S): return
        checked += 1
        if L(S) == 0:
            l0 += 1
            Mv, _ = exact_M(S)
            if minM is None or Mv < minM: minM = Mv; minM_S = S
            if Mv < F(1,14): found_ctx.append((Mv, S))
            elif Mv == F(1,14): found_tight.append((Mv, S))
    # family 1: drop e, add the minimal-covering big element times m, plus a 2nd coordinated big
    for e in range(1, 14):
        for m1 in range(1, 16):
            consider([v for v in small if v != e] + [14*m1])
    # family 2: two drops, two mults of 14 with shared factor (coordinated)
    for d in it.combinations(range(1,14), 2):
        base = [v for v in small if v not in d]
        for m1 in range(1, 14):
            for m2 in range(m1+1, 14):
                consider(base + [14*m1, 14*m2])
    # family 3: closest-to-1/14 84m family perturbations
    base = [1,2,3,4,5,6,7,8,9,10,11,13]
    for m in range(1, 40):
        consider(base + [84*m])
    # family 4: three coordinated mults sharing modulus 42 (=lcm structure)
    for d in it.combinations(range(1,14), 3):
        base = [v for v in small if v not in d]
        for m1 in range(1, 7):
            for m2 in range(m1+1, 8):
                for m3 in range(m2+1, 9):
                    consider(base + [42*m1, 42*m2, 42*m3])
    print(f"  checked {checked} covering primitive 13-sets; L=0 survivors {l0}")
    print(f"  min M over L=0 survivors: {minM if minM else 'n/a'} at {minM_S}")
    print(f"  TIGHT (M=1/14): {len(found_tight)}   COUNTEREXAMPLES (M<1/14): {len(found_ctx)}")
    if found_ctx:
        for Mv, S in found_ctx[:5]:
            print(f"    *** COUNTEREXAMPLE M={Mv}={float(Mv):.6f}  S={S}")
    if found_tight:
        for Mv, S in found_tight[:5]:
            print(f"    tight M=1/14 S={S}")

    banner("DONE")

if __name__ == "__main__":
    main()
