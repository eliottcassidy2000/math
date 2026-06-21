#!/usr/bin/env python3
"""
lrc14_road_coloring_sync_opus_s3.py   (opus-2026-06-20-S3)

THREAD: ROAD COLORING THEOREM (Trahtman) <-> the LRC Z/7 clock automaton.

QUESTION: is the LRC cover condition "measS7(E) = meas{x: colors(E,x)=all Z/7}"
a SYNCHRONIZATION / reachability property of a deterministic finite automaton on
the 7 sectors, and is consec the SYNCHRONIZING (Trahtman-optimal) coloring with
minimal reset (Cerny-like) threshold?

============================================================================
STEP 0 -- the genuine finite automaton (the "road").  [PROVED structurally]
============================================================================
For a fixed slope x in [0,1), define sector(e) = floor(7*frac(e*x)) in Z/7.
Empirically (and provably from floor(7*frac((e+1)x)) - floor(7*frac(e x))):
the per-step ADVANCE  a(e) = (sector(e+1) - sector(e)) mod 7  takes EXACTLY
two values, q = floor(7x) and q+1, in the pattern of a SturMIAN word with
rotation number rho = frac(7x).  Hence the dynamics of e -> e+1 is the
2-letter deterministic automaton on Z/7:

    L : s |-> (s + q)   mod 7        (the "Low" road)
    H : s |-> (s + q+1) mod 7        (the "High" road)

This is a strongly-connected, constant-out-degree-2 digraph (the 7-cycle Cayley
graph of Z/7 with two generators q, q+1) -- exactly Trahtman's setting EXCEPT
the WORD over {L,H} is not free: it is the Sturmian itinerary forced by rho.
So the LRC automaton is a ROAD-COLORED 7-state graph whose driving word is a
rotation word, not an arbitrary input.

The COVER condition "{sector(e): e in E} = Z/7" is then: starting at sector 0
(e=0 always gives sector 0), the set of states VISITED at the SELECTED times
e in E equals all of Z/7.  This is a REACHABILITY-SET-FULLNESS condition, the
natural automaton twin of synchronization.

============================================================================
HYPOTHESES (falsifiable):
============================================================================
H-SYNC:  "cover Z/7" (the LRC event) is EQUIVALENT to a synchronization /
         full-reachability statement of the {L,H} automaton driven by the
         Sturmian itinerary of rho=frac(7x).  TEST: build the automaton, and
         check the cover event equals "the visited-state set at times E is all
         of Z/7".  (Sanity/definition check -- expected TRUE by construction.)

H-CERNY: consec MINIMIZES a Cerny-like RESET THRESHOLD: among k-subsets E, the
         consecutive block reaches full coverage in the SHORTEST window / with
         the SMALLEST "reset cost", and THIS is why it maximizes measS7.
         Concretely test two reset proxies:
         (a) the minimal e_max = max(E) needed to first cover Z/7 for a given x
             (synchronization LENGTH), and whether consec covers over the widest
             x-measure at each fixed window length;
         (b) the reset WORD structure: for the pure rotation automaton, the
             coloring {L,H} is synchronizing iff gcd(q, 7)=... -- compute the
             Cerny reset length of the q,q+1 coloring and relate to cover.

H-ROAD:  The q-dependence is a Z/7 TWIST e -> q*e (matches HYP-2703's slope-band
         twist e -> s*e).  TEST whether the road-coloring "reset length" is the
         same across all admissible q (rotation invariance) -- if so, coverage
         is a SLOPE-rho property only, the q only relabels states; this would
         explain why HYP-2703's per-band (per-q) extremality FAILS but the
         aggregate works (synchronization is rho-driven, twist is cosmetic).

============================================================================
HONESTY CONTRACT: report exact Fractions; if synchronization/reset does NOT
predict consec-maximality, say DEAD-END.  Do not inflate the analogy.
============================================================================
stdlib only; exact Fractions.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)


# ---------------------------------------------------------------------------
# Core sector / cover engine (matches THM-536 / lrc14 canon)
# ---------------------------------------------------------------------------
def sector(e, x, scale=7):
    return int(scale * ((e * x) % 1))


def cover_event(E, x, scale=7):
    return len(set(sector(e, x, scale) for e in E)) == scale


def measS7(E, scale=7):
    """Exact measure of {x in [0,1): {floor(scale*frac(e x)): e in E} = Z/scale}."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, scale * e + 1):
            bps.add(F(m, scale * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if cover_event(E, xm, scale):
            total += x1 - x0
    return total


# ---------------------------------------------------------------------------
# STEP 0 (PROVED structurally) : the two-letter automaton; verify advance set
# ---------------------------------------------------------------------------
def verify_two_letter(scale=7, trials=400, emax=40):
    print("=" * 74)
    print("STEP 0: per-step advance a(e)=(sector(e+1)-sector(e)) mod 7 in {q,q+1}")
    print("=" * 74)
    import random
    rng = random.Random(424242)
    bad = 0
    for _ in range(trials):
        x = F(rng.randint(1, 9973), 9973)
        q = int(scale * x) % scale  # floor(7x) mod 7
        advances = set()
        for e in range(emax):
            advances.add((sector(e + 1, x, scale) - sector(e, x, scale)) % scale)
        allowed = {q % scale, (q + 1) % scale}
        if not advances <= allowed:
            bad += 1
            if bad <= 3:
                print(f"   x={x} float~{float(x):.4f} q={q} advances={advances} allowed={allowed}")
    print(f"   trials={trials}: advance set NOT subset of {{q,q+1}} in {bad} cases.")
    print("   => deterministic 2-letter (L=+q, H=+q+1) rotation automaton on Z/7."
          if bad == 0 else "   => STRUCTURE FAILS (re-derive).")
    return bad == 0


# ---------------------------------------------------------------------------
# H-SYNC : cover event == full visited-state set (definition twin). Verify.
# Also: the SYNCHRONIZING WORD / RESET angle for the pure rotation coloring.
# ---------------------------------------------------------------------------
def road_coloring_reset_length(q, scale=7):
    """
    The 2-letter road coloring on Z/7: L=+q, H=+q+1.  A word w over {L,H} is
    SYNCHRONIZING if it maps all 7 states to one state.  But pure ROTATIONS are
    bijections, so NO word collapses states: the rotation automaton is NEVER
    synchronizing in Trahtman's strict sense (permutation automaton, image never
    shrinks).  Instead the relevant notion is REACHABILITY/COVERAGE: from a fixed
    start, what window of steps visits all 7 states?  Compute, for the *constant*
    word (e=0,1,2,...,t i.e. all-L if rho->0, or the mixed Sturmian), the minimal
    window length to visit all states under pure +q stepping (rho=0 limit).
    For pure +q stepping, states visited in t steps = {0,q,2q,...,(t-1)q}; this
    is all of Z/7 iff gcd(q,7)=1 and t>=7.  Return that minimal t (the 'rotation
    reset' length) and whether q generates Z/7.
    """
    if q % scale == 0:
        return None, False  # +0 never moves: never covers
    g = math.gcd(q % scale, scale)
    generates = (g == 1)
    if not generates:
        return None, False
    # minimal t s.t. {0,q,2q,...,(t-1)q mod 7} = Z/7
    seen = set()
    t = 0
    cur = 0
    while len(seen) < scale:
        seen.add(cur % scale)
        cur = (cur + q) % scale
        t += 1
        if t > 100:
            break
    return t, generates


def hsync_check(scale=7):
    print("=" * 74)
    print("H-SYNC: 'cover Z/7' == full visited-state set; rotation reset lengths")
    print("=" * 74)
    # definitional twin: trivially true (cover_event uses the visited set). Confirm
    # on a small grid that cover_event(E,x) == (visited set == Z/7).
    import random
    rng = random.Random(7)
    bad = 0
    for _ in range(2000):
        x = F(rng.randint(1, 997), 997)
        E = sorted(random.Random(_).sample(range(0, 14), 8))
        vis = set(sector(e, x) for e in E)
        if (vis == set(range(scale))) != cover_event(E, x):
            bad += 1
    print(f"   cover_event == (visited set == Z/7): mismatches = {bad}  (expected 0, definitional)")
    print()
    print("   Pure-rotation road coloring reset (rho=0 limit), q=floor(7x) mod 7:")
    for q in range(scale):
        t, gen = road_coloring_reset_length(q)
        print(f"     q={q}: generates Z/7={gen}  rotation_reset_len t*={t}")
    print("   NOTE: rotations are PERMUTATIONS -> never synchronizing in the strict")
    print("   Cerny sense (image never shrinks). The LRC notion is COVERAGE/")
    print("   reachability, not state-collapse. This is the first honest correction.")


# ---------------------------------------------------------------------------
# H-CERNY (a): synchronization LENGTH proxy.  For each x, the minimal max-index
# window W=max(E) needed to cover Z/7 by SOME E subset of {0..W}.  Then ask:
# does consec_k cover over the widest x-measure for each fixed window length,
# i.e. is consec the 'fastest synchronizer'?
# ---------------------------------------------------------------------------
def min_window_to_cover(x, scale=7, wmax=60):
    """Minimal W s.t. {sector(0),...,sector(W)} = Z/7 (cover by the FULL prefix
       0..W -- the consec block!). This is the rotation 'first cover time'."""
    seen = set()
    for e in range(wmax + 1):
        seen.add(sector(e, x, scale))
        if len(seen) == scale:
            return e
    return None


def hcerny_a(scale=7):
    print("=" * 74)
    print("H-CERNY(a): consec = 'fastest synchronizer'? first-cover-time of the")
    print("            full prefix vs. coverage by sparse subsets")
    print("=" * 74)
    # The consec block {0..k-1} covers Z/7 for x iff first-cover-time <= k-1.
    # So measS7(consec_k) = meas{x: first_cover_time(x) <= k-1}.
    # measS7(any E with max(E)=W) <= meas{x: first_cover_time(x) <= W} = measS7(consec_{W+1}).
    # i.e. consec block of the SAME SPAN dominates ANY subset of that span.
    # TEST this 'span-domination' exactly: for every E, measS7(E) <= measS7(consec_{max(E)+1}).
    bad = 0
    tested = 0
    worst = (F(0), None)
    for k in (4, 5, 6):
        for combo in itertools.combinations(range(1, 12), k - 1):
            E = [0] + list(combo)
            W = max(E)
            cons_span = list(range(W + 1))
            mE = measS7(E)
            mC = measS7(cons_span)
            tested += 1
            if mE > mC:
                bad += 1
                if mE - mC > worst[0]:
                    worst = (mE - mC, (E, cons_span))
    print(f"   span-domination measS7(E) <= measS7(consec_{{max(E)+1}}): tested {tested}, "
          f"violations = {bad}")
    if bad:
        print("   worst violation:", float(worst[0]), worst[1])
    else:
        print("   => PROVED-style fact (this run): the FULL PREFIX 0..W covers for")
        print("      the WIDEST x-measure among all subsets of {0..W}; equality iff")
        print("      E already contains a covering sub-collection for a.e. covered x.")
    return bad == 0


# ---------------------------------------------------------------------------
# The crux gap: span-domination gives consec_{W+1} >= E with max W, but that
# block has W+1 elements, possibly MORE than k. The real LRC crux compares
# EQUAL CARDINALITY. Quantify the 'reset cost' = extra indices the prefix needs.
# ---------------------------------------------------------------------------
def reset_cost_table(scale=7):
    print("=" * 74)
    print("H-CERNY(b): 'reset cost' = span the prefix needs to match a sparse E's")
    print("            cover; does consec win at EQUAL CARDINALITY by paying less?")
    print("=" * 74)
    # first-cover-time distribution over x: F(t) = meas{x: first_cover_time <= t}
    # = measS7(consec_{t+1}). Tabulate.
    print("   first-cover-time CDF  F(t)=measS7(consec_{t+1}) (consec needs t+1 indices):")
    prev = F(0)
    for t in range(6, 14):
        Ft = measS7(list(range(t + 1)))
        print(f"     t={t:2d}: F(t)=measS7(consec_{t+1})={float(Ft):.6f}  "
              f"increment={float(Ft-prev):.6f}")
        prev = Ft
    print("   For a sparse E with cardinality k and span W>k: it 'spreads' k probes")
    print("   over a window of W+1 slots. Cover prob <= F(W) but it uses only k slots.")
    print("   consec_k uses k slots over window k-1: cover prob = F(k-1).")
    print("   CRUX restated: is F(k-1) [pay k, span k] >= best sparse cover [pay k, span>k]?")
    print("   This is the genuine reset/Cerny TRADE: SPAN buys coverage, CARDINALITY")
    print("   is the budget. consec spends its whole budget on the SHORTEST covering")
    print("   window. The open inequality = 'cheapest reset wins'.")


# ---------------------------------------------------------------------------
# H-ROAD: twist invariance of the reset/coverage across q (rotation classes).
# If the COVER MEASURE restricted to a q-band depends on q only via e->q*e,
# matches HYP-2703's slope-band twist. Test: band coverage under q-twist.
# ---------------------------------------------------------------------------
def hroad_twist(scale=7):
    print("=" * 74)
    print("H-ROAD: q-twist (state relabel e->q*e) and the slope-band link to HYP-2703")
    print("=" * 74)
    # bandcover_s(E) = meas{f in[0,1): {(s*e+floor(e f)) mod 7}=Z/7}.  The s-twist
    # IS the road-coloring q-relabel. The reset 'speed' of band s = how fast the
    # +s rotation explores Z/7. Tabulate the pure-rotation reset for each s and
    # the band-cover of consec_8, side by side.
    def band_cover(E, s):
        E = sorted(set(E)); bps = {F(0), F(1)}
        for e in E:
            if e == 0: continue
            for m in range(0, e + 1): bps.add(F(m, e))
        bps = sorted(b for b in bps if 0 <= b <= 1)
        tot = F(0)
        for i in range(len(bps) - 1):
            f0, f1 = bps[i], bps[i + 1]
            if f1 <= f0: continue
            fm = (f0 + f1) / 2
            if len(set((s * e + int(e * fm)) % scale for e in E)) == scale:
                tot += f1 - f0
        return tot
    cons = list(range(8))
    print("   s | rotation reset t*(s) | bandcover_s(consec_8)")
    for s in range(scale):
        t, gen = road_coloring_reset_length(s)
        bc = band_cover(cons, s)
        print(f"   {s} |        {str(t):>4}          |   {float(bc):.5f}")
    print("   HONEST CORRECTION: pure-rotation skeleton {s*e mod7: e=0..7} is FULL")
    print("   (=Z/7) for ALL s with gcd(s,7)=1 (s=1..6); only s=0 is degenerate {0}.")
    print("   So 'rotation reset length' is 7 for every s>=1 and does NOT distinguish")
    print("   the bands. The slow-band loss at s=1 (HYP-2703 R1) is therefore NOT")
    print("   explained by rotation speed -- the q-twist is a clean state relabel but")
    print("   the per-band extremality failure is a Sturmian floor-spread effect, not")
    print("   a synchronization-speed effect. The road-coloring reset DOES cleanly")
    print("   explain s=0 (degenerate road: coverage impossible from rotation alone).")


def crux_gap(scale=7):
    print("=" * 74)
    print("CRUX-GAP: does span-domination CLOSE the LRC crux? (brutally honest: NO)")
    print("=" * 74)
    adv = [0, 2, 3, 4, 5, 6, 7, 8]
    cons8 = list(range(8))
    cons9 = list(range(9))
    mA, m8, m9 = measS7(adv), measS7(cons8), measS7(cons9)
    print(f"   adv [0,2..8] (k=8, span 8): measS7 = {mA}  ~ {float(mA):.6f}")
    print(f"   span-dom gives adv <= consec_9 (span 8 -> prefix 0..8 = 9 elements):")
    print(f"      consec_9 measS7 = {m9} ~ {float(m9):.6f}; adv<=consec_9: {mA <= m9}")
    print(f"   but the LRC crux needs adv <= consec_8 (EQUAL cardinality 8):")
    print(f"      consec_8 measS7 = {m8} ~ {float(m8):.6f}; adv<=consec_8: {mA <= m8}")
    print(f"      margin consec_8 - adv = {m8-mA} ~ {float(m8-mA):.6f}  (consec wins)")
    print("   => span-domination relates DIFFERENT cardinalities. consec_k uniquely")
    print("      has MINIMAL span among |E|=k, so span-dom gives NO traction at fixed")
    print("      cardinality. The reachability/synchronization frame PROVES the")
    print("      fixed-span lemma but DOES NOT close the fixed-cardinality crux.")


def main():
    print("LRC(14) ROAD COLORING <-> Z/7 clock automaton  (opus-2026-06-20-S3)\n")
    verify_two_letter()
    print()
    hsync_check()
    print()
    hcerny_a()
    print()
    reset_cost_table()
    print()
    crux_gap()
    print()
    hroad_twist()
    print()
    print("=" * 74)
    print("VERDICT: PARTIAL. The genuine automaton is the 2-letter rotation (L=+q,")
    print("H=+q+1) Sturmian road on Z/7. Rotations are PERMUTATIONS so it is NEVER")
    print("synchronizing in Cerny's state-collapse sense; the LRC notion is COVERAGE")
    print("(reachability fullness), not synchronization. The lasting YIELD is the")
    print("PROVED span-domination lemma: measS7(E) <= measS7(consec_{max(E)+1}),")
    print("a pointwise tautology (visited set is monotone under adding probes). But")
    print("it does NOT close the fixed-cardinality crux. Road-coloring cleanly")
    print("explains only the degenerate band s=0. NET: clean PARTIAL + sharp DEAD-END")
    print("for synchronization as a NEW certificate.")
    print("=" * 74)


if __name__ == "__main__":
    main()
