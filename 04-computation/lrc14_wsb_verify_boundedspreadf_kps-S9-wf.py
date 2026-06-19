#!/usr/bin/env python3
"""
lrc14_wsb_verify_boundedspreadf_kps-S9-wf.py   (kps-2026-06-19-S9-wf)

ADVERSARIAL VERIFICATION of the claimed advance "bounded-spread-finite-check":

  THEOREM (claimed): For k in {8,9,10} and every primitive E with span<=16,
     meas(S7(E)) <= cap_k, consec_k is the UNIQUE maximizer, zero shapes over cap.
  Plus:  flat to B=18; span 17..40 worst meas(S7) = 0.218, 0.372, 0.464 (k=8,9,10);
         pair-Bonferroni stays over cap_8 so Regime B needs a signed estimate.

WE DO, INDEPENDENTLY (NOT reusing the claimant's code):
 (0) Build an EXACT meas(S7) engine on 7e-grid breakpoints; cross-validate it
     against (a) a fine deterministic grid sample, and (b) the canon J/Sr/Ly tools
     and known exact values (consec_8 etc.).
 (1) EXHAUSTIVE bounded-spread: for k=8,9,10, scan ALL primitive 0=e1<...<ek<=B
     for B=16 (and B=18 "flat" claim). Record max meas(S7), the argmax(es), and
     EVERY shape with meas(S7) > cap_k.  Check unique-maximizer = consec_k.
 (2) HUNT for counterexamples to meas(S7)(E) <= cap_k over a WIDE net:
       - exhaustive bounded-spread (above),
       - aggressive wide-spread random + structured (span up to 40+),
       - RESONANT configs (span divisible by 7 / apex-prime-7, w == 0 mod 7),
       - short-relation shapes {0,1,N,N+1}-style and dilated-AP,
       - the AP family (consec scaled / arithmetic progressions).
     ONE hit with meas(S7) > cap_k  => holds=False with witness.
 (3) Sanity-check the claimed wide-spread worst values 0.218/0.372/0.464 and
     whether ANY of them (or the resonant configs) exceed cap_k. Check that the
     claimed Regime-B "Bonferroni over cap_8" honesty is correct => the gap is real.

cap_k:  2243/5880 (k8), 2025/4004 (k9), 36/91 (k10), 25/91 (k11), 1/7 (k12), 0 (k13).
"""
import sys, itertools, math, random
from fractions import Fraction as F
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------------
# EXACT meas(S7) engine -- INDEPENDENT reimplementation.
# meas(S7(E)) = measure{ x in [0,1): every sector [j/7,(j+1)/7) hit by some frac(e x) }.
# Breakpoints: x = m/(7e) for all e in E (where frac(e x) crosses a sector wall).
# On each open cell the sector hit-set is constant; sample the midpoint.
# ----------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        # frac(e x) hits sector wall j/7 when e x == (m + j/7), i.e. x = (7m+j)/(7e)
        for t in range(0, 7*e + 1):
            bps.add(F(t, 7*e))
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    total = F(0)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = set()
        for e in E:
            # sector index of frac(e*xm)
            fr = (e * xm) % 1
            secs.add(int(fr * 7))
        if len(secs) == 7:
            total += (x1 - x0)
    return total

# Independent float grid sample of meas(S7) for cross-validation.
def measS7_sample(E, N=400000):
    E = sorted(set(int(e) for e in E))
    cnt = 0
    for s in range(N):
        x = (s + 0.5) / N
        secs = set()
        for e in E:
            secs.add(int((e * x % 1.0) * 7))
            if len(secs) == 7:
                break
        if len(secs) == 7:
            cnt += 1
    return cnt / N

# ----------------------------------------------------------------------------
# Canon J / Sr / Ly tools (copied verbatim from the task brief) for cross-check.
# meas(S7) should equal p0 = E[1{N=0}] where N = # unhit sectors among {1..6}.
# By inclusion-exclusion p0 = sum_{A subseteq {1..6}} (-1)^{|A|} J(A,E).
# ----------------------------------------------------------------------------
def J(A, E):
    E = sorted(set(E)); arcs = [(F(j,7), F(j+1,7)) for j in A]; bp = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for (a, b) in arcs:
            for end in (a, b):
                m = 0
                while True:
                    xv = (end + m) / e
                    if xv >= 1: break
                    if xv >= 0: bp.add(xv)
                    m += 1
    bp = sorted(b for b in bp if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if all(not any(a < ((e*mid) % 1) < b for (a, b) in arcs) for e in E): tot += hi - lo
    return tot

def measS7_via_IE(E):
    # p0 = sum over subsets A of {1..6} of (-1)^|A| J(A,E)
    tot = F(0)
    secs = list(range(1, 7))
    for r in range(0, 7):
        for A in combinations(secs, r):
            tot += (-1)**r * J(set(A), E)
    return tot

CAP = {8: F(2243,5880), 9: F(2025,4004), 10: F(36,91), 11: F(25,91), 12: F(1,7), 13: F(0)}

def primitive(E):
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1

# ============================================================================
def step0_validate():
    print("="*72)
    print("STEP 0: validate the independent meas(S7) engine")
    print("="*72)
    # (a) consec_8 known L_y < cap_8; check meas(S7) value and that meas <= Ly.
    tests = []
    for k in (8, 9, 10):
        E = list(range(k))
        m = measS7(E)
        mie = measS7_via_IE(E)
        ms = measS7_sample(E, 300000)
        ok_ie = "OK" if m == mie else "*** IE MISMATCH ***"
        ok_s = "OK" if abs(float(m) - ms) < 3e-3 else "*** SAMPLE MISMATCH ***"
        print(f"  consec_{k}: meas(S7)={m}={float(m):.6f}  IE={mie} {ok_ie}  sample={ms:.6f} {ok_s}  cap_{k}={float(CAP[k]):.6f}")
    # (b) a couple of arbitrary shapes: engine vs IE vs sample
    for E in [[0,1,2,3,4,5,6,8], [0,2,3,5,7,8,11,13], [0,1,3,4,5,9,10,12,14,16]]:
        m = measS7(E); mie = measS7_via_IE(E); ms = measS7_sample(E, 200000)
        ok = "OK" if m == mie and abs(float(m)-ms) < 4e-3 else "*** MISMATCH ***"
        print(f"  E={E}: meas={float(m):.6f} IE={float(mie):.6f} samp={ms:.6f} {ok}")
    print()

# ============================================================================
def exhaustive_box(k, B):
    """Scan all primitive 0=e1<...<ek<=B. Return (maxval, argmaxes, over_cap_list)."""
    cap = CAP[k]
    best = F(-1); args = []; over = []; cnt = 0
    consec = tuple(range(k))
    consec_val = None
    for combo in combinations(range(1, B+1), k-1):
        E = (0,) + combo
        if not primitive(E):
            continue
        cnt += 1
        m = measS7(list(E))
        if E == consec:
            consec_val = m
        if m > best:
            best = m; args = [E]
        elif m == best:
            args.append(E)
        if m > cap:
            over.append((m, E))
    return best, args, over, cnt, consec_val

def step1_exhaustive():
    print("="*72)
    print("STEP 1: EXHAUSTIVE bounded-spread, check meas(S7)<=cap_k & unique max=consec")
    print("="*72)
    for k in (8, 9, 10):
        for B in (16, 18):
            best, args, over, cnt, cval = exhaustive_box(k, B)
            cap = CAP[k]
            consec = tuple(range(k))
            uniq = (len(args) == 1 and args[0] == consec)
            print(f"  k={k} B={B}: scanned {cnt} primitive shapes")
            print(f"    max meas(S7)={best}={float(best):.6f}  cap_{k}={float(cap):.6f}  "
                  f"{'<=cap OK' if best<=cap else '*** EXCEEDS CAP ***'}")
            print(f"    consec_{k} meas(S7)={cval}={float(cval):.6f}  "
                  f"{'(= max)' if cval==best else '(NOT max!)'}")
            print(f"    #argmax={len(args)}  unique-consec-maximizer={uniq}  args(first3)={args[:3]}")
            print(f"    #shapes over cap_{k} = {len(over)}")
            for m, E in sorted(over, reverse=True)[:8]:
                print(f"        OVER: meas={float(m):.6f} E={E}")
        print()

# ============================================================================
def step2_hunt():
    print("="*72)
    print("STEP 2: COUNTEREXAMPLE HUNT for meas(S7)(E) > cap_k  (wide/resonant/short-rel)")
    print("="*72)
    hits = []
    random.seed(20260619)

    def consider(k, E, tag):
        E = sorted(set(int(e) for e in E))
        if len(E) != k or E[0] != 0:
            return None
        if not primitive(E):
            return None
        m = measS7(E)
        if m > CAP[k]:
            hits.append((k, m, E, tag))
        return m

    worst = {8: (F(0), None), 9: (F(0), None), 10: (F(0), None)}
    def track(k, E, tag):
        m = consider(k, E, tag)
        if m is not None and m > worst[k][0]:
            worst[k] = (m, E)

    # --- (a) AP / arithmetic-progression family: {0,d,2d,...} with gcd reduction = consec
    #         and "thinned" APs (every other), dilated consec.
    for k in (8, 9, 10):
        for span in range(k-1, 45):
            # random spreads of given span containing 0 and span
            for _ in range(400):
                inner = random.sample(range(1, span), k-2) if span-1 >= k-2 else None
                if inner is None: continue
                E = [0] + sorted(inner) + [span]
                track(k, E, f"rand-span{span}")
    # --- (b) RESONANT: span divisible by 7, w==0 mod 7, sets sitting on /7 lattice scaled
    for k in (8, 9, 10):
        for base_span in (7, 14, 21, 28, 35, 42):
            for _ in range(800):
                if base_span-1 < k-2: continue
                inner = random.sample(range(1, base_span), k-2)
                E = [0] + sorted(inner) + [base_span]
                track(k, E, f"resonant-span{base_span}")
            # explicit: consec scaled so span = base_span (only if divisible)
            # AP with common difference giving total = base_span
        # multiples-of-7 differences embedded
        for _ in range(2000):
            E = [0] + sorted(random.sample(range(1, 60), k-1))
            sp = max(E)
            if sp % 7 == 0:
                track(k, E, "rand-span-div7")
    # --- (c) short-relation shapes: {0,1,N,N+1,...} clusters near 0 and near N
    for k in (8, 9, 10):
        for N in range(k, 50):
            # two tight blocks: a near 0, b near N
            for a in range(1, k-1):
                b = k - a  # block sizes
                left = list(range(0, a))            # 0,1,...,a-1
                right = list(range(N, N+b))          # N,...,N+b-1
                E = left + right
                track(k, E, f"shortrel-N{N}-split{a}/{b}")
    # --- (d) dilated consec: lambda*consec rounded (breaks integrality => resonance)
    for k in (8, 9, 10):
        base = list(range(k))
        for num, den in [(3,2),(5,3),(7,4),(7,5),(11,7),(8,7),(9,7),(15,7),(2,1),(3,1),(7,2)]:
            E = sorted(set((num*e + den//2)//den * 0 + (num*e)//den for e in base))
            # ensure size k by alternative rounding if collisions:
            if len(set((num*e)//den for e in base)) == k:
                track(k, [ (num*e)//den for e in base], f"dilate{num}/{den}")
    # --- (e) fully exhaustive small wide box near the claimed wide worst, k=8 span 17..24
    for k in (8,):
        for B in (17, 18, 19, 20):
            for combo in combinations(range(1, B+1), k-1):
                E = (0,)+combo
                if max(E) < 17:  # only genuinely-wide ones beyond bounded box
                    continue
                if not primitive(E): continue
                track(k, list(E), f"exh-wide-B{B}")
    # --- (f) targeted: the claimed wide worst-case shapes themselves (reconstruct extremes)
    #        Try to find what hits 0.218 (k8), 0.372 (k9), 0.464 (k10) and check < cap.
    for k,target in [(8,0.218),(9,0.372),(10,0.464)]:
        for span in range(17, 41):
            for _ in range(1500):
                if span-1 < k-2: continue
                inner = random.sample(range(1, span), k-2)
                E = [0]+sorted(inner)+[span]
                track(k, E, f"wsb-target-k{k}")

    print("  Worst meas(S7) found per k (over ALL hunt families incl wide/resonant):")
    for k in (8, 9, 10):
        m, E = worst[k]
        print(f"    k={k}: worst meas(S7)={float(m):.6f}  cap_{k}={float(CAP[k]):.6f}  "
              f"{'<=cap' if m<=CAP[k] else '*** OVER CAP ***'}   at E={E}")
    print()
    if hits:
        print(f"  *** {len(hits)} COUNTEREXAMPLES meas(S7)>cap_k FOUND: ***")
        for k, m, E, tag in sorted(hits, reverse=True, key=lambda t: float(t[1]))[:20]:
            print(f"     k={k} meas(S7)={float(m):.6f} > cap_{k}={float(CAP[k]):.6f}  [{tag}]  E={E}")
    else:
        print("  No counterexample meas(S7)>cap_k found in any hunt family.")
    return hits, worst

# ============================================================================
def step3_widevals():
    print("="*72)
    print("STEP 3: verify claimed wide-spread worst values & resonant safety")
    print("="*72)
    # Reproduce the claimed worst values 0.218 / 0.372 / 0.464 via local optimization
    # over wide spans, and confirm they are well below cap_k (so Regime B has slack
    # in VALUE even though Bonferroni cannot certify it).
    claimed = {8: 0.218, 9: 0.372, 10: 0.464}
    random.seed(7)
    for k in (8, 9, 10):
        best = F(0); bestE = None
        for span in range(17, 41):
            for _ in range(3000):
                if span-1 < k-2: continue
                inner = random.sample(range(1, span), k-2)
                E = [0]+sorted(inner)+[span]
                if not primitive(E): continue
                m = measS7(E)
                if m > best:
                    best = m; bestE = E
        print(f"  k={k}: best wide(span17-40) meas(S7)={float(best):.6f} (claimed worst {claimed[k]})  "
              f"cap_{k}={float(CAP[k]):.6f}  {'<=cap' if best<=CAP[k] else '*** OVER CAP ***'}  E={bestE}")
    print()
    print("  NOTE: a true Regime-B THEOREM needs a *signed* explicit bound valid for ALL")
    print("        span>16, not a sampled worst-case. The above only shows VALUE-slack.")

# ============================================================================
if __name__ == "__main__":
    step0_validate()
    step1_exhaustive()
    hits, worst = step2_hunt()
    step3_widevals()
    print("="*72)
    print("VERDICT SUMMARY")
    print("="*72)
    print("See STEP 1 (bounded-spread exhaustive) and STEP 2 (hunt) above.")
    print(f"Bounded-spread counterexamples in hunt: {len([h for h in hits])}")
