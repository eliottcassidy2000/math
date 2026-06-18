"""
Adversarial verification of the 'lrc-lever-overlap' claimed advance on LRC(14) case S3.

We verify EVERY load-bearing claim from scratch with exact Fraction arithmetic:
  (A) The arc-width lemma implication C(S) => M(S)>=1/14  (sanity check the criterion machinery).
  (B) LEMMA L_collapse (single-gap cluster band) -- re-derive the geometry exactly.
  (C) LEMMA L_ruler (periodic ruler) -- re-derive exactly.
  (D) WIDTH-INEQ claim: 11*w0 > wT => width(J_0) > thresh, and width(J_0)>thresh on tight clusters.
  (E) The headline claim: C(S) holds for EVERY covering S3 set (via SOME removed runner & SOME ruler member).
  (F) The 'v=Vmax rule' (G2 counterexample), and whether 'remove the max' ever genuinely fails.
  (G) HUNT for a covering S3 set with M(S) < 1/14 (a true LRC(14) counterexample) -- the real target.
  (H) HUNT for a covering S3 set where C(S) FAILS for ALL choices of removed runner (criterion gap).

Author: kind-pasteur-S3-wf  (adversarial verifier)
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

# ---------------- exact tools (copied verbatim from prompt) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

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

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

# ---------------- criterion C(S) (the PROVED implication target) ----------------
def thresh(v): return F(1, 7*v)

def C_via(S, v):
    """Does removing runner v give widest arc > 1/(7v)?  (the arc-width criterion piece)"""
    A = [x for x in S if x != v]
    return Wwidth(A) > thresh(v)

def C_holds(S):
    """C(S): exists v in S with W(S\{v}) > 1/(7v). Returns (bool, list_of_witness_v)."""
    wits = [v for v in S if C_via(S, v)]
    return (len(wits) > 0, wits)

# ---------------- case classification ----------------
def case_of(S):
    S = sorted(set(S))
    Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'  # k>=2 and Vmax>=13*Vmin

# ============================================================================
# (A) sanity: criterion machinery + Mval agree that C(S)=>M>=1/14 on samples
# ============================================================================
def test_A_criterion_implies_M():
    print("="*70)
    print("(A) sanity: when C(S) holds, is M(S) >= 1/14?  (PROVED implication)")
    fails = 0; checked = 0
    rng = random.Random(1)
    for _ in range(400):
        S = gen_random_S3(rng)
        if S is None: continue
        ok, wits = C_holds(S)
        if ok:
            checked += 1
            m = Mval(S)
            if m < F(1, 14):
                fails += 1
                print("  VIOLATION of implication:", S, "M=", m)
    print(f"  checked {checked} S3 sets with C(S) true; implication failures: {fails}")
    return fails == 0

# ============================================================================
# generators for covering S3 sets
# ============================================================================
def gen_random_S3(rng, max_attempts=200):
    """Generate a random primitive covering 13-set in case S3."""
    for _ in range(max_attempts):
        # small part: pick a covering-friendly base, then a large cluster
        # strategy: pick 11-12 small speeds from 1..13 plus a cluster of large speeds
        nsmall = rng.randint(8, 12)
        small = sorted(rng.sample(range(1, 14), nsmall))
        nlarge = 13 - nsmall
        if nlarge < 2: continue
        # cluster around V0 with small spread
        V0 = rng.randint(14, 400)
        spread = rng.randint(14, 60)
        large = set()
        tries = 0
        while len(large) < nlarge and tries < 100:
            large.add(V0 + rng.randint(0, spread)); tries += 1
        large = sorted(large)
        if len(large) < nlarge: continue
        S = sorted(set(small) | set(large))
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if not is_covering(S): continue
        if case_of(S) != 'S3': continue
        return S
    return None

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def make_covering_small(rng):
    """Return a covering-compatible small part (subset of 1..13) that covers as much of 2..14 as possible."""
    # we need full set to cover 2..14; 14 covered by a multiple of 14 (in the large cluster usually) or 7&2
    pass

# ============================================================================
# (B) re-derive L_collapse geometry EXACTLY
#     J_j = [(14j+1)/(14 w0), (14j+13)/(14 wT)]; claim: tau in J_j => all u in L' safe
# ============================================================================
def Lcollapse_window(w0, wT, j):
    lo = F(14*j+1, 14*w0)
    hi = F(14*j+13, 14*wT)
    return (lo, hi)

def test_B_Lcollapse():
    print("="*70)
    print("(B) re-derive L_collapse: tau in J_j  =>  every u in [w0,wT] is 1/14-safe?")
    # Claim in proof: for tau in J_j and w0<=u<=wT, u*tau in [w0*lo, wT*hi] subset [j+1/14, j+13/14].
    # Check: w0*lo = w0*(14j+1)/(14 w0) = (14j+1)/14 = j + 1/14.  GOOD.
    #        wT*hi = wT*(14j+13)/(14 wT) = (14j+13)/14 = j + 13/14. GOOD.
    # But the proof needs: for EVERY u in [w0,wT], u*tau in [j+1/14, j+13/14].
    #   u*tau >= w0*lo only if tau>=0 (yes) and u>=w0 -> u*tau>=w0*tau>=w0*lo? NO: u>=w0 gives u*tau>=w0*tau, and tau>=lo gives w0*tau>=w0*lo. OK u*tau>=w0*lo=j+1/14. GOOD lower.
    #   u*tau <= wT*tau <= wT*hi = j+13/14.  GOOD upper.
    # So frac(u*tau) lands in [j+1/14,j+13/14] mod 1 ONLY if the whole interval [u*lo, u*hi] stays within ONE period [j, j+1].
    #   The geometry above shows u*tau in [j+1/14, j+13/14] as a REAL number (not mod 1) for all u in [w0,wT]. Since that real interval is within (j, j+1), frac is well-defined and ||u tau||>=1/14. GOOD.
    # CAVEAT: this requires the window nonempty (lo<hi) AND that intermediate u keep u*tau in the SAME [j+1/14,j+13/14]; verified above the bound is uniform. So L_collapse is geometrically SOUND.
    # Now numerically confirm on adversarial w0,wT,j.
    fails = 0; checks = 0
    rng = random.Random(7)
    for _ in range(3000):
        w0 = rng.randint(14, 500); wT = w0 + rng.randint(0, 80)
        j = rng.randint(0, 5)
        lo, hi = Lcollapse_window(w0, wT, j)
        if lo >= hi: continue
        checks += 1
        # sample many u in [w0,wT] and tau in (lo,hi); check safety
        for u in range(w0, wT+1):
            # test endpoints and midpoint of window
            for tau in (lo, (lo+hi)/2, hi):
                if F(0) < tau:  # tau>0
                    if nrm(u*tau) < F(1,14) - F(1, 10**9):  # strict, with tiny slack
                        # exact check
                        if nrm(u*tau) < F(1,14):
                            fails += 1
    print(f"  checked {checks} windows; safety violations: {fails}")
    print("  ANALYTIC: w0*lo=j+1/14, wT*hi=j+13/14 exactly; uniform bound for all u in [w0,wT]. SOUND.")
    return fails == 0

# ============================================================================
# (C) re-derive L_ruler EXACTLY: I_j^r = ((14j+1)/(14r),(14j+13)/(14r)), period 1/r
# ============================================================================
def test_C_Lruler():
    print("="*70)
    print("(C) re-derive L_ruler: safe arcs of single runner r are I_j^r, period 1/r")
    # ||r tau|| >= 1/14  <=>  frac(r tau) in [1/14, 13/14]  <=>  r tau mod 1 in [1/14,13/14]
    #   <=> tau in union_j [(j+1/14)/r, (j+13/14)/r] = ((14j+1)/(14r),(14j+13)/(14r)). period 1/r. CORRECT.
    fails = 0
    rng = random.Random(11)
    for _ in range(2000):
        r = rng.randint(14, 300); j = rng.randint(0, r-1)
        lo = F(14*j+1, 14*r); hi = F(14*j+13, 14*r)
        for tau in (lo, (lo+hi)/2, hi):
            if nrm(r*tau) < F(1,14):
                fails += 1
    print(f"  single-runner safe-arc formula violations: {fails}")
    print("  ANALYTIC: exact iff. SOUND.")
    return fails == 0

# ============================================================================
# (D) WIDTH-INEQ:  claim 11*w0 > wT  =>  width(J_0) > thresh(Vmax)
#     width(J_0) = 13/(14 wT) - 1/(14 w0) = (13 w0 - wT)/(14 w0 wT)
#     thresh = 1/(7 Vmax)
#     width>thresh  <=>  (13 w0 - wT)/(14 w0 wT) > 1/(7 Vmax)
#                   <=>  (13 w0 - wT) * Vmax > 2 w0 wT      (the stated WIDTH-INEQ)
#     Since Vmax>=wT (Vmax is the max, wT max of L'), suffices (13 w0 - wT)*wT > 2 w0 wT
#                   <=>  13 w0 - wT > 2 w0  <=>  11 w0 > wT.   CHECK the algebra & hunt failures.
# ============================================================================
def test_D_widthineq():
    print("="*70)
    print("(D) WIDTH-INEQ algebra + adversarial hunt")
    # algebra check on randoms
    bad_algebra = 0
    rng = random.Random(13)
    for _ in range(20000):
        w0 = rng.randint(14, 1000); wT = w0 + rng.randint(0, 200); Vmax = wT + rng.randint(0, 500)
        width = F(13,14)/wT - F(1,14)/w0  # = 13/(14 wT) - 1/(14 w0)
        th = thresh(Vmax)
        lhs_gt = width > th
        ineq = (13*w0 - wT)*Vmax > 2*w0*wT
        if lhs_gt != ineq:
            bad_algebra += 1
    print(f"  width(J_0)>thresh  <=>  (13w0-wT)Vmax>2w0wT  : mismatches = {bad_algebra}")
    # sufficiency: 11 w0 > wT and Vmax>=wT  => WIDTH-INEQ
    bad_suff = 0
    for _ in range(20000):
        w0 = rng.randint(14, 1000); wT = rng.randint(14, 11*w0 - 1); Vmax = wT + rng.randint(0, 1000)
        if 11*w0 > wT and Vmax >= wT:
            if not ((13*w0 - wT)*Vmax > 2*w0*wT):
                bad_suff += 1
    print(f"  '11 w0 > wT and Vmax>=wT => WIDTH-INEQ' counterexamples = {bad_suff}")
    # But: is width(J_0)>thresh ENOUGH for L_collapse? No -- J_0 must ALSO lie in a P-safe arc.
    print("  NOTE: width(J_0)>thresh is necessary-not-sufficient for L_collapse (needs P-safety too).")
    return bad_algebra == 0 and bad_suff == 0

# ============================================================================
# (E) HEADLINE: does C(S) hold for every covering S3 set (some removed runner)?
#     We search broadly and report any S3 set where C(S) FAILS entirely.
#     Plus check the known exception set.
# ============================================================================
KNOWN_EXC = [1,2,3,5,7,8,9,10,11,12,13,27,28]

def test_E_headline_C(num=2500):
    print("="*70)
    print("(E) HEADLINE: C(S) for every covering S3 set? hunting failures of C(S) (ALL v)")
    rng = random.Random(99)
    Cfail = []; checked = 0
    for _ in range(num):
        S = gen_random_S3(rng)
        if S is None: continue
        checked += 1
        ok, wits = C_holds(S)
        if not ok:
            Cfail.append(S)
    # check the known exception explicitly
    okx, witx = C_holds(KNOWN_EXC)
    print(f"  checked {checked} S3 sets; C(S)-total-failures: {len(Cfail)}")
    if Cfail:
        for S in Cfail[:10]:
            print("   C FAILS:", S, " M=", Mval(S))
    print(f"  known-exception {KNOWN_EXC}: C holds? {okx}, witnesses {witx}")
    return Cfail

# ============================================================================
# (F) the 'remove the max' (v=Vmax) rule: how often does it FAIL while C still holds?
# ============================================================================
def test_F_vmax_rule(num=2000):
    print("="*70)
    print("(F) v=Vmax rule: how often does C_via(S,Vmax) FAIL (but C holds via other v)?")
    rng = random.Random(123)
    vmax_fail = 0; total = 0; both_fail = []
    margins = []
    for _ in range(num):
        S = gen_random_S3(rng)
        if S is None: continue
        total += 1
        Vmax = max(S)
        A = [x for x in S if x != Vmax]
        ratio = Wwidth(A) * 7 * Vmax  # >1 means C_via(Vmax) holds
        margins.append(ratio)
        if ratio <= 1:
            vmax_fail += 1
            ok, wits = C_holds(S)
            if not ok:
                both_fail.append(S)
    mn = min(margins) if margins else None
    mx = max(margins) if margins else None
    print(f"  checked {total}; v=Vmax FAILS on {vmax_fail} sets; margin range [{mn},{mx}]")
    print(f"  sets where BOTH Vmax fails AND C fails entirely: {len(both_fail)}")
    for S in both_fail[:10]:
        print("   ", S, "M=", Mval(S))
    # known exc
    Sx = KNOWN_EXC; Vx = max(Sx); Ax=[x for x in Sx if x!=Vx]
    print(f"  known-exc Vmax={Vx}: ratio = {float(Wwidth(Ax)*7*Vx):.4f} (claim ~0.94, <1)")
    return both_fail

# ============================================================================
# (G) THE REAL TARGET: hunt a covering S3 set with M(S) < 1/14 (true counterexample to LRC14)
#     Broad search + adversarial near-tight / structured sets.
# ============================================================================
def test_G_hunt_M(num=4000):
    print("="*70)
    print("(G) HUNT: covering S3 set with M(S) < 1/14  (would REFUTE LRC(14))")
    rng = random.Random(2024)
    worst = None; worstM = F(1)
    below = []
    checked = 0
    for _ in range(num):
        S = gen_random_S3(rng)
        if S is None: continue
        checked += 1
        m = Mval(S)
        if m < worstM:
            worstM = m; worst = S
        if m < F(1, 14):
            below.append((S, m))
    print(f"  checked {checked} S3 sets; min M = {worstM} = {float(worstM):.5f} at {worst}")
    print(f"  sets with M < 1/14: {len(below)}")
    for S, m in below[:10]:
        print("    COUNTEREXAMPLE?", S, "M=", m, float(m))
    return below, worst, worstM

# ============================================================================
# (H) adversarial structured S3 families: APs, harmonic small parts, boundary clusters,
#     wide clusters (wT/w0>=11), degenerate L'={single}, large Vmax.
# ============================================================================
def covering_check_and_M(S):
    S = sorted(set(S))
    if len(S) != 13: return None
    if not is_covering(S): return None
    return S

def test_H_structured():
    print("="*70)
    print("(H) adversarial structured S3 families")
    results = []
    # family 1: small={1..11}+two larges 14m, 14m' tight
    bases_small = [
        [1,2,3,4,5,6,7,8,9,10,11],
        [1,2,3,4,5,6,7,8,9,11,13],
        [1,2,3,5,6,7,9,10,11,12,13],
        [2,3,4,5,6,7,9,10,11,12,13],
        [1,2,3,4,6,7,8,9,11,12,13],
    ]
    cluster_templates = [
        # (anchor, offsets) -- ensure 14|something for covering
        (14, [0,1]), (28,[0,1]), (14,[0,3]), (98,[0,5]),
        (14,[0,14]),  # wide
        (140,[0,1]), (140,[0,40]), (280,[0,1]),
    ]
    rng = random.Random(5)
    for small in bases_small:
        for anchor, offs in cluster_templates:
            large = [anchor+o for o in offs]
            S = covering_check_and_M(small+large)
            if S is None: continue
            if case_of(S) != 'S3': continue
            m = Mval(S)
            ok, wits = C_holds(S)
            results.append((S, m, ok, wits))
    # family 2: arithmetic-progression clusters with various V0
    for V0 in [14, 30, 70, 140, 300, 700, 1400]:
        for step in [1, 2, 3, 14]:
            for nlarge in [2,3,4,5]:
                large = [V0 + step*i for i in range(nlarge)]
                # need a multiple of 14 for covering -> include if V0 multiple of 14 else add small 7&2
                nsmall = 13 - nlarge
                small = list(range(1, 1+nsmall))
                S = covering_check_and_M(sorted(set(small+large)))
                if S is None: continue
                if len(S)!=13: continue
                if case_of(S) != 'S3': continue
                m = Mval(S)
                ok, wits = C_holds(S)
                results.append((S, m, ok, wits))
    # family 3: very large Vmax stress (persistence at large w0)
    # NOTE: Wwidth cost ~ sum of speeds; cap V0 so safe_components stays tractable.
    for V0 in [1000, 3000, 7000]:
        large = [V0, V0+3, V0+9, V0+12, V0+15]
        small = [1,2,3,5,6,7,9,11]  # 8 small + 5 large = 13
        S = covering_check_and_M(sorted(set(small+large)))
        if S is None:
            # adjust to ensure covering & 13
            small = [1,2,3,4,6,7,9,11]
            S = covering_check_and_M(sorted(set(small+large)))
        if S is None: continue
        if len(S)!=13: continue
        if case_of(S)!='S3': continue
        m = Mval(S)
        ok, wits = C_holds(S)
        results.append((S, m, ok, wits))

    minM = F(1); worst=None; Cfails=[]
    for S, m, ok, wits in results:
        if m < minM: minM=m; worst=S
        if not ok: Cfails.append((S,m))
    print(f"  {len(results)} structured S3 sets tested")
    print(f"  min M among structured = {minM} = {float(minM):.5f} at {worst}")
    print(f"  C(S) total-failures among structured: {len(Cfails)}")
    for S,m in Cfails[:20]:
        print("    C FAILS:", S, "M=", m, float(m))
    # also show any with M<1/14
    below=[(S,m) for S,m,ok,wits in results if m<F(1,14)]
    print(f"  structured sets with M<1/14: {len(below)}")
    for S,m in below[:20]:
        print("    M<1/14:", S, "M=", m)
    return results, minM, worst, Cfails

# ============================================================================
# (I) the known exception: full exact audit
# ============================================================================
def test_I_known_exc():
    print("="*70)
    print("(I) full exact audit of the claimed single exception")
    S = KNOWN_EXC
    print("  S =", S, " covering:", is_covering(S), " case:", case_of(S))
    print("  Vmin,Vmax =", min(S), max(S), " 13*Vmin =", 13*min(S))
    m = Mval(S)
    print("  M(S) =", m, "=", float(m), "  >=1/14?", m>=F(1,14), " (claim 2/19)")
    print("  2/19 =", float(F(2,19)))
    # per-v criterion
    for v in S:
        A=[x for x in S if x!=v]
        w=Wwidth(A); r=w*7*v
        print(f"   remove {v:3d}: W={float(w):.5f} ratio={float(r):.4f} C={'YES' if r>1 else 'no '}")
    return m

# ============================================================================
def main():
    print("ADVERSARIAL VERIFICATION: lrc-lever-overlap advance on LRC(14) S3")
    print()
    A = test_A_criterion_implies_M()
    B = test_B_Lcollapse()
    C = test_C_Lruler()
    D = test_D_widthineq()
    Cfail_E = test_E_headline_C()
    bothfail_F = test_F_vmax_rule()
    below_G, worst_G, worstM_G = test_G_hunt_M()
    res_H, minM_H, worst_H, Cfails_H = test_H_structured()
    mI = test_I_known_exc()
    print()
    print("="*70)
    print("SUMMARY")
    print(f"  (A) C=>M>=1/14 implication holds on samples: {A}")
    print(f"  (B) L_collapse geometry sound: {B}")
    print(f"  (C) L_ruler geometry sound: {C}")
    print(f"  (D) WIDTH-INEQ algebra+sufficiency sound: {D}")
    print(f"  (E) #covering S3 sets where C(S) fails entirely (random): {len(Cfail_E)}")
    print(f"  (F) #sets where Vmax-rule fails AND C fails: {len(bothfail_F)}")
    print(f"  (G) min M over random S3 = {worstM_G}={float(worstM_G):.5f}; #below 1/14: {len(below_G)}")
    print(f"  (H) min M over structured S3 = {minM_H}={float(minM_H):.5f}; #C-fails: {len(Cfails_H)}")
    print(f"  (I) known-exc M = {mI} = {float(mI):.5f} (>=1/14: {mI>=F(1,14)})")

if __name__ == '__main__':
    main()
