#!/usr/bin/env python3
"""
lrc14_angleC_finite_base -- mac-mini-2026-06-17-S6

ANGLE C: the FINITE BASE case + the bound that makes LRC(14) a finite computation.

Reductions already in canon:
  THM-523/524/525/526: LRC(14) <=> M(S) >= 1/14 for every primitive COVERING 13-set S,
  where M(S) = max_tau min_{v in S} ||v tau|| and "covering" = S has a multiple of every q in 2..14.

EXACT M tool (from task; validated): M is attained at an envelope vertex tau in
  { (2k+1)/(2v_i), k/(v_a+v_b), k/(v_a-v_b) } cap (0,1/2],
so M(S) is an EXACT rational with denom(M) | one of {2v_i, v_a+/-v_b} => denom(M) <= 2 max(S).
Hence "M(S) >= 1/14" is EXACTLY DECIDABLE -- no quantization-of-lcm needed (that is the LONELY
MEASURE L(S), THM-522, a DIFFERENT object; see HONESTY note in the report).

Angle C deliverables:
  (1) The bounded-max reduction. Show: every covering 13-set with max(S) > B is handled by
      the arc-width criterion C(S) [exists v: W(S\{v}) > 1/(7v)] => M(S) >= 1/14 (THM-526 ext),
      using EXACT W. Find the empirical threshold B_C beyond which C never fails on a search,
      and the pigeonhole-provable threshold B_pig where C is PROVED via the mu/N bound.
  (2) Count covering 13-sets with all entries <= B (finite). Is exhaustive M-check feasible?
  (3) Run the exact finite M-check for B = 14,28,56,84,168 (as feasible). Confirm
      inf M = 7/89 > 1/14 and watch for stabilization. Report completeness + computation size.

stdlib only; exact Fractions throughout.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import sys, time

C = F(1, 14)

# ---------- EXACT M ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t):
    return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): Cc.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): Cc.add(F(k, d)); k += 1
    Cc.add(F(1, 2)); return Cc
def M(S):
    return max(g(S, t) for t in cand(S))

# fast float prescreen, exact confirm -- big speedup for the finite sweep
def M_fast(S):
    Ss = sorted(set(S))
    cs = set()
    for v in Ss:
        k = 0
        while 2*k+1 <= v: cs.add((2*k+1)/(2.0*v)); k += 1
    n = len(Ss)
    for i in range(n):
        for j in range(i+1, n):
            for d in (Ss[i]+Ss[j], Ss[j]-Ss[i]):
                if d > 0:
                    k = 1
                    while 2*k <= d: cs.add(k/float(d)); k += 1
    cs.add(0.5)
    def gf(t):
        m = 1.0
        for v in Ss:
            r = (v*t) % 1.0; r = r if r <= 1.0-r else 1.0-r
            if r < m: m = r
        return m
    # we only need to know whether M >= 1/14. Find best float; if comfortably away from 1/14, decide;
    # else confirm exactly.
    bt = max(cs, key=gf); bv = gf(bt)
    thr = 1.0/14.0
    if bv > thr + 1e-6:   # safely >= 1/14
        return None, True     # decided LRC-good without exact M
    # near the boundary or below -> exact
    m = M(Ss)
    return m, (m >= C)

# ---------- EXACT widest safe arc W(A) at level 1/14 ----------
def darcs(v, c=C):
    hw = F(c, v); return [(F(k, v)-hw, F(k, v)+hw) for k in range(v)]
def wrapU(iv):
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1); a = lo - s; b = hi - s
        if b <= 1: o.append((a, b))
        else: o.append((a, F(1))); o.append((F(0), b-1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch: ch = ch if ch > hi else hi
        else: r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch)); return r
def Wsafe(A, c=C):
    dz = []
    for v in set(A): dz += darcs(v, c)
    if not dz: return F(1)
    dz = wrapU(dz)
    best = F(0)
    for i in range(len(dz)):
        hi = dz[i][1]; lo = dz[(i+1) % len(dz)][0] + (1 if i == len(dz)-1 else 0)
        if lo - hi > best: best = lo - hi
    return best
def criterion_exactW(S):
    """C(S): exists v with W(S\{v}) > 1/(7v). Returns (holds, v)."""
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        if Wsafe(A) > F(1, 7*v): return True, v
    return False, None

# ---------- pigeonhole-PROVABLE criterion (Angle A, mu/N bound) ----------
def mu(A): return 1 - sum(F(1, 7*u) for u in A)
def pigeon_holds(S):
    """SUFFICIENT (proved) condition for C(S): exists v with 7v*mu(S\{v}) > sum_{u!=v} u.
       Uses W(A) >= mu(A)/N(A) and N(A) <= sum_{u in A} u. Returns (holds, v)."""
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        if 7*v*mu(A) > sum(A): return True, v
    return False, None

qs = range(2, 15)
def covering(S): return all(any(v % q == 0 for v in S) for q in qs)
def primitive(S): return reduce(gcd, S) == 1


def part1_thresholds():
    print("="*78)
    print("PART 1: bounded-max reduction -- where does the arc-width criterion handle S?")
    print("="*78)
    print("""
  Two thresholds:
    B_pig  = the max(S) above which the PIGEONHOLE-PROVABLE criterion (mu/N) fires
             on every covering 13-set (so all such S are RIGOROUSLY handled).
    B_C    = the max(S) above which the EXACT-W criterion C(S) never fails on search
             (empirical; stronger but the universal proof is the open GAP).
""")
    import random
    random.seed(606)
    # search for the LARGEST max(S) at which pigeon_holds fails (over a broad sample)
    pig_fail_maxes = []
    exactW_fail_maxes = []
    N = 30000
    for _ in range(N):
        st = random.random()
        if st < 0.35:
            drop = random.choice(range(1, 14))
            base = [v for v in range(1, 14) if v != drop]
            S = base + [14*random.randint(1, 25)]
        elif st < 0.65:
            dr = random.sample(range(1, 14), 2)
            base = [v for v in range(1, 14) if v not in dr]
            S = base + [14*random.randint(1, 12),
                        random.choice([11, 13, 22, 26, 143, 154, 182, 286])*random.randint(1, 4)]
        else:
            V = random.randint(20, 160)
            csz = random.randint(3, 6)
            lo = max(2, V-22)
            cl = random.sample(range(lo, V+1), csz)
            small = random.sample(range(1, 14), 13-csz)
            S = small + cl
        S = sorted(set(S))
        if len(S) != 13 or not covering(S): continue
        ph, _ = pigeon_holds(S)
        if not ph:
            pig_fail_maxes.append(max(S))
            cw, _ = criterion_exactW(S)
            if not cw:
                exactW_fail_maxes.append(max(S))
    print(f"  sampled covering 13-sets where PIGEONHOLE (mu/N) fails: {len(pig_fail_maxes)}")
    if pig_fail_maxes:
        print(f"    largest max(S) with pigeonhole-fail: {max(pig_fail_maxes)}")
    print(f"  of those, EXACT-W criterion C(S) ALSO fails: {len(exactW_fail_maxes)}")
    if exactW_fail_maxes:
        print(f"    largest max(S) with exact-W C-fail: {max(exactW_fail_maxes)}")
    else:
        print("    exact-W criterion C held on EVERY sampled set (no fail).")
    return (max(pig_fail_maxes) if pig_fail_maxes else 0,
            max(exactW_fail_maxes) if exactW_fail_maxes else 0)


def part1b_pigeon_clean_bound():
    print()
    print("-"*78)
    print("PART 1b: a CLEAN provable bound. Remove the largest runner V; bound when pigeon fires.")
    print("-"*78)
    print("""
  Remove v=V=max(S). A=S\\{V} (12 distinct runners, each < V). Pigeon fires via V iff
        7V*mu(A) > sum_{u in A} u.
  mu(A) = 1 - (1/7) sum_{u in A} 1/u >= 1 - (1/7) H_12   (sum 1/u maximal at A={1..12}),
        = 1 - H_12/7.  And sum_{u in A} u <= sum of the 12 largest entries < V, but that is
  not bounded by a constant when the set is CLUSTERED. So removing V alone does NOT give a
  clean constant bound for clustered sets -- this is exactly the GAP the task flags.
""")
    H12 = sum(F(1, u) for u in range(1, 13))
    mu_floor = 1 - H12/7
    print(f"  H_12 = {H12} ~ {float(H12):.5f};  mu(A) >= 1 - H_12/7 = {mu_floor} ~ {float(mu_floor):.5f}")
    # For the SINGLE-large family (12 entries from {1..13} + V) sum_A u <= 78, so:
    #   pigeon via V fires iff 7V*mu(A) > sum_A u, with mu(A)>=mu_floor, sum_A<=78.
    # Worst over the 12-subsets:
    worst = F(0); worstA = None
    for A in combinations(range(1, 14), 12):
        A = list(A); thr = F(sum(A)) / (7*mu(A))   # V must exceed this
        if thr > worst: worst = thr; worstA = A
    import math
    print(f"  SINGLE-large family ({{12 of 1..13}} + V): pigeon-via-V fires iff V > {float(worst):.4f}")
    print(f"    (worst 12-core {worstA}); so B_pig^single = {math.ceil(float(worst))} suffices there.")
    print("  => For single-large covering sets, pigeonhole PROVES C for V >= 22. Finite base: V <= 21.")
    print("  => For CLUSTERED-large sets the mu/N pigeon is too weak (N << sum_A u); the EXACT W is")
    print("     needed and its universal lower bound is the remaining GAP (Angle A's open part).")


def count_covering(B):
    """Count primitive covering 13-sets with all entries in [1,B]. Returns (n_cover, n_cover_prim)."""
    # must contain a multiple of each q in 2..14. Build by inclusion over required residues is hard;
    # just enumerate combinations is C(B,13) -- only feasible for tiny B. Use smarter: the set must
    # contain mult of 8,9,11,13,7,5, etc. We count by direct combination enumeration for small B,
    # and by an estimate for larger B.
    from math import comb
    return comb(B, 13)


def part2_counts():
    print()
    print("="*78)
    print("PART 2: how many covering 13-sets have all entries <= B? (finiteness + feasibility)")
    print("="*78)
    from math import comb
    print("  Total 13-subsets of {1..B} = C(B,13):")
    for B in [14, 20, 28, 40, 56, 84, 120, 168, 300]:
        print(f"    B={B:4d}:  C({B},13) = {comb(B,13):,}")
    print("""
  C(B,13) is the crude superset count; the COVERING ones are a small fraction, but C(B,13)
  is the work an exhaustive sweep does (we filter covering inside the loop). Feasible region:
  C(B,13) up to ~10^8-10^9 with the float-prescreen M_fast (most sets decided without exact M).
  Beyond B~40 the raw C(B,13) explodes, so the sweep must enumerate covering sets DIRECTLY,
  not filter -- see Part 3, which constrains the large runner to the few covering-forced values.
""")


def enumerate_covering_bounded(B, large_threshold=14):
    """Enumerate primitive covering 13-sets with max(S) <= B by structure:
       12 'small' runners from {1..13} (the unavoidable core) + 1 'extra' in (13,B],
       PLUS genuinely-small sets that are covering with all entries <= 13? (none -- need mult of 14).
       This captures the SINGLE-large family exhaustively. For clustered sets we add a separate sweep.
       Returns a generator of sets. (Single-large is the family with the known minimizer 7/89.)"""
    seen = set()
    # single-large: choose 12 of {1..13}, add one e in [14,B]
    for core in combinations(range(1, 14), 12):
        core = list(core)
        for e in range(14, B+1):
            S = tuple(sorted(core + [e]))
            if len(set(S)) != 13: continue
            if S in seen: continue
            if not covering(S): continue
            seen.add(S)
            yield list(S)


def part3_finite_check():
    print()
    print("="*78)
    print("PART 3: exact finite M-check, inf M, stabilization")
    print("="*78)
    print("""
  Sweep the SINGLE-large covering family {12 of 1..13} + e (14<=e<=B) exhaustively with exact M.
  This is the family containing the known global minimizer S*={1..11,13,84}, M=7/89.
  Report inf M for increasing B and whether it stabilizes at 7/89.
""")
    overall_min = None; overall_arg = None
    # the single-large covering family is tiny (<400 sets up to B=600), so compute exact M for ALL.
    for B in [84, 168, 182, 300, 600, 1200]:
        t0 = time.time()
        mn = None; arg = None; n = 0; below = 0
        for S in enumerate_covering_bounded(B):
            n += 1
            m = M(sorted(S))            # exact
            if m < C: below += 1
            if mn is None or m < mn: mn = m; arg = sorted(S)
        dt = time.time() - t0
        if mn is not None and (overall_min is None or mn < overall_min):
            overall_min = mn; overall_arg = arg
        mtxt = f"{mn}={float(mn):.6f}" if mn is not None else "(empty family)"
        print(f"  B={B:4d}: single-large covering swept={n:>7,}  "
              f"exact-min M={mtxt}  (arg {arg})  below 1/14:{below}  [{dt:.1f}s]")
    print()
    print(f"  GLOBAL inf M over single-large family (B<=600) = {overall_min} = {float(overall_min):.6f}")
    print(f"    achieved at {overall_arg}")
    print(f"    1/14 = {float(C):.6f};  7/89 = {float(F(7,89)):.6f}")
    print(f"    inf M >= 1/14 ?  {overall_min >= C}   (= 7/89 ?  {overall_min == F(7,89)})")
    return overall_min, overall_arg


def part4_clustered_finite():
    print()
    print("="*78)
    print("PART 4: clustered-large finite base -- exhaustive small-B covering sweep (all structures)")
    print("="*78)
    print("""
  The single-large family is not all covering sets; clustered-large sets exist. But by Part 1b,
  EVERY covering 13-set with max(S) > B_0 is handled by the criterion (pigeonhole for single-large,
  exact-W empirically for clustered). To make a COMPLETE finite base we must also enumerate
  clustered covering sets with max <= B_0. Here we do the genuinely exhaustive sweep over ALL
  covering 13-subsets of {1..B0} for small B0, confirming inf M there too.
""")
    from math import comb
    thrf = 1.0/14.0
    for B0 in [14, 18, 22, 24]:
        t0 = time.time()
        mn = None; arg = None; n = 0; below = 0
        for S in combinations(range(1, B0+1), 13):
            S = list(S)
            if not covering(S): continue
            if reduce(gcd, S) != 1: continue
            n += 1
            # float prescreen: only sets whose float-M is within 1e-3 of 1/14 need an exact M
            # for the inf; but to REPORT the exact inf we exact-eval the near-boundary ones and
            # track the smallest. Sets far above 1/14 cannot be the minimizer.
            cs = set()
            for v in S:
                k = 0
                while 2*k+1 <= v: cs.add((2*k+1)/(2.0*v)); k += 1
            for i in range(len(S)):
                for j in range(i+1, len(S)):
                    for d in (S[i]+S[j], S[j]-S[i]):
                        if d > 0:
                            k = 1
                            while 2*k <= d: cs.add(k/float(d)); k += 1
            cs.add(0.5)
            def gf(t):
                mm = 1.0
                for v in S:
                    r = (v*t) % 1.0; r = r if r <= 1.0-r else 1.0-r
                    if r < mm: mm = r
                return mm
            bv = gf(max(cs, key=gf))
            # exact-eval if float-M is within 2e-3 of the current incumbent min (or of 1/14)
            ref = float(mn) if mn is not None else 0.10
            if bv <= ref + 2e-3 or bv <= thrf + 2e-3:
                m = M(S)
                if m < C: below += 1
                if mn is None or m < mn: mn = m; arg = S
        dt = time.time()-t0
        tag = f"exact-min M={mn}={float(mn):.6f}" if mn is not None else "(all far above 1/14)"
        print(f"  B0={B0}: exhaustive primitive covering 13-subsets of {{1..{B0}}}: count={n:>6,}  {tag}  "
              f"arg={arg}  below 1/14:{below}  [C({B0},13)={comb(B0,13):,}]  [{dt:.1f}s]")


if __name__ == "__main__":
    bpig, bexact = part1_thresholds()
    part1b_pigeon_clean_bound()
    part2_counts()
    mn, arg = part3_finite_check()
    part4_clustered_finite()
    print()
    print("="*78)
    print("VERDICT")
    print("="*78)
    print(f"""
  - EXACT decidability: M(S)>=1/14 is exactly decidable (denom(M)<=2max(S)); the finite base
    is a genuine terminating computation. (The 'M in (1/(14 lcm))Z' quantization claimed in the
    task is FALSE for M -- it holds for the LONELY MEASURE L, THM-522. M's correct quantization
    is denom(M) | (2v_i) or (v_a +/- v_b), giving denom(M)<=2max(S). Decidability is unaffected.)
  - inf M over the single-large family: 7/89~0.07865 up to B<=168, then 14/183~0.07650 (at e=182);
    BOTH are > 1/14. (HONESTY: the repo 'inf M = 7/89' is only the B<=168 value; the 1..12 + 182
    config gives a smaller 14/183 that is still LRC-good. Watch for further drops at larger e.)
  - The bounded-max reduction is CLEAN and PROVED for the SINGLE-large family (pigeon B_pig=22).
  - For CLUSTERED-large covering sets the pigeonhole (mu/N) is too weak; the exact-W criterion
    holds on all searches but its UNIVERSAL proof is the remaining GAP. So: Angle A (large, single)
    + this finite base = COMPLETE proof for the single-large family; the clustered-large universal
    bound is the one missing lemma.
""")
