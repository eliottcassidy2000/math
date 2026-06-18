#!/usr/bin/env python3
"""
REFRAME 7 — the others-clear condition as a simultaneous Diophantine
congruence-interval (Bohr-set) system.

THM-524 reduces LRC(N) at S to: SOME binding-pair crossing tau* = num/D
(D = v_a +/- v_b) clears all runners at level >= 1/N, i.e.

    for every v in S:   ||v * num/D|| >= 1/N
    <=>  v*num mod D  lies in the closed interval [ceil(D/N), D - ceil(D/N)]
         (the "central band" of width D - 2*ceil(D/N) + 1 residues).

So fixing D, the question "does crossing tau=num/D clear all 13 at level 1/N"
is a SIMULTANEOUS congruence-interval problem in the single unknown num mod D:
    num must avoid, for each v, the two "forbidden caps" {r : v*r mod D in [0,c-1] U [D-c+1, D-1]}
    where c = ceil(D/N).

This script:
 (1) Reformulates M(S) >= 1/N exactly as "exists (D, num) with the band property",
     and verifies the reformulation is EQUIVALENT to the M-tool (no gap, no slop).
 (2) Counting / pigeonhole: at the optimal D for the covering minimizer, how many
     num in [1, D-1] clear all 13?  Is it always >= 1?  How does the clearing-count
     behave (is there ever exactly 0 for the true M-attaining D -- which would be a
     contradiction, a useful internal check)?
 (3) CRT / density test: the covering condition (S has a multiple of every q in 2..N)
     -- does it FORCE a clearing num at SOME D via a counting lower bound?  We test the
     naive union-bound: each runner v forbids a fraction ~ 2*(c-1)/D... no, ~ 2*c/D of
     residues but the forbidden set for v is a union of gcd(v,D) arithmetic-progression
     caps. Compute the EXACT forbidden-measure and union-bound slack; report whether the
     union bound is ever positive (would PROVE clearance at that D).
"""
from fractions import Fraction as F
from math import gcd, ceil
import itertools, random

N = 14  # LRC(14): want M(S) >= 1/14

# ---------- exact M tool (verbatim from task) ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t):
    return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

# ---------- REFRAME 7 primitives ----------
def clears_at_level_int(S, D, num, k):
    """Does crossing tau = num/D clear every runner at level >= k/D, i.e.
       v*num mod D in [k, D-k] for all v?  (level k/D <-> band [k, D-k].)"""
    for v in S:
        r = (v * num) % D
        if not (k <= r <= D - k):
            return False
    return True

def central_band_for_level(D, level):
    """Smallest integer k with k/D >= level; the band is [k, D-k]."""
    # k = ceil(level * D)
    k = ceil(level * D)  # level is a Fraction
    return k

def count_clearing_num(S, D, k):
    """Number of num in [1, D-1] with v*num mod D in [k, D-k] for all v.
       (num=0 always fails: gives r=0 < k.)"""
    cnt = 0
    cl = []
    for num in range(1, D):
        if clears_at_level_int(S, D, num, k):
            cnt += 1; cl.append(num)
    return cnt, cl

def reframe7_M_via_bands(S, Nval=N):
    """Compute M(S) purely through the band/congruence reformulation, and the
       attaining (D, num). For each candidate denominator D (from pairwise sums/diffs
       and 2*v_i for the all-odd peak), and each num in [1,D-1] coprime-or-not,
       the cleared level is (min over v of (v*num mod D treated as central distance))/D.
       Return best level and where attained. This must equal M(S)."""
    S = sorted(set(S))
    Ds = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            Ds.add(S[i] + S[j])
            if S[j] - S[i] > 0: Ds.add(S[j] - S[i])
    for v in S:
        Ds.add(2 * v)  # peak denominators
    best = F(0); bestat = None
    for D in Ds:
        for num in range(1, D):
            # central distance of v*num mod D from 0, as fraction of D:
            # dist = min(r, D-r) where r = v*num mod D
            lvl = min(min((v * num) % D, D - ((v * num) % D)) for v in S)
            val = F(lvl, D)
            if val > best:
                best = val; bestat = (D, num)
    return best, bestat

# ============================================================
print("=" * 78)
print("REFRAME 7: others-clear as a simultaneous congruence-interval system")
print("=" * 78)

# ----- (0) named covering cases -----
named = {
    "generic core-core {1..5,7..13,84}": {1,2,3,4,5,7,8,9,10,11,12,13,84},
    "resonant flank-w {1..11,13,84}":    set(range(1,12))|{13,84},
}
print("\n[0] Named covering cases: M-tool vs band-reformulation")
for name, S in named.items():
    m, at = M(S)
    bm, bat = reframe7_M_via_bands(S)
    D = at.denominator; num = at.numerator
    ok = (m == bm)
    print(f"  {name}")
    print(f"     M-tool      M={m}  tau*={at}  (D={D}, num={num})")
    print(f"     band-reform M={bm} at (D,num)={bat}  AGREE={ok}")
    # the band picture at the attaining D
    k = central_band_for_level(D, F(1, N))
    cnt, cl = count_clearing_num(S, D, k)
    # the actual attained level band
    klvl = m.numerator if m.denominator == D else None
    print(f"     band [k={k},D-k={D-k}] at level 1/14 (c=ceil(D/14)={k}): "
          f"#num in [1,{D-1}] clearing ALL 13 at level>=1/14 = {cnt}")
    if cnt and len(cl) <= 12: print(f"        clearing num: {cl}")

# ----- (1) EQUIVALENCE check: band-reformulation == M-tool on random covering sets -----
print("\n[1] Equivalence of band-reformulation and M-tool (random + covering sets)")
random.seed(7)
def make_covering_13set(maxspeed=120):
    """Random primitive 13-set that is COVERING: has a multiple of every q in 2..14."""
    for _ in range(4000):
        S = set()
        # seed with small speeds to help coverage
        base = random.sample(range(1, 30), random.randint(6, 11))
        S |= set(base)
        while len(S) < 13:
            S.add(random.randint(1, maxspeed))
        S = set(list(S)[:13])
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if is_covering(S):
            return S
    return None
def gcd_all(S):
    from functools import reduce
    return reduce(gcd, S)
def is_covering(S, Nval=N):
    for q in range(2, Nval + 1):
        if not any(v % q == 0 for v in S):
            return False
    return True

agree = 0; tot = 0; disagree_ex = []
trials = 0
while tot < 400 and trials < 60000:
    trials += 1
    S = make_covering_13set()
    if S is None: continue
    m, _ = M(S)
    bm, _ = reframe7_M_via_bands(S)
    tot += 1
    if m == bm: agree += 1
    else:
        if len(disagree_ex) < 5: disagree_ex.append((sorted(S), str(m), str(bm)))
print(f"  band-reform == M-tool on {tot} covering 13-sets: agree={agree}/{tot}")
if disagree_ex:
    print("  DISAGREEMENTS:", disagree_ex)
else:
    print("  -> band reformulation is EXACTLY equivalent to M (as expected; tautological).")

# ----- (2) Counting: at the M-attaining D, how many num clear all 13 at level 1/14? -----
print("\n[2] Clearing-num COUNT at the M-attaining denominator D (covering sets)")
print("    For S with M(S) >= 1/14, the attaining (D,num) certifies clearance.")
print("    Question: how many num in [1,D-1] clear all 13 at level>=1/14? min over sets?")
random.seed(11)
counts = []
min_count_case = None; min_count = 10**9
nsets = 0; trials = 0
while nsets < 250 and trials < 60000:
    trials += 1
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N):  # would be an LRC(14) counterexample -- record separately
        print(f"  !! COUNTEREXAMPLE candidate S={sorted(S)} M={m} < 1/14")
        continue
    D = at.denominator
    k = central_band_for_level(D, F(1, N))
    cnt, cl = count_clearing_num(S, D, k)
    nsets += 1
    counts.append(cnt)
    if cnt < min_count:
        min_count = cnt; min_count_case = (sorted(S), str(m), D, cnt)
from collections import Counter
print(f"  scanned {nsets} covering sets with M>=1/14.")
print(f"  clearing-count histogram (at M-attaining D, level 1/14): {dict(sorted(Counter(counts).items()))}")
print(f"  MIN clearing-count = {min_count}; case: S={min_count_case[0]}")
print(f"     M={min_count_case[1]} D={min_count_case[2]} clears={min_count_case[3]} num")
print("  NOTE: min clearing-count >= 1 is GUARANTEED tautologically (the attaining num clears).")
print("        The real question is whether a COUNTING bound (independent of knowing M)")
print("        forces count >= 1. See [3].")

# ----- (3) CRT / union-bound: does covering force a clearing num at the M-attaining D? -----
print("\n[3] Union-bound / CRT density test at the M-attaining D")
print("    Forbidden set for runner v at level 1/14, modulus D: ")
print("       F_v = {num in [1,D-1] : v*num mod D in [0,k-1] U [D-k+1,D-1]}, k=ceil(D/14).")
print("    |F_v| (exact) and union-bound slack  (D-1) - sum_v |F_v|.")
print("    If slack > 0 at some D, that PROVES a clearing num exists at that D (no M needed).")

def forbidden_count(v, D, k):
    """ |{num in [1,D-1] : v*num mod D in caps}|.
        v*num mod D takes each value in the coset {multiples of gcd(v,D)} equally often.
        The caps are [0,k-1] U [D-k+1, D-1] (2k-1 residues incl 0). num ranges over 1..D-1.
        Count exactly. """
    cnt = 0
    for num in range(1, D):
        r = (v * num) % D
        if r < k or r > D - k:  # r in [0,k-1] or [D-k+1, D-1]
            cnt += 1
    return cnt

random.seed(13)
slack_positive = 0; slack_tested = 0
best_slack = -10**9; best_slack_case = None
worst_slack = 10**9
slack_at_Mstar = []
nsets = 0; trials = 0
while nsets < 200 and trials < 60000:
    trials += 1
    S = make_covering_13set()
    if S is None: continue
    m, at = M(S)
    if m < F(1, N): continue
    D = at.denominator
    k = central_band_for_level(D, F(1, N))
    tot_forb = sum(forbidden_count(v, D, k) for v in S)
    slack = (D - 1) - tot_forb
    nsets += 1; slack_tested += 1
    slack_at_Mstar.append(slack)
    if slack > 0: slack_positive += 1
    if slack > best_slack: best_slack = slack; best_slack_case = (sorted(S), D, k, slack, str(m))
    if slack < worst_slack: worst_slack = slack
print(f"  at M-attaining D: union-bound slack positive in {slack_positive}/{slack_tested} sets.")
print(f"  best slack = {best_slack}; worst slack = {worst_slack}")
if best_slack_case:
    print(f"  best-slack case S={best_slack_case[0]} D={best_slack_case[1]} k={best_slack_case[2]} slack={best_slack_case[3]} M={best_slack_case[4]}")
print("  INTERPRETATION: union bound at the M-attaining D is essentially always NEGATIVE")
print("  (13 runners each forbid ~2k/D ~ 2/14 ~ 14% of residues -> total ~ 13*2/14 ~ 1.86 > 1).")
print("  So a naive union bound CANNOT prove clearance: the forbidden caps overlap heavily,")
print("  and that overlap is exactly the content LRC needs (it is NOT given by counting).")

# ----- (3b) but maybe a BETTER D exists with positive slack? sweep all candidate D -----
print("\n[3b] Sweep ALL candidate D for the worst covering set; is slack EVER positive?")
worst = sorted(named["resonant flank-w {1..11,13,84}"])  # the unique minimizer 7/89
S = set(worst)
Ds = set()
for i in range(len(worst)):
    for j in range(i+1, len(worst)):
        Ds.add(worst[i]+worst[j])
        if worst[j]-worst[i] > 0: Ds.add(worst[j]-worst[i])
posD = []
for D in sorted(Ds):
    if D < N: continue
    k = central_band_for_level(D, F(1, N))
    if 2*k - 1 >= D: continue  # band empty
    tot_forb = sum(forbidden_count(v, D, k) for v in S)
    slack = (D-1) - tot_forb
    if slack > 0:
        posD.append((D, slack))
print(f"  S={worst} (the unique LRC minimizer, M=7/89)")
print(f"  candidate D with POSITIVE union-bound slack at level 1/14: {posD[:20]}")
print(f"  (#positive-slack D = {len(posD)} out of {len([D for D in Ds if D>=N])} candidate D >= 14)")
if posD:
    D0 = posD[0][0]; k0 = central_band_for_level(D0, F(1,N))
    cnt, cl = count_clearing_num(S, D0, k0)
    print(f"  At first positive-slack D={D0}: actual clearing-num count = {cnt}  (clears={cl[:8]})")
    print("  -> positive slack at SOME D *would* give a counting proof of clearance there,")
    print("     IF that D's crossing also achieves a binding pair. But note: a clearing num")
    print("     at an arbitrary D gives M >= 1/14 regardless (||v num/D|| >= 1/14 for all v).")

print("\nDONE.")
