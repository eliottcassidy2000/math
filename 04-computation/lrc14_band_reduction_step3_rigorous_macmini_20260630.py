#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S61
=======================
TWO rigorous deliverables for LRC14 (n=14, covering 13-sets, M(S)=max_t min_v ||v t||):

TASK A -- THE BAND-PRIME REDUCTION (a minor related statement, fully proven):
  LRC14  <=>  M(S) >= 1/14 for every covering 13-set that is a +-transversal-or-multiple
              mod EACH of the band primes {17,19,23}.
  Because: if S fails the transversal/multiple condition at some p in {17,19,23}, the radius-1
  witness theorem (klein-S42, PROVED) gives M(S) >= 2/p >= 2/23 > 1/14, so LRC14 holds for S
  with room to spare. The band primes are EXACTLY the primes in (n, 2n/ (14/183)] = (14, 26.1] .
  For primes p<=13 the condition is FREE (covering => a multiple of p).

TASK B -- MAKE klein-S45 STEP 3 (the "budget leaves ~1 slot") RIGOROUS, no exhaustive search:
  Step 3 reformulated rigorously as: (i) the (T_p) condition above; (ii) a missing core speed k
  BREAKS the pair {k,p-k} mod p for p in (12+k,25], whose only patch is a speed = +-k mod p that
  is LARGE (>= p-k >= 13); (iii) the CARDINALITY FLOOR: no multiple of 23 => >= 11 speeds in
  distinct pairs mod 23 (each speed lies in exactly one pair), so |S| >= 11, <= 2 spare.
  THEN we PROVE rigorously that cardinality CANNOT finish: a single (huge) CRT speed
  w = 0 mod 182, w = +-k mod 17*19*23 satisfies ALL large-speed obligations at once, so the
  cardinality lower bound on |S| is exactly 11 (achieved by the construction) -- the residual is
  NOT a counting gap but the residue/Diophantine claim that this CRT patch digs an M-hole
  (Step 4 = HYP-3745, perturbation-proved; general case = multi-family inexhaustibility = LRC14).

This script VERIFIES every rigorous claim numerically/exactly and prints a certificate table.
"""
from fractions import Fraction as F
from itertools import combinations
import random

# ---------- exact M(S) via complete breakpoint set (MISTAKE-86 guard) ----------
def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for d in (Sg[i]-Sg[j], Sg[i]+Sg[j]):
                if d > 0:
                    for kk in range(1, d): cand.add(F(kk, d))
    best = F(0); at = None
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in Sg)
        if g > best: best = g; at = t
    return best, at

def covers(S, n=14):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))

def is_pm_transversal(S, p):
    """True iff {+-s mod p : s in S, p nmid s} hits every nonzero pair of (Z/p)*."""
    pairs_needed = set(range(1, (p-1)//2 + 1))          # pair reps 1..(p-1)/2
    have = set()
    for s in S:
        r = s % p
        if r == 0: continue
        have.add(min(r, p-r))
    return pairs_needed <= have

def has_mult(S, p):
    return any(s % p == 0 for s in S)

n = 14
CONSTR = list(range(1, 13)) + [182]          # {1..12,182}, the believed covering-min set
target = F(n, n*n-n+1)                         # 14/183
print(f"n={n}  target = n/Phi6 = {target} = {float(target):.6f}   (1/14 = {float(F(1,14)):.6f})\n")

# ---------------------------------------------------------------------------
# 0. The band primes are exactly the primes p with 2/p > 14/183, i.e. p <= 23, AND p > 14.
# ---------------------------------------------------------------------------
print("="*78)
print("0. BAND PRIMES = primes p in (n, 26.14], i.e. {17,19,23}.  (2/p > 14/183 <=> p<=26)")
print("="*78)
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))
band = [p for p in range(2, 40) if isprime(p) and F(2,p) > target]
print(f"  primes with 2/p > 14/183 :  {band}   (all p<=23)")
print(f"  of these, p>n=14         :  {[p for p in band if p>14]}  = the BAND PRIMES")
print(f"  p<=13 are covered FREE   :  covering => multiple of p for every prime p<=13")
for p in band:
    print(f"    2/{p} = {float(F(2,p)):.4f}  >1/14? {F(2,p)>F(1,14)}   >14/183? {F(2,p)>target}")
print()

# ---------------------------------------------------------------------------
# 1. (T_p) REFORMULATION (rigorous, from the radius-1 witness theorem):
#    M(S) < 2/p  =>  (multiple of p)  OR  (+-transversal mod p).
#    Verify the CONTRAPOSITIVE numerically: not(T_p) => M(S) >= 2/p, for p in {17,19,23},
#    over many covering sets (random + structured). (The theorem itself is klein-S42 HYP-3741.)
# ---------------------------------------------------------------------------
print("="*78)
print("1. (T_p) REFORMULATION  M(S)<2/p => mult-of-p or +-transversal mod p   [p=17,19,23]")
print("   verifying contrapositive  not(T_p) => M(S) >= 2/p  on covering sets")
print("="*78)
random.seed(1)
viol = 0; tested = 0
def rand_covering(maxspeed=60):
    while True:
        S = sorted(random.sample(range(1, maxspeed), 13))
        if covers(S, n): return S
samples = [CONSTR, list(range(2,15)), list(range(1,13))+[13,182],
           list(range(1,13))+[26,182], [1,2,3,4,5,6,7,8,9,10,11,12,182]]
for _ in range(140): samples.append(rand_covering())
for S in samples:
    M,_ = Mexact(S)
    for p in (17,19,23):
        Tp = has_mult(S,p) or is_pm_transversal(S,p)
        if not Tp:
            tested += 1
            if M < F(2,p):
                viol += 1
                print(f"  !! VIOLATION p={p}: S={S} M={M} < 2/{p}")
print(f"  checked not(T_p) cases: {tested};  violations of M>=2/p: {viol}   "
      f"=> reformulation holds on all samples" if viol==0 else "  !!! REFORMULATION FAILED")
print()

# ---------------------------------------------------------------------------
# 2. TASK A -- THE REDUCTION: covering set failing some (T_p) has M >= 2/23 > 1/14.
#    So LRC14 only needs the triple-band-transversal covering sets.  Show the construction
#    IS such a set (so the subclass is nonempty and contains the extremal candidate).
# ---------------------------------------------------------------------------
print("="*78)
print("2. TASK A REDUCTION: a covering 13-set that fails (T_p) for some p in {17,19,23}")
print("   has M >= 2/p >= 2/23 = %.4f > 1/14 = %.4f  -> LRC14 holds for it with room." %
      (float(F(2,23)), float(F(1,14))))
print("   => LRC14  <=>  M>=1/14 on covering sets that are transversal-or-mult mod 17,19,AND 23")
print("="*78)
print("   the construction {1..12,182}:")
for p in (17,19,23):
    print(f"     mod {p}: mult? {has_mult(CONSTR,p)}   +-transversal? {is_pm_transversal(CONSTR,p)}"
          f"   (182 mod {p} = {182%p})")
Mc,tc = Mexact(CONSTR)
print(f"   => construction is a TRIPLE band-transversal covering set; M={Mc}={float(Mc):.5f} at t={tc}")
print(f"   so the hard subclass is non-empty and the extremal candidate lives in it.\n")

# ---------------------------------------------------------------------------
# 3. TASK B -- STEP 3 RIGOROUS:  missing k breaks pair {k,p-k} mod p at p in (12+k,25];
#    the only patch is a LARGE speed = +-k mod p (>= p-k >= 13).  Tabulate, all k.
# ---------------------------------------------------------------------------
print("="*78)
print("3. TASK B (Step 3 rigorous): missing core k breaks pair {k,p-k} mod p, p in (12+k,25]")
print("   patch = speed =+-k mod p; smallest non-missing rep is p-k (>=13) -> LARGE.  Per k:")
print("="*78)
print(f"   {'k':>2} | {'broken band primes p in (12+k,25]':35} | {'patch values (=+-k mod p, in 13..p-1)'}")
for k in range(1, 13):
    bp = [p for p in (17,19,23) if 12+k < p <= 25]
    patches = {p: [v for v in range(13, p) if v % p in (k % p, (-k) % p)] for p in bp}
    pat = "; ".join(f"{p}:{patches[p]}" for p in bp) if bp else "(none -- k=11,12: large-multiple subcase)"
    print(f"   {k:>2} | {str(bp):35} | {pat}")
print("   (matches klein-S45 Step1: k=1,3->{17,19,23}; k=6->{19,23}; k=10->{23}; k=12->{})\n")

# ---------------------------------------------------------------------------
# 4. THE CARDINALITY FLOOR (rigorous) AND ITS TIGHTNESS (the honest residual):
#    no mult of 23 => >=11 speeds in distinct pairs mod 23 => |S|>=11, <=2 spare.
#    AND: a single CRT speed satisfies ALL large-speed obligations -> cardinality min = 11,
#    so counting CANNOT reach 14; the residual is Step 4 (the patch's hole), not budget.
# ---------------------------------------------------------------------------
print("="*78)
print("4. CARDINALITY FLOOR (rigorous) + TIGHTNESS (why counting can't finish Step 3)")
print("="*78)
# 4a. floor: each speed lies in exactly one pair mod 23 (if 23 nmid s); 11 pairs => >=11 speeds.
print("  4a. no mult of 23 => +-transversal mod 23 needs all 11 pairs hit; each speed hits")
print("      exactly ONE pair mod 23 => >= 11 speeds.  |S|=13 => at most 2 spare.")
# construction tightness:
res23 = sorted(s % 23 for s in CONSTR)
pairs_hit = {}
for s in CONSTR:
    pr = min(s%23, 23-(s%23))
    pairs_hit.setdefault(pr, []).append(s)
redundant = sum(len(v)-1 for v in pairs_hit.values())
print(f"      construction residues mod23: {res23}")
print(f"      pairs covered: {len(pairs_hit)}/11 ; redundant speeds: {redundant}  (=> exactly 2 spare, TIGHT)")
# 4b. tightness: one CRT speed kills all large-speed obligations.
M3 = 17*19*23
print(f"  4b. large-speed obligations for missing k<=4 (no mult of band primes):")
print(f"      O1:13|L  O2:14|L  O3:L=+-k mod23  O4:L=+-k mod19  O5:L=+-k mod17")
# exhibit a single L satisfying all (CRT): L = 0 mod 182, = k mod (17*19*23)
def crt(rem_mod):
    # rem_mod: list of (r, m), pairwise coprime
    R, M = 0, 1
    for r, m in rem_mod:
        # solve R + M*x = r (mod m)
        inv = pow(M % m, -1, m)
        x = ((r - R) * inv) % m
        R += M * x; M *= m
    return R % M, M
for k in (1,2,3,4):
    L, MM = crt([(0,182),(k % M3, M3)])
    ok = (L%13==0 and L%14==0 and L%23 in (k%23,(-k)%23) and L%19 in (k%19,(-k)%19) and L%17 in (k%17,(-k)%17))
    print(f"      k={k}: single CRT speed L={L} (mod {MM})  satisfies O1..O5 simultaneously? {ok}")
print("      => the minimum number of speeds for O1..O5 is 1 (one giant CRT speed).")
print("      => CARDINALITY lower bound on |S| stays at 11 (the mod-23 floor); it never reaches 14.")
print("      => Step 3 CANNOT be closed by counting; the entire residual is whether that CRT")
print("         patch digs an M-hole (Step 4 = HYP-3745 / multi-family inexhaustibility).")
print()

# ---------------------------------------------------------------------------
# 5. SANITY: the all-but-one-core perturbations (missing exactly k, core else present + 182)
#    realize the patch=p-k and verify M>target -- consistent with klein-S45's table.
# ---------------------------------------------------------------------------
print("="*78)
print("5. SANITY: near-construction missing exactly core k (rest of core + 182 + patch p-k)")
print("="*78)
for k in range(1, 13):
    core = [v for v in range(1,13) if v != k]
    bp = [p for p in (17,19,23) if 12+k < p <= 25]
    # canonical escape: add smallest patches p-k for each broken prime, then 182, trim/pad to 13
    extra = sorted({p-k for p in bp} | {182})
    S = sorted(set(core) | set(extra))
    # pad with a benign large speed if <13, trim if >13 (keep covering)
    while len(S) < 13: S.append(max(S)+1)
    S = S[:13]
    if not covers(S, n):  # ensure covering by forcing 182 and a mult of 11 etc. (core has 11)
        S = sorted(set(core) | {182} | {p-k for p in bp}); S = (S+[max(S)+1]*13)[:13]
    M,t = Mexact(S)
    print(f"   k={k:>2}: S={S}  cover={covers(S,n)!s:5} M={str(M):8}={float(M):.4f} {'>target' if M>target else '<=TARGET!!'}")
print("\nDONE.")
