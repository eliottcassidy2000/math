"""
opus-2026-07-19-S400, part 2.

(A) BRANCH LEMMA (proof + exact verification): any 12-family with >= 11 even
speeds and one odd speed has M >= 1/12.
  Proof (two-sheet dodge, THM-760 mechanism at c=2): write evens 2u_i, odd b.
  Let tau be a maximizer for U = {u_1..u_11}: min_i ||u_i tau|| >= 1/12 by
  LRC(12) (11 speeds, settled).  The two lifts t = tau/2 and t = tau/2 + 1/2
  give identical even-runner distances ||2u_i t|| = ||u_i tau||, while the odd
  runner's two positions differ by exactly b/2 == 1/2 (mod 1), so one lift has
  ||b t|| >= 1/4.  Hence M >= min(1/12, 1/4) = 1/12.  QED (elementary).
  This closes the gcd-branch of the S-T cage: 11-even families are NOT in the
  micro-gap [1/13, 113/1466) since 1/12 > 113/1466.
  Verification: exact sheet-dodge witness on random 11-even families.

(B) PROBE of S-T Conjecture 7.1 (universal witness denominator) at k=13:
  "exists D s.t. for any d >= D every non-tight coprime v has a witness in
  (1/d)Z" -- probed against the repo's ladder {1..11,13,12m} (m=3..8,
  M = m/(12m+5) > 1/14, non-tight) and controls.  For each d, GOOD(d) iff
  some j/d has min_v ||v j/d|| >= 1/14 exactly.  Reports BAD d's; persistent
  BAD d at large d = evidence the conjecture is delicate at k=13.
"""
from fractions import Fraction
from math import gcd
import random

def dist(v, num, den):
    r = (v*num) % den
    return Fraction(min(r, den-r), den)

# ---------- (A) branch lemma verification ----------
print("=== (A) 11-even branch: exact sheet-dodge verification ===")
random.seed(4001)
ok = True
for trial in range(300):
    U = random.sample(range(1, 200), 11)
    b = random.choice(range(1, 400, 2))
    W = [2*u for u in U] + [b]
    # find tau = maximizer for U among q<=90 (LRC(12) guarantees >= 1/12 exists
    # at SOME tau; the scan finds a >=1/12 point -- always succeeds)
    best = (Fraction(0), None)
    for q in range(2, 91):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(dist(u, a, q) for u in U)
            if m > best[0]:
                best = (m, Fraction(a, q))
    mU, tau = best
    if mU < Fraction(1, 12):
        print(f"  scan failed to reach 1/12 for U={U} (got {mU}) -- enlarge scan")
        ok = False
        continue
    # two sheets: t = tau/2 and tau/2 + 1/2
    t1 = tau/2
    t2 = tau/2 + Fraction(1, 2)
    def fam_min(t):
        return min(dist(w, t.numerator, t.denominator) for w in W)
    got = max(fam_min(t1), fam_min(t2))
    if got < Fraction(1, 12):
        print(f"  VIOLATION: W={W} sheet-dodge gives only {got}")
        ok = False
print(f"  300 random 11-even families: sheet-dodge >= 1/12 in all cases: {ok}")

# ---------- (B) Conjecture 7.1 probe ----------
print("\n=== (B) S-T Conjecture 7.1 probe at k=13 (threshold 1/14) ===")
thr_n, thr_d = 1, 14   # min >= 1/14  <=>  14*r >= den at each check (exact)
def good_d(V, d):
    for j in range(0, d):
        okj = True
        for v in V:
            r = (v*j) % d
            r = min(r, d - r)
            if 14*r < d:   # ||v j/d|| = r/d < 1/14
                okj = False
                break
        if okj:
            return True
    return False

fams = {f"ladder m={m} {{1..11,13,{12*m}}} M={m}/{12*m+5}":
            list(range(1,12))+[13,12*m] for m in range(3, 9)}
fams["control {1..12,15} (loose)"] = list(range(1,13))+[15]
fams["tight AP {1..13} (M=1/14, tight -- witnesses only at s/14)"] = list(range(1,14))
for name, V in fams.items():
    bad = [d for d in range(20, 1501) if not good_d(V, d)]
    tail_bad = [d for d in bad if d >= 1000]
    print(f"  {name}:")
    print(f"    #BAD d in [20,1500] = {len(bad)}; largest BAD = {max(bad) if bad else None}; "
          f"#BAD >= 1000: {len(tail_bad)}")
    if bad and len(bad) < 40:
        print(f"    BAD list: {bad}")
print("\n  Reading: if non-tight families keep BAD d's arbitrarily high, a")
print("  uniform D in Conjecture 7.1 must exceed them; a persistent-arithmetic")
print("  BAD pattern (e.g. all d in a residue class BAD) would refute it.")
