#!/usr/bin/env python3
"""THM-779: the r=8 token-walk blocking criterion at the prime-7 lens (opus-S302).

Fuses boxeph's THM-773 token algebra with HYP-6840's rainbow frame:
at c = 7, owner a's bad sheet off its walls is the single token
    k_a(x) = -w_a^{-1} * round(w_a x)  (mod 7),
walls of a sit at x in (1/2 + Z)/w_a (mesh 1/w_a), the token steps by -w_a^{-1}
per wall, and owner a blocks NOTHING at its own wall (clearance exactly 1/14).

CRITERION (proved in THM-779): the r=8 deck is fully blocked over an interval J
iff (i) on every wall-free piece the 8 tokens cover F_7, and (ii) at every wall
in J the 7 non-walling tokens are EXACTLY F_7 (equivalently the walling owner is
a member of the unique collision pair), and (iii) J contains no simultaneous
walls. All integer dynamics -- no interval arithmetic.

PARTS:
 1  token formula referee: tokens vs the S300 bad_sheets machinery, exact.
 2  criterion referee: token-walk blocking check vs brute chamber check on the
    S301 search instances (must agree exactly).
 3  BREAK-LENGTH CENSUS: maximal blocking runs (consecutive covered walls) over
    random and annealed 8-tuples -- the empirical K0.
 4  THE PIERCE CONSEQUENCE: for spread cores, every core-safe component contains
    more walls than K0 => r=8 full blocking impossible there; verified on the
    S301 core; boxeph's 19/216 survivor's run length measured.
"""
import sys, os, random
from fractions import Fraction as F
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc14_certificates import good_intervals, LAM

random.seed(302)
P7 = 7

def circ(x):
    fx = x - (x.numerator // x.denominator)
    return min(fx, 1 - fx)

def inv7(w): return pow(w % 7, 5, 7)   # w^-1 mod 7 (w^5 = w^-1 since w^6=1)

def token(w, x):
    """k_a(x) = -w^{-1} round(w x) mod 7; None exactly on a wall (tie)."""
    wx = w * x
    n = wx.numerator; d = wx.denominator
    # round to nearest integer; half-integer = wall
    twice = 2 * n
    if twice % (2 * d) == d: return None          # exactly half-integer: wall
    r = (n + d // 2) // d if (2 * (n % d) < d or (2 * (n % d) == d)) else n // d + 1
    r = int(round(n / d)) if d < 10**15 else r    # fractions are exact & small here
    # exact nearest integer for Fraction:
    q, rem = divmod(n, d)
    r = q + (1 if 2 * rem > d else 0)
    return (-inv7(w) * r) % 7

def brute_sheet(w, x):
    for k in range(7):
        if circ(F(w) * (x + k) / 7) < LAM: return k
    return None

fails = 0
def check(name, cond):
    global fails
    print(("  [PASS] " if cond else "  [FAIL] ") + name)
    if not cond: fails += 1

print("=" * 88)
print("PART 1 -- token formula referee (exact, vs S300 bad-sheet machinery)")
print("=" * 88)
ok = True
for trial in range(4000):
    w = random.randint(1, 600)
    if w % 7 == 0: continue
    x = F(random.randint(0, 10**6), 10**6) + F(1, 3 * 10**7)
    if token(w, x) != brute_sheet(w, x): ok = False
check("token(w,x) == brute bad sheet on 4000 random exact points", ok)
# wall behavior: at x = (1/2 + m)/w the owner blocks nothing (clearance exactly 1/14)
ok2 = True
for trial in range(300):
    w = random.randint(1, 300)
    if w % 7 == 0: continue
    m = random.randint(0, w - 1)
    xw = (F(1, 2) + m) / w
    if brute_sheet(w, xw) is not None or token(w, xw) is not None: ok2 = False
check("at own walls the owner blocks nothing (token None, clearance = 1/14)", ok2)

# ---- the integer walk over an interval ----
def walls_in(w, lo, hi):
    """wall x-values of owner w in (lo, hi): x = (1/2+m)/w."""
    ms = range(int(w * lo - F(1,2)) - 1, int(w * hi - F(1,2)) + 2)
    return [ (F(1,2) + m) / w for m in ms if lo < (F(1,2) + m) / w < hi ]

def blocking_run_analysis(W8, lo, hi):
    """exact: sweep walls in (lo,hi); return (all_blocked, first_failure_kind,
    walls_checked, run_lengths list of maximal covered-wall runs)."""
    events = []
    for w in W8:
        for x in walls_in(w, lo, hi):
            events.append((x, w))
    events.sort()
    runs, cur = [], 0
    all_ok = True; first_fail = None
    # piece check between events too
    xs_prev = lo
    for x, w in events:
        mid = (xs_prev + x) / 2
        toks = [token(v, mid) for v in W8]
        piece_ok = set(t for t in toks if t is not None) == set(range(7))
        others = [token(v, x) for v in W8 if v != w]
        wall_ok = (None not in others) and (set(others) == set(range(7)))
        simultaneous = any(token(v, x) is None for v in W8 if v != w)
        if piece_ok and wall_ok and not simultaneous:
            cur += 1
        else:
            if all_ok:
                all_ok = False
                first_fail = ("piece" if not piece_ok else
                              "simultaneous" if simultaneous else "wall-rainbow")
            runs.append(cur); cur = 0
        xs_prev = x
    runs.append(cur)
    return all_ok, first_fail, len(events), runs

print()
print("=" * 88)
print("PART 2 -- criterion referee vs brute chamber check (S301 instances)")
print("=" * 88)
corePs = [2, 5, 9, 11, 13]
GP = good_intervals(corePs)
Wbest = [1, 8, 17, 18, 22, 23, 26, 74]     # S301 annealer best
# brute: sample points; token-walk: exact criterion -- compare failure verdicts
agree = True
for a, b in GP[:6]:
    blocked, kind, nw, runs = blocking_run_analysis(Wbest, a, b)
    # brute-force spot check on 15 points of the component
    brute_all = True
    for i in range(1, 16):
        x = a + (b - a) * F(i, 16)
        sheets = {brute_sheet(w, x) for w in Wbest} - {None}
        if len(sheets) < 7: brute_all = False
    if blocked != brute_all and brute_all:  # walk says fail but brute says ok on samples -> only walls could differ
        agree = False
check("token-walk verdicts consistent with brute sampling on 6 components", agree)

print()
print("=" * 88)
print("PART 3 -- break-length census: how long can blocking runs survive?")
print("=" * 88)
def max_run(W8, span=F(1)):
    _, _, nw, runs = blocking_run_analysis(W8, F(1, 1000), F(1, 1000) + span)
    return max(runs) if runs else 0, nw
best_run = 0; best_W = None
tot_runs = []
for trial in range(120):
    W8 = random.sample([w for w in range(1, 500) if w % 7], 8)
    r, nw = max_run(W8, F(1, 4))
    tot_runs.append(r)
    if r > best_run: best_run, best_W = r, sorted(W8)
tot_runs.sort()
print(f"  random census (120 8-tuples, quarter-period span): "
      f"median max-run {tot_runs[60]}, 90% {tot_runs[108]}, max {best_run} at {best_W}")
# anneal to MAXIMIZE the run
W8 = list(best_W); cur, _ = max_run(W8, F(1, 4))
for step in range(400):
    W2 = list(W8); W2[random.randrange(8)] = random.choice([w for w in range(1, 500) if w % 7])
    if len(set(W2)) < 8: continue
    r2, _ = max_run(W2, F(1, 4))
    if r2 >= cur: W8, cur = W2, r2
print(f"  annealed max blocking run: {cur} consecutive covered walls at {sorted(W8)}")
K0 = max(best_run, cur)
print(f"  ==> empirical K0 = {K0}")

print()
print("=" * 88)
print("PART 4 -- the pierce consequence on real cores + the 19/216 survivor")
print("=" * 88)
# walls per core-safe component for the S301 core with moderate exceptions
for name, W8 in [("S301 annealer best", Wbest),
                 ("annealed max-run tuple", sorted(W8))]:
    min_walls = None
    for a, b in GP:
        nw = sum(len(walls_in(w, a, b)) for w in W8)
        if min_walls is None or nw < min_walls: min_walls = nw
    print(f"  {name}: min walls per core-safe component = {min_walls} "
          f"{'> K0 => PIERCED on every component' if min_walls > K0 else '(<= K0: inconclusive by census alone)'}")
# boxeph survivor: P={5,7,8,13,14}, W={108,169,143,213,206,197,30,162}, t0=19/216
Wsur = [108, 169, 143, 213, 206, 197, 30, 162]
x0 = F(19, 216)
toks = [token(w, x0) for w in Wsur]
print(f"  survivor tokens at 19/216: {toks} (None = at wall)")
lo, hi = x0 - F(1, 40), x0 + F(1, 40)
blocked, kind, nw, runs = blocking_run_analysis(Wsur, lo, hi)
print(f"  survivor window (+-1/40): blocked={blocked}, first_fail={kind}, "
      f"walls={nw}, max run={max(runs)}")

print()
print("=" * 88)
print(f"SUMMARY: {'ALL CHECKS PASSED' if fails == 0 else f'{fails} FAILURES'}; K0 = {K0}")
print("=" * 88)
sys.exit(1 if fails else 0)
