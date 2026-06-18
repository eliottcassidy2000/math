#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the "easy-dominates-hard" LRC(14) reduction.

LENS: re-derive from scratch
  (1) the danger-comb width/spacing of a parked runner w == 0 (mod 14), and
  (2) the meas(G_C) lower bound (claimed 7/858 at the drop-6 core).
Then test the heart of the reduction: "thin comb cannot cover a fat G_C."

We work with EXACT rationals (fractions.Fraction) throughout. No floats for decisions.

Definitions (matching canon):
  ||x|| = distance from x to nearest integer.
  Danger set of speed v (gap target 1/14):  D_v = { tau in [0,1) : ||v tau|| < 1/14 }.
  G_C = [0,1) \ union_{v in C} D_v   (the CLOSED-or-OPEN lonely set; we use the
        OPEN danger sets, so G_C = { tau : ||v tau|| >= 1/14 for all v in C }).
  meas(G_C) = Lebesgue measure of G_C.

We represent every set as a finite union of rational intervals and compute exact measures.
"""

from fractions import Fraction as F
from math import gcd
from functools import reduce

# ----------------------------------------------------------------------
# Exact interval arithmetic on [0,1): list of (a,b) with 0<=a<b<=1, disjoint, sorted.
# ----------------------------------------------------------------------

def normalize(intervals):
    """Merge a list of (a,b) Fraction intervals into disjoint sorted form, clipped to [0,1]."""
    cl = []
    for a, b in intervals:
        if a < 0: a = F(0)
        if b > 1: b = F(1)
        if a < b:
            cl.append((a, b))
    if not cl:
        return []
    cl.sort()
    out = [list(cl[0])]
    for a, b in cl[1:]:
        if a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1][1] = b
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]

def measure(intervals):
    return sum((b - a for a, b in intervals), F(0))

def complement(intervals):
    """Complement within [0,1) of a normalized interval list."""
    intervals = normalize(intervals)
    out = []
    prev = F(0)
    for a, b in intervals:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        out.append((prev, F(1)))
    return out

def intersect(A, B):
    """Intersection of two normalized interval lists."""
    A = normalize(A); B = normalize(B)
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out

# ----------------------------------------------------------------------
# Danger set D_v = { tau in [0,1) : ||v tau|| < 1/14 } as exact intervals.
#   ||v tau|| < 1/14  <=>  v tau is within 1/14 of an integer k
#   <=>  tau in ( (k - 1/14)/v , (k + 1/14)/v )  for k = 0,1,...,v-1 (and wrap).
# Each tooth has width 2/(14 v) = 1/(7 v). There are v teeth. Total = v/(7v) = 1/7.
# We use a general half-width h (default 1/14) so we can test variants.
# ----------------------------------------------------------------------

def danger(v, h=F(1, 14)):
    """Danger set of speed v with half-width h on the v*tau axis. Returns intervals in [0,1)."""
    out = []
    for k in range(0, v + 1):
        a = (F(k) - h) / v
        b = (F(k) + h) / v
        out.append((a, b))  # normalize clips to [0,1]; k=0 left part and k=v right part wrap
    return normalize(out)

def lonely_set(C, h=F(1, 14)):
    """G_C = complement of union of danger sets of speeds in C."""
    U = []
    for v in C:
        U += danger(v, h)
    U = normalize(U)
    return complement(U)

def lonely_measure(C, h=F(1, 14)):
    return measure(lonely_set(C, h))

# ----------------------------------------------------------------------
# Exact gap M(S) tool (verbatim from the prompt, validated).
# ----------------------------------------------------------------------

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
        if v > b:
            b = v; at = t
    return b, at

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 1: Danger-comb geometry of a parked runner w == 0 (mod 14)")
print("=" * 78)

for w in [14, 28, 84, 98, 168]:
    Dw = danger(w)
    m = measure(Dw)
    # number of teeth, tooth width, tooth-center spacing
    tooth_w = F(1, 7 * w)
    spacing = F(1, w)
    print(f"w={w:4d}: meas(D_w)={m}  (={float(m):.6f})  teeth={w}  "
          f"tooth_width=1/{7*w}={float(tooth_w):.6g}  spacing=1/{w}={float(spacing):.6g}")
    # sanity: measure should be exactly 1/7
    assert m == F(1, 7), f"meas(D_w) != 1/7 for w={w}: {m}"
print("CHECK: meas(D_w) = 1/7 exactly for all w (each tooth width 1/(7w), w teeth).")
print()
print("So the comb D_w has w teeth at centers k/w (k=0..w-1), each of HALF-width")
print("  (1/14)/w = 1/(14w). The GAP between consecutive teeth (safe between teeth)")
print("  has length  1/w - 2/(14w) = 1/w - 1/(7w) = (6/7)/w = 6/(7w).")
print("  As w -> infinity, the comb becomes arbitrarily FINE: the largest safe")
print("  interval BETWEEN teeth has length 6/(7w) -> 0.  <-- 'thin comb' intuition.")
print()

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 2: meas(G_C) lower bound -- re-derive the drop-6 extremal 7/858")
print("=" * 78)

# The conjectured extremal 12-core (AP {1..13} with 6 dropped):
C_drop6 = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13]
mG = lonely_measure(C_drop6)
print(f"C_drop6 = {C_drop6}")
print(f"meas(G_C) = {mG}  (= {float(mG):.8f})")
print(f"claimed extremal 7/858 = {F(7,858)} (= {float(F(7,858)):.8f})")
print(f"MATCH 7/858 ? {mG == F(7,858)}")
print()

# Show the actual safe intervals (the 'mouths').
G = lonely_set(C_drop6)
print(f"G_C has {len(G)} safe components (mouths):")
for a, b in G:
    print(f"   [{a}, {b}]   length {b-a} = {float(b-a):.8f}")
print(f"   total = {measure(G)}")
print()

# Sweep ALL 12-subsets of {1..13} (drop one element) -- find the minimum meas(G_C).
print("Sweep: meas(G_C) over all 12-subsets obtained by dropping one of {1..13}:")
rows = []
full = list(range(1, 14))
for drop in full:
    C = [v for v in full if v != drop]
    m = lonely_measure(C)
    rows.append((m, drop))
rows.sort()
for m, drop in rows:
    print(f"   drop {drop:2d}: meas(G_C) = {str(m):>14s}  = {float(m):.8f}")
print(f"MIN over single-drops = {rows[0][0]} at drop {rows[0][1]}")
print()

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 3: THE CORE TEST -- can meas(G_C) be SMALLER than the comb spacing?")
print("=" * 78)
print("""
The reduction's intuition: a parked runner w==0(mod14) has danger comb D_w with
teeth spaced 1/w and the safe-between-teeth gaps shrink like 6/(7w). For the
covering set S = C u {w} to be a counterexample, D_w would have to ERASE all of
G_C, i.e. G_C must be entirely covered by D_w's teeth.

If meas(G_C) >= some c > 0 and the comb is 'thin' (fine teeth), the intuition says
the teeth can't swallow a fat G_C. BUT the real question for the LOWER bound on
L(S) = meas(G_C \\ D_w) is NOT 'spacing vs measure' -- it is how much of G_C the
teeth overlap. We test directly: L(S) = meas(G_C cap complement(D_w)).
""")

def L_exact(C, w, h=F(1, 14)):
    """L(S) = meas{tau : ||v tau|| >= 1/14 for all v in C u {w}} = meas(G_C \\ D_w)."""
    G = lonely_set(C, h)
    Dw = danger(w, h)
    safe = intersect(G, complement(Dw))
    return measure(safe)

# Verify against the canon numbers: drop-6 + 84 -> the famous covering set, M=7/89.
S_famous = sorted(C_drop6 + [84])  # wait: canon says {1..11,13,84}. Let's check both.
print("Canon baseline covering sets:")
for C, w, name in [
    ([1,2,3,4,5,6,7,8,9,10,11,13], 84, "{1..11,13} u {84}"),
    (C_drop6, 84, "drop-6 core {1..5,7..13} u {84}"),
]:
    S = sorted(C + [w])
    Ms, ts = M(S)
    Ls = L_exact(C, w)
    mGc = lonely_measure(C)
    print(f"  {name}: S={S}")
    print(f"     M(S) = {Ms} = {float(Ms):.6f}   (>1/14={float(F(1,14)):.6f}? {Ms>=F(1,14)})  at tau={ts}")
    print(f"     meas(G_C) = {mGc} = {float(mGc):.8f}")
    print(f"     L(S)=meas(G_C minus D_w) = {Ls} = {float(Ls):.8f}")
print()

# Now the adversarial sweep: for the drop-6 core, sweep w over multiples of 14
# and ALSO arbitrary w (parked or not) and find the MINIMUM L(S).  Does L stay > 0?
print("Drop-6 core: sweep w = 14*m, m=1..40, compute L(S) and M(S) exactly.")
worst = None
for m in range(1, 41):
    w = 14 * m
    C = C_drop6
    L = L_exact(C, w)
    S = sorted(C + [w])
    # primitivity: gcd of S
    gg = reduce(gcd, S)
    Sp = [x // gg for x in S]
    Ms, _ = M(S)
    flag = ""
    if L <= 0:
        flag = "  <-- L=0 !!"
    if Ms < F(1, 14):
        flag += "  <-- M<1/14 COUNTEREXAMPLE !!"
    if worst is None or L < worst[0]:
        worst = (L, w, Ms)
    if m <= 12 or L <= F(1, 5000) or Ms < F(1, 14):
        print(f"   m={m:2d} w={w:4d}: L(S)={str(L):>16s}={float(L):.8f}  M(S)={float(Ms):.6f}{flag}")
print(f"   MIN L over m=1..40: L={worst[0]} = {float(worst[0]):.8f} at w={worst[1]}, M={float(worst[2]):.6f}")
print()

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 4: ADVERSARIAL HUNT -- is meas(G_C) >= 7/858 for ALL 12-subsets?")
print("=" * 78)
print("""
The decoupling lever needs meas(G_C) bounded below UNIFORMLY (OPEN-Q-108).
Canon's sharpened target: meas(G_C) >= 7/858 for every 12-subset C of distinct
positive integers. The danger-comb-vs-fat-G_C intuition is MEANINGLESS unless
this lower bound actually holds. We hunt for a 12-subset with meas(G_C) < 7/858.

NOTE: meas(G_C) = 0 is EASY to achieve if C itself is a 'covering-like' 12-set
(e.g. if C already contains a multiple of every q in 2..12 -- then by LRC the
12-runner gap can still be >= 1/13, but the gap-1/14 SAFE MEASURE could be tiny).
This is exactly the regime the reduction must survive. Let's probe.
""")

import itertools, random

TARGET = F(7, 858)
print(f"TARGET threshold 7/858 = {float(TARGET):.8f}")
print()

# 4a. meas(G_C) can it be 0?  A 12-set C with meas(G_C)=0 would mean the 12 combs
#     already cover [0,1) -- i.e. C is itself a 'measure-covering' 12-set.
#     By PROVEN LRC(12) the *gap* M(C) >= 1/13 > 1/14, so the CLOSED lonely set
#     (>= 1/14) is nonempty, but its MEASURE (open >1/14) can still be 0?
#     Test: does meas(G_C)=0 ever happen for a primitive 12-set?
print("4a. Hunt for meas(G_C) very small or zero among structured 12-sets.")
min_found = (TARGET, list(C_drop6))
# structured families: drop-2-from-{1..14}, dilated APs, AP+stranger, near-covering
fams = []
# all 2-subsets dropped from {1..14}  -> 12-subsets of {1..14}
base14 = list(range(1, 15))
for combo in itertools.combinations(base14, 12):
    fams.append(list(combo))
# AP {1..12} dilated/shifted and AP {1..11,13}-style + a large stranger replacing one
for w in range(14, 200):
    fams.append(sorted(C_drop6[:-1] + [w]))         # replace 13 by w
    fams.append(sorted([1,2,3,4,5,7,8,9,10,11,13] + [w]))  # 12 -> w variants
print(f"   testing {len(fams)} structured 12-sets ...")
zero_meas = []
for C in fams:
    if len(set(C)) != 12: 
        continue
    m = lonely_measure(C)
    if m < min_found[0]:
        min_found = (m, C)
    if m == 0:
        zero_meas.append(C)
print(f"   MIN meas(G_C) over structured = {min_found[0]} = {float(min_found[0]):.8f}")
print(f"     at C = {min_found[1]}")
print(f"   below 7/858? {min_found[0] < TARGET}")
print(f"   #(meas=0) found = {len(zero_meas)}")
if zero_meas[:3]:
    for C in zero_meas[:3]:
        print(f"     meas=0 example: C={C}  M(C)={float(M(C)[0]):.6f}")
print()

# 4b. Random primitive 12-subsets (bounded speeds) -- broad net.
print("4b. Random primitive 12-subsets, speeds in [1, V], hunt meas(G_C) < 7/858.")
random.seed(20260617)
def primitive(S):
    return reduce(gcd, S) == 1
best_rand = (TARGET, None)
n_below = 0
n_tested = 0
n_zero = 0
for V in [20, 30, 50, 90, 150]:
    for _ in range(4000):
        C = sorted(random.sample(range(1, V + 1), 12))
        if not primitive(C):
            continue
        n_tested += 1
        m = lonely_measure(C)
        if m == 0:
            n_zero += 1
        if m < best_rand[0]:
            best_rand = (m, C)
        if m < TARGET:
            n_below += 1
print(f"   tested {n_tested} primitive random 12-sets")
print(f"   #(meas < 7/858) = {n_below}   #(meas = 0) = {n_zero}")
print(f"   MIN meas(G_C) = {best_rand[0]} = {float(best_rand[0]):.8f} at C={best_rand[1]}")
print(f"   below 7/858? {best_rand[0] < TARGET}")
print()

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 5: 'Thin comb cannot cover fat G_C' -- does the geometry hold?")
print("=" * 78)
print("""
CLAIM under test: a single comb D_w (w==0 mod 14) cannot erase all of G_C when
meas(G_C) is bounded below. Quantitatively the decoupling bound is
   L(S) = meas(G_C \ D_w) >= (6/7) meas(G_C) - r/(7w)
where r = #components(G_C) (the boundary cost: each safe component can lose at
most one tooth's worth of its two ends; the canon's 'r/(7w)' term).

Re-derive r/(7w): each tooth has width 1/(7w). A connected safe component of G_C
can intersect many teeth in its interior (those are full teeth, removing a
(6/7)-fraction is wrong -- inside a component the comb removes the FULL tooth
measure, fraction (1/(7w))/(1/w) = 1/7 of length on average). The (6/7) factor
is meas-density of the comb COMPLEMENT = 1 - 1/7 = 6/7. So the heuristic:
   L(S) ~ (6/7) meas(G_C)  minus  boundary effects (partial teeth at component ends).
Boundary effect per component end <= half a tooth = 1/(14w); 2 ends, r components
=> total boundary loss <= 2r/(14w) = r/(7w).  THIS IS THE CANON BOUND. Verify it
holds EXACTLY (L >= (6/7)meas - r/(7w)) and find when it goes NEGATIVE (the gap).
""")

def components(C, h=F(1,14)):
    return len(lonely_set(C, h))

print("Drop-6 core (r=4 components): test decoupling bound L >= (6/7)meas - r/(7w).")
mGc = lonely_measure(C_drop6)
r = components(C_drop6)
print(f"   meas(G_C)={mGc}={float(mGc):.8f}, r={r}, (6/7)meas={float(F(6,7)*mGc):.8f}")
print(f"   {'w':>5} {'L(S) exact':>16} {'(6/7)meas-r/(7w)':>18} {'bound holds?':>12} {'bound>0?':>9}")
viol = 0
for m in range(1, 30):
    w = 14*m
    L = L_exact(C_drop6, w)
    bound = F(6,7)*mGc - F(r, 7*w)
    holds = L >= bound
    if not holds: viol += 1
    if m <= 10 or not holds or bound <= 0:
        print(f"   {w:>5} {float(L):>16.8f} {float(bound):>18.8f} {str(holds):>12} {str(bound>0):>9}")
print(f"   decoupling-bound violations (L < (6/7)meas - r/(7w)): {viol}")
print(f"   bound positive for w >= r/(6 meas) = {F(r,1)/(F(6)*mGc)} = {float(F(r,1)/(F(6)*mGc)):.3f}")
print()

# ----------------------------------------------------------------------
print("=" * 78)
print("PART 6: THE REAL GAP -- 3+ coordinated growing speeds (canon's GAP A)")
print("=" * 78)
print("""
The decoupling bound iterates ONE large speed at a time. The OPEN regime: the
12-core C ITSELF contains several large coordinated speeds, so meas(G_C) is NOT
a fixed positive constant -- it could in principle shrink toward 0 as the
coordinated speeds grow. THIS is what 'thin comb vs fat G_C' must survive: if
meas(G_C) itself -> 0, the whole reduction collapses.

We hunt: 12-sets with MULTIPLE large coordinated speeds, minimize meas(G_C).
The conjecture says it stays >= 7/858. Try to BREAK it.
""")

import itertools, random
# Coordinated families: small AP core + several multiples of 14 (or of 7, 12, etc.)
print("6a. {1..k} + several coordinated large multiples of 14, minimize meas(G_C).")
best = (F(1), None)
tested = 0
for ncore in [9, 10, 11]:
    core = list(range(1, ncore+1))
    nlarge = 12 - ncore
    for combo in itertools.combinations([14*j for j in range(1, 16)], nlarge):
        C = sorted(core + list(combo))
        if len(set(C)) != 12: continue
        tested += 1
        m = lonely_measure(C)
        if m < best[0]:
            best = (m, C)
print(f"   tested {tested} coordinated-14 cores; MIN meas(G_C)={best[0]}={float(best[0]):.8f}")
print(f"     at C={best[1]}; below 7/858? {best[0] < F(7,858)}")
print()

print("6b. Coordinated multiples of 7,12,84 (lcm structure) -- try to shrink G_C.")
best2 = (F(1), None)
tested2 = 0
mults = []
for base in [7, 12, 13, 84, 91, 168]:
    mults += [base*j for j in range(1, 8)]
mults = sorted(set(m for m in mults if m <= 600))
for ncore in [8, 9, 10]:
    core = list(range(1, ncore+1))
    nlarge = 12 - ncore
    cnt = 0
    for combo in itertools.combinations(mults, nlarge):
        C = sorted(set(core + list(combo)))
        if len(C) != 12: continue
        cnt += 1
        if cnt > 60000: break
        tested2 += 1
        m = lonely_measure(C)
        if m < best2[0]:
            best2 = (m, C)
print(f"   tested {tested2} mixed-coordinated cores; MIN meas(G_C)={best2[0]}={float(best2[0]):.8f}")
print(f"     at C={best2[1]}; below 7/858? {best2[0] < F(7,858)}")
print()

# 6c. The decisive structural question: can meas(G_C) -> 0 along a family?
print("6c. DECISIVE: drop-6 core with the SMALL speeds also scaled -- does any")
print("    primitive 12-set drive meas(G_C) strictly below 7/858? Greedy descent.")
random.seed(11111)
cur = list(C_drop6)
cur_m = lonely_measure(cur)
improved = True
rounds = 0
while improved and rounds < 200:
    improved = False
    rounds += 1
    # try replacing one speed by a nearby value
    for idx in range(12):
        for newv in range(1, 200):
            if newv in cur: continue
            cand_set = sorted(cur[:idx] + [newv] + cur[idx+1:])
            if len(set(cand_set)) != 12: continue
            m = lonely_measure(cand_set)
            if m < cur_m:
                cur, cur_m = cand_set, m
                improved = True
                break
        if improved: break
print(f"   greedy-descent final meas(G_C)={cur_m}={float(cur_m):.8f}")
print(f"     at C={cur}; below 7/858? {cur_m < F(7,858)}")
print()
