# opus-2026-07-17-S362 -- HYP-7450: TAXONOMY OF THE SURVIVING 11%.
# THM-1040: the seven-moduli sieve kills ~89% of comparable 13-families.
# What do the survivors look like, and which existing tool covers each?
from fractions import Fraction as F
from math import gcd, floor, comb
import random, itertools
LAM = F(1,14)
MODULI = [8,9,10,11,12,13,14]
def blocks_all(V): return all(any(v % q == 0 for v in V) for q in MODULI)

def teeth(x, lo, hi):
    w = LAM / x; out = []
    for j in range(floor((lo - w)*x), floor((hi + w)*x) + 2):
        a, b = max(F(j,x) - w, lo), min(F(j,x) + w, hi)
        if a < b: out.append((a,b))
    return out
def inter(u, v):
    out, i, j = [], 0, 0
    while i < len(u) and j < len(v):
        a, b = max(u[i][0], v[j][0]), min(u[i][1], v[j][1])
        if a < b: out.append((a,b))
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out
def mu(S):
    cur = teeth(S[0], F(0), F(1))
    for x in S[1:]: cur = inter(cur, teeth(x, F(0), F(1)))
    return sum(b-a for a,b in cur)

def bonf5(V):
    """1 - S1 + S2 - S3 + S4 - S5, exact."""
    tot = F(1)
    for k in range(1, 6):
        Sk = F(0)
        for C in itertools.combinations(V, k): Sk += mu(list(C))
        tot += (-1)**k * Sk
    return tot

def near_ap(V, tol=2):
    d = [V[i+1]-V[i] for i in range(len(V)-1)]
    return max(d) - min(d) <= tol
def is_dilate(V):
    g = 0
    for v in V: g = gcd(g, v)
    return g > 1

print("TAXONOMY OF FAMILIES BLOCKING ALL SEVEN MODULI (the surviving ~11%)")
random.seed(362)
blockers = []
while len(blockers) < 18:
    m = random.randint(23, 120)
    V = sorted(random.sample(range(m, 13*m), 13))
    if blocks_all(V): blockers.append(V)
print(f"  collected {len(blockers)} blockers; classifying each against existing tools\n")
cnt = {'dilate':0, 'nearAP':0, 'bonf5':0, 'survives':0}
for V in blockers:
    d = is_dilate(V); a = near_ap(V)
    b = bonf5(V)
    covered = d or a or (b > 0)
    if d: cnt['dilate'] += 1
    if a: cnt['nearAP'] += 1
    if b > 0: cnt['bonf5'] += 1
    if not covered: cnt['survives'] += 1
    tag = []
    if d: tag.append('dilate')
    if a: tag.append('nearAP')
    tag.append(f"BONF5={float(b):+.4f}")
    print(f"    {str(V[:4])[:-1]}...]  {'  '.join(tag)}"
          f"{'   <-- SURVIVES ALL' if not covered else ''}")
n = len(blockers)
print()
print(f"  of {n} blockers:  dilate {cnt['dilate']}, near-AP {cnt['nearAP']},"
      f" BONF5>0 {cnt['bonf5']}")
print(f"  SURVIVING ALL THREE TOOLS: {cnt['survives']}/{n}")
print()
print(f"  net: the sieve removes ~89%; of the ~11% remaining, BONF5 alone")
print(f"  certifies {100*cnt['bonf5']/n:.0f}% -- the true residual is the rest.")

print()
print("=" * 66)
print("THE THRESHOLD TEST: does BONF5 turn positive above a speed floor?")
print("  (survivors above clustered at SMALL speeds -- if there is a threshold")
print("   V0 with BONF5 > 0 for all blockers of min-speed >= V0, the TRUE")
print("   RESIDUAL IS BOUNDED, and a bounded residual is a finite census.)")
random.seed(3620)
for lo, hi in [(23, 70), (150, 300), (600, 900)]:
    got = []
    tries = 0
    while len(got) < 4 and tries < 4000:
        tries += 1
        m = random.randint(lo, hi)
        V = sorted(random.sample(range(m, 13*m), 13))
        if V[0] < lo: continue
        if blocks_all(V): got.append(V)
    vals = [float(bonf5(V)) for V in got]
    pos = sum(1 for v in vals if v > 0)
    vals.sort()
    print(f"    min-speed in [{lo:4d},{hi:4d}): {pos}/{len(vals)} have BONF5>0 ;"
          f" range [{vals[0]:+.4f}, {vals[-1]:+.4f}]")
print()
print("  READ: if the positive fraction rises to 1 with the speed floor, the")
print("  residual is a BOUNDED family set and closable by finite census;")
print("  if it plateaus below 1, the residual is unbounded and needs theory.")
