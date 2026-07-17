# opus-2026-07-17-S340 -- HYP-7230: THE 7-WALL CRUMB LAW (empirical, exact).
# At k = 7 comparable speeds the union bound saturates; loneliness in a
# window needs pair-overlap mass. This measures, over random 7-blocks and
# fee-scale windows: (a) how often the window is genuinely dead (uncovered
# empty) vs rescued; (b) in rescued windows, the total pair-overlap mass
# vs the union-bound deficit; (c) WHICH pair supplies the largest in-window
# overlap (gcd? ratio proximity? position?) -- the empirical law for the
# last wall (mac-mini's lane; my sawtooth/beat structure is the instrument).
from fractions import Fraction
from math import floor, gcd
import random, itertools
from collections import Counter

F = Fraction

def teeth_in(x, a, b):
    w = F(1, 14 * x)
    out = []
    for j in range(floor((a - w) * x), floor((b + w) * x) + 2):
        lo, hi = max(F(j, x) - w, a), min(F(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def union_len(ivs):
    ivs = sorted(ivs)
    tot = F(0); cur = None
    for lo, hi in ivs:
        if cur is None or lo > cur[1]:
            if cur: tot += cur[1] - cur[0]
            cur = [lo, hi]
        else:
            cur[1] = max(cur[1], hi)
    if cur: tot += cur[1] - cur[0]
    return tot

def inter_len(u, v):
    out, i, j = F(0), 0, 0
    u, v = sorted(u), sorted(v)
    while i < len(u) and j < len(v):
        a, b = max(u[i][0], v[j][0]), min(u[i][1], v[j][1])
        if a < b: out += b - a
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out

random.seed(340)
print("THE 7-WALL CRUMB LAW (exact, fee-scale windows):")
dead = live = 0
rescue_gcd = Counter()
rescue_adj = Counter()
ratios = []
for trial in range(300):
    m = random.randint(20, 200)
    B = sorted(random.sample(range(m, 13 * m), 7))
    V = B[-1]
    L = F(3, 2 * V)                     # the fee scale (BlockSplitLift eps)
    a0 = F(random.randint(0, 10**6), 10**6 + 3)
    b0 = a0 + L
    per = [teeth_in(x, a0, b0) for x in B]
    allt = [iv for t in per for iv in t]
    covered = union_len(allt)
    if covered >= L:
        dead += 1
        continue
    live += 1
    # the crumb: S1 - covered = total overlap mass reclaimed
    S1 = sum(union_len(t) for t in per)
    crumb = S1 - covered
    # which pair supplies the max in-window overlap
    best = None
    for (i, j) in itertools.combinations(range(7), 2):
        ov = inter_len(per[i], per[j])
        if best is None or ov > best[0]: best = (ov, i, j)
    ov, i, j = best
    if ov > 0:
        g = gcd(B[i], B[j])
        rescue_gcd[min(g, 5) if g < 5 else '5+'] += 1
        rescue_adj['adjacent' if j == i + 1 else 'non-adj'] += 1
        ratios.append(float(B[j]) / B[i])
print(f"  300 windows at L = 3/(2V): dead (fully covered) = {dead}, "
      f"live = {live}  ({100*live/300:.0f}% live)")
print(f"  rescuing-pair gcd distribution: {dict(rescue_gcd)}")
print(f"  rescuing-pair adjacency: {dict(rescue_adj)}")
if ratios:
    ratios.sort()
    print(f"  rescuing-pair ratio: median {ratios[len(ratios)//2]:.2f}, "
          f"p90 {ratios[int(len(ratios)*0.9)]:.2f}")
print()
print("  scan: dead-window FRACTION per block (how much of the circle is dead")
print("  at fee scale -- the window-choice freedom upstream must beat this):")
fracs = []
for trial in range(25):
    m = random.randint(20, 120)
    B = sorted(random.sample(range(m, 13 * m), 7))
    V = B[-1]
    L = F(3, 2 * V)
    deadc = 0; N = 140
    for s in range(N):
        a0 = F(s, N)
        per = [teeth_in(x, a0, a0 + L) for x in B]
        if union_len([iv for t in per for iv in t]) >= L: deadc += 1
    fracs.append(deadc / N)
fracs.sort()
print(f"  25 blocks x 140 positions: dead fraction median "
      f"{fracs[len(fracs)//2]:.2f}, max {fracs[-1]:.2f}, min {fracs[0]:.2f}")
print()
print("  -> if max dead-fraction < 1: window choice ALWAYS exists at fee")
print("     scale; the 7-wall crumb is a POSITIONAL pigeonhole, not a")
print("     per-window analytic bound.")
