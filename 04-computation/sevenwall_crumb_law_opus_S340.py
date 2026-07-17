# opus-2026-07-17-S340 -- HYP-7230: THE 7-WALL CRUMB LAW.
# At k = 7 comparable speeds the union bound saturates; loneliness in a
# window needs pair-overlap mass.  The first section is sampled telemetry on
# random 7-blocks and fee-scale windows.  The second section replaces the old
# 140-start proxy by the exact circular positional statistic: if the closed
# danger union has proper circular component lengths ell_i, the measure of
# starts a whose oriented length-L window lies wholly in the union is
# sum_i max(ell_i-L,0).  The whole-circle case has measure one.
# Tournament analysis is deliberately not imposed on this positional step:
# pair orientation does not preserve circular component lengths or erosion.
# The faithful alternate vertices are the merged danger-union components.
from fractions import Fraction
from math import floor, gcd
import random, itertools
from collections import Counter

F = Fraction

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

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

def danger_circle_intervals(x):
    """Closed danger teeth for speed x, split at the 0/1 chart boundary."""
    w = F(1, 14 * x)
    out = []
    for j in range(x):
        lo, hi = F(j, x) - w, F(j, x) + w
        if lo < 0:
            out.append((lo + 1, F(1)))
            out.append((F(0), hi))
        else:
            out.append((lo, hi))
    return out

def circular_component_lengths(B):
    """Exact lengths of the closed circular components of the danger union.

    The first and last linear pieces are merged when they are the two chart
    pieces of one component through zero.  Return (whole_circle, lengths).
    """
    ivs = sorted(iv for x in B for iv in danger_circle_intervals(x))
    merged = []
    for lo, hi in ivs:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    linear_measure = sum((hi - lo for lo, hi in merged), F(0))
    if len(merged) == 1 and merged[0] == [F(0), F(1)]:
        return True, (F(1),)
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        wrap = merged[0][1] + (1 - merged[-1][0])
        lengths = [hi - lo for lo, hi in merged[1:-1]] + [wrap]
    else:
        lengths = [hi - lo for lo, hi in merged]
    require(all(0 < ell < 1 for ell in lengths), "bad proper circular component")
    require(sum(lengths, F(0)) == linear_measure, "wrap merge changed union measure")
    return False, tuple(sorted(lengths))

def exact_dead_start_measure(B, L):
    """Measure of starts a with the circular window [a,a+L] fully covered."""
    require(F(0) <= L <= 1, "window length must lie in [0,1]")
    whole, lengths = circular_component_lengths(B)
    if whole:
        return F(1), lengths
    measure = sum((max(ell - L, F(0)) for ell in lengths), F(0))
    require(F(0) <= measure <= 1, "dead-start measure left [0,1]")
    return measure, lengths

def fmt_fraction(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"

random.seed(340)
print("THE 7-WALL CRUMB LAW (sampled window telemetry; exact arithmetic per sample):")
dead = live = 0
best_gcd = Counter()
best_adj = Counter()
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
    # Which of all 21 pairs supplies the largest in-window overlap.  This is
    # descriptive telemetry, not an identification with the guaranteed close
    # adjacent pair in the sorted block.
    best = None
    for (i, j) in itertools.combinations(range(7), 2):
        ov = inter_len(per[i], per[j])
        if best is None or ov > best[0]: best = (ov, i, j)
    ov, i, j = best
    if ov > 0:
        g = gcd(B[i], B[j])
        best_gcd[min(g, 5) if g < 5 else '5+'] += 1
        best_adj['adjacent' if j == i + 1 else 'non-adj'] += 1
        ratios.append(float(B[j]) / B[i])
print(f"  sampled 300 windows at L = 3/(2V): dead (fully covered) = {dead}, "
      f"live = {live}  ({100*live/300:.0f}% live)")
best_observations = sum(best_gcd.values())
print(f"  best positive-overlap pair observations = {best_observations}/{live}")
print(f"  best-pair gcd distribution: {dict(best_gcd)}")
print(f"  best-pair adjacency: {dict(best_adj)}")
if ratios:
    ratios.sort()
    print(f"  best-pair ratio: median {ratios[len(ratios)//2]:.2f}, "
          f"p90 {ratios[int(len(ratios)*0.9)]:.2f}")
print(f"  best pair is nonadjacent in {best_adj['non-adj']}/{best_observations} "
      "positive-overlap samples;")
print("  it is not the guaranteed close adjacent pair.")
print()
print("EXACT POSITIONAL EROSION ON THE SAME 25 SAMPLED BLOCKS:")
print("  for proper circular components of lengths ell_i,")
print("  dead-start measure = sum_i max(ell_i - L, 0);")
print("  the whole-circle danger union is handled separately with measure 1.")
grid_fracs = []
exact_measures = []
for trial in range(25):
    m = random.randint(20, 120)
    B = sorted(random.sample(range(m, 13 * m), 7))
    V = B[-1]
    L = F(3, 2 * V)
    exact_dead, component_lengths = exact_dead_start_measure(B, L)
    exact_measures.append(exact_dead)
    deadc = 0; N = 140
    for s in range(N):
        a0 = F(s, N)
        per = [teeth_in(x, a0, a0 + L) for x in B]
        if union_len([iv for t in per for iv in t]) >= L: deadc += 1
    grid_fracs.append(F(deadc, N))
exact_measures.sort()
grid_fracs.sort()
print("  exact dead-start measures (sorted):")
print("  " + ",".join(fmt_fraction(value) for value in exact_measures))
print(f"  exact sampled-block summary: positive={sum(value > 0 for value in exact_measures)}/25, "
      f"median={fmt_fraction(exact_measures[len(exact_measures)//2])} "
      f"({float(exact_measures[len(exact_measures)//2]):.6f}), "
      f"max={fmt_fraction(exact_measures[-1])} ({float(exact_measures[-1]):.6f}), "
      f"min={fmt_fraction(exact_measures[0])} ({float(exact_measures[0]):.6f})")
print("  sampled 140-start grid frequencies (not circle measures or bounds):")
print(f"  median={float(grid_fracs[len(grid_fracs)//2]):.2f}, "
      f"max={float(grid_fracs[-1]):.2f}, min={float(grid_fracs[0]):.2f}")
print()
print("  -> each of these 25 sampled blocks has exact dead-start measure < 1,")
print("     so each sampled block admits a fee-scale live start.")
print("  -> this is not a uniform 7% theorem and does not yet provide the")
print("     upstream window-choice freedom needed to close the seven-wall.")
print("  -> carrier audit: merged circular components preserve positional")
print("     erosion; a runner-pair tournament does not preserve their lengths.")
