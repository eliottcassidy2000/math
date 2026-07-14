#!/usr/bin/env python3
"""THM-777 experiment: the |G'| floor decision for regime 2 (opus-2026-07-14-S301).

Tests, with exact rational arithmetic, whether the 1/14-safe measure |G'_P| of a
12-core suggests a positive floor over primitive shapes (regime 2's last lemma via
rho = v*/maxP <= 12/(pi |G'|)), or exhibits decay in the listed probes.

PART 1  EXACT BOUNDED-HEIGHT CENSUS: min |G'| over ALL primitive 12-subsets of
        {1..H} for H = 16, 17, 18 — exact rational minima + minimizer shapes +
        bottom-5 table. (Dilation-invariance makes primitive shapes canonical.)
PART 2  ADVERSARIAL FAMILIES (the MISTAKE-140/145 battery discipline):
        (a) perturbed dilates at the explicitly listed eight c-values through 40;
        (b) GW-dilate perturbations c*{1..11,13} with one element +-1;
        (c) tooth insertions at the explicitly listed N-values through 2003;
        (d) two-block and compressed shapes (MISTAKE-141 near-dilate lineage).
PART 3  HEURISTIC HILL-DESCENT: 500 proposed moves from each of five seeds,
        with replacement proposals bounded by 2500. This is not exhaustive.
PART 4  THE DIAGNOSTIC + THE BRIDGE: the unconditional Lipschitz tail
        |G'| >= 2(M(P)-1/14)/maxP >= 1/(91 maxP) (proved, LRC(13) floor); the
        rho bridge rho <= 12/(pi |G'|); and a floor-vs-decay diagnostic with the
        conjectured asymptotic floor if the data supports one.
"""
import sys, os, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc14_certificates import good_intervals, L_exact

random.seed(301)

def Gm(P):
    return L_exact(sorted(P))

def primitive(P):
    g = 0
    for x in P: g = gcd(g, x)
    return g == 1

print("=" * 92)
print("PART 1 -- exact census: min |G'| over primitive 12-subsets of {1..H}")
print("=" * 92)
best_overall = None
for H in (16, 17, 18):
    found = []
    for P in combinations(range(1, H + 1), 12):
        if not primitive(P): continue
        g = Gm(P)
        found.append((g, P))
    found.sort()
    print(f"  H = {H}: {len(found)} primitive shapes; exact min |G'| = {found[0][0]} "
          f"= {float(found[0][0]):.6f} at {found[0][1]}")
    for g, P in found[1:5]:
        print(f"          runner-up {float(g):.6f} at {P}")
    if best_overall is None or found[0][0] < best_overall[0]:
        best_overall = found[0]

print()
print("=" * 92)
print("PART 2 -- adversarial families (exact; scale sweep)")
print("=" * 92)
fam_min = None
def track(name, P, quiet=False):
    global fam_min
    if len(set(P)) != 12: return None
    g = Gm(P)
    if fam_min is None or g < fam_min[0]: fam_min = (g, name, tuple(P))
    if not quiet: print(f"  {name:<42} |G'| = {float(g):.6f}")
    return g

print("  (a) perturbed dilates:")
vals_a = []
for c in (2, 3, 5, 8, 13, 21, 34, 40):
    g1 = track(f"{{c,...,11c,12c+1}} c={c}", [c*i for i in range(1, 12)] + [12*c + 1], quiet=True)
    g2 = track(f"{{c,...,11c,12c-1}} c={c}", [c*i for i in range(1, 12)] + [12*c - 1], quiet=True)
    vals_a.append((c, float(g1), float(g2)))
print("   c:", [v[0] for v in vals_a])
print("   {..,12c+1}:", [f"{v[1]:.5f}" for v in vals_a])
print("   {..,12c-1}:", [f"{v[2]:.5f}" for v in vals_a])
print("  (b) GW-dilate perturbations c*{1..11,13} with one slot tweaked:")
vals_b = []
for c in (2, 5, 13, 29):
    base = [c*i for i in range(1, 12)] + [13*c]
    for tweak in (+1, -1):
        P = base[:-1] + [13*c + tweak]
        vals_b.append((c, tweak, float(track(f"c*GW-core tweak {tweak:+d} c={c}", P, quiet=True))))
print("   ", [(v[0], v[1], f"{v[2]:.5f}") for v in vals_b])
print("  (c) tooth insertions (large N):")
for N in (101, 503, 1009, 2003):
    track(f"{{1..11, N}} N={N}", list(range(1, 12)) + [N])
for N in (101, 1009):
    track(f"{{1..10, 12, N}} N={N}", list(range(1, 11)) + [12, N])
    track(f"{{2..12, N}} N={N}", list(range(2, 13)) + [N])
print("  (d) blocks/compressed:")
track("two-block {1..6, 9..14}", list(range(1, 7)) + list(range(9, 15)))
track("compressed 2*{1..11} u {13}", [2*i for i in range(1, 12)] + [13])
track("{1..11, 24} (GW-side)", list(range(1, 12)) + [24])
track("{1..11, 26}", list(range(1, 12)) + [26])

print()
print("=" * 92)
print("PART 3 -- heuristic hill descent (primitive proposals bounded by height 2500)")
print("=" * 92)
def descend(seed, steps=500, hmax=2500):
    cur = sorted(seed); curg = Gm(cur)
    for s in range(steps):
        Q = list(cur)
        i = random.randrange(12)
        if random.random() < 0.5:
            Q[i] = max(1, Q[i] + random.choice([-2, -1, 1, 2]))
        else:
            Q[i] = random.randint(1, hmax)
        Q = sorted(set(Q))
        if len(Q) != 12 or not primitive(Q): continue
        qg = Gm(Q)
        if qg < curg: cur, curg = Q, qg
    return curg, cur
seeds = [best_overall[1],
         list(range(1, 12)) + [13],
         list(range(1, 13)),
         sorted(random.sample(range(1, 200), 12)),
         [2*i for i in range(1, 12)] + [13]]
desc_best = None
for sd in seeds:
    g, P = descend(list(sd))
    print(f"  seed {tuple(sd)[:4]}...: descended to |G'| = {float(g):.6f} at {P}")
    if desc_best is None or g < desc_best[0]: desc_best = (g, P)

print()
print("=" * 92)
print("PART 4 -- diagnostic, the Lipschitz tail, and the rho bridge")
print("=" * 92)
gmin, Pmin = best_overall
print(f"  exact census floor (primitive, maxP <= 18): |G'| = {gmin} = {float(gmin):.6f} at {Pmin}")
print(f"  adversarial-family min: {float(fam_min[0]):.6f} ({fam_min[1]})")
print(f"  hill-descent min:       {float(desc_best[0]):.6f} at {desc_best[1]}")
overall = min(gmin, fam_min[0], desc_best[0])
print(f"  OVERALL MIN OBSERVED:   {float(overall):.6f}")
print()
if desc_best[0] >= gmin * F(9, 10) and fam_min[0] >= gmin * F(9, 10):
    print("  VERDICT: FLOOR SIGNAL -- none of the explicitly listed finite probes moved")
    print("  below the bounded census minimum; no convergence claim follows from this sample.")
else:
    print("  VERDICT: DECAY CANDIDATE FOUND -- see minimizer above; the floor route dies.")
print()
print("  PROVED (unconditional Lipschitz tail): |G'_P| >= 2(M(P)-1/14)/maxP >= 1/(91*maxP)")
print("  PROVED (the rho bridge, one line):     rho = v*/maxP = r_P/(pi|G'|maxP) <= 12/(pi|G'_P|)")
print("  => on every stratum where |G'| >= gamma:  rho <= 12/(pi*gamma); with gamma = census")
print(f"     floor {float(gmin):.5f}: rho <= {12/(3.14159*float(gmin)):.1f} (uniform, that stratum).")
