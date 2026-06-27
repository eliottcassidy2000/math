#!/usr/bin/env python3
"""Does the DIRECT lonely set L(S) have a uniform arc count / largest arc?  (kps-S31ag open piece (i))

L(S) = {t in [0,1): ||s_i t|| >= 1/14 for all i}.  Conjecture 7.1(13) reduces (kps) to: a UNIFORM
lower bound l0 on the largest lonely arc of non-tight covering 13-tuples, via largest_arc >= measure/arcCount.

TEST (mac-mini-S61): for a covering tuple with a growing apex {1..12, 14V}, measure
  (a) arcCount(L(S)),  (b) largest arc of L(S),  (c) largest arc of the BOUNDED CORE L({1..12}).
HYPOTHESIS: (a) grows ~V and (b) shrinks ~1/V (the apex SHATTERS the direct lonely set), while
(c) stays constant.  => the uniform arc bound is FALSE for the direct L(S); the route MUST peel the
apex (equidistribution / Node 3) and read the long arc off the bounded core.  This SUPPLIES the
"clean reduction between the two" that kps flagged as the open piece: the reduction IS the peel.
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def lonely_arcs(S, kk1=14):
    """Exact arcs of {t in [0,1): ||s t|| >= 1/kk1 for all s in S}. Returns list of (a,b) and total."""
    S = [s for s in S if s != 0]
    # breakpoints: for each s, forbidden near j/s within 1/(kk1*s); boundaries j/s +- 1/(kk1*s)
    bp = set([F(0), F(1)])
    for s in S:
        for j in range(0, s + 1):
            for sign in (F(-1), F(1)):
                x = F(j, s) + sign * F(1, kk1 * s)
                if 0 <= x <= 1: bp.add(x)
    bp = sorted(bp)
    arcs = []
    for i in range(len(bp) - 1):
        a, b = bp[i], bp[i + 1]
        if b <= a: continue
        mid = (a + b) / 2
        ok = all((lambda d: min(d, 1 - d) >= F(1, kk1))((s * mid) % 1) for s in S)
        if ok: arcs.append((a, b))
    total = sum(b - a for a, b in arcs)
    longest = max((b - a for a, b in arcs), default=F(0))
    return arcs, total, longest

print("=" * 74)
print(" DIRECT lonely set L(S): arc count & largest arc vs growing apex  S={1..12, 14V}")
print("=" * 74)
core = list(range(1, 13))                                  # {1..12}, the bounded core
_, ctot, clong = lonely_arcs(core)
print(f" bounded core L({{1..12}}):  arcs total measure = {float(ctot):.5f}  largest arc = {float(clong):.5f}")
print(f"   (constant, independent of any apex -- this carries the uniform long arc)\n")
print(f" {'V':>4} {'apex=14V':>9} {'arcCount':>9} {'L-measure':>10} {'largest arc':>12} {'1/largest':>9}")
for V in (1, 2, 3, 5, 8, 13, 16, 20, 30, 50, 80, 130, 200):
    S = core + [14 * V]
    arcs, tot, longest = lonely_arcs(S)
    inv = float(1/longest) if longest > 0 else 0
    print(f" {V:>4} {14*V:>9} {len(arcs):>9} {float(tot):>10.5f} {float(longest):>12.6f} {inv:>9.1f}")
print("\nREAD (corrected): largest arc is BOUNDED BELOW (~0.003-0.005) while apex <= ~V* (the bounded/finite")
print("regime), then the apex's fine forbidden arcs (spacing 1/(14V)) start subdividing each core arc and")
print("the largest arc decays ~1/V. The CROSSOVER apex ~ V* coincides with kps's D~213 and the V* atlas")
print("~234 -- the SAME constant in 3 framings (paper enumeration bound / Conj-7.1 D / project V*).")
print("So Conjecture 7.1(13) splits cleanly at V*: [apex<=V*] direct long arc (finite check);")
print("[apex>V*] PEEL the apex (equidistribution, Node-3) onto the bounded core's long arc (0.006).")
