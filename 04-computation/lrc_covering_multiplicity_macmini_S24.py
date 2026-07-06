#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S24 -- HYP-4532: FRESH CONNECTION -- the COVERING-MULTIPLICITY
second-moment structure of the danger arcs.

Per opus-S114 (safe routes are reformulations; need a real tool) + codex-S243
(Markov is indirect: LRC = anti-approximation), a fresh quantitative angle:
mu(t) = #{i : ||v_i t|| < beta} is the COVERING MULTIPLICITY.  safe = Leb{mu=0}.
Int(mu) = 12*2beta = 24beta = 1.92 (fixed).  The danger arcs are a
COVERING-SYSTEM (arcs at m/v_i); safe=0 <=> they cover.

TEST whether the second moment gives a handle:
 (1) Var(mu) = Int(mu^2) - 1.92^2; Int(mu^2) = sum_{i,j} pairwise-overlap =
     the theta second-moment (THM-594, gcd-controlled).  Does the AP MINIMIZE
     Var(mu) (evenest cover = equioscillation)?
 (2) The covering-system view: is the AP the 'exact-cover-like' (min overlap)
     config?  Mirsky-Newman: no exact cover with distinct moduli.
 (3) Does a second-moment / Paley-Zygmund bound relate Var(mu) to safe>0?
     (honest: 2nd moment gives UPPER bounds on safe; test if the AP is the
     unique Var-minimizer, a fresh rigidity handle even if not a lower bound.)
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

BETA = 2 / 25

def dist(x):
    return abs(x - round(x))

def multiplicity_stats(W, beta=BETA, grid=60000):
    """Int(mu), Int(mu^2), Var(mu), safe=Leb{mu=0}, via a fine grid."""
    s0 = s1 = s2 = 0
    for j in range(grid):
        t = (j + 0.5) / grid
        mu = sum(1 for v in W if dist(v * t) < beta)
        s1 += mu; s2 += mu * mu
        if mu == 0:
            s0 += 1
    n = grid
    Emu = s1 / n; Emu2 = s2 / n
    return dict(safe=s0 / n, Emu=Emu, Emu2=Emu2, Var=Emu2 - Emu * Emu)

def pairwise_overlap_sum(W, beta=BETA):
    """sum_{i!=j} Leb{both danger} + diagonal.  Int(mu^2) exactly-ish via THM-594:
    overlap(v,w) at radius beta ~ (2beta)^2 * (period structure) + gcd resonance.
    Here: numeric estimate of the pairwise overlap sum (the theta 2nd moment)."""
    # Int(mu^2) = sum_i Leb(danger_i) + sum_{i!=j} Leb(danger_i cap danger_j)
    # = 12*2beta + pairwise.  Estimate pairwise via gcd-controlled overlap.
    # exact-ish: overlap(v,w) = Leb{||vt||<b and ||wt||<b}.
    def overlap(v, w):
        g = gcd(v, w)
        # coprime part: v/g, w/g; overlap ~ (2b)^2 for large lcm, + resonance
        # numeric on a modest grid scaled to lcm
        L = (v * w) // g
        grid = min(4000, max(400, 4 * L))
        c = 0
        for j in range(grid):
            t = (j + 0.5) / grid
            if dist(v * t) < BETA and dist(w * t) < BETA:
                c += 1
        return c / grid
    tot = 0.0
    for i in range(len(W)):
        for j in range(len(W)):
            if i == j:
                tot += 2 * BETA
            else:
                tot += overlap(W[i], W[j])
    return tot

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def main():
    print("=" * 80)
    print("COVERING MULTIPLICITY mu(t)=#{i:||v_i t||<2/25}; safe=Leb{mu=0}")
    print("Int(mu) = 24*beta = 1.92 fixed.  Does the AP MINIMIZE Var(mu)?")
    print("=" * 80)
    AP = tuple(range(1, 13))
    st = multiplicity_stats(AP)
    print(f"  AP {{1..12}}: safe={st['safe']:.5f}, E[mu]={st['Emu']:.4f}, "
          f"E[mu^2]={st['Emu2']:.4f}, Var(mu)={st['Var']:.4f}")
    random.seed(24)
    rows = []
    for _ in range(400):
        W = primitive(tuple(sorted(random.sample(range(1, 40), 12))))
        if len(set(W)) != 12:
            continue
        st = multiplicity_stats(W, grid=30000)
        rows.append((st['Var'], st['safe'], W))
    rows.sort()  # by Var ascending
    print(f"\n  {len(rows)} non-AP families, sorted by Var(mu) ASC:")
    print(f"  {'Var(mu)':>10} {'safe':>9}   family")
    for V, sf, W in rows[:5]:
        print(f"  {V:>10.4f} {sf:>9.5f}   {list(W)}")
    print("  ...")
    for V, sf, W in rows[-3:]:
        print(f"  {V:>10.4f} {sf:>9.5f}   {list(W)}")
    apVar = multiplicity_stats(AP)['Var']
    below = sum(1 for V, sf, W in rows if V < apVar)
    print(f"\n  AP Var(mu) = {apVar:.4f}; non-AP families with SMALLER Var: {below}/{len(rows)}")
    print(f"  ({'AP is the Var-MINIMIZER (evenest cover = equioscillation)' if below==0 else 'AP is NOT the Var-minimizer'})")
    # correlation Var vs safe
    import statistics
    Vs = [V for V, sf, W in rows]; sfs = [sf for V, sf, W in rows]
    mV, ms = statistics.mean(Vs), statistics.mean(sfs)
    num = sum((a-mV)*(b-ms) for a,b in zip(Vs,sfs))
    den = (sum((a-mV)**2 for a in Vs)*sum((b-ms)**2 for b in sfs))**0.5
    print(f"  corr(Var(mu), safe) = {num/den:+.3f} "
          f"(does more variance => more safe/gap? the covering handle)")
    print()
    print("  ASSESSMENT: 2nd moment gives UPPER bounds on safe (Cauchy-Schwarz:")
    print("  Leb{mu>=1} >= E[mu]^2/E[mu^2]); it does NOT lower-bound safe (the")
    print("  density floor needs the full signed theta-series / Selberg minorant).")
    print("  BUT if the AP uniquely minimizes Var(mu), that is a fresh RIGIDITY")
    print("  handle (evenest-cover = equioscillation = AP), independent of safe.")

if __name__ == "__main__":
    main()
