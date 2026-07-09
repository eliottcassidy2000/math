#!/usr/bin/env python3
r"""
lrc14_composition_adversary_monad_S2.py   (monad-explorer-2026-07-09-S2, HYP-5717 part 2)

ADVERSARIAL HARDENING of the low-anchor composition (MISTAKE-102 discipline: random
sampling misses structured extremizers).  The structured adversaries:

  1. NEAR-DILATED-AP primitive covering families: 14*{1..13} is covering with M = 1/14
     EXACTLY (dilation invariance) but NON-PRIMITIVE (gcd 14; the reduction divides it
     out to the non-covering {1..13}).  The honest adversary is its PRIMITIVE
     perturbations: one entry +-1/+-2 (gcd 1, covering kept, M near 1/14 with a genuine
     small margin).  If the composition fails there, THAT is the true hard class.
  2. The @91 7-structured covering cluster (mac-mini's hembed counterexample family).
  3. Dilated-AP patches at c = 7, 14 with 2 entries perturbed.
  4. Hill-climb MINIMIZING the composition's best margin over primitive covering sets.

Also records per-failure WHICH piece fails: no drift-valid gap at any j (existence),
or gap exists but G_P conflict (intersection).
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import numpy as np

_s = open('/home/bigo/math/04-computation/lrc14_clamp_port_composition_monad_S2.py').read()
exec(_s[:_s.rfind('if __name__')])

def is_primitive(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    return g == 1

def best_composition_margin(v, S_L_list=None, detail=False):
    """max over (j, gap, split) of the drift-condition margin g(1-r) - 1/7 - 2ra,
    over j with tau in G_P.  Returns (margin, fired, whichfail)."""
    v = sorted(v)
    Vmax = v[-1]
    if S_L_list is None:
        S_L_list = sorted(set([Vmax // 4, Vmax // 3, Vmax // 2, 2 * Vmax // 3,
                               int(0.8 * Vmax), Vmax - 1]))
    best = None
    sawgap = False
    for S_L in S_L_list:
        E_L = [Vmax - vi for vi in v if Vmax - vi <= S_L]
        P = [vi for vi in v if Vmax - vi > S_L]
        if not E_L:
            continue
        r = F(S_L, Vmax)
        if r >= 1:
            continue
        for j in range(Vmax):
            teeth = sorted(set(F((e * j) % Vmax, Vmax) for e in E_L))
            nt = len(teeth)
            for t in range(nt):
                lo = teeth[t]
                hi = teeth[t + 1] if t + 1 < nt else teeth[0] + 1
                g = hi - lo
                a = lo
                if a + g > 1:
                    continue
                m = g * (1 - r) - F(1, 7) - 2 * r * a
                if m > 0:
                    sawgap = True
                    phi = a + g / 2
                    tau = (j + phi) / Vmax
                    if in_G_P(tau, P):
                        if best is None or m > best:
                            best = m
    if best is not None:
        return float(best), True, "fires"
    return -1.0, False, ("gap-but-GP-conflict" if sawgap else "no-valid-gap")

if __name__ == "__main__":
    print("=" * 100)
    print("ADVERSARY 1 -- primitive perturbations of 14*{1..13} (covering, near-tight)")
    print("=" * 100)
    base = [14 * i for i in range(1, 14)]
    rows = []
    for idx in range(13):
        for d in (-2, -1, 1, 2):
            v = sorted(set(base[:idx] + [base[idx] + d] + base[idx + 1:]))
            if len(v) != 13 or not is_primitive(v) or not is_covering(v):
                continue
            m, fired, why = best_composition_margin(v)
            rows.append((fired, why, v, m))
    fails = [r for r in rows if not r[0]]
    print(f"  tested {len(rows)} primitive covering perturbations: "
          f"fired {len(rows)-len(fails)}, FAILED {len(fails)}")
    for fired, why, v, m in fails[:10]:
        print(f"    FAIL ({why}): {v}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("ADVERSARY 2 -- the @91 covering cluster + 7-structured relatives")
    print("=" * 100)
    for name, v in [("@91 7-struct", [9,16,24,33,40,47,54,62,65,70,77,84,91]),
                    ("7*{1..12}+91", sorted([7*i for i in range(1,13)] + [91])),
                    ("7*{1..13} patched", sorted([7*i for i in range(1,13)] + [92]))]:
        v = sorted(set(v))
        if len(v) != 13:
            print(f"  [{name}] skipped (|v| = {len(v)})")
            continue
        pc = is_primitive(v) and is_covering(v)
        m, fired, why = best_composition_margin(v)
        print(f"  [{name}] primitive-covering={pc}: "
              f"{'FIRES margin ' + format(m, '.5f') if fired else 'FAIL: ' + why}")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("ADVERSARY 3 -- hill-climb MINIMIZING the best composition margin (primitive covering)")
    print("=" * 100)
    rng = random.Random(2026070902)
    cur = [14 * i for i in range(1, 13)] + [183]
    cur = sorted(set(cur))
    def fitness(v):
        if len(v) != 13 or not is_primitive(v) or not is_covering(v):
            return None
        m, fired, why = best_composition_margin(v, S_L_list=[max(v)//3, max(v)//2, 2*max(v)//3, max(v)-1])
        return (m if fired else -1.0, fired, why)
    fit = fitness(cur)
    if fit is None:
        cur = sorted(set([14*i for i in range(1,14)] ) - {14*13} | {181, 183} )[:13]
        fit = fitness(cur)
    worst = (fit[0] if fit else 1e9, list(cur), fit[2] if fit else "?")
    tested = 0
    for step in range(150):
        v2 = list(cur)
        i = rng.randrange(13)
        v2[i] = max(2, v2[i] + rng.choice([-3, -2, -1, 1, 2, 3, 7, -7, 14, -14]))
        v2 = sorted(set(v2))
        f2 = fitness(v2)
        if f2 is None:
            continue
        tested += 1
        if f2[0] < worst[0]:
            worst = (f2[0], v2, f2[2])
            print(f"  step {step}: worst margin {f2[0]:.5f} ({f2[2]}) at {v2}")
            sys.stdout.flush()
        if f2[0] < (fit[0] if fit else 1e9) or rng.random() < 0.3:
            cur, fit = v2, f2
    print(f"  hill-climb ({tested} valid): worst = {worst[0]:.5f} ({worst[2]}) at {worst[1]}")
    m, fired, why = best_composition_margin(worst[1])
    print(f"  full-split recheck of worst: {'FIRES margin ' + format(m, '.5f') if fired else 'FAIL: ' + why}")
    sys.stdout.flush()
