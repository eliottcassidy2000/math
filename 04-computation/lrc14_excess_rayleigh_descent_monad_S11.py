#!/usr/bin/env python3
r"""
lrc14_excess_rayleigh_descent_monad_S11.py   (monad-explorer-2026-07-07-S11, HYP-5117)

ADVERSARIAL DESCENT on the excess-mass Rayleigh certificate (Part 5 of the cubic-moment
gate session).  The certificate: mu_{1/7}(E) >= RAY(E) = a^T M^{-1} a where
a_i = E[V_{theta_i}], M_ij = E[V_{theta_i} V_{theta_j}], V_t = sum_r (g_r - t)_+
(vanishes exactly on Bad).  Exact rationals per shape.

HUNT: minimize RAY(E) over 8-element shapes with JUMP moves (MISTAKE-119 discipline)
+ all structured families from the fleet's adversary corpus.  Also: k=13 spot check
(bar m_P = 0.0565), and the mu-vs-RAY scatter to see how much of the truth the
certificate captures.

Question to answer: is inf_E RAY(E) attained at the AP (a moment-level R2), and is
the infimum above the k=8 bar 0.675 (klein) / 0.6185 (kps)?
"""
import random, sys
from fractions import Fraction as F
from math import gcd

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])

WIDTHS8 = [F(1,7), F(31,210), F(8,49), F(6,35), F(19,105), F(1,5), F(23,105), F(1,4), F(2,7)]

def normalize(E):
    E = sorted(set(E))
    m0 = E[0]
    E = [e - m0 for e in E]
    g = 0
    for e in E[1:]:
        g = gcd(g, e)
    if g > 1:
        E = [e // g for e in E]
    return E

def ray_floor(E, widths=WIDTHS8):
    aV, MV, m3V, vmaxV = excess_moments(E, widths)
    pz = aV[0] * aV[0] / MV[0][0] if MV[0][0] else F(0)
    at3 = 1 - atom_zero_bound_3mom(aV[0], MV[0][0], m3V, vmaxV[0])
    ray, coef = rayleigh_floor(aV, MV)
    return ray, pz, at3

if __name__ == "__main__":
    rng = random.Random(51120260707)
    k = 8
    print("=" * 100)
    print("PART 5 -- DESCENT ON RAY (k=8). bars: klein 0.675 / kps T_8 0.6185. AP_8 RAY = 0.8214")
    print("=" * 100)
    cache = {}
    def RAY(E):
        t = tuple(E)
        if t not in cache:
            try:
                cache[t] = ray_floor(E)
            except Exception:
                cache[t] = (None, None, None)
        return cache[t]

    # seed pool: structured adversaries
    seeds = [
        list(range(8)),
        [0,1,2,3,4,5,6,8], [0,1,2,3,4,5,6,13], [0,1,2,3,5,6,7,9],
        [0,2,4,6,8,10,11,12], [0,2,4,6,7,8,10,12],
        [0,1,2,3,40,41,42,43], [0,1,2,3,10,11,12,13],
        [0,1,2,3,4,5,6,40], [0,1,2,4,8,16,32,64],
        [0,1,3,7,12,20,30,44], [0,1,2,3,4,5,7,14], [0,3,6,9,12,15,18,22],
        [0,1,2,3,4,5,12,13],   # AP + double far pair
        [0,1,1+7,2+7,3+7,4+7,5+7,6+7],  # shifted block
        [0,7,14,21,28,35,42,49],        # dilated AP d=7 (== AP by dilation inv)
        [0,1,2,4,5,6,8,9],              # near-AP with holes
        [0,1,3,4,6,7,9,10],             # d=3 pairs
        [0,2,3,5,7,8,10,12],
        [0,1,2,3,5,7,11,13],            # primes-ish
    ]
    best = (None, None)
    results = []
    for E in seeds:
        E = normalize(E)
        r, pz, at3 = RAY(E)
        if r is None: continue
        results.append((float(r), E))
        if best[0] is None or r < best[0]:
            best = (r, E)
        print(f"  seed {str(E):42s} RAY={float(r):.4f} PZ={float(pz):.4f} AT3={float(at3):.4f}")
        sys.stdout.flush()

    # random + jump descent
    cur = normalize(list(range(8)))
    curf = float(RAY(cur)[0])
    for step in range(200):
        move = rng.random()
        E2 = list(cur)
        if move < 0.5:
            i = rng.randrange(k)
            E2[i] = rng.randint(0, 60)
        elif move < 0.8:
            i, j = rng.randrange(k), rng.randrange(k)
            E2[i] = E2[j] + rng.choice([1, 2, -1, -2, 7, 14])
        else:
            E2 = [rng.randint(0, 60) for _ in range(k)]
        E2 = normalize(E2)
        if len(E2) != k:
            continue
        r = RAY(E2)[0]
        if r is None:
            continue
        rf = float(r)
        results.append((rf, E2))
        if best[0] is None or r < best[0]:
            best = (r, E2)
            print(f"  step {step}: NEW WORST RAY={rf:.4f} at {E2}")
            sys.stdout.flush()
        if rf < curf or rng.random() < 0.3:
            cur, curf = E2, rf
    print(f"\nWORST RAY over {len(cache)} shapes: {float(best[0]):.4f} at {best[1]}")
    print(f"  exact worst RAY = {best[0]}")
    ranked = sorted(results)[:10]
    print("  10 lowest RAY shapes:")
    for rf, E in ranked:
        print(f"    {rf:.4f}  {E}")

    print()
    print("=" * 100)
    print("PART 6 -- k=13 SPOT CHECK (bar m_P = 0.0565; klein ledger 0.056; AP_13 mu = 0.4425)")
    print("=" * 100)
    widths13 = [F(1,7), F(31,210), F(8,49), F(6,35), F(1,5)]
    for name, E in [("AP_13", list(range(13))),
                    ("GW (12->24)", [0,1,2,3,4,5,6,7,8,9,10,12,23]),
                    ("2AP interlace", [0,2,4,6,8,10,12,14,16,18,20,22,23]),
                    ("two-block 7+6", [0,1,2,3,4,5,6,30,31,32,33,34,35]),
                    ("crux 2:1 mix", [0,2,4,6,8,10,12,14,16,18,20,21,22])]:
        E = normalize(E)
        aV, MV, m3V, vmaxV = excess_moments(E, widths13)
        pz = aV[0] * aV[0] / MV[0][0]
        at3 = 1 - atom_zero_bound_3mom(aV[0], MV[0][0], m3V, vmaxV[0])
        ray, coef = rayleigh_floor(aV, MV)
        sd = ShapeData(E, max_mixed=1)
        print(f"{name:18s} mu={1-float(sd.pbad):.4f} | PZ={float(pz):.4f} "
              f"AT3={float(at3):.4f} RAY={float(ray):.4f}  (bar 0.0565)")
        sys.stdout.flush()
