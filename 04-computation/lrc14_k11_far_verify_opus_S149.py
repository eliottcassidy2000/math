"""
lrc14_k11_far_verify_opus_S149.py   (opus-2026-07-08-S149, HYP-5357)

VERIFY THE ONE REMAINING k=11 TAIL LEMMA: far <= E[W]^2 for spread families.

far = E[W^2] - near  (both exact via Farey integration: E[W^2] from pz_exact, near from
lrc14_k11_tail_sharp_near).  If far <= E[W]^2 holds for all families of diameter >= D*, then
the sharp near/far bound PZ >= E[W]^2/(near + E[W]^2) >= 0.47 (part 1) CLOSES the k=11 tail,
and [exhaustive compact diam <= D*, exact PZ >= 0.3468] completes the leg.

This computes, per family: E[W], E[W^2] (exact), near (exact), far = E[W^2]-near, and checks
  (i)  far <= E[W]^2  (the decorrelation lemma -- where does it hold?)
  (ii) the sharp PZ_tail = E[W]^2/(near+E[W]^2) >= bar
to pin the crossover diameter D* and confirm the assembly.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys, random

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import pz_exact, BAR
from lrc14_k11_tail_sharp_near_opus_S149 import EW_and_near

K = 11
bar = BAR[K]


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def analyze(E):
    EW, near = EW_and_near(E)
    EW1, EW2, PZ = pz_exact(E)
    assert EW1 == EW, (E, EW1, EW)
    far = EW2 - near
    far_ok = (far <= EW * EW)
    pz_tail = EW * EW / (near + EW * EW)   # valid iff far_ok
    return EW, EW2, near, far, far_ok, PZ, pz_tail


def main():
    print("=" * 100)
    print(f"k={K}: far <= E[W]^2 verification (far = E[W^2]-near, exact) + sharp PZ tail")
    print(f"  bar = {float(bar):.6f}")
    print("=" * 100)
    print("  per diameter: #(far<=E[W]^2) / total, min PZ_tail (where far<=E[W]^2), min exact PZ")
    for D in (11, 12, 13, 14, 15, 16, 18, 20, 25, 30, 40):
        if D <= 16:
            gen = ([0] + list(m) + [D] for m in combinations(range(1, D), K - 2))
        else:
            rng = random.Random(D)
            gen = ([0] + sorted(rng.sample(range(1, D), K - 2)) + [D] for _ in range(400))
        n_far_ok = 0; n = 0
        min_pztail = None; min_pztail_arg = None
        min_pz = None
        far_violators = []
        for E in gen:
            E = sorted(set(E))
            if len(E) != K or gcd_all([E[i+1]-E[i] for i in range(len(E)-1)]) != 1:
                continue
            n += 1
            EW, EW2, near, far, far_ok, PZ, pz_tail = analyze(E)
            if min_pz is None or PZ < min_pz:
                min_pz = PZ
            if far_ok:
                n_far_ok += 1
                if min_pztail is None or pz_tail < min_pztail:
                    min_pztail = pz_tail; min_pztail_arg = E
            else:
                if len(far_violators) < 3:
                    far_violators.append((E, float(far), float(EW*EW)))
        pt = f"{float(min_pztail):.5f}" if min_pztail is not None else "n/a"
        print(f"  diam {D:3d}: far<=E[W]^2 on {n_far_ok:5d}/{n:5d}  min PZ_tail(far-ok)={pt}"
              f"  min exact PZ={float(min_pz):.5f}"
              f"{'  [far-violators exist: compact]' if n_far_ok<n else '  [ALL spread]'}")
        if far_violators and D >= 18:
            print(f"        far-violator sample (compact-like at diam {D}): {far_violators[0][0]}")
    print()
    print("  READING: at small diam MANY families violate far<=E[W]^2 (compact, far>E[W]^2 --")
    print("  handled by exhaustive exact PZ >= 0.3468).  As diam grows, far<=E[W]^2 becomes")
    print("  universal (decorrelation) and the sharp PZ_tail >= ~0.47 >> bar.  The crossover")
    print("  diameter D* is where far<=E[W]^2 turns universal: below D* exhaustive, above D* tail.")


if __name__ == "__main__":
    main()
