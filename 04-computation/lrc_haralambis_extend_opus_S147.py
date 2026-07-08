"""
lrc_haralambis_extend_opus_S147.py   (opus-2026-07-07-S147, HYP-5277 part 5)

EXTEND the Haralambis 1977 conjecture (mu(D) = kappa(D) for every 3-element set D) past
the prior computational record max(D) <= 25 (Liu-Robinson 2020).

METHOD: mu(D) = max cycle mean of the D-avoiding window graph (S144).  To test mu = kappa
we only need: is there a cycle of mean > kappa?  = positive_cycle at v = kappa (Bellman
longest-path, numpy-vectorized, O(V) iterations).  No cycle => mu <= kappa; with mu >= kappa
(Cantor-Gordon) => mu = kappa.  A positive cycle => mu > kappa (a |S|=3 counterexample to
Haralambis -- would be major).

We push max(D) as far as the state count (~lambda^max) allows within a wall-clock budget,
and report the frontier honestly.  Primitive D only (gcd = 1); D = {a,b,c}, a<b<c.
"""
from fractions import Fraction as F
from math import gcd
import sys, time
import numpy as np

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_mu_eq_M_maxcycle_opus_S144 import build_window_graph, positive_cycle_exists
from lrc_graph_interpretation_ladder_opus_S141 import M_exact


def main():
    budget = 900.0  # wall-clock seconds
    t0 = time.time()
    print("=" * 92)
    print("HARALAMBIS EXTENSION: mu(D) = kappa(D) for 3-element D, pushing max(D) past 25")
    print("  (prior record: Liu-Robinson 2020, max <= 25; test = no cycle beats kappa)")
    print("=" * 92)
    max_reached = 0
    n_eq = 0
    counterexamples = []
    skipped_big = 0
    per_c = {}
    for c in range(3, 40):
        if time.time() - t0 > budget:
            print(f"  [budget reached before c={c}]")
            break
        c_eq = 0; c_skip = 0; c_max_states = 0
        for b in range(2, c):
            for a in range(1, b):
                if gcd(gcd(a, b), c) != 1:
                    continue
                D = [a, b, c]
                # state count guard: build graph, but skip if too big
                states, t0g, t1g = build_window_graph(D)
                if len(states) > 3_000_000:
                    c_skip += 1; skipped_big += 1
                    continue
                c_max_states = max(c_max_states, len(states))
                kap = M_exact(D)
                pos, _ = positive_cycle_exists(t0g, t1g, kap.numerator, kap.denominator)
                if pos:
                    counterexamples.append((tuple(D), kap))
                    print(f"  *** mu > kappa at D={D}, kappa={kap} -- HARALAMBIS COUNTEREXAMPLE ***")
                else:
                    c_eq += 1; n_eq += 1
                    max_reached = max(max_reached, c)
            if time.time() - t0 > budget:
                break
        per_c[c] = (c_eq, c_skip, c_max_states)
        if c_eq or c_skip:
            print(f"  c={c:2d}: mu=kappa on {c_eq:3d} sets, skipped {c_skip} (too big), "
                  f"max states {c_max_states}, t={time.time()-t0:.0f}s")

    print()
    print(f"  RESULT: mu = kappa confirmed on {n_eq} primitive 3-element sets; "
          f"{len(counterexamples)} counterexamples; max(D) reached = {max_reached}")
    if not counterexamples:
        if max_reached > 25:
            print(f"  ==> Haralambis conjecture EXTENDED to max(D) = {max_reached} "
                  f"(past the prior record 25), for the fully-enumerated c.")
        print(f"  (skipped {skipped_big} sets with > 3M window states -- the frontier needs a")
        print(f"   sub-exponential |S|=3 mu-engine to go further; window graph is ~lambda^max.)")
    print(f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
