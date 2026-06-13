#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: FAIR matched-B exhaustive comparison of inf_S Gstar(S)
across n, to read the n-asymptotic of  n * inf.   monad-compute-2026-06-06.

The prior table (signed_lrc_inf_highB) computed inf at DIFFERENT B per n
(n=6@B<=18, n=7@B<=13, n=8@B<=11), so the n*inf values are not directly
comparable -- larger n was cut off at smaller B.  Here we run the SAME proven
exhaustive search at a COMMON B for n=5..8, then also push each n to its
feasible ceiling, so the n-trend of n*inf is read on equal footing.

Imports the verified exhaustive routine from signed_lrc_inf_highB_monad_s3.py
(THM-426 cut maximization, exact rational maximin, set-level pruning).

Open question (THM-429/HYP-2296): does n*inf -> 0 (decay ~1/n^2) or stay
bounded away from 0 (true Theta(1/n) floor)?  Lower bound: Gstar >= 2/((n-1)(n-2))
~ 2/n^2 (THM-429 measure bound), so n*inf >= ~2/n -> 0 is ALLOWED; the question
is whether the TRUE inf realizes it.
"""
import importlib.util, sys, time
from fractions import Fraction as F

spec = importlib.util.spec_from_file_location(
    "highb", "04-computation/signed_lrc_inf_highB_monad_s3.py")
hb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hb)

def run(n, B):
    t = time.time()
    g, sets, ns = hb.search(n, B)
    dt = time.time() - t
    return g, sets, ns, dt

def main():
    print("=" * 80)
    print("SIGNED-LRC matched-B exhaustive inf_S Gstar(S)   monad-compute")
    print("=" * 80)

    # Phase 1: common B across n=5..8 (fair n-comparison)
    COMMON_B = 13
    print(f"\n--- Phase 1: COMMON B={COMMON_B} (fair n-comparison) ---")
    print(f"{'n':>3} {'#sets':>7} {'inf Gstar':>12} {'float':>9} {'n*inf':>8} "
          f"{'1/n':>8} {'2/n^2':>8} {'sec':>6}")
    print("-" * 80)
    phase1 = []
    for n in (5, 6, 7, 8):
        g, sets, ns, dt = run(n, COMMON_B)
        phase1.append((n, g, sets[0] if sets else None))
        print(f"{n:>3} {ns:>7} {str(g):>12} {float(g):>9.5f} {float(n*g):>8.4f} "
              f"{float(F(1,n)):>8.5f} {2.0/n**2:>8.5f} {dt:>6.1f}")
        sys.stdout.flush()

    print("\nPhase 1 minimizers:")
    for n, g, V in phase1:
        print(f"   n={n}: inf={g}  V={V}")

    # Phase 2: per-n feasible ceiling (best known upper bound on inf)
    print(f"\n--- Phase 2: per-n ceiling (best feasible B) ---")
    ceilings = [(5, 22), (6, 17), (7, 14), (8, 12)]
    print(f"{'n':>3} {'B':>3} {'#sets':>7} {'inf Gstar':>12} {'float':>9} "
          f"{'n*inf':>8} {'sec':>6}")
    print("-" * 80)
    for n, B in ceilings:
        g, sets, ns, dt = run(n, B)
        brk = "BREAK" if g < F(1, n) else ""
        print(f"{n:>3} {B:>3} {ns:>7} {str(g):>12} {float(g):>9.5f} "
              f"{float(n*g):>8.4f} {dt:>6.1f}  {brk}  {sets[0] if sets else None}")
        sys.stdout.flush()

    print("\nDONE.")

if __name__ == "__main__":
    main()
