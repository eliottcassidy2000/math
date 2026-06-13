#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: EXTEND the n*inf_S Gstar(S) asymptotic to n=9 and
push n=6,7,8 to higher speed-bound B, to decide HYP-2296 Part C:
   does n*inf -> 0 (decay ~1/n^2)  or  stay bounded away from 0 (Theta(1/n))?

monad-compute-2026-06-07.  Continues a34958b (matched-B n=5..8).

Key correction over the matched-B run: those inf values are B-TRUNCATED UPPER
BOUNDS.  The true inf drops as B grows (e.g. n=6: 3/20 @B=13 -> 6/41 @B=18).
So for the asymptotic we must read each n's inf at its highest feasible B and
watch whether n*inf has STABILIZED (B-converged) before comparing across n.

This driver runs a per-(n,B) plan, prints/saves one line as soon as each
(n,B) finishes (incremental: killable without losing earlier rows), and flags
B-convergence (inf unchanged from the previous B at that n).

Uses the verified exhaustive routine search() from
signed_lrc_inf_highB_monad_s3.py (THM-426 cut maximization, exact rational
maximin, set-level pruning).
"""
import importlib.util, sys, time
from fractions import Fraction as F

spec = importlib.util.spec_from_file_location(
    "highb", "04-computation/signed_lrc_inf_highB_monad_s3.py")
hb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hb)

# Per-n plan: ascending B so we can watch B-convergence of the inf.
# Ordered so cheaper rows land first (incremental safety).
PLAN = [
    (5, [16, 20, 24]),
    (6, [14, 16, 18, 20]),
    (7, [13, 14, 15, 16]),
    (8, [12, 13, 14]),
    (9, [11, 12, 13, 14]),
]

def main():
    print("=" * 90)
    print("SIGNED-LRC inf_S Gstar(S): n=9 extension + high-B convergence   monad-compute-2026-06-07")
    print("=" * 90)
    print("Question (HYP-2296 C): n*inf -> 0 (1/n^2) or bounded>0 (Theta(1/n))?")
    print(f"{'n':>3} {'B':>3} {'#sets':>7} {'inf Gstar':>12} {'float':>9} "
          f"{'n*inf':>8} {'1/n':>8} {'2/n^2':>8} {'conv?':>6} {'sec':>7}")
    print("-" * 90)
    sys.stdout.flush()
    for n, Bs in PLAN:
        prev = None
        for B in Bs:
            t = time.time()
            g, sets, ns = hb.search(n, B)
            dt = time.time() - t
            conv = "YES" if (prev is not None and g == prev) else ""
            mini = sets[0] if sets else None
            print(f"{n:>3} {B:>3} {ns:>7} {str(g):>12} {float(g):>9.5f} "
                  f"{float(n*g):>8.4f} {float(F(1,n)):>8.5f} {2.0/n**2:>8.5f} "
                  f"{conv:>6} {dt:>7.1f}   min={mini}")
            sys.stdout.flush()
            prev = g
        print("-" * 90)
        sys.stdout.flush()
    print("DONE.")

if __name__ == "__main__":
    main()
