#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-470 lab, part 2: THE t=7 WALL (kind-pasteur-2026-06-11-S1, HYP-2396).

Trigger: the dyadic 1-jet algebra F2J is feature-UNSAT PER-SIZE at (3,7) — and
F2J refines F2, so the bi-dyadic algebra also dies per-size at t=7 (one size past
THM-465 A's frontier).  Gap-determined feature algebras are ANTITONE in t (a
witness at t+1 restricts to one at t: same classes, fewer demands), so each such
algebra's SAT region is an interval [3, t_dead).

Master question: the FINEST gap-determined algebra — the full translation-
invariant one, Finv(x,y) = the gap vector itself — at (3,7).
  * If feature-UNSAT: NO translation-invariant strong witness exists on [7]^3,
    subsuming every gap-determined rung (F2, F2J, F2X, leading-digit, ...) at all
    t >= 7; together with R(1,2)=3 and R(2,2)=5 (where invariant = free, THM-453 G)
    this is the first support for R(n,2) = 2n+1 — which would REFUTE HYP-2363
    (strong-witness frontier) and redirect the constructive-Specker program toward
    non-invariant or non-strong witnesses.
  * If SAT: the wall is algebra-specific; the invariant ladder continues above the
    jet rung and the table is a new witness one size past the known frontier.

Control first: Finv at (3,6) must be SAT (F2 SAT there already implies it).
"""

import time
from erdos592_bidyadic_rule_macmini_s3 import solve_featured, independent_verify


def Finv(x, y):
    return tuple(y[i] - x[i] for i in range(len(x)))


def main():
    t0 = time.time()
    print("=== control: full invariant algebra at (3,6) — must be SAT ===", flush=True)
    r, w = solve_featured(3, 6, Finv, "INV", tlimit=3600)
    if r:
        ok = independent_verify(3, 6, w[0])
        print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)

    print("\n=== THE WALL: full invariant algebra at (3,7) ===", flush=True)
    r7, w7 = solve_featured(3, 7, Finv, "INV", tlimit=7200)
    if r7:
        ok = independent_verify(3, 7, w7[0])
        print(f"   independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
        print("   -> wall is algebra-specific; invariant ladder lives at t=7", flush=True)
    elif r7 is False:
        print("   -> NO translation-invariant strong witness on [7]^3;", flush=True)
        print("      all gap-determined algebras dead for ALL t >= 7 (antitone);", flush=True)
        print("      R(n,2) = 2n+1 pattern: 3, 5, 7?  (HYP-2396)", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
