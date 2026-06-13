#!/usr/bin/env python3
"""Finish the n=8 slice of the C' exhaustive box certificate
(monad-compute-2026-06-03-S597).

The combined run lrc_Cprime_exhaustive_box_monad_s597.py exhausted n=4,5,6,7
(522,496 primitive multiple-of-8 ... multiple-of-7 configs, 0 tight) but timed
out partway through n=8. This driver completes the n=8 box (B=24, K=3) on its
own using the SAME verified routines, so the certificate covers n=4..8.

Reuses is_loose / enumerate_box from the main script unchanged (early-exit exact
positivity == open_safe_measure>0, validated 0 mismatches). A single tight
config would refute C'.
"""
import sys
sys.path.insert(0, "04-computation")
from lrc_Cprime_exhaustive_box_monad_s597 import enumerate_box, is_loose  # noqa: E402

N, B = 8, 24


def main():
    print("=" * 72)
    print("C' EXHAUSTIVE BOX CERTIFICATE -- n=8 completion (monad-compute-S597)")
    print(f"n={N}, B={B} (K={B // N}): every primitive (n-1)-subset of {{1..{B}}}")
    print("containing a multiple of 8, tested loose EXACTLY (positivity).")
    print("=" * 72)
    total = loose = tight = 0
    tights = []
    for V in enumerate_box(N, B):
        total += 1
        if is_loose(V, N):
            loose += 1
        else:
            tight += 1
            tights.append(V)
        if total % 50000 == 0:
            print(f"    ... scanned {total} ({loose} loose, {tight} tight)",
                  flush=True)
    print()
    print(f"n={N} B={B}: {total} configs, {loose} loose, {tight} tight  -> "
          f"{'PASS' if tight == 0 else 'FAIL'}")
    if tight:
        print("*** C' COUNTEREXAMPLE(S) ***")
        for t in tights[:50]:
            print("   TIGHT (M=1/8):", t)
    else:
        print("CERTIFICATE: every primitive multiple-of-8 speed set in the box is")
        print("LOOSE -- 0 exceptions. Combined with the main run, C' is now")
        print("EXHAUSTIVELY certified for n=4,5,6,7,8 within the listed boxes.")


if __name__ == "__main__":
    main()
