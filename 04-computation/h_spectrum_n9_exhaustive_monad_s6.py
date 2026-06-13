#!/usr/bin/env python3
"""h_spectrum_n9_exhaustive_monad_s6.py — EXHAUSTIVE n=9 H-spectrum (monad-compute-2026-06-04-S6).

Extends the exhaustive n=8 H-spectrum (monad-compute-2026-06-04-S1,
`h_spectrum_n8_exhaustive_monad.py`, complete census over 2^28 labeled tournaments,
gaps in [1,609] = {7,21}) UP ONE LEVEL to n=9.

n=9 has A000568(9) = 191,536 isomorphism classes. A labeled census (2^36 ~ 6.9e10)
is infeasible in Python, but H(T) is an isomorphism invariant, so enumerating the
191,536 ISO CLASSES yields the COMPLETE n=9 H-spectrum (every achieved H value).

Engine: reuse the validated isomorph-free orderly generator from
`h21_finite_check_v2_monad_s4.py` (color-refinement canonical form; validated
no-cap against A000568 through n=7, and here through n=9). We run it with an
effectively infinite c3 cap so NOTHING is pruned -> all iso classes are produced.

For each iso class we compute H(T) = number of directed Hamiltonian paths
(Redei count, always odd) via the standard 2^m bitmask DP (H_count from the engine).

OUTPUT:
  * per-level iso-class counts n=2..9, checked against A000568,
  * distinct H values at n=9, min/max,
  * full GAP set: odd integers in [1, max_H] NOT achieved by any n=9 tournament,
  * explicit confirmation of the {7,21} low-gap status and what unlocks at n=9
    that was absent at n=8 (e.g. did the high-end n=8 gaps {611,...} fill in?),
  * H histogram (value -> #iso classes) saved alongside.

Session: monad-compute-2026-06-04-S6.
"""
import sys, os, time
sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from h21_finite_check_v2_monad_s4 import (
    extend, beats_from_canon, H_count, c3_count, is_strong,
)

A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}
NOCAP = 10 ** 9  # effectively infinite c3 cap -> no pruning -> all iso classes

# n=8 exhaustive spectrum facts (from h_spectrum_n8_exhaustive_monad, OPEN-Q line ~1162):
N8_GAPS_LOW = {7, 21}                       # the only gaps in [1,609]
N8_MAXH = 661
N8_HIGH_GAPS = {611, 615, 617, 619, 623, 625, 635, 647, 655}  # high-end sparseness < 661


def generate_level(maxn):
    """Orderly-generate all iso classes, returning {n: set_of_canon_codes} only for n=maxn,
    printing per-level counts and checking against A000568."""
    print(f"  {'n':>3} {'classes':>10} {'expected':>10} {'mark':>9} {'cum_s':>9}")
    t0 = time.time()
    R = {1: {0}}
    for k in range(1, maxn):
        nxt = set()
        for code in R[k]:
            extend(beats_from_canon(code, k), k, NOCAP, nxt)
        R[k + 1] = nxt
        del R[k]
        n = k + 1
        exp = A000568.get(n)
        mark = "OK" if exp == len(nxt) else "MISMATCH"
        print(f"  {n:>3} {len(nxt):>10,} {str(exp):>10} {mark:>9} {time.time()-t0:>8.2f}s")
        if exp is not None and exp != len(nxt):
            print(f"  ABORT: iso-class count mismatch at n={n}")
            return None
    return R[maxn]


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    print("=" * 78)
    print(f"  EXHAUSTIVE H-SPECTRUM at n={maxn} via isomorph-free enumeration")
    print("=" * 78)
    print()
    classes = generate_level(maxn)
    if classes is None:
        sys.exit(1)

    print(f"\n  [H] computing H over all {len(classes):,} iso classes at n={maxn} ...")
    t0 = time.time()
    hist = {}        # H value -> #iso classes
    strong_hist = {} # H value -> #strong iso classes
    for code in classes:
        beats = beats_from_canon(code, maxn)
        h = H_count(beats, maxn)
        hist[h] = hist.get(h, 0) + 1
        if is_strong(beats, maxn):
            strong_hist[h] = strong_hist.get(h, 0) + 1
    print(f"      done in {time.time()-t0:.1f}s")

    achieved = set(hist)
    # all H are odd (Redei); sanity check
    all_odd = all(h % 2 == 1 for h in achieved)
    minH, maxH = min(achieved), max(achieved)
    n_distinct = len(achieved)

    # gap set: odd integers in [1, maxH] not achieved
    gaps = [h for h in range(1, maxH + 1, 2) if h not in achieved]
    low_gaps = [g for g in gaps if g <= 609]

    print("\n" + "=" * 78)
    print(f"  RESULTS — n={maxn} EXHAUSTIVE H-SPECTRUM")
    print("=" * 78)
    print(f"  iso classes            : {len(classes):,}  (A000568({maxn})={A000568.get(maxn)})")
    print(f"  all H odd              : {all_odd}")
    print(f"  distinct H values      : {n_distinct}")
    print(f"  H range                : [{minH}, {maxH}]")
    print(f"  # strong-achieved H vals: {len(strong_hist)}")
    print()
    print(f"  GAP SET (odd, not achieved, in [1,{maxH}]): {len(gaps)} values")
    print(f"    low gaps (<=609)     : {low_gaps}")
    high_gaps = [g for g in gaps if g > 609]
    print(f"    high gaps (>609)     : {len(high_gaps)} values"
          + (f" (first/last: {high_gaps[0]}..{high_gaps[-1]})" if high_gaps else ""))
    print()
    print("  --- comparison to exhaustive n=8 (h_spectrum_n8_exhaustive_monad) ---")
    print(f"    n=8 low gaps [1,609] : {sorted(N8_GAPS_LOW)}  (permanent {{7,21}})")
    print(f"    n=9 low gaps [1,609] : {sorted(low_gaps)}")
    if set(low_gaps) == N8_GAPS_LOW:
        print("    => {7,21} REMAIN the only low gaps at n=9 (consistent with permanent gap set).")
    else:
        newly_filled = N8_GAPS_LOW - set(low_gaps)
        newly_opened = set(low_gaps) - N8_GAPS_LOW
        print(f"    => CHANGED: newly filled {sorted(newly_filled)}, newly opened {sorted(newly_opened)}")
    print(f"    n=8 high gaps (>609) : {sorted(N8_HIGH_GAPS)} (max H=661)")
    filled_n8_high = sorted(g for g in N8_HIGH_GAPS if g in achieved)
    print(f"    of those, achieved at n=9: {filled_n8_high}")
    print()
    # what new H values appear at n=9 (top of range that n=8 could not reach)
    print(f"  n=8 max H = {N8_MAXH}; n=9 max H = {maxH}  (delta {maxH-N8_MAXH})")

    # save full histogram
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "05-knowledge", "results")
    histfile = os.path.join(outdir, f"h_spectrum_n{maxn}_histogram_monad_s6.tsv")
    with open(histfile, "w") as f:
        f.write("H\tcount_all\tcount_strong\n")
        for h in sorted(hist):
            f.write(f"{h}\t{hist[h]}\t{strong_hist.get(h,0)}\n")
    print(f"\n  histogram saved: {histfile}")

    # also save gap list
    gapfile = os.path.join(outdir, f"h_spectrum_n{maxn}_gaps_monad_s6.txt")
    with open(gapfile, "w") as f:
        f.write(f"# n={maxn} exhaustive H-spectrum gap set (odd, not achieved, in [1,{maxH}])\n")
        f.write(f"# iso classes={len(classes)}, distinct H={n_distinct}, range=[{minH},{maxH}]\n")
        f.write(f"low_gaps(<=609): {low_gaps}\n")
        f.write(f"high_gaps(>609): {high_gaps}\n")
    print(f"  gap list saved : {gapfile}")
    print("\n  DONE.")


if __name__ == "__main__":
    main()
