#!/usr/bin/env python3
"""h21_window_closure_monad_s5.py — CLOSE the H=21 finite window.

Builds directly on the (bug-fixed, MISTAKE-054) v2 engine
`h21_finite_check_v2_monad_s4.py`, reusing its validated isomorph-free
c3<=CAP generator. New question answered here:

  Does the population of *strong* tournaments with c3 <= 10 die out at finite m,
  thereby closing the H=21 finite window WITHOUT relying on Busch's recurrence?

Reduction (HYP-2193, recapped):
  * H is multiplicative over strong components, and H=7 is a permanent gap
    (THM-029), so 21 = 3*7 is impossible as a product => any H(T)=21 needs a
    SINGLE strong component with H=21.
  * For a strong tournament, alpha_1 (all directed odd cycles) >= c3, and
    H = I(Omega,2) >= 1 + 2*alpha_1 >= 1 + 2*c3. So H=21 forces c3 <= 10.
  Hence: "no strong tournament has c3 <= 10 AND H=21" closes the gap.

The decisive structural fact this script verifies computationally:
  min number of 3-cycles in a strong tournament on n vertices == n-2 (Moon).
  => at n >= 13, every strong tournament has c3 >= 11 > 10, so the strong
     c3<=10 population is EMPTY and the window is closed at n=13.

Per level m we report:
    classes (c3<=10), #strong, min c3 over strong, min c3 over all c3<=10,
    min H over strong, and #(H==21).  Closure is declared at the first m with
    #strong == 0.

Each level is printed + flushed immediately so the critical m=13 result is
saved even if a later level is cut off by the timeout.

Session: monad-compute-2026-06-04-S5.
"""
import sys, os, time
sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from h21_finite_check_v2_monad_s4 import (
    extend, canon, beats_from_canon, is_strong, c3_count, H_count, validate,
)

CAP = 10


def run(maxm):
    print("=" * 78)
    print(f"  H=21 WINDOW CLOSURE: strong c3<={CAP} population vs m  (target m={maxm})")
    print("=" * 78)
    if not validate():
        print("ABORT: A000568 validation failed")
        return
    print(f"\n  [main] generate c3<={CAP} classes, scan strong population per level:\n")
    print(f"    {'m':>3} {'classes':>12} {'strong':>9} {'minc3_str':>10} "
          f"{'minc3_all':>10} {'minH_str':>9} {'H=21':>5}  {'gen/scan(s)':>12}")
    R = {1: {0}}
    closed_at = None
    minc3_strong_seq = []
    for k in range(1, maxm):
        cur = R[k]
        nxt = set()
        t0 = time.time()
        for code in cur:
            extend(beats_from_canon(code, k), k, CAP, nxt)
        R[k + 1] = nxt
        del R[k]
        gen_t = time.time() - t0
        m = k + 1
        t1 = time.time()
        strong = 0
        minc3_str = None
        minc3_all = None
        minH = None
        h21 = 0
        for ccode in nxt:
            beats = beats_from_canon(ccode, m)
            c3 = c3_count(beats, m)
            if minc3_all is None or c3 < minc3_all:
                minc3_all = c3
            if is_strong(beats, m):
                strong += 1
                if minc3_str is None or c3 < minc3_str:
                    minc3_str = c3
                H = H_count(beats, m)
                if minH is None or H < minH:
                    minH = H
                if H == 21:
                    h21 += 1
        scan_t = time.time() - t1
        minc3_strong_seq.append((m, minc3_str))
        print(f"    {m:>3} {len(nxt):>12,} {strong:>9,} "
              f"{str(minc3_str):>10} {str(minc3_all):>10} "
              f"{str(minH):>9} {h21:>5}  {gen_t:>6.1f}/{scan_t:<5.1f}")
        if h21:
            print(f"    *** H=21 FOUND at m={m} (c3<=10 strong) — investigate! ***")
        if strong == 0 and closed_at is None and m >= 4:
            closed_at = m
            print(f"    *** CLOSURE: no strong tournament with c3<={CAP} at m={m}. "
                  f"H=21 window CLOSED for all m>={m}. ***")
            # one more level for margin, then stop generating (sets shrink, but
            # if strong is 0 it stays 0 since c3 only grows under extension)
            break
    print("\n  min c3 over strong tournaments per m (expect == n-2, Moon):")
    print("    m      :", " ".join(f"{m:>3}" for m, _ in minc3_strong_seq))
    print("    minc3  :", " ".join(f"{str(c):>3}" for _, c in minc3_strong_seq))
    print("    n-2    :", " ".join(f"{m-2:>3}" for m, _ in minc3_strong_seq))
    matches = all((c is None) or (c == m - 2) for m, c in minc3_strong_seq if c is not None)
    print(f"    min c3(strong) == n-2 for all computed m with strong>0: {matches}")
    if closed_at:
        print(f"\n  RESULT: H=21 finite window CLOSED at m={closed_at} "
              f"(strong c3<={CAP} population empty).")
    else:
        print(f"\n  RESULT: window not yet closed by m={maxm}; strong c3<={CAP} "
              f"population still nonempty.")
    print("  DONE.")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxm", type=int, default=14)
    args = ap.parse_args()
    run(args.maxm)
