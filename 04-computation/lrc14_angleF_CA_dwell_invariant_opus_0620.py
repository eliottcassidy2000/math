#!/usr/bin/env python3
"""ANGLE F -- Cellular-automaton conserved/monotone quantity for measS7.

MODEL.  Offset set E (0 in E, |E|=k).  Color of clock e at slow time x in [0,1):
    color(e,x) = floor(7*frac(e*x))  in  Z/7.
As x increases, clock e's color advances at "speed" e: a +1 step (mod 7) happens
at every x = j/(7e), j integer.  At each such breakpoint EXACTLY ONE clock jumps
to the next cell -- a CA update (events almost surely distinct for primitive E).
Occupancy o(x) in {0,1}^7: o_c=1 iff some clock has color c.  N(x)=#occupied.
    measS7(E) = meas{ x : o(x) = all-ones } = p0.
e=0 pins color 0 forever (sector-0 anchor); so the inner cover {1..6} is what
varies.  We work the FULL Z/7 (all-7) cover to match measS7 exactly.

ANGLE F asks: is there a CA invariant (Lyapunov / winding / charge) that
EXPLAINS or BOUNDS why consec = {0,1,...,k-1} has the largest all-ones dwell?

This script:
  (1) Builds the exact event-driven CA trace (Fraction breakpoints).
  (2) Computes measS7 exactly = total dwell in all-ones state; cross-checks
      vs the independent p0_exact engine.
  (3) Tests CA invariants:
      (a) total winding W(x)=sum_e floor(7*e*x) -- a strictly monotone "clock"
          (the CA's intrinsic time).  Is the all-ones set an interval of W?
      (b) dwell-time STATISTICS of the all-ones state: number of all-ones
          intervals (return count), longest dwell, mean dwell.
      (c) the "deficit charge" D(x) = sum over missing colors -- does consec
          minimize integrated time-at-deficit?
      (d) per-event analysis: which events ENTER and which LEAVE all-ones,
          and a candidate conserved quantity (enter-count = leave-count, a
          topological/winding identity).
  (4) consec vs spread comparison at k=8,9,10 on dwell statistics.

stdlib only; exact Fractions.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations


# ---------------------------------------------------------------------------
# CA event-driven trace
# ---------------------------------------------------------------------------

def breakpoints(E):
    """All x in [0,1] where some clock changes color (or endpoints)."""
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for j in range(0, 7 * e + 1):
            bps.add(F(j, 7 * e))
    return sorted(bps)


def color_of(e, x):
    """floor(7*frac(e*x)) in 0..6."""
    y = (F(e) * x) % 1
    return (y.numerator * 7) // y.denominator


def occ_vector(E, x):
    occ = [0] * 7
    for e in E:
        occ[color_of(e, x)] = 1
    return tuple(occ)


def ca_trace(E):
    """List of atoms (lo, hi, occ_tuple). occ over full Z/7."""
    bps = breakpoints(E)
    atoms = []
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        atoms.append((lo, hi, occ_vector(E, mid)))
    return atoms


# ---------------------------------------------------------------------------
# measS7 from the trace, plus dwell statistics
# ---------------------------------------------------------------------------

def measS7_and_stats(E):
    atoms = ca_trace(E)
    all_ones = (1,) * 7
    total = F(0)
    # dwell runs: consecutive all-ones atoms merge into one interval
    runs = []  # list of (length, n_atoms)
    cur_len = F(0)
    cur_atoms = 0
    in_run = False
    for lo, hi, occ in atoms:
        w = hi - lo
        if occ == all_ones:
            total += w
            cur_len += w
            cur_atoms += 1
            in_run = True
        else:
            if in_run:
                runs.append((cur_len, cur_atoms))
            cur_len = F(0)
            cur_atoms = 0
            in_run = False
    if in_run:
        runs.append((cur_len, cur_atoms))
    n_runs = len(runs)
    longest = max((r[0] for r in runs), default=F(0))
    mean = total / n_runs if n_runs else F(0)
    return {
        "measS7": total,
        "n_runs": n_runs,
        "longest": longest,
        "mean_dwell": mean,
        "n_atoms": len(atoms),
    }


# ---------------------------------------------------------------------------
# CA invariants
# ---------------------------------------------------------------------------

def winding_total(E, x):
    """W(x) = sum_e floor(7*e*x) -- intrinsic CA clock, strictly increasing."""
    return sum((F(7 * e) * x).__floor__() for e in E)


def is_all_ones_an_interval_in_W(E):
    """Is the all-ones set {x : o(x)=1^7} an INTERVAL in the winding clock W?

    i.e. does the all-ones state, ordered by W, occupy a contiguous block?
    Returns (is_interval, n_runs_in_x, n_runs_in_W).
    Since W is strictly increasing in x (for primitive E), x-order = W-order,
    so this is just whether there's a single x-run.  We report it to confirm.
    """
    stats = measS7_and_stats(E)
    return stats["n_runs"] == 1, stats["n_runs"]


def deficit_charge_integral(E):
    """Integrated deficit charge: int_0^1 (7 - N(x)) dx = sum over colors of
    time-uncovered.  By linearity = sum_c meas{color c missed}.  Lower = better.
    measS7 = meas{deficit 0}.  Returns (int_deficit, mean_N)."""
    atoms = ca_trace(E)
    total_def = F(0)
    for lo, hi, occ in atoms:
        N = sum(occ)
        total_def += (hi - lo) * (7 - N)
    return total_def, 7 - total_def  # mean_N = 7 - int_deficit


def enter_leave_events(E):
    """At each breakpoint, classify the transition o(x-) -> o(x+):
       ENTER  : becomes all-ones
       LEAVE  : was all-ones, no longer
       INTERNAL: neither.
    Conserved quantity to test: #ENTER == #LEAVE (topological winding of the
    all-ones state on the circle [0,1) with periodic identification)."""
    atoms = ca_trace(E)
    all_ones = (1,) * 7
    flags = [a[2] == all_ones for a in atoms]
    n = len(flags)
    enters = leaves = 0
    for i in range(n):
        prev = flags[(i - 1) % n]  # periodic
        cur = flags[i]
        if cur and not prev:
            enters += 1
        if prev and not cur:
            leaves += 1
    return enters, leaves


# ---------------------------------------------------------------------------
# independent cross-check engine (full Z/7 cover)
# ---------------------------------------------------------------------------

def p0_exact_full(E):
    """measS7 via independent breakpoint sweep, full Z/7 cover."""
    bps = breakpoints(E)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        seen = {color_of(e, mid) for e in E}
        if len(seen) == 7:
            total += hi - lo
    return total


# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------

# measS7(consec_k) for k=7..12 (prompt tabulation, k=7 first)
CONSEC_TAB = {
    7: F(31, 210), 8: F(481, 1470), 9: F(2447, 5880),
    10: F(8899, 17640), 11: F(3419, 5880), 12: F(121103, 194040),
}


def consec(k):
    return tuple(range(k))


def fmt(q):
    return f"{q} ({float(q):.6f})"


def main():
    print("=" * 78)
    print("ANGLE F -- CA dwell-time invariant for measS7 (full Z/7 cover)")
    print("=" * 78)

    # --- sanity: consec values match tabulated, engine cross-check ---
    print("\n[1] consec measS7 vs tabulated + cross-check + dwell stats")
    print(f"{'k':>3} {'measS7':>14} {'tab?':>5} {'xcheck?':>8} "
          f"{'n_runs':>7} {'longest':>12} {'mean_dwell':>12} {'atoms':>6}")
    for k in (8, 9, 10, 11, 12):
        E = consec(k)
        st = measS7_and_stats(E)
        cross = p0_exact_full(E)
        tab = CONSEC_TAB.get(k)
        ok_tab = "OK" if tab == st["measS7"] else "FAIL"
        ok_x = "OK" if cross == st["measS7"] else "FAIL"
        print(f"{k:>3} {fmt(st['measS7']):>14} {ok_tab:>5} {ok_x:>8} "
              f"{st['n_runs']:>7} {float(st['longest']):>12.6f} "
              f"{float(st['mean_dwell']):>12.6f} {st['n_atoms']:>6}")

    # --- invariant (a): is all-ones a single interval in W? ---
    print("\n[2] CA invariant (a): is the all-ones set a single interval?")
    print("    (W strictly increasing => x-order = winding-order)")
    for k in (8, 9, 10):
        E = consec(k)
        single, nr = is_all_ones_an_interval_in_W(E)
        print(f"  consec_{k}: single-interval={single}  (n_runs={nr})")

    # --- invariant (d): enter==leave conservation ---
    print("\n[3] CA invariant (d): #ENTER == #LEAVE (winding conservation)")
    for k in (8, 9, 10):
        E = consec(k)
        en, le = enter_leave_events(E)
        print(f"  consec_{k}: enters={en} leaves={le} balanced={en == le}")

    # --- invariant (c): deficit charge ---
    print("\n[4] deficit-charge integral int(7-N)dx and mean N (full Z/7)")
    print(f"{'k':>3} {'int_deficit':>16} {'mean_N':>12}")
    for k in (8, 9, 10):
        E = consec(k)
        d, mN = deficit_charge_integral(E)
        print(f"{k:>3} {fmt(d):>16} {float(mN):>12.6f}")


if __name__ == "__main__":
    main()
