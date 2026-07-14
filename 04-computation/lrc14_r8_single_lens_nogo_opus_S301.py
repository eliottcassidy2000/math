#!/usr/bin/env python3
"""The r=8 single-lens deck NO-GO certificate (opus-2026-07-14-S301, HYP-6840).

GOAL: exhibit a covering 13-speed family V = 7P u W, |W| = 8 (all 7∤w), with an
EXACT interval-arithmetic certificate that the open exception-danger set D(W)
CONTAINS the closed preimage T = 7^{-1}(Gbar_P) of the core's closed 1/14-safe set.
Consequence: for EVERY core-safe t0 and EVERY sheet k, some exception is strictly
unsafe -- the lens-7 deck pierce (THM-767/771 style) can NEVER fire on V. This is
a rigorous method boundary: at r = 8 the slack (sum of counts = 8c/7 > c) is
realizable as a TOTAL cover of the core-safe deck, so no pointwise count/event
argument at this lens decides V.

SIGNAL COMPLEMENT: the same V is then closed by other routes (exact M >= 1/14 via
the library; capped envelope / band protocol on the top peel) -- separating the
method's boundary (real) from any doubt about the conjecture (none).

Construction (v2 after an honest negative with an overlapping core): the SPREAD
core P = {2,5,9,11,13} has small safe measure (|Gbar| ~ 2/7, the right target
size) and carries q = 2,5,9,11,13,7,14 itself; one structured exception 360
carries q = 3,4,6,8,10,12; SEVEN free slots are chosen by greedy exact set-cover
of T over candidates up to height 1200, with randomized pair repair. The first
attempt (core {1,2,3,4,6}, |Gbar| = 4/7, six free slots) failed and is recorded
in the .out history as the honest negative that motivated the target-size fix.
"""
import sys, os, random
from fractions import Fraction as F
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc14_certificates import good_intervals, M_exact, is_covering, h_band_protocol, LAM

random.seed(3010)
C = 7

def closed_minus_open(target, open_arcs):
    """exact: portions of closed intervals `target` NOT covered by the OPEN arcs.
    Closed semantics throughout: an open arc (lo,hi) never covers its own
    endpoints, so those points survive as (possibly degenerate) closed pieces --
    exactly the strict-danger semantics the no-go needs. [a,b] \\ (lo,hi) =
    [a, lo] (if a <= lo < b) u [hi, b] (if a < hi <= b); a degenerate [x,x] is
    removed only when lo < x < hi."""
    rem = [(a, b) for a, b in target]
    for lo, hi in open_arcs:
        if hi <= lo: continue
        nxt = []
        for a, b in rem:
            if hi <= a or lo >= b:
                if a == b and lo < a < hi:
                    continue                  # isolated point swallowed
                nxt.append((a, b)); continue
            # now lo < b and hi > a: real interior overlap
            if lo >= a: nxt.append((a, lo))   # keeps the endpoint lo
            if hi <= b: nxt.append((hi, b))   # keeps the endpoint hi
        rem = nxt
    return rem

def danger_arcs(w, lam=LAM):
    """open danger arcs of runner w on [0,1): ((m-lam)/w, (m+lam)/w)."""
    arcs = []
    for m in range(w + 1):
        lo, hi = (m - lam) / w, (m + lam) / w
        arcs.append((max(lo, F(0)), min(hi, F(1))))
    return arcs

def subtract_runner(rem, w, lam=LAM):
    """fast exact: subtract runner w's OPEN danger arcs from closed pieces `rem`,
    enumerating per piece only the m-values whose arc can touch it."""
    out = []
    for a, b in rem:
        pieces = [(a, b)]
        m_lo = int((w * a - lam)) - 1
        m_hi = int((w * b + lam)) + 1
        for m in range(m_lo, m_hi + 1):
            lo, hi = (m - lam) / w, (m + lam) / w
            nxt = []
            for x, y in pieces:
                if hi <= x or lo >= y:
                    if x == y and lo < x < hi: continue
                    nxt.append((x, y)); continue
                if lo >= x: nxt.append((x, lo))
                if hi <= y: nxt.append((hi, y))
            pieces = nxt
            if not pieces: break
        out.extend(pieces)
    return out

def measure(iv):
    return sum(b - a for a, b in iv)

# --- the target: T = 7^{-1}(closed core-safe set) ---
P = [2, 5, 9, 11, 13]
GP = good_intervals(P)                       # closed 1/14-safe intervals of the core
T = []
for j in range(C):
    for a, b in GP:
        T.append(((a + j) / C, (b + j) / C))
T.sort()
print(f"core P = {P}: |Gbar_P| = {measure(GP)} = {float(measure(GP)):.5f}, "
      f"{len(GP)} components; target T = 7^-1(Gbar_P): {len(T)} closed intervals, "
      f"measure {float(measure(T)):.5f}")

# --- structured exceptions for covering ---
W = [360]                                    # 8*45: q=3,4,6,8,10,12
# q=10 check: 360 = 10*36 ok. distinct, 7∤w ok.
rem = T
for w in W: rem = subtract_runner(rem, w)
print(f"after structured W={W}: uncovered measure {float(measure(rem)):.5f} "
      f"in {len(rem)} pieces", flush=True)

# --- greedy exact set-cover for the remaining six slots ---
cands = [w for w in range(2, 1200) if w % 7 != 0 and w not in W
         and w not in [7 * p for p in P]]
NFREE = 8
while len(W) < NFREE and rem:
    best_w, best_gain, best_rem = None, F(0), None
    for w in cands:
        if w in W: continue
        r2 = subtract_runner(rem, w)
        gain = measure(rem) - measure(r2)
        # prefer strictly fewer leftover pieces on ties
        if gain > best_gain or (gain == best_gain and best_rem is not None
                                and r2 and len(r2) < len(best_rem)):
            best_w, best_gain, best_rem = w, gain, r2
    if best_w is None: break
    W.append(best_w); rem = best_rem
    print(f"  greedy pick w={best_w}: uncovered now {float(measure(rem)):.6f} "
          f"in {len(rem)} pieces", flush=True)

covered = (not rem)
if not covered and len(W) == NFREE:
    # local repair: try swapping the last picks
    print("  greedy incomplete -- attempting randomized repair on last two slots...")
    base6 = W[:NFREE - 2]
    rem6 = T
    for w in base6: rem6 = subtract_runner(rem6, w)
    for trial in range(4000):
        pair = random.sample(cands, 2)
        if any(w in base6 for w in pair): continue
        r2 = subtract_runner(subtract_runner(rem6, pair[0]), pair[1])
        if not r2:
            W = base6 + pair; covered = True
            print(f"  repair found: {pair}")
            break

print()
V = sorted([C * p for p in P] + W)
print(f"V = {V}")
print(f"distinct 13 speeds: {len(set(V)) == 13}; covering: {is_covering(V)}")
if covered:
    # FINAL EXACT CERTIFICATE: recheck from scratch
    final = T
    for w in W: final = subtract_runner(final, w)
    print(f"NO-GO CERTIFICATE: D(W) contains T exactly -> {not final} "
          f"(recheck uncovered = {float(measure(final)):.9f}, pieces = {len(final)})")
    print("=> the lens-7 deck pierce can NEVER fire on V: every core-safe t0 has")
    print("   every sheet strictly blocked by some exception. Method boundary at r=8.")
else:
    print("NO CERTIFICATE FOUND in this construction -- record as honest negative.")

print()
print("--- SIGNAL COMPLEMENT: V is still lonely, by other routes ---")
MV = M_exact(V)
print(f"M_exact(V) = {MV} = {float(MV):.5f} >= 1/14: {MV >= F(1, 14)}")
lay, det = h_band_protocol(V)
print(f"band protocol on V: layer {lay}, detail {det}")
# other lenses: does any other scale c' admit a THM-761/767 route?
for c2 in (2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 143, 360):
    Pc = [v // c2 for v in V if v % c2 == 0]
    r2 = 13 - len(Pc)
    if len(Pc) >= 1:
        print(f"  lens c={c2}: |P|={len(Pc)}, r={r2}" +
              ("  <- THM-761 regime (r<=6)" if r2 <= 6 else ""))
