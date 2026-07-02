#!/usr/bin/env python3
"""klein-2026-07-02-S97 -- HYP-4005(a): the W-BAND SWEEPS (rate lemma's finite leg).
For each row k, the binding damped shape E* (from HYP-4003/4004): explicit K(E*) =
sum_r |y_r| sum_{|A|=r} 2 comp(A,E*) r, threshold W0 = K/margin; exact sweep of
L_y(E* u {w}) for w in (spread, W0] confirming <= cap_{k+1}."""
from fractions import Fraction as F
import itertools
exec(open("04-computation/ly_far_element_rate_lemma_klein.py").read().split('print("="*92)')[0])

BINDING = {8: [0,1,2,3,4,5,6,7], 9: [0,1,2,3,4,5,6,7,8], 10: [0,1,2,3,4,5,6,7,8,9],
           11: [0,1,2,3,4,5,6,7,8,10,12]}
for k, E in BINDING.items():
    row = k + 1
    rmax = max(r for r,_ in DUALS[row])
    K = F(0)
    for r, y in DUALS[row]:
        K += abs(y) * sum(2 * len(avoid_set(E, list(A))) * r
                          for A in itertools.combinations(range(1,7), r))
    margin = CAPS[row] - Ly_damped(E, row)
    W0 = K / margin
    # sweep w: exact L_y(E u {w}) vs cap_row for w up to min(W0, cap at 400 for runtime)
    Wcap = min(int(W0) + 1, 400)
    over = [w for w in range(max(E)+1, Wcap+1) if Ly(E+[w], row) > CAPS[row]]
    print(f"row {row} (E*={E}): K = {K} = {float(K):.1f}; margin = {float(margin):.5f}; "
          f"W0 = {float(W0):.0f}; swept w <= {Wcap}: OVER-CAP = {len(over)} {over[:6]}")
print("NOTE: sweep truncated at 400 where W0 exceeds it (runtime); remainder = rate-lemma regime")
print("only if W0 <= 400, else the [400, W0] band remains queued (report per row above).")
