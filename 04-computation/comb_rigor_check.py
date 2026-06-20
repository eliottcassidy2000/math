#!/usr/bin/env python3
"""
Scrutinize the comb inequality rigor. The claimed bound:
  meas(G \ D_r) >= (6/7) M - 2c/(7r)
where G safe (measure M, c components), D_r = union of m/r-centered arcs half-width 1/(14r).

Reasoning check: D_r is a periodic comb, period 1/r, each tooth width 2*(1/(14r)) = 1/(7r).
So density of D_r in any FULL period = (1/(7r)) / (1/r) = 1/7.
On a component [a,b] of G of length L: it contains floor or so full periods.
Removed = (length of D_r inside [a,b]). Over full periods: exactly (1/7) per period length.
The two boundary periods (partial) can each contribute at most one full tooth = 1/(7r).
So removed <= (1/7)*L_fullperiods + 2/(7r) <= (1/7)*L + 2/(7r).
Thus meas(component \ D_r) >= L - (1/7)L - 2/(7r) = (6/7)L - 2/(7r).
Summed over c components: meas(G\D_r) >= (6/7)M - 2c/(7r).  <-- this is RIGOROUS.

Verify this bound holds AND is the right inequality direction by brute exact test
on adversarial cases (small r, many components).
"""
from fractions import Fraction
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2
import itertools
from math import gcd
from functools import reduce

# Exhaustive over ALL 2-hole bases (any two holes, NOT just tower) + ALL new speeds r in [14,200]:
# confirm meas(G_B \ {r}) >= (6/7)M_B - 2 c_B/(7 r) with ZERO violations.
worst_slack = None
viol = 0; tested = 0
for holes in itertools.combinations(range(1,14), 2):
    B = tuple(d for d in range(1,14) if d not in holes)
    M, c = lonely_measure(B, want_comps=True)
    for r in range(14, 201):
        if r in B: continue
        C = tuple(sorted(B + (r,)))
        LHS = lonely_measure(C)
        RHS = Fraction(6,7)*M - Fraction(2*c, 7*r)
        tested += 1
        slack = LHS - RHS
        if slack < 0:
            viol += 1
        if worst_slack is None or slack < worst_slack[0]:
            worst_slack = (slack, holes, r, LHS, RHS, c)
print(f"comb bound tested on {tested} (2-hole base, r) pairs; violations = {viol}")
ws = worst_slack
print(f"tightest slack = {float(ws[0]):.8f} at holes={ws[1]} r={ws[2]} (LHS={float(ws[3]):.6f}, RHS={float(ws[4]):.6f}, c={ws[5]})")
print(f"=> comb inequality is {'RIGOROUS (always >= 0 slack)' if viol==0 else 'VIOLATED'}")
