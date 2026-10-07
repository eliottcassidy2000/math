#!/usr/bin/env python3
"""Audit A: debt/D ratios of the limit chains (HYP-9240 'in-phase debts within 0.75D-1.1D'), landing-time ratio s0/D,
and the maximal odd-step excess (odd steps - Terras steps/2 until {1,2}) of small N."""
from a2_barrier import limit_chain_scaled, collatz_terras_to_12
recs = [limit_chain_scaled(D) for D in range(3, 3001)]
inph = [(r['k']/r['D'], r['D']) for r in recs if r['kind'] == 'cycle' and r['inphase']]
outp = [(r['k']/r['D'], r['D']) for r in recs if r['kind'] == 'cycle' and not r['inphase']]
for lo in (3, 100, 163, 500):
    a = [x for x in inph if x[1] >= lo]
    print(f"in-phase D >= {lo}: k/D in [{min(a)[0]:.3f}, {max(a)[0]:.3f}]  argmin D={min(a)[1]}, argmax D={max(a)[1]}")
print(f"out-of-phase k/D range: [{min(outp)[0]:.3f}, {max(outp)[0]:.3f}]")
print("violations of 0.75D-1.1D (in-phase):", [D for x, D in inph if not (0.75 <= x <= 1.1)])
print(f"landing time s0/D for D >= 100: [{min(r['s0']/r['D'] for r in recs if r['D'] >= 100):.3f}, "
      f"{max(r['s0']/r['D'] for r in recs if r['D'] >= 100):.3f}]")
for NM in (880, 2527):
    best = max((o - s/2, N) for N in range(1, NM + 1) for (s, o) in [collatz_terras_to_12(N)])
    print(f"max odd-step excess over N <= {NM}: {best[0]} at N = {best[1]}")
