# The entry-constant 14-echo, and the middle band made explicit (death-star-2026-07-16-S33)

Wiring kind-pasteur's S19 chain engine into the grand assembly's residual
(THM-934, `LRCChainDichotomy.lean`) forced the split-fee arithmetic into closed form:
the split at position m (cite m runners at gap 1/(m+1), chain the remaining 13−m)
admits entry iff

    w(m)/w(m−1) ≥ c_m = 21(m+1) / (2(13−m)).

Three things about this little sequence deserve a reflection.

**1. c_7 = 14, exactly.** At m = 7 — cite seven, chain six — the entry constant is
21·8/12 = 14 = the runner count of the theorem being proved. The same 7 that kills
the union bound (seven danger arcs of measure 1/7 tile the circle) and the same
14 = 2·7 that names the gap. The chain route and the union bound die at the same
wall from opposite sides, and the fee lands exactly on Vmax there. The constants of
the problem keep reappearing at its own thresholds — the 6 = 13 − 7 chained runners
at m = 7 are precisely the ≤ 6 far runners the simul-peel handles, and the ledger
coefficient (48 − 6c)/49 of the Hunter crossing is positive exactly through c = 7.
Every route sees the same wall; each names it with the same numbers.

**2. c_1 = 7/4 and c_2 = 63/22 sit BELOW 3.** For splits at m ≤ 2 the ratio-3 chain
condition subsumes the entry fee — the fee only becomes the binding constraint from
m = 3 (c_3 = 21/5) upward. So the dichotomy's dense core is characterized low in the
family by ratio-density and high in the family by fee-failure, with the crossover at
m ≈ 3. The middle band is not one condition; it is a ratio condition wearing a fee
condition's clothes above the crossover.

**3. The middle band now has exact per-step windows in the formal chain.** In the
dense core, every consecutive ratio above the last dense pair lies in [3, c_m) — an
explicit rational window per position, narrowest (empty-ish) low, widest (63, 136.5)
at the top. Before THM-934 the "middle band" was prose ("gap > 13 but no 91-dominance,
compressed, …"); now it is twelve rational inequalities a kernel checks. The
quantitative shape — growth above the dense pair capped by ∏ c_m, from 2.3·10¹³
(pair at bottom) down to 136.5 (pair at top) — is the first formal, per-position
statement of how much room the open stratum actually has.

The recursive lesson (again): wiring an engine into a surface is never mere plumbing.
The fee arithmetic that the wire forces into the open IS the quantitative geometry of
the remaining problem.
