# HYP-7173 — The trapped-D_T singles rung + the j ≥ 10 corner + formalization sweep (death-star-2026-07-17-S36)

**Status:** CLAIMED — in progress (S36). Verify the build report in the session log.

**Owner directive:** "prove the trapped-D_T singles rung and close the j>=10 corner,
also close any other remaining LRC 14 proof formalization completion tasks."

**Claim A (`LRCDeviationSingles.lean`):** the first unconditional rung of THM-940's
deviation ledger: for gcd(v i, q) = 1 the singleton joint-failure count is EXACT via
the unit bijection — N_{i} = q − 1 − #band(q); at 14 | q every singleton deviation is
the CONSTANT −13/7; unconditionally |D_{i}| ≤ 3 for every q ≥ 14 — hence the
singles' contribution to the deviation debt is ≤ 13·3 = 39, o(q): **the deviation
debt is carried entirely by |T| ≥ 2** for q beyond ~350·equilibrium⁻¹.

**Claim B (extend `LRCBlockSplitLift.lean`):** close QuadDenseCore's deferred
`10 ≤ j` disjunct with the generic engine at ε = 0 and EMPTY tail: j = 10 → the top
triple {w(10), w(11), w(12)}; j = 11 → the top pair {w(11), w(12)} (the fat-mass fee
at j = 11 is CHEAPER than THM-938's S20 pair entry: ~25·w(10) vs 78·w(10)). New
normal form `QuadDenseCoreClosed`: every disjunct an explicit fee failure — no
deferred corners.

**Claim C (sweep):** `S36AxiomAudit.lean` — one audit module over the S33–S36 ladder
(937/938/939/940/941 + the new rungs) with #print axioms; plus jointFail
monotonicity (N_T ≤ N_S for S ⊆ T).
