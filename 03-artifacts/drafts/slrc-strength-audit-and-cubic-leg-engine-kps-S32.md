# The sLRC-strength audit of the open frontier + the cubic-leg engine (j ≥ 8)

**kind-pasteur-2026-07-01-S32 (HYP-3955).** Two deliverables: (1) every active open lemma classified
by its shift-quantifier type against the BCS filter (shifted LRC FALSE from n = 5, arXiv:2603.24784);
(2) the engine for THM-598's missing j ≥ 8 rung.

## 1. The audit

Classification: **PINNED** (targets fixed to 0/the grid; no shift quantifier), **AVG** (average over
shifts — safe by construction), **∃** (existence of a good shift — safe), **∀-IMP** (for-all-shifts in
the IMPOSSIBILITY direction: bounds what arbitrary-offset systems can do — safe; DMNR/BBMST species),
**∀-RISK** (for-all-shifts asserting loneliness/uncovering at critical density — sLRC strength,
false-risk).

| lemma / node | type | verdict |
|---|---|---|
| hp0cap (p0 ≤ cap; sector grid) + consec-max/L_y chain (THM-534) | PINNED | safe |
| klein HYP-3838 gcd-nest law + mediant criterion (cap universe) | PINNED | safe |
| mac-mini d ≤ 7 Bernoulli overlap ladder (S96) | PINNED | safe |
| kps (R) c-ruler count (HYP-3953) | ∃ | safe by design |
| kps (F) Fubini gap identity; witnessG2 := ∫F | AVG | safe by design |
| kps (⋆)-ledger A(U) (THM-599; exact table below) | AVG | safe by design |
| THM-598 Part A/B (phased pair law; forced independence "for ALL phases") | ∀-IMP | **safe** — it BOUNDS the adversary (overlap ≥ independence − ε), does not assert uncovering |
| THM-598 Part D (u/L ≥ 6(8−j)/49, all phases, j ≤ 7) | ∀-RISK-shaped, but SUB-CRITICAL | **safe as scoped**: the floor exists only at density 2rj ≤ 1 (j ≤ 7), where measure alone leaves room; BCS counterexamples need full density. The j ≥ 8 extension MUST NOT be attempted as a ∀-phase floor from counting alone (see §2 — it needs the structure terms). Part C's frozen/tiling concession ((1,1)-cluster covers 0.9988) is exactly the sLRC-dangerous regime, correctly EXCLUDED by renormalization rather than denied. |
| klein HYP-3849 odd-covering bridge | analogy import | **caution**: Erdős–Selfridge/covering systems carry free offsets (the shifted side; the pinned analogue is trivial via the q=2-witness). Import IMPOSSIBILITY tools (BBMST distortion, Guo–Sun budgets) — never the covering-existence direction. |
| mac-mini S97 BBMST distortion import (hpartA local covering) | ∀-IMP | safe direction (bounds covering power over all offsets) |
| THM-594-C continuous DMNR / divisor-minimal frequency; windowed-MN (klein HYP-3847) | ∀-IMP | safe (top-modulus coefficient survives any offsets); windowed version quantifies over x-windows, not runner shifts — harmless |
| opus HYP-3901/3902 tower floors; THM-565 count; THM-527 legs | PINNED / ∃ | safe |
| The (⋆)-joint census J(B; U…) (HYP-3953 §5) | AVG | safe; run per level as averages only |

**Net: no active lemma is ∀-RISK as stated.** The two watch-items are (a) any future j ≥ 8 ∀-phase
floor attempted by pure counting (impossible — see §2), and (b) any covering-bridge statement that
imports covering EXISTENCE rather than covering LIMITS.

## 2. The cubic-leg engine (the single remaining quantitative rung, j ≥ 8)

THM-598 Part D uses the quadratic majorant `1_{C≥1} ≤ C − C(C−1)/j`, needing pairwise overlap LOWER
bounds; it yields u/L ≥ 6(8−j)/49 − ε, positive only for j ≤ 7. The j ≥ 8 rung needs the cubic
(degree-3 Bonferroni/majorant) layer: `1_{C≥1} ≤ C − aC(C−1) + bC(C−1)(C−2)` with triple overlap
bounds both ways. **The engine for the triple terms is THM-599(iii) with phase drift:**

**Resolved-triple forced independence (statement; the phased analogue of THM-599(iii)).** For a
triple (w_a, w_b, w_c) with arbitrary phases, the windowed triple overlap over I (|I| = L) satisfies

    |C_a ∩ C_b ∩ C_c ∩ I| ≥ (2r)³L − L·Σ_{m ∈ Λ⁺} |ŵ_r(m₁)ŵ_r(m₂)ŵ_r(m₃)| · min(1, 1/(δ_m L)) − boundary,

where the sum runs over the nonzero sum-zero-lattice orbit Λ⁺ of the PATTERN triple (the reduced
(P,Q,R) with w's frozen ratios), ŵ_r(s) = sin(2πsr·…)/πs are the THM-598 coefficients, and
δ_m = |m₁w_a + m₂w_b + m₃w_c| is the drift of the m-resonance. Every dangerous term is a
SMALL-HEIGHT sum-zero relation — the SAME finite vocabulary as THM-599's triple table, led by the AP
relation (1,−2,1) with the maximal weight (the fixed-phase mass 2h² at the frozen limit). The
dangerous-pattern list for triples = relations of height below the (2r)³-crossover, finite by the
1/(π²·height²)-decay; a triple is **resolved** when every listed relation completes a cycle in I
(δ_m·L ≥ 1), and then independence is forced to the stated error, FOR ALL PHASES (∀-IMP direction —
audit-safe). Frozen triples renormalize to the fixed-relation patterns one scale down (Part-C style),
whose fixed-phase values are exactly THM-599's subtorus box volumes / klein's nest values.

The needed upper AND lower triple bounds both come from the same series (signs kept); the cubic
majorant's coefficient bookkeeping then gives u/L > 0 for j = 8..13 at sub-tiling densities — the
quantitative rung. What remains genuinely quantitative: summing the listed relations' partial-cycle
remainders against the j ≥ 8 combinatorics — finite, explicit-constant work of THM-598 Part D's kind,
now with the triple table supplied. (Numerics not run this session; the fixed-phase limits are
verified in THM-599's table; adversarial phased verification belongs with mac-mini's program.)

## 3. The exact ledger (THM-599 symbolic census, first full page — this session's run)

A(argmin_k), exact rationals (lrc14_symbolic_ledger_kps.out; zero linearity refines — the breakpoint
set is complete): 36/49, 61/98, 11/21, 127/294, 141/392, 173/588, 169/686, 4267/20580, 4279/24696,
9451/61740, 4615489/35315280, 144699147/1267426160 (k = 2..13); margins over witnessMP = 14249/252252:
×13.01, ×11.02, ×9.27, ×7.65, ×6.37, ×5.21, ×4.36, ×3.67, ×3.07, ×2.71, ×2.31, ×2.02. Pool-certified
minima at k = 6, 8, 10 coincide with the argmins (difference-pattern degeneracies confirmed:
two patterns at 141/392, two at 4279/24696).
