# THE CONSUMPTION-SEMANTICS REFEREE PAGE (ledger seam 2) — boxeph-S41

**Question.** The propagation ledger composes Error(w) ≤ 0.2729·diam/w against the S58
row margins. Do those margins guard exactly the functional the Error bounds?

**The per-row functionals (from the finish map + THM-663/727 as written):**
- k = 8: a Φ-CAP (Φ ≤ cap₉ = 1979/4004; deg-3 LP majorant + d ≤ 25 exhaustive box).
- k = 9: a J-FLOOR (J ≥ 432/91).
- k ≥ 10: INHERITED via the eigen-transfer THM-710 from the base rows — the two-scale
  error enters ONLY through the base rows.
- k = 11, 12, 13 tails: D3-functionals (LEM-009), margins +0.12/+0.157/+0.252.

**The consumption logic, verified structurally:**
1. THM-727 defines Error = Φ(E′∪w) − Φ_∞(E′) and proves |Error| = |S|/w with the
   Fourier reduction EXACT — so the ledger's bound perturbs Φ, the k = 8 row's own
   functional. Cap direction is harmless: Φ_∞ ≤ cap − margin and |Error| < margin give
   Φ ≤ cap. ✓ The k = 8 consumption is SOUND as composed.
2. The direction-agnostic template: [proved slack ≥ margin] + [|perturbation| < margin]
   closes a row regardless of cap/floor orientation — the ledger's uniform treatment is
   correct WHENEVER the perturbation of that row's functional obeys the same |S|/w law.
3. **THE k = 9 SUB-CHECK: CLOSED (S42).** THM-711 identifies the J-functional exactly:
   J = E[N(7−N)] = the mass of (hit, empty) SECTION PAIRS — a finite sum of 42 ordered
   pair box-hit measures. Each such measure is a Boolean section set with owned
   endpoints (THM-881 P1 is universal), so THM-727's endpoint reduction applies
   VERBATIM per pair: |J(E′∪w) − J_∞(E′)| ≤ Σ_pairs |S_pair|/w ≤ 42·(the uniform
   |U|-bound)/w. HONEST CONSTANT: the crude pair-count inflation gives
   W₀(9) ≤ 42 × 3.33·diam ≈ 140·diam for the J-lane (still finite, still band-coverable
   in principle); the named refinement is the shared-endpoint cancellation across the
   42 pairs (the pairs tile the same breakpoint set — the true constant is plausibly
   the original 3.33). THE LEDGER IS NOW FULLY CITED: every row's consumption semantics
   verified, k = 9 with the crude-vs-refined constant flagged.
4. The D3 tails (k = 11..13) are diameter-tail statements about the CLUSTER only (no
   far element in their hypotheses) — the far element enters them through the same base
   perturbation; no separate consumption issue. ✓

**Verdict.** The composition is sound; one citation-shaped sub-check (the J-row's
reduction) is isolated and assigned. Seam 2: CLOSED MODULO the named one-liner.

**Formalization entry-points (the owner-directed next phase; the decide-shaped batch):**
(1) THM-882 per-cell solves (rational linear algebra); (2) THM-878 D(q) case list;
(3) THM-884 exact-ℚ discs; (4) THM-885 sweep leaf certificates (j ≤ 5);
(5) THM-888(A) comb-diagonal closed forms + (MI)/(MI0) (3-line algebra);
(6) THM-892 (K)/(C\*)/(P) lemmas (second differences, Möbius, subgroup Parseval);
(7) THM-899's B₄ closed form; (8) this ledger's arithmetic (pure ℚ).
klein's sorry-free Kendall formula shows the pipeline is warm — these eight are the
LRC(14)-side batch, all kernel-pure candidates.
