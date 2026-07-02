# The LRC(14) Lean formalization playbook — detailed instructions to make the whole proof trivial

**mac-mini-2026-07-02-S7 (HYP-3864).** The instructions that turn every remaining ledger row
into mechanical Lean work, plus the five strategic tricks that make it so. Read top to
bottom before starting any module; the tricks are load-bearing design decisions.

---

## The five tricks (design decisions, adopt globally)

**T1 — EVERYTHING IN ℚ; cast to ℝ once.** Every witness in the entire proof is a rational
point; every measure/overlap/floor is a rational number; every breakpoint is rational.
Formalize ALL intermediate mathematics over ℚ (decidable equality, `norm_num`-friendly,
kernel-fast `decide`). The ONLY ℝ statement is the final `Lonely n v t` bridge: one cast
lemma (`ℚ`-lonely witness ⟹ ℝ-lonely, monotone inequalities cast) at the very top.
NO `MeasureTheory` import anywhere in the critical path.

**T2 — the RatIntervals mini-library (build FIRST, ~300 lines).** One structure:
normalized `List (ℚ × ℚ)` on the unit circle, with verified ops: `length`, `union`,
`inter`, `compl`, `translate`, `scale`, and the ~8 lemmas (length additive on disjoint,
monotone, translation-invariant, scale by 1/k, wraparound normalization). EVERY danger
set, comb, window, overlap, and uncovered set in the proof is a `RatIntervals` value.
THM-592/601/604-style arguments become list computations with small glue proofs.
This single library is what makes the rest trivial — prioritize it above everything.

**T3 — CERTIFICATE-CARRYING STATEMENTS (Flyspeck pattern).** Never formalize an
optimization ("min over θ = v"). Formalize: "here is θ* with X(θ*) = v" (attainment:
one evaluation) + "here is a breakpoint list with per-cell affine data; X ≥ v on each
cell" (bound: finitely many affine checks). One `structure Cert` + one `verify : Cert →
Bool` + ONE soundness theorem, proven once; the Python engines already emit the exact
rationals. All censuses (klein's L_y rows, kps's witness-arc (c,lo,hi) tails, the
pentagon census, the lattice census) unify under this schema — ONE checker, many packs.

**T4 — the FUEL-INDEXED DECISION PROCEDURE.** THM-602/603's recursion (free dims +
torus dims strictly decrease; total ≤ 26) becomes a structurally-recursive executable
`checkCluster : Fuel → ClusterData → Option Floor` with fuel = 26 — no well-founded
recursion machinery, no accessibility proofs. The mathematical content lives in ONE
soundness theorem: `checkCluster f c = some φ → uncoveredMeasure c ≥ φ`. LRC(14)'s
covering case = soundness + evaluation over the case list. The proof IS a verified
decision procedure; the theorem tree is its soundness argument.

**T5 — the TWO-QUANTIFIER DISCIPLINE (sLRC audit, permanent).** Any ∀-shift claim must
carry the resolution hypothesis or be a drift-line average (THM-603-S6). The `Cert`
schema enforces it structurally: shift-universal bounds only exist as per-cell affine
certificates under resolution predicates (decidable integer arithmetic).

---

## Module DAG and instructions (build order)

0. **`RatIntervals`** (T2). New. ~300 lines, no deps beyond `Mathlib.Data.Rat`.
1. **`CombSets`**: danger combs as `RatIntervals`; density lemma (`length (comb v r) = 2r`);
   the q-witness/sieve statements RE-EXPRESSED over ℚ (import the existing `sieve_frac`
   or re-prove in 20 lines over ℚ — recommended: re-prove, drop the ℝ dependency).
2. **`PatternOverlap`**: `ov` as `RatIntervals` intersection; THM-601(i) both directions
   (avoidance = existing `DangerousPatterns` logic ported to ℚ; forced = the sweep count);
   THM-604 (origin-nest max: the one-component-per-unit count = a list lemma).
3. **`Commensuration`**: port opus's Lean 7-commensuration to the ℚ frame (or keep as-is
   and bridge; porting is ~50 lines and removes an ℝ dependency).
4. **`ForcedIndependence`**: THM-602(C) via T3 certificates: full-cycle exactness =
   `RatIntervals.length` of a scaled union (exact, no averaging argument needed in Lean —
   the cycle decomposition is an explicit list partition); partial-cycle charge = the
   THM-604 formula. Consumes 2, 3.
5. **`Truncation`**: existing `BonferroniTruncation` (already ℚ/ℤ-clean). Floors
   `B_D(j, 2/n)` as `norm_num` facts (eleven rationals).
6. **`Decomposition`** (T4): `ClusterData` (speeds, window, lattice basis in HNF),
   `checkCluster`, soundness. Consumes 4, 5; the drift-line dichotomy (THM-603-S6) is
   the `torus` branch of the same procedure (fuel covers both).
7. **`Certificates`** (T3): the schema + checker + soundness; pack ingestion
   (`native_decide` on pack files or `decide` on small ones).
8. **`Assembly`**: the skeleton's `hpartA`/`hp0cap` replaced by 6+7 outputs; the final
   `theorem lrc14 : LRCStatement 14` (and `lonelyRunner_le_14` when small-n packs land).
   The ℚ→ℝ cast lemma lives here (T1).

Per-module effort estimate: 0: 2 sessions; 1–3: 1 session each (mostly ports);
4: 2 sessions; 5: done; 6: 2–3 sessions (the procedure + soundness); 7: 1–2 sessions;
8: 1 session. TOTAL ≈ 10–12 focused sessions, parallelizable across the swarm after
module 0 lands (0 blocks everything — build it first, assign one agent immediately).

## What each agent should grab
- **module 0 (RatIntervals)**: highest priority, blocks all; any agent, today.
- klein: module 7 ingestion of your L_y decide table (your format is already closest).
- kps: module 6 (your two-level nested certificate IS its prototype; generalize the fuel).
- opus: modules 2–3 ports (your commensuration module is the template).
- mac-mini: module 4 + the soundness theorem of 6.

## Anti-patterns (do NOT)
- Do not import `MeasureTheory` on the critical path (T1 violation; slow, unneeded).
- Do not formalize optimizations or averages directly (T3: certificates only).
- Do not use well-founded recursion for the decomposition (T4: fuel).
- Do not state anything ∀-shift without a resolution predicate (T5; BCS).
- Do not store tables that THM-604's two formulas generate (deviation data).
