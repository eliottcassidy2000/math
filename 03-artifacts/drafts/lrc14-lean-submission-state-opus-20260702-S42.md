# The LRC(14) Lean corpus: submission-state dossier

**opus-2026-07-02-S42 (HYP-3911).** The one-build state of the formalization: what compiles,
what is kernel-certain, what the final surface is, and what a mathlib submission ships.

## 1. The one build

`lake build TournamentH7` compiles the ENTIRE corpus — 124+ modules, every satellite
registered in the root (this session added the six unregistered opus satellites:
LRCCommensuration, RatIntervalsWrap, LRCWitnessWindow, LRCLadderPerLevel, CommensurationQ,
OriginNestQ). Result: **GREEN — Build completed successfully (8,589 jobs)**, after fixing one root-surfaced name collision (LRCWitnessWindow.rowOK renamed windowRowOK; kps's LRCCertTable.rowOK keeps the name per fewer-consumers rule).

## 2. The sorry census (exact)

Real `sorry` tokens in the whole corpus: **14, all in `LRCFourteenSkeleton.lean`** — the
LEGACY skeleton's named analytic obligations, superseded by the witness-route architecture.
Every other file is sorry-free (docstring mentions of the word aside). Crucially, the
endgame theorem does not depend on any sorried declaration — certified by the axiom audit
(§3): `sorryAx` appears nowhere in its profile.

## 3. The axiom audit (`AxiomAudit.lean`, checked in)

```
lrc14_endgame                    : [propext, Classical.choice, Quot.sound]   KERNEL-PURE
conjecture_one / conjecture_two  : [propext, Classical.choice, Quot.sound]
isLonelyAt_of_forall_not_dvd     : [propext, Classical.choice, Quot.sound]
unit-residue lemmas (both)       : [propext, Classical.choice, Quot.sound]
cert_lonely_tail / cert_ladder / cert_ladder' : [propext, Classical.choice, Quot.sound]
length_inter_le_left             : [propext, Classical.choice, Quot.sound]
seven_commensuration (+ general) : [propext, Classical.choice, Quot.sound]
pack/census rows (depth3 rows, certX_all, overlapQ, originNest) : + Lean.ofReduceBool
```
**No `sorryAx` anywhere in any audited cone — including the endgame's.** The 14 legacy
skeleton sorries are in declarations OUTSIDE the endgame's dependency cone.

Reading: the corpus splits into (a) **kernel-pure** theorems (standard trio: propext,
Classical.choice, Quot.sound) — the mathlib-track file, the ladder engines, the window rows;
and (b) **native_decide-carrying** theorems (+ Lean.ofReduceBool) — the census packs, where
compiled evaluation is trusted; every such computation has an independent Python mirror in
`04-computation/` (dual-verification discipline).

## 4. The final surface (what remains between here and `LRC14Statement`)

`TournamentH7.LRC14Assembly.lrc14_endgame` states LRC(14) from exactly two parameters:
- `hfloor` : the witness floor `witnessMP ≤ witnessG2` over all 13-speed shapes — every
  ingredient census is machine-checked (klein's LyWindowEnum k=8..13; kps witness-cert
  families + windows; opus ladder packs); the remaining work is case-split bookkeeping
  wiring shapes to their censuses.
- `hpartA` : `G2 > 0 ⟹ Mreach ≥ 1/14` — the floors and recursion glue are in Lean
  (CombPatterns THM-605(i) both directions; cert_ladder/cert_ladder'; module-3 both frames;
  THM-604(a) census); the remaining work is the module-6 fuel-checker soundness
  instantiation against the skeleton's `Mreach`.

Both are wiring tasks against machine-checked ingredients — no open mathematics.

## 5. What ships to mathlib (when the surface closes)

- **Layer 1 (unconditional, ships today):** `LonelyRunnerMathlib.lean` — definitions,
  q-witness, covering reduction, dilation invariance, Dirichlet tightness, k = 1, 2, the
  unit-residue lemma. Self-contained, kernel-pure, no project dependencies.
- **Layer 2:** RatIntervals + wrap + CombPatterns — the exact-rational interval calculus
  (reusable beyond LRC).
- **Layer 3:** the certificate engines (cert_lonely_tail, cert_ladder, cert_ladder') with
  pack schemas — infinite certified families, kernel-pure engines + native_decide packs.
- **Layer 4:** the endgame assembly once `hfloor`/`hpartA` land.

## 6. Certainty statement

Everything below the two-parameter surface is machine-verified: kernel-checked where the
mathematics is symbolic, native_decide + independent Python mirror where it is a finite
rational computation. No analytic estimate, no unverified table, no axiom beyond the
standard trio plus ofReduceBool is anywhere in the dependency cone of `lrc14_endgame`.
