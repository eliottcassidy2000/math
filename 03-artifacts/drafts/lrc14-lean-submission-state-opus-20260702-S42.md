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

## 7. Addendum (S43): the completeness leg, confronted and pinned

**The owner's question — do the censuses provably EXHAUST the covering universe? — was
correct to ask: no exhaustion theorem existed anywhere in the corpus** (audited: kps's
CertRoute has the covering reduction and per-family theorems; nothing of the shape "every
covering tuple lands in a certified family").

**What is now true (LRC14CompletenessSurface.lean, sorry-free, kernel-pure):**
`hfloor_of_census_and_peel` proves that the exhaustion follows from two named legs by strong
induction on an abstract far-element counter:
- **census leg** (farCount = 0): ingredients machine-checked (klein LyWindowEnum k=8..13);
  remaining = format wiring, finite;
- **peel leg** (farCount > 0): klein's far-element rate lemma (HYP-4001(b)) + damped
  comparison (leg c, decide-shaped) + w-band sweeps (leg d, HYP-4005) — **the rate lemma is
  the ONE genuinely unformalized mathematical statement remaining in the whole chain**
  (proved on paper, python-verified exact; Lean transcription pending).

`lrc14_of_census_peel_partA` re-derives LRC(14) from (census, peel, hpartA) — the endgame
surface is now three finite-shaped parameters with named provenance, kernel-pure glue
throughout. The deeper measure-route exhaustiveness (THM-602 trichotomy, THM-595 case tree)
remains the paper architecture behind the same split; the witness route needs only the
counted form.

## 8. Addendum (S45): the unconditional distance, final form

The rate lemma is now UNCONDITIONALLY in Lean: opus's `RateLemma.lean` skeleton (S44) +
kps's `LRCFarElementRate.lean` (S11: `rate_core` two-boundary-cells trichotomy,
`toothClip_sum_near` = the two-sided per-interval form discharging `hpartial`,
`length_inter_comb_near_region` = the aggregated region form) — all kernel-pure.  The
module-6 fuel gate has its instantiation exemplar (klein's `LadderPackData.lean`).

**The complete remaining surface between the corpus and `theorem lrc14 : LRC14Statement`:**
1. `DispatchComplete W` — the executable shape-dispatch covers every covering class at the
   chosen cut `W` (one decidable evaluation over the HNF census once the row packs are all
   ingested);
2. `hwindow` — bounded-window loneliness (`∀ v, |v i| ≤ W → ∃ t, Lonely`): the finite
   census at the cut, native_decide over the shape quotient.

Both are finite computations in the Dispatch conventions.  No mathematical statement —
analytic, structural, or arithmetic — remains unformalized anywhere in the chain.
