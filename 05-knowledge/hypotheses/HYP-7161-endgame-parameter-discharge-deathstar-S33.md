# HYP-7161 — The endgame parameter discharge: maximal unconditional composition (death-star-2026-07-16-S33)

**Status:** CLAIMED — in progress (S33)

**Claim (verify-first: read the build, not this stub, for the final state):**
The two parameters of `lrc14_grand_assembly` are `cite : LRCUpTo13` (the sanctioned
citation node — stays named per owner directive) and `hresidual : ResidualObligation`.
S33 works the second: compose the existing sorry-free reducers
(`residualObligation_of_ledger` / `_of_B5` / the dissociated chain
`LRCDissociatedAssembly` → `LRCDissociatedRigidity`, the detuned-dispatch reduce,
kind-pasteur's HYP-7145 stratification: generic bulk + opus THM-928(A) lacunary ratio ≥ 19
+ structured sliver) into the sharpest machine-checked surface, prove new glue where a
finite gap remains, and report the exact residual honestly.

**Known at claim time:**
- census leg CLOSED: `hdistinct20_from_data` (native_decide, 6084 rows) + `hwindowW_closed 22`
  kernel-pure (LEM-024 route).
- peel/far leg reduced: `CoveringFarLonely` ⟸ top-ratio-91 (cite) + dominant peel (cite)
  + spread13 + detuned + coarse + common-residue — remaining class = ResidualObligation
  (covering ∧ gap>13 ∧ compressed ∧ distinct |speeds| ∧ max ≥ 23 ∧ no detuned ∧ no coarse
  ∧ no common residue).
- Discrete reducers exist: HasLiveRuler (pair-sum ledger), B5 > 0 (discrete Bonferroni).

**Needs evidence:** whether the proved strata (lacunary ≥ 19, rigidity, block-gluing
THM-933) jointly cover the residual class, or what exact sub-class remains.

**S33 mid-session state (checkpoint):** frontier mapped — the fleet's minimal analytic
target is `SafeMeasureFloorPrimitive` (opus S207; "prove the floor, not the rigidity",
kps S127); the 7-wall crossing has Hunter (klein) + telescoping engine (mac-mini)
proved with the pair-floor per-cell obligation open; the paper program's remaining
items are codex's K6 formula + the manifest's ten Lean items (boxeph S49/S50).
NEW THIS SESSION (→ THM-937): `LRCChainDichotomy.lean` — the split-search dichotomy
wiring kps-S19's `cite_chain_lonely` into `ResidualObligation` (never wired before):
LRC14Statement ⟸ cite + DenseCoreObligation, where the dense core adds to the
residual the explicit certificate {dense pair ratio < 3 at j; ratios ≥ 3 above j;
entry fee 2(12−k)w(k+1) < 21(k+2)w(k) failing ∀ k ≥ j}. Referee: 200k-tuple
exhaustiveness PASS, 50k planted citable tails all closed
(`chain_dichotomy_referee_deathstar_S33.out`). Lean build in progress at checkpoint
time — VERIFY the build report in the session log before consuming.

## Collaborative machine-checked discharge — codex-S23

`TournamentH7.LRCEndgameParameterDischarge` now gives the sharpest composition
available from the proved detuned reducers.  It closes two pieces of the exceptional
dispatch without adding hypotheses:

- if `nonMultCard v g = 2` and every `g`-multiple is already a `2*g`-multiple, the
  existing two-detuned lift theorem supplies the lonely time;
- if `nonMultCard v g = 3` and every exceptional coordinate has reduced denominator
  at least `8`, the existing three-detuned coarse theorem supplies the lonely time.

Consequently `lrc14_from_deep_detuned_and_dissociated_B5` reduces the complete sharp
assembly to the sanctioned `cite : LRCUpTo13` plus exactly two named suppliers:

1. `DeepExceptionalDetunedDispatch`: the nonterminating two-adic tower and the
   three-detuned branch with a reduced denominator below `8`;
2. `DissociatedB5Supply`: positive depth-five Bonferroni mass on the primitive,
   covering, gap-family, compressed, distinct-absolute-speed, nonstructured,
   nondetuned residual.

No generic branch, primitive peel, witness-attainment step, or assembly implication
remains implicit.  The subsequent live pull adds two orthogonal refinements.  THM-937
closes every residual family with a citable ratio-`3` tail, leaving its explicit
`ChainDenseCore`.  THM-935 proves that `B5` only sees exact-support `2`--`5` relations
and reduces universal algorithmic exhaustion to the support-`3,4,5` lattice tails
`T_s(H)`.  Boxeph S53 also refutes the unqualified S52 diagonal-dominance conjecture:
cross terms are negative off resonance but positive on six resonant frames, where a
mode-coherence correction is required.  These three facts sharpen the second supplier,
not the now-closed detuned dispatch glue.

**Assumption challenge.**  The quotient vertices here are divisibility levels and
detuned coordinates, not runners, arcs, or fixed circle sections.  This quotient
preserves `nonMultCard`, reduced denominators, and immediate lift termination, while it
forgets phase chronology and circle-component incidence; those data must remain in the
actual tuple `v` and in the downstream lonely-time consumer.

## Dense-core composition and verified THM-937 repair — codex-S25

The clean build found that the THM-937 checkpoint had not completed elaboration; its
four local compatibility failures are now repaired, and `lonely_or_denseCore`,
`residualObligation_of_denseCore`, and `lrc14_of_denseCore` audit with foundational
axioms only.

`TournamentH7.LRCDenseCoreEndgame` then composes the two reductions at the correct
logical level.  Inside `ResidualObligationPrimitive`, it first branches on a `d=2,3`
detuning and consumes `MultiDetunedDispatch`; on the complementary dissociated branch,
it sorts absolute speeds and consumes `lonely_or_denseCore`.  A successful chain split
is transported directly back to the original tuple.  Only the other branch calls the
new `DenseCoreDissociatedB5Supply`.  Thus no invalid implication “lonely implies
positive B5” is used.

The current narrowest checked theorem is
`lrc14_from_deep_detuned_and_denseCore_dissociated_B5`.  Besides the sanctioned
`LRCUpTo13` citation, its exact suppliers are:

1. `DeepExceptionalDetunedDispatch` (nonterminating two-adic lift or a three-detuned
   small reduced denominator);
2. positive `B5` only on the primitive, dissociated residual whose sorted absolute
   speeds satisfy `ChainDenseCore`.

THM-935 identifies the likely mathematical content of supplier 2—universal
support-`3,4,5` relation-lattice tails after the proved support-2 tail—but that
relation-mass theorem is not yet the Lean supplier.
