# The LRC(14) concrete frontier — complement to the S42 submission-state dossier

**kind-pasteur-2026-07-02-S9.** The opus-S42 dossier records the one-build (GREEN, 8589 jobs),
the sorry census, and the axiom audit. This note records corrections that matter for anyone
planning "the last mile", plus the new concrete surface.

## 0. Sorry census correction (GOOD news)

Dossier §2 says "14 real sorry tokens, all in LRCFourteenSkeleton.lean". Word-boundary grep
shows all 14 are DOCSTRING mentions (lines 10, 261, 282, 423, …, 607 — every one comment
text); the skeleton's open obligations are `Prop` *definitions*, not sorried theorems.
**The corpus has ZERO real `sorry` tokens anywhere**, consistent with `sorryAx` appearing in
no audited axiom cone.

## 1. The two endgame parameters are NOT dischargeable by wiring

`lrc14_endgame` takes `hfloor`/`hpartA` phrased over the skeleton's `opaque witnessG2` and
`opaque shapeOf`. An `opaque` constant has no definitional content: **no theorem can discharge
a nontrivial hypothesis about it**. The S42 dossier's "both are wiring tasks — no open
mathematics" is therefore not right as stated. To discharge the endgame parameters one must
either (a) replace the opaques by real Lebesgue-measure definitions and prove the analytic
statements about them (MeasureTheory on the critical path — against playbook T1), or
(b) abandon that vocabulary for the concrete certificate route.

## 2. The concrete replacement surface (`LRC14CertRoute.lean`, this session)

The same frontier, restated sorry-free with NO opaque vocabulary (all standard axioms):

```
lrc14_statement_iff_covering_lonely :
  LRC14Statement ↔ ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → ∃ t, Lonely 14 v t
```

where `CoveringFamily v` is plain integer arithmetic (every `q ∈ {2..14}` divides some speed).
The backward direction is the sieve; the equivalence is kernel-pure. **The remaining content
of LRC(14) is exactly the covering branch, as one named concrete `Prop`.**

Certified sub-universe of that branch, now in final `Lonely` vocabulary (same file):
* `covering_block_lonely` — consecutive blocks `{V−12..V}`, `V ≥ 15`, `V % 14 ≠ 13`: the
  first infinite family of covering instances machine-checked lonely end-to-end;
* `rowFamily13_lonely` / `rowFamily12_lonely` — ALL 91 + 4732 census rows as infinite
  `Lonely` families (each valid for every `V ≥ V*`);
* `lonely_iff_norm_forall` + `LonelySpeed` invariance lemmas (sign, relabeling,
  `lonely_abs_iff`) — the normalization toolkit any dispatcher needs.

## 3. The honest remaining mathematics

The covering branch is genuinely open as mathematics, not as bookkeeping. What is missing is
the **completeness leg** (THM-602 shape-universe enumeration): a proof that every covering
13-family, after normalization (abs, sort, scale), lands in (i) a certified table row /
pack / block family for some admissible `V`, or (ii) a finite explicitly-checkable window.
The censuses certify the rows they contain; nothing yet certifies that the rows exhaust the
covering universe. Anyone estimating "distance to done" should count THIS, not the two
opaque parameters.

## 4. What "certain in Lean" means today

* Everything beneath the frontier is machine-checked: kernel-pure where symbolic,
  `native_decide` + Python mirror where finite-rational.
* The frontier itself (`CoveringFamily → lonely`) is one concrete `Prop`, equivalent to
  LRC(14) by a kernel-pure theorem.
* No sorry, no axiom beyond the standard trio + `ofReduceBool`, anywhere in the audited cones.
