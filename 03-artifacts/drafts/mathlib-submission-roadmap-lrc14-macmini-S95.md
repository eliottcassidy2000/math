# Mathlib submission roadmap — the LRC(14) formalization (mac-mini-2026-07-01-S95)

## What is now machine-checked (sorry-free, lake-verified, toolchain v4.30.0)

Project: `04-computation/lean/TournamentH7/` (mathlib dependency, all builds green).

1. **`TournamentH7.LonelyRunner`** (oracle-S18, pre-existing): the definitions
   (`Lonely`, observer tie-freeness) and the **q-witness sieve** `sieve_frac`
   (THM-523 core: a reduced fraction `a/q` with no speed divisible by `q` is lonely)
   + `counterexample_needs_all_divisors` (counterexamples are covering sets).
2. **`TournamentH7.LRCWitnessAttainment`** (mac-mini-S25, pre-existing): the clearance
   `min_i dist(v_i t, ℤ)` is continuous and periodic; sup attained; `reach ≥ 1/14 → ∃ lonely t`.
3. **`TournamentH7.LRCUnitResidue`** (NEW, S95): the **unit-residue lemma**
   (THM-593 Part A + the stability addendum): missing unit residue ⟹ the explicit
   shifted witness `a/q − 1/(q(V+1))` is lonely with margin `1/q + 1/(q(V+1))`.
   Both the raw inequality (`unit_residue_improvement`) and the `Lonely` packaging.
4. **`TournamentH7.PolygonMirskyNewman`** (NEW, S95): the **regular-polygon coloring
   theorem** (discrete Mirsky–Newman/Davenport–Rado): a k ≥ 2 coloring of a regular
   n-gon's vertices into regular-polygon classes has two congruent classes
   (`two_congruent_classes`), via the root-of-unity pole argument
   (`cls_eq_image`, `cls_charSum`, geometric-sum dichotomy).
   This is the discrete twin of THM-594(C) (continuous Mirsky–Newman).

## Suggested PR sequence (mathlib-general content first)

- **PR 1 (`Mathlib.Combinatorics.???.MirskyNewman`)**: the polygon theorem.
  Self-contained, classical, referee-friendly; check for prior art
  (covering-systems formalizations exist outside mathlib; DMNR itself appears absent).
  Generalization worth doing for the PR: state over `ZMod n` cosets and derive the
  `range n` version; add the ℤ-version (no exact cover of ℤ by APs with distinct moduli)
  as a corollary via periodization.
- **PR 2 (`Mathlib.NumberTheory.LonelyRunner.Basic`)**: `Lonely`, dilation invariance,
  `sieve_frac`, `counterexample_needs_all_divisors`, witness attainment.
  These are the standard reductions cited across the LRC literature
  (Bohman–Holzman–Kleitman; Tao's blog; Malikiosis–Santos–Schymura) — good
  foundational value independent of LRC(14).
- **PR 3**: the unit-residue lemma + tight-set corollaries (needs `M(S)` as a def:
  `iSup` of the clearance; the sup-attainment from module 2 makes the tight-set
  reading precise).
- Repo-specific (NOT for mathlib): the radius-derivative structure (THM-592) needs
  piecewise-linear/BV infrastructure — park until PRs 1–3 land.

## Style debts before submission
- rename to mathlib conventions (`Nat.` namespace hygiene, `theorem` vs `lemma`),
  docstring format, `#lint` pass, remove repo-specific commentary from docstrings,
  minimize imports (`Mathlib.Tactic` → targeted imports).
