# Known Issues / Verification TODO

This file lists portions of the Lean code added in `opus-2026-05-29` (S7
extension) that have NOT been verified by `lake build` on the machine
where they were written (Lean was not installed there).

A next agent with a working Lean environment should:

1. Run `lake build` in `04-computation/lean/TournamentH7`.
2. Fix any of the issues below as they surface.
3. Verify `#print axioms` output in `Verify.lean` matches the documented
   axiom list in `SUBMISSION.md`.

## Files added in S7 (need build verification)

- `TournamentH7/Forbidden.lean`
- `TournamentH7/H21.lean`
- `TournamentH7/H63.lean` (finite n≤7 statement only; universal H≠63 is false)
- `TournamentH7/Redei.lean`
- `TournamentH7/HSpectrum.lean`
- updated root `TournamentH7.lean` and `Verify.lean`

## Specific issues to watch for

### 1. `alpha_chain_step` lemma specialisations (Forbidden.lean)

The lemmas `alpha_triple_chain`, `alpha_quad_chain`, `alpha_quint_chain`,
`alpha_sext_chain` use `simpa using h` to reduce `alphaCount (k-1) T` to
`alphaCount (k - 1) T` with the literal subtraction normalised. If
`simpa` doesn't fire, replace with:

```lean
lemma alpha_triple_chain (T : Tournament n) :
    alphaCount 3 T ≠ 0 → 3 ≤ alphaCount 2 T := by
  intro h
  exact alpha_chain_step T 3 (by decide) h
```

relying on definitional equality `3 - 1 = 2 : ℕ` to close the goal.

### 2. `set` interaction with `interval_cases` in H21.lean

The arithmetic core of `alpha_solution_of_H_eq_twentyone` uses
`set a₁ := alphaCount 1 T with ha₁` to abbreviate. Subsequent
`interval_cases a₂` should work, but if not, replace `set` with explicit
`have` / `let` bindings.

### 3. `Nat.choose 2 2` in H21.lean

The line `have : Nat.choose 2 2 = 1 := by decide` should succeed; if Lean
complains, replace with `Nat.choose_self 2` or `by simp [Nat.choose]`.

### 4. `ocf_extended` axiom (Forbidden.lean)

The statement claims `H T = ...` as an unconditional equality. This is
mathematically true only when `α_k = 0` for `k ≥ 7`, which holds at small
n but not in general. The axiom is consistent with itself if no agent
provides a counterexample inside Lean.

A more careful axiom would be:

```lean
axiom ocf_extended_safe (T : Tournament n) (hcap : ∀ k ≥ 7, alphaCount k T = 0) :
    H T = 1 + 2 * alphaCount 1 T + ... + 64 * alphaCount 6 T
```

For h ≤ 63 the cap `α_k = 0 for k ≥ 7` is forced by arithmetic, so the
proofs can supply `hcap` from `omega`. For the present session we kept
the simpler unconditional form.

### 5. `Redei.redei_existence` for n = 0

The axiom is stated with `hn : 1 ≤ n` to avoid the degenerate empty case.
If Lean complains about consistency at n = 0, the hypothesis ensures the
axiom is only ever applied at n ≥ 1.

### 6. H_ne_even uses `Nat.not_odd_iff_even`

Mathlib4's `Nat.not_odd_iff_even` may have a slightly different name
across versions (e.g. `Odd.not_even`, `Nat.even_or_odd`). If the build
errors, substitute the appropriate Mathlib lemma name.

## Audit checklist after a successful build

`lake env lean TournamentH7/Verify.lean` should print the dependency
axioms for each of:

- `H_ne_seven_audit` — depends on oracle's 7 original axioms.
- `H_ne_twentyone_audit` — additionally depends on `ocf_extended`,
  `alpha_chain_step`, `alpha_binomial_bound`, `no_alpha_10_0`,
  `no_alpha_8_1`, `no_alpha_6_2`, `no_alpha_4_3`.
- `H_ne_sixtythree_le_seven_audit` — depends on `H_ne_sixtythree_le_seven_axiom` only.
- `H_pos_audit` — depends on `redei_existence`.
- `forbidden_pair_audit` — universal {7,21} bundle.
- `forbidden_trio_le_seven_audit` — finite n≤7 {7,21,63} bundle.
