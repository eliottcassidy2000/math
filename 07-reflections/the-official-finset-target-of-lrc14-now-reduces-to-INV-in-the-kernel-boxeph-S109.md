# The official Finset bridges to `INVcov` and `ResidualINV`

> **SECOND CORRECTION (THM-1158 / MISTAKE-170):** literal `INVcov` is refuted
> by `2*{1,...,13}`.  Its Finset consumer is kernel-valid but has a false
> premise.  The `ResidualINV` bridge and equivalence remain exact; they are the
> honest formal interface.  Descriptions below of `INVcov` as a useful open
> supplier are retained only as historical record.

*boxeph-2026-07-18-S109; corrected codex-2026-07-18 (MISTAKE-166).*

The representation bridge is sound, but its premise must be stated precisely.  The
formalization now exports two related theorems:

```lean
LRC14_finset_of_INVcov :
  LRCUpTo13 → INVcov → LRC14.LRC14

LRC14_finset_of_residual_INV :
  LRCUpTo13 → ResidualINV → LRC14.LRC14
```

They serve different purposes.  `INVcov` is the useful, noncircular sufficient open
target.  `ResidualINV` is the exact counterexample target and records the shape of a
hypothetical obstruction, but—under the cited AP-core bridge—it is logically
equivalent to the working LRC(14) statement.  The residual theorem therefore closes a
formal representation gap; it does not turn LRC(14) into a smaller open problem.

## The two inverse interfaces

The stronger interface is

```lean
def INVcov : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v →
    (¬ ∃ t, Lonely 13 v t) →
    ∃ vstar, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar
```

Here `Covering v` says that every modulus `2..14` divides at least one speed.
`INVcov` is sufficient for LRC(14): on the no-`Lonely14` branch, the divisor sieve
forces Covering and band monotonicity forces no `Lonely13`; dominance then feeds
`ap_core_bridge`.

The exact residual interface is

```lean
def ResidualINV : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → Covering v →
    (¬ ∃ t, Lonely 14 v t) →
    ∃ vstar, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar
```

The file proves both directions needed to calibrate it:

```lean
residualINV_of_INVcov : INVcov → ResidualINV

residualINV_iff_LRC14 (cite : LRCUpTo13) :
  ResidualINV ↔
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t, Lonely 14 v t
```

Thus `INVcov` is the genuine sufficient supplier one may try to prove independently;
`ResidualINV` is a faithful language for counterexample elimination.  Neither is
currently identified with Tao's `n=12` inverse/AP-uniqueness theorem.  Such an
identification would require a separate structural bridge.

The old unqualified implication

```text
(¬ ∃ t, Lonely 13 v t)  →  13-fold dominance
```

is false.  The canonical family `{1,…,13}` has exact maximum `M = 1/14`, hence no
`Lonely 13` time, but has no 13-fold dominant speed.  Adding the covering premise is
not cosmetic: it is the information contributed by the no-`Lonely14` branch.

## The representation gap bridged

The working theorem has shape

```lean
∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t : ℝ, Lonely 14 v t
```

whereas the ledger's official target is

```lean
def LRC14.LRC14 : Prop :=
  ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
    ∃ t ∈ Icc (0 : ℝ) 1,
      ∀ w ∈ W, ∀ a : ℤ, 1 / 14 ≤ |(w : ℝ) * t - a|
```

Two elementary transfers close this mismatch.

1. **Enumerate the Finset.**  `W.equivFinOfCardEq hWcard : W ≃ Fin 13` gives an
   equivalence whose inverse defines `v i := ((e i : ℕ) : ℤ)`.  Its range is exactly
   `W`, and positivity transfers from the ledger hypothesis.

2. **Normalize the time.**  The proved lemma `lonely_fract` replaces a witness `t`
   by `Int.fract t ∈ [0,1)`.  Indeed,
   `vᵢ · fract(t) - m = vᵢ · t - (m + vᵢ⌊t⌋)`, and the parenthesized term is still an
   integer quantified by `Lonely`.

For `w ∈ W`, choose `i = e.symm ⟨w, hw⟩`.  Then `(v i : ℝ) = (w : ℝ)`, so the
pointwise Lonely inequality is exactly the ledger predicate.

## Kernel-checked formalization ladder

All the following rungs audit to Lean's standard axiom trio
`[propext, Classical.choice, Quot.sound]`:

| rung | theorem | precise content |
|---|---|---|
| S105 | `ap_core_bridge` | 13-fold dominance + LRC(≤13) implies `Lonely 14` |
| S106 | `sieve_dispatch` | not Covering implies `Lonely 14` |
| S107 | `density_far_extension` | a sufficiently separated far speed extends a lonely frame |
| S108 | `M_split` | no `Lonely14` forces Covering and no `Lonely13` |
| S108 | `LRC14_of_INVcov` | LRC(≤13) + `INVcov` implies working LRC(14) |
| S108 | `residualINV_iff_LRC14` | under LRC(≤13), the residual target is exact/equivalent |
| S109 | `LRC14_finset_of_INVcov` | the sufficient route lands on `LRC14.LRC14` |
| S109 | `LRC14_finset_of_residual_INV` | the exact residual route lands on `LRC14.LRC14` |

## Net

- **Delivered:** `LRCFinsetBridge.lean` exports the time-normalization lemma and both
  official-target bridges.
- **Useful open theorem:** proving `INVcov` would prove the official LRC(14) target,
  given the explicit `LRCUpTo13` citation.
- **Exact diagnostic language:** `ResidualINV` describes precisely what dominance
  would have to do on a counterexample, and the formal equivalence theorem prevents
  us from mistaking that restatement for progress on the open mathematics.
- **Remaining structural question:** connect a genuinely stronger theorem—possibly a
  covering inverse, compact-class classification, or tournament/quotient invariant—to
  this exact residual domain without assuming no `Lonely14` in the supplier itself.

Cross-links:
[[the-M-split-and-the-complete-kernel-checked-reduction-of-lrc14-to-INV-boxeph-S108]],
[[the-density-route-discharge-is-now-kernel-pure-lean-boxeph-S107]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
[[the-elementary-half-of-lrc14-is-now-kernel-pure-lean-the-ap-core-bridge-boxeph-S105]].
