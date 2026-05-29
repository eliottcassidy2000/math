# Iso-Invariance as a Cleaner Axiom Base

**Session:** oracle-2026-05-29-S3
**Files:** `04-computation/lean/TournamentH7/TournamentH7/{Iso, IsoProperties}.lean`

## Observation

A natural way to organise tournament theorems is by what they actually
prove about *isomorphism classes* rather than about individual labelled
tournaments.  Every meaningful tournament invariant (Hamiltonian-path
count, regularity, SC-ness, …) is preserved under the natural notion
of tournament isomorphism:

> An iso `f : T₁ ≅ T₂` is an `Equiv.Perm (Fin n)` whose
> `arc_eq : ∀ i j, T₁.arc i j = T₂.arc (f.perm i) (f.perm j)`.

This session's Lean experience showed that recasting theorems through
this lens has surprising organisational benefits.

## Concrete observations

### (1) `outDegree` is iso-invariant — provable in Lean *without any axioms*.

The proof is just `Finset.card_bij` with the iso's underlying
permutation as the bijection.  No `ocf`, no `moonMoser`, no
`tilde_score_*` needed — only the definitions plus mathlib's
extensionality.

### (2) `IsRegular` is iso-invariant — provable, again axiom-free.

Direct corollary of (1): if every score equals (n−1)/2 in T₁, every
score equals (n−1)/2 in T₂ since the iso just permutes scores.

### (3) `IsSelfComplementary` has a clean iff:
> `IsSelfComplementary T  ↔  T ≅ op T`.

The `≅` formulation factors out the existential quantifier over
permutations and makes the structure visible.

## The pattern

Many theorems in the project canon are stated for "tilings" (labelled
tournaments with a base path), but the underlying mathematical content
is about iso classes.  When you formalise these in Lean, the canonical
iso framework lets you separate:

* the **labelling-dependent statement** (the one in the project canon,
  often about specific tilings or specific vertex roles); from
* the **labelling-invariant content** (an iso class fact, the actual
  theorem one would publish).

Concrete examples this session:

| Project statement | Iso-class version |
|-------------------|--------------------|
| Tile-complement T̃ score at vertex 0 ≠ T score at vertex 0 | The score *sequence* of T̃ differs from that of T at the base-path sink |
| Paley(p) regular ⟹ Paley(p) ∉ SF (one tiling rep) | regular iso class ⟹ iso class ∉ SF (any tiling rep) |
| THM-280: grid-symmetric *tiling* gives SC class | THM-280 lifts to iso classes via vertex reversal |

The "iso-class version" is usually the more mathematical statement;
the project canon's "tiling statement" is the *computational* version
used to verify or implement.

## Recommendation

The project's `IsSelfComplementary` and `IsSelfFlip` definitions
should standardise on the iso-class form.  This Lean session introduced:

```lean
def IsSelfComplementary (T : Tournament n) : Prop := T ≅ op T
def IsSelfFlip (T : Tournament n) : Prop :=
  ∃ σ : Equiv.Perm (Fin n), ∀ i j, T.arc i j = (tilde T).arc (σ i) (σ j)
```

Both factor out the existential cleanly.  Future canon definitions
should follow this pattern.

## A formalisation-driven canon principle (sequel)

In `formalization-driven-decompositions.md` we noted that the
formalisation revealed a sharper `alpha_descent` axiom.  This
session reveals an analogous pattern:

> **Iso-class statements are cleaner than tiling statements** —
> the Lean proof of an iso-class statement is often shorter, depends
> on fewer axioms, and exposes the underlying mathematical structure.

This is the "iso-class principle".  Combined with the
"arithmetic/structural decomposition" of forbidden-H, we have:

1. *Sharpen* the axioms (use the strongest natural form).
2. *Decompose* by content type (arithmetic vs structural).
3. *Lift* to iso classes wherever possible.

These three meta-principles emerged from the Lean formalisation
process — not from informal reasoning.  They are concrete,
checkable, and applicable to the rest of the project canon.

## Concrete TODO from this session

1. **Audit project canon definitions** for iso-class vs tiling
   distinctions.  Where canon defines a property via a specific
   labelling, add an iso-class form alongside.
2. **Prove H_iso_invariant** in Lean (sketched in `IsoProperties.lean`
   but not yet completed).
3. **Define `Tournament.iso_class T`** as a quotient type, formalising
   the iso-class graph G_n's vertex set in Lean.
