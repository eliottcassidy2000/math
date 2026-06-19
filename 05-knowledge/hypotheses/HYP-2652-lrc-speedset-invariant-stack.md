# HYP-2652: LRC speed-set invariant stack and near-tight tail laws

Status: OPEN / structural addendum.

Source: codex-2026-06-19-S33.

This complements HYP-2650 and the incoming HYP-2651 core-gap atlas.  HYP-2650 says the determinant is an addressed
wall/crossing sheaf, not a scalar.  HYP-2652 tests the same principle on exact
13-speed rows at the LRC14 threshold and isolates a proof-facing stack:

```text
CRT obligation
  -> additive anti-coset / relation density
  -> safe-component boundary-owner geometry
  -> binding denominator readout
```

## Evidence

Script:

```text
04-computation/lrc_invariant_atlas_codex_20260619.py
```

Output:

```text
05-knowledge/results/lrc_invariant_atlas_codex_20260619.out
```

Reflection:

```text
07-reflections/lrc-invariant-stack-codex-20260619.md
```

The atlas profiles 640 primitive 13-speed rows: AP/Goddyn-Wong/tail families,
S3 cluster rows, covering-biased random rows, and generic primitive rows.  It
computes exact `M(S)`, exact level-`1/14` safe components, residue/q-covering
profiles, additive relation counts, height-1 equal-subset anti-coset counts,
safe-component boundary-owner data, owner-pair tournament fingerprints, and
binding-denominator readouts.

Best raw correlations with exact `M(S)`:

```text
corr(M, three_term_count)        = -0.919
corr(M, safe_measure)            = +0.902
corr(M, longest_consecutive_run) = -0.896
corr(M, difference_energy)       = -0.890
corr(M, four_sum_collisions)     = -0.889
corr(M, support6plus_subsets)    = -0.722
```

Interpretation: the hard rows are relation-dense.  This is consistent with
HYP-2606/HYP-2646 and with the anti-coset direction of HYP-2612, with the
post-S13 caveat that CASE-thm538 disputes the full zero-padded `K` support-six
floor.  The `support6plus_subsets` invariant here is an active equal-subset /
additive-relation proxy, not a claim that short zero-padded Fourier relations
vanish.  Relation density creates the obstruction, but the address/sign/boundary
data decides whether the obstruction survives valuation.

## Tail Laws

The atlas found two exact near-tight laws in the tested bank:

```text
S = {1,2,...,12,13a}       has M(S) = a/(13a+1)   for tested a=2..9.
S = {1,2,...,11,13,12a}    has M(S) = a/(12a+5)   for tested a>=3.
```

The Goddyn-Wong row `{1,2,...,11,13,24}` is the exceptional tight seed at
`a=2`, with `M=1/14`.

These are binding-denominator invariants.  In the first family the denominator
is `13a+1`; in the second it is `12a+5` after the tight seed.  The formulas
should be proved directly from THM-524's active crossing language.

## Tournament Analysis

The session tried safe-component boundary owner pairs as tournament vertices.
The owner-tournament shadow compresses the endpoint geometry but does not
determine `M` by itself; its fibers still mix rows substantially.  This
supports the recurring lesson: tournament fingerprints are useful quotients,
but the faithful invariant must retain boundary addresses and component
lengths before scalarization.

## Proposed Proof Use

1. Treat CRT/q-coverage as a gate, not a determinant.
2. Use additive anti-coset counts to identify AP/Goddyn-Wong/tail-like hard
   normal forms.
3. Use safe-component boundary owners as the geometric carrier of the witness
   reservoir.
4. Use binding denominators as readouts to be forced by the first three layers,
   not as primary proof inputs.

The open claim is not that these invariants prove LRC14.  It is that they form
the correct resolution ladder for the gap-side speed-set structure, parallel to
HYP-2650's addressed wall/crossing sheaf on the broader invariant-separation
side.
