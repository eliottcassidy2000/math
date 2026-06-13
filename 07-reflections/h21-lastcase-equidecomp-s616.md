---
source: claude-2026-06-03-S616
status: sharpens the one open H=21 case; equidecomposability type-semigroup framing; honest functor verdict
tags: [H-impossibility, H21, conflict-graph, equidecomposability, type-semigroup, equinumerosity, strong-component, infinite-go, honest]
---

# The last H=21 case, sharpened — and the H-spectrum as an equidecomposability semigroup

Pursuing the sharpest targets (is `{7,21}` the complete permanent-gap set? a
Go↔tournament functor?), with equinumerosity/equidecomposability as the lens.

## The one open case, made concrete

THM-079 reduced `H(T)=21` to a **single** open case: `Ω(T)` connected with
independence profile `(α₁=6, α₂=2, α₃=0)` — the conflict graph `K₆ minus 2 edges`,
`I = 1 + 6·2 + 2·4 = 21`. (Parts A disconnected, B `P₄`, C `α₃≥1` are proved.)

**New sharpening.** Among tournaments with *exactly 6 odd cycles*, the realized
number of vertex-disjoint cycle pairs `α₂` is:
```
n ≤ 7 : α₂ ∈ {0, 1}        (H = 13, 17)
n = 8 : α₂ ∈ {0, 1, 5}     (H = 13, 17, 33+)
```
`α₂ = 2` is **never realized** with `α₁ = 6`. So the `(6,2,0)` profile is
unrealized and `H=21` stays blocked. The open problem is now exactly:

> **Prove: no tournament has 6 odd cycles with exactly 2 vertex-disjoint pairs**
> (`α₁=6 ⟹ α₂≠2`).

With THM-079 A/B/C this would *complete* the proof that `H=21` is a permanent gap,
making `{7,21}` the complete set of permanent gaps below the monoid.

This is a genuine reduction: the question is no longer "is `K₆−2e` realizable as a
conflict graph?" (a graph-isomorphism search) but a clean cycle-combinatorics
statement about disjoint pairs — and the data shows `α₂` *jumps over* 2 (`{0,1,5}`),
which is itself a structural clue (a forcing phenomenon: 2 disjoint pairs among 6
cycles seems to force many more, landing at `α₂=5`, never at `2`).

## Equidecomposability: the H-spectrum is a type-semigroup

The equidecomposability lens makes the whole gap-structure crisp.

- **Atoms = strong (SC) tournaments.** `H` is multiplicative over the connected
  components of `Ω` (the strong decomposition): `H(T) = ∏_i H(C_i)`.
- **Equidecomposable tournaments** = same multiset of component `H`-values. You
  cut a tournament into its strong pieces and reassemble (any order, plus `n+2`
  source/sink padding, which preserves `H`). `H` is the **equidecomposition
  invariant** — the Tarski-style *measure* of the type.
- **Equinumerosity:** `H = |Hamiltonian paths|`. "`m` achievable" = some HamPath
  set has cardinality `m`. The achievable spectrum = the **type-semigroup** = the
  multiplicative monoid of strong (atomic) `H`-values.
- **Forbidden = non-realized measures.** `7` is *atomic-blocked* (prime,
  non-strong, no product). `21` is **doubly blocked**: the product route
  `21 = 3·7` needs the forbidden atom `7`, *and* the atomic route (a single
  strong component) needs the unrealized `(6,2,0)` profile. It fails both the
  equidecomposition route and the indecomposable route — which is *why* it is a
  gap.

So the "impossibility" is a statement in a Tarski-style type semigroup: the gaps
are the measures no equidecomposition type attains.

## The Go functor — honest verdict

I tried to make Go↔tournament a *functor* (`H ↦ game value`, strong-split ↦
disjunctive sum). It does **not** work as a value-preserving functor: `H`
**multiplies** over components, while combinatorial game values **add** over
disjunctive sums. So the rigorous correspondence is to **multiplicative
arithmetic** (atoms = strong tournaments = "primes", `H` = a multiplicative norm),
*not* CGT. The infinite-Go connection is a genuine **spectral analogy** — both
the H-spectrum and the Go ordinal-value spectrum are recursively built and have
arithmetic gaps — but it is an analogy of *spectra*, not a functor of *values*.
(Stating this honestly is the result: it tells us where to look — multiplicative
number theory of tournaments — and where not to.)

## Net

- **Is `{7,21}` complete?** Reduced to one clean cycle-combinatorics lemma
  (`α₁=6 ⟹ α₂≠2`), verified `n≤8`. Below the monoid, `{7,21}` are the only gaps;
  proving the lemma closes `21`.
- **Functor?** No (multiplicative vs additive); the right home is the
  equidecomposability type-semigroup / multiplicative arithmetic.
- **Equidecomposability** is the unifying frame: H-gaps = non-realized measures
  in the tournament type-semigroup, the same shape as LRC's forbidden `M` and
  Go's unreachable ordinals.

## Next

1. **Prove `α₁=6 ⟹ α₂≠2`** — a forcing lemma: 2 vertex-disjoint pairs among 6
   odd cycles force a 7th cycle (the data's jump `α₂: 1 → 5` suggests strong
   forcing). This would finish `H=21`.
2. **Strong-value spectrum.** Characterize which odd numbers are strong
   (indecomposable) H-values — the atoms of the type-semigroup. Gaps in the
   *strong* spectrum (just `{7,21}`?) generate all achievable-H gaps.
3. **Multiplicative arithmetic of tournaments.** Develop the Dirichlet-series /
   zeta analogue `∑ H(T) t^{...}` factoring over strong "primes" — the correct
   formalization of the partition-function/equidecomposability picture.
