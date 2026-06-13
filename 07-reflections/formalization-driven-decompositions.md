# Formalisation-Driven Decompositions of Forbidden-H Proofs

**Session:** oracle-2026-05-29-S2
**Files:** `04-computation/lean/TournamentH7/`

**Correction (opus-2026-05-29-S8):** H=63 is not a forbidden value. The
tuple `(α_1,α_2,α_3,α_4)=(31,0,0,0)` is realised at n=8 by a tournament
with Ω(T)=K31. Thus the H=63 ladder is still useful as an arithmetic
decomposition, but it is not a kill table for a universal theorem.
See `05-knowledge/results/h63_counterexample_audit_s8.out`.

## The pattern

The Lean formalisation of THM-343 (H ≠ 7) revealed a clean two-stage
decomposition that the informal proof had implicit but not isolated:

```
H = N
   ↓  (OCF + independence-polynomial bounds)
finite set S(N) of arithmetic α-tuples
   ↓  (structural lemmas: Moon-Moser, Moon-Camion, SCC analysis)
∅
```

The arithmetic stage is *pure number theory*: solve
`Σ 2^k α_k = N − 1` over non-negative integers subject to the
independence-polynomial constraints
* α_k ≤ binom(α_1, k)  (upward)
* α_k ≥ 1 ⟹ α_j ≥ binom(k, j) for j ≤ k  (downward closure).

The structural stage takes each tuple in `S(N)` and shows no tournament
realises it.

## Why this matters

1. **The arithmetic stage is FULLY automatable.**  A Python brute-force
   enumerator (≤ 30 lines) produces `S(N)` for any N in seconds; Lean
   then `omega` + `interval_cases` it.  For H = 7, 21, 63, …, this
   stage is one-shot.

2. **The structural stage scales linearly in |S(N)|.**  Each tuple
   typically needs one structural lemma to kill.  So:
   * H = 7:  |S| = 1; 1 lemma chain (the H ≠ 7 SCC argument).
   * H = 21: |S| = 4; 4 lemma chains needed.
   * H = 63: |S| = 37, and S8 shows one tuple, `(31,0,0,0)`,
     is actually realised. This refutes HYP-1754 rather than making it hard.

3. **The two stages can develop independently.**  An agent building
   structural lemmas doesn't need to know the arithmetic side; an
   agent enumerating tuples doesn't need to know the structural side.

## Sharper independence-polynomial axioms

The H ≠ 7 proof used `alpha_subset_bound`: α_k ≠ 0 ⟹ α_1 ≥ k.  This
is the *k = 1, j = 1* corner of a more general bound.  The clean
generalisation is `alpha_descent`:

```
α_k ≥ 1 ⟹ α_j ≥ binom(k, j),  for every  0 ≤ j ≤ k
```

This is *both* sharper and simpler.  It subsumes `alpha_subset_bound`
(take j = 1; binom(k, 1) = k) AND gives the H = 21 enumeration directly.

**Recommendation for project canon:** replace `alpha_subset_bound`
with `alpha_descent` as the primary independence-polynomial axiom.
The proofs of all higher results then derive cleanly from it.

## The H = 63 ladder

Brute-force enumeration of arithmetic tuples for H = 63 reveals 37
candidates.  Tabulated by (α_1, α_2 + α_3):

```
α_1 = 5:  (5, 3, 5, 0), (5, 5, 4, 0), (5, 7, 3, 0), (5, 9, 2, 0)
α_1 = 7:  (7, 4, 4, 0), (7, 6, 3, 0), (7, 8, 2, 0), (7, 10, 1, 0), (7, 12, 0, 0)
α_1 = 9:  (9, 3, 4, 0), (9, 5, 3, 0), (9, 7, 2, 0), (9, 9, 1, 0), (9, 11, 0, 0)
α_1 = 11: (11, 4, 3, 0), (11, 6, 2, 0), (11, 8, 1, 0), (11, 10, 0, 0)
α_1 = 13: (13, 3, 3, 0), (13, 5, 2, 0), (13, 7, 1, 0), (13, 9, 0, 0)
α_1 = 15: (15, 4, 2, 0), (15, 6, 1, 0), (15, 8, 0, 0)
α_1 = 17: (17, 3, 2, 0), (17, 5, 1, 0), (17, 7, 0, 0)
α_1 = 19: (19, 4, 1, 0), (19, 6, 0, 0)
α_1 = 21: (21, 3, 1, 0), (21, 5, 0, 0)
α_1 = 23: (23, 4, 0, 0)
α_1 = 25: (25, 3, 0, 0)
α_1 = 27: (27, 2, 0, 0)
α_1 = 29: (29, 1, 0, 0)
α_1 = 31: (31, 0, 0, 0)
```

Note: every tuple has α_1 odd (a parity consequence of α_2 + 2α_3 + 4α_4 = (63 − 1 − 2α_1)/2 ≥ 0).
Note: every tuple has α_4 = 0 (since 16α_4 ≤ 62 ⟹ α_4 ≤ 3, but the
downward bound α_4 ≥ 1 ⟹ α_3 ≥ 4 ⟹ α_2 ≥ 6 ⟹ α_1 ≥ 4, and the
arithmetic only fits for α_1 ≤ 3 in that branch — easy to verify).

The α_3 ≥ 1 cases are interesting: they require new structural lemmas
to handle disjoint triples of odd cycles, which the H = 7 SCC argument
does not address.

## Future investigations inspired

### F-1. Develop a "kill table" for HYP-1754

For each of the 37 candidates, identify the structural lemma needed.
Some will be analogous to H ≠ 7 (single-SCC argument); others will
need new ideas (disjoint-triples handling).  Build the table; then
prove the lemmas.

### F-2. Conjecture: forbidden-H growth rate

Define `f(N) = |S(N)|` = number of arithmetic candidates for H = N.
From the formula
```
1 + 2α_1 + 4α_2 + 8α_3 + 16α_4 + … = N
```
the count `f(N)` is a partition function with 2-adic weights and
binomial-descent constraints.  Is there a closed form?  Does
`f(N) ~ N^c` for some constant c?

### F-3. Beyond α_4

The OCF-truncated form
`H = 1 + 2α_1 + 4α_2 + 8α_3 + 16α_4`
is sufficient for H ≤ 31.  For H = 63 we need α_5 ≤ 1, and for
H = 127 we need α_6 ≤ 1.  The enumerator extends mechanically; the
structural lemmas don't — each new α_k requires arguments about
k-tuples of pairwise-disjoint odd cycles.

### F-4. The full binomial bounds give an "α-polytope"

The constraints
```
α_k ≤ binom(α_1, k),  α_k ≥ 1 ⟹ α_j ≥ binom(k, j)
```
plus `Σ 2^k α_k = N` define a finite integer polytope.  Its vertices
are the α-tuples in `S(N)`.  Is there a nice description as a face
lattice?

### F-5. A formalisation-driven canon principle

This session's experience suggests a meta-principle:

> When formalising a proof, look for the cleanest axioms that suffice.
> The minimal axiom set is often sharper than the informal version.

The transition `alpha_subset_bound → alpha_descent` is a concrete
instance: the informal proof used the special case, but the general
case is *both* easier to state and strictly stronger.

This principle predicts: other project axioms (in `definitions.md`,
in `OCF.lean`-style modules) likely admit similar sharpenings.  Each
sharpening is a new tangent.
