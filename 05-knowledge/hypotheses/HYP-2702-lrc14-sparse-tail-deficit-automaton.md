---
id: HYP-2702
title: LRC14 sparse residual tails should be absorbed by a generated-context deficit automaton
status: OPEN; exact sparse-tail scout supports singleton/coherent neutralization
source: codex-2026-06-20-S64
tangent: T937
depends_on:
  - HYP-2697
  - HYP-2698
  - HYP-2700
related:
  - HYP-2701
  - HYP-2699
  - HYP-2694
  - THM-557
  - THM-558
  - THM-555
---

# HYP-2702 - Sparse Tail Deficit Automaton

## Claim

After HYP-2700, consecutive blocks already dominate arbitrary bounded cluster
shapes on the top residual layers `|R| >= 5` and on the true full-cover
functional `p_0`.  The remaining HYP-2697/HYP-2698 obstruction is the sparse
tail `|R| <= 3`, where nonconsecutive shapes can beat the consecutive block on
named residual packets.

The concurrent HYP-2701 survival middle-mass gate is signal for this route:
the right coordinate changes when the basis changes.  HYP-2699 diagonalizes
Bonferroni `U4` in a survival basis; HYP-2702 should analogously diagonalize
the sparse residual tail in a generated-context / miss-zeta basis before
scalarizing.

After the S64 scout, the actual HYP-2701 computation landed as
`lrc14_truewide_survival_middle_mass_codex_s64.py`: it verifies the exact
survival currency `p1+p2+p3+p4-4p6`, finds three k=8 floor failures but no cap
failures in the true-wide box, and finds zero floor/cap failures for k=9 and
above in the audited boxes.  This strengthens the analogy: the true-wide gate
is controlled only after the survival basis change, while the sparse-tail gate
is controlled only after the generated-context basis change.

The expected repair is that those sparse-tail wins are absorbed by the
generated context word:

```text
x -> w_x(R),  R subset {1,...,6},
```

where `w_x` comes from the finite sector-mask OR/deletion automaton and is best
read in miss-zeta product coordinates.  The finite proof object should be a
**deficit automaton** whose states record

```text
(|R|, sector geometry of R, hit-count deficit, miss-zeta context weight).
```

The automaton vertices are not runners or sectors alone.  They are proof-state
deficits: a residual packet that wins for an arbitrary shape, together with the
minimum generated context pressure needed to neutralize that win.  Edges compare
state transitions under adding singleton context carriers, merging singleton
context carriers into coherent blocks, or lifting `|R|` toward the already-safe
large-residual layers.

## Evidence From The S64 Scout

The exact scout
`04-computation/lrc14_sparse_tail_deficit_automaton_codex_s64.py` stores its
run in
`05-knowledge/results/lrc14_sparse_tail_deficit_automaton_codex_s64.out`.

It uses exact Fraction arithmetic over the same `Z/7` sector-mask breakpoints
as HYP-2698.  To keep the run finite and directly tied to the known
obstruction, it takes a full bounded census for size `3` (`max(E)<=13`) and the
HYP-2698 near-consecutive frontier for sizes `4..6`.

Sparse coordinate wins are real and geometrically concentrated:

```text
size 3 full span<=13: 56 coordinate-violator shapes
  |R|=2: max gain 5/588 at C=(0,1,4), R=(5,6)
  |R|=3: max gain 5/294 at C=(0,1,3), R=(3,5,6)

size 4 frontier: 7 violator shapes, including |R|=4 max gain 3/196
size 5 frontier: 8 violator shapes, |R|=3,4 only
size 6 frontier: 1 violator shape, C=(0,1,2,3,4,6)
```

The dominant residual geometries are low-layer and spread/collar packets, not
full cover:

```text
|R|=3 gaps (1,2,4): 464 coordinates
|R|=2 gaps (1,6):   295 coordinates
|R|=3 gaps (1,1,5):  64 coordinates
```

Every tested sparse-tail coordinate violator is neutralized by the singleton
product context at total size `7`:

```text
size 3, r=4: 56 tested, 0 failures, min delta 20/16807
size 4, r=3:  7 tested, 0 failures, min delta 37/16807
size 5, r=2:  8 tested, 0 failures, min delta 199/24010
size 6, r=1:  1 tested, 0 failures, min delta 69/6860
```

The exact coherent generated-context scan also has zero failures:

```text
size 3, r=4, contexts=5: 280 tests, min 20/16807 at [1+1+1+1]
size 4, r=3, contexts=3:  21 tests, min 37/16807 at [1+1+1]
size 5, r=2, contexts=2:  16 tests, min 199/24010 at [1+1]
size 6, r=1, contexts=1:   1 test,  min 69/6860 at [1]
```

Thus the weakest generated context in this scout is always the all-singleton
context.  This sharpens HYP-2698's proof route: first prove the symmetric
hit-count kernel inequality, then prove that merging singleton context carriers
into coherent blocks cannot decrease the margin.

The all-singleton context has a useful cellular-automaton interpretation.  On
a residual packet of size `j`, one singleton carrier applies the monotone death
rule

```text
j -> j-1 with probability j/7,
j -> j   with probability (7-j)/7.
```

Equivalently, for a fixed `t`-sector requirement, the probability that `r`
singletons cover it is the HYP-2698 kernel

```text
g_r(t)=7^-r sum_{a=0}^t (-1)^a binom(t,a)(7-a)^r.
```

So the first proof target is a finite one-dimensional death-chain inequality
for hit-count laws.  The second proof target is a context-merging monotonicity
lemma: coherent blocks should dominate the singleton death chain for the
compression margin, as all exact minima in the scout occur at all-singletons.

Tournament Analysis on proof carriers is transitive:

```text
miss_zeta_product_word
> coherent_OR_contexts
> singleton_hit_count_kernel
> survival_basis_signal
> large_R_stratification
> raw_coordinate_weights
```

The fingerprint is `score_hist={0:1,1:1,2:1,3:1,4:1,5:1}`,
`directed_3cycles=0`, singleton SCCs, and one Hamiltonian path.  The quotient
preserves generated-context pressure and rejects raw residual coordinates as
the losing scalarization.

## Assumption Challenge

The session explicitly rejects the easy vertex choices:

- **not runners:** raw runner vertices forget the generated residual language;
- **not residual coordinates:** arbitrary positive weights make the theorem
  false;
- **not binary tournament arcs:** the `1/2`-scale tournament is only a correlate
  of the exact `Z/7` coloring.

The quotient preserves the LRC predicate "context-generated residual pressure
times cluster capacity" and destroys the individual carrier phases.  The
challenged assumption is that sparse residual packets must be handled by a
separate ad hoc finite check; HYP-2702 tries to make them a finite automaton
frontier attached to HYP-2698's miss-zeta product word.

## Status

No LRC14 proof is claimed.  The scout supports the generated-context basis
route and identifies the singleton death-chain kernel as the first proof target.

## S65 Death-Chain Band Addendum

The follow-up scout
`04-computation/lrc14_death_chain_band_automaton_codex_s65.py` stores its run
in
`05-knowledge/results/lrc14_death_chain_band_automaton_codex_s65.out`.
It uses the S64 engine and decomposes the singleton death-chain margin into the
seven HYP-2703 slope bands.

The exact frontier result is sharper than the original target:

```text
72 sparse-coordinate frontier rows tested
aggregate singleton death-chain failures: 0
first-order hit-count stochastic-dominance failures: 72
negative per-band cells: 192 / 504
minimum aggregate margin: 20/16807 at size 3, C=(0,1,3)
coherent context tests, including bounded-worst extras: 329, failures 0
```

So the first proof target is not first-order stochastic dominance of hit-count
laws and not per-band monotonicity.  Both are false on every audited frontier
row or on many band cells.  The correct statement is a signed seven-band
death-chain inequality:

```text
sum_s band_delta_s(C) > 0,
```

followed by generated-context monotonicity.  The weakest coherent context is
still the all-singleton partition, with equality only at `[1+1+...+1]`.

The S65 bounded singleton slices also found no failures beyond the coordinate
frontier:

```text
size 3 span<=13: 56 shapes, 0 failures
size 4 span<=13: 261 shapes, 0 failures
size 5 span<=9:  124 shapes, 0 failures
size 6 span<=9:  125 shapes, 0 failures
```

The cellular-automaton sign examples make the obstruction legible.  The worst
aggregate row has band word `-0+0+0-`: two slow negative bands are outweighed by
two positive bands.  The most negative single band row has word `-++-++-`.
Thus HYP-2702 now matches HYP-2703/HYP-2704: local proof routes fail, while the
completed signed aggregate survives.

## S65 Cross-Domain Quotient Atlas

The second S65 scout
`04-computation/lrc14_gibbs_quantum_roadcoloring_bridge_codex_s65.py` stores
its run in
`05-knowledge/results/lrc14_gibbs_quantum_roadcoloring_bridge_codex_s65.out`.
It treats the user's Gibbs, Arnold cat map, road-coloring, Hopfield, double-slit
amplitude, Fubini-Study, and Clifford/T-gate suggestions as quotient tests.

The productive exact lemma is the road-coloring shadow.  Let states be residual
masks `R subset {1,...,6}` and let letter `a` delete `a` from `R` when present.
This monotone deletion automaton is not strongly connected, but its reset-word
counts are exactly the singleton death-chain kernel:

```text
#{words of length r resetting a fixed t-mask}
  = sum_{j=0}^t (-1)^j binom(t,j)(7-j)^r
g_r(t) = 7^-r times this count.
```

This gives a short proof sketch for why the all-singleton context is the right
first object: it is the random road-coloring reset probability of the deficit
automaton.

The other external lenses are useful mostly as guardrails:

```text
Gibbs which-band reweighting:
  every row has negative bands; smallest beta* = 183.619629.
  Any adversarial bias toward observed losing bands can kill positivity.

Arnold cat map over F_7:
  torus cycles {1:1, 8:6}; projective cycles [[0,4,2,6],[1,3,5,inf]].
  31/72 rows have a negative cat-cycle sum, and 72/72 have a negative prefix.

Hopfield/Hebbian memory:
  only two binary band-sign patterns occur, `-++-++-` and `-+++++-`.
  This is a descriptive attractor alphabet, not a proof.

Fubini-Study hit-count geometry:
  projective distance correlates with margin (corr 0.675197), but does not
  encode the signed kernel order.

Double-slit / Clifford-magic reading:
  measuring individual bands destroys interference; HYP-2707's Clifford/magic
  split says coarse stabilizer-like quotients are tractable but too lossy.
  HYP-2702 needs the generated signed context, the analogue of retaining the
  non-Clifford phase resource before scalarization.
```

Tournament Analysis after the cross-domain atlas again ranks proof quotients,
not runners:

```text
miss_zeta_generated_magic
> road_coloring_death_chain
> Gibbs_uniform_band_sum
> Hopfield_band_attractor
> hitcount_Fubini_Study_geometry
> cat_map_prefix_orbit
> which_path_abs_band_measure
```

Fingerprint: `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`,
`directed_3cycles=0`, singleton SCCs, and one Hamiltonian path.

This addendum also resolves the "what are edges and vertices?" assumption
challenge for this phase: vertices are not runners, arcs, bands alone, Gibbs
states, or neural classes.  They are proof quotients of the generated deficit
language.  The quotient must preserve the LRC predicate "generated context
resets residual deficits before final scalar comparison" and must not measure
which band carried the intermediate debt.

## Post-Rebase HYP-2705/HYP-2707 Signal

After the S65 atlas checkpoint, `origin/main` added HYP-2705 and strengthened
HYP-2707.  HYP-2705 owns the true-wide projective synchronizer formulation:

```text
actual survival currency
  = death-chain boundary + signed projective deviation.
```

It proves the Fubini-Study projection identity for the missed-depth
death-chain and records that the depth quotient has no negative Gibbs
obstruction after two decorrelated hits.  HYP-2702 should not duplicate that
true-wide route.  Its distinct contribution is the sparse-tail frontier:
band-level/FOSD quotients fail on raw sparse residuals, while singleton
road-coloring resets and generated miss-zeta context still neutralize them.

HYP-2707 now has a proved tournament-side Clifford core: `c3 mod 2` is a
GF(2) quadratic form, and its bilinear rank gives THM-555's parity expectation
by the Boolean Gauss-sum formula.  The transfer lesson for HYP-2702 is that
raw hit-count and band quotients are stabilizer-like coarse shadows; the
generated signed residual language is the sparse-tail analogue of the magic
phase resource.

The later HYP-2707 magic-spread update adds an important caution: magic is not
the same as being extremal or large.  On the tournament side, the H-maximizer
can live in a stabilizer score class while the largest magic spread occurs in
an unbalanced class.  The LRC analogue should therefore measure how much the
generated residual profile is unpredictable from the coarse death-depth /
hit-count quotient, not merely how large `p0` or a sparse coordinate is.
