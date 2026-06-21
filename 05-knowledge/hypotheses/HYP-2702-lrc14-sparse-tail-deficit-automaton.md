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
