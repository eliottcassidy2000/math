---
source: codex-2026-06-03-S604
status: SYNTHESIS / atlas; anti-Poisson read of LRC, coimage fundamentals, strong tournament residuals, and adjacent repo threads
tags: [LRC, anti-Poisson, coimage, p0-collapse, additive-chains, strong-tournaments, tournament-analysis, Helly, Vitali-wall, Collatz, partition-functions]
---

# Everything Through the Anti-Poisson Frame

The repo now has enough pieces to make "anti-Poisson" mean something sharper
than "weirdly not random."  The Poisson side is the free rare-event baseline:
independent arcs, independent residue hits, independent conflict constraints,
many small causes superposed.  The anti-Poisson side is the arithmetic
correlation that empties the ground cell anyway.

For LRC, THM-406 makes this exact:

```text
p_0 = sum_j (-1)^j S_j,
S_j = total j-fold danger-arc overlap,
{p_k} = spectral measure / coimage of the coverage-depth observable.
```

So the collapse rows are not merely AP imitations.  They are rows where the
whole overlap tower cancels, while the boundary witness skeleton still exists.
That is why the user's sporadic additive chains matter.  `(1,3,4,7)` and
`(1,3,4,5,9)` are not just extra examples; they say the anti-Poisson locus is
larger than the arithmetic progression and is controlled by additive/shell
relations rather than by AP membership.

## The Basic Picture

The free picture says:

```text
many sparse danger events -> Poisson-like depth -> p_0 positive.
```

The anti-Poisson picture says:

```text
structured arithmetic correlations -> all-orders overlap cancellation -> p_0 = 0.
```

The decisive object is the coimage: forget all distinctions that the observable
cannot see, keep exactly the pushforward/spectral law that still decides the
question.  In LRC the coimage is `{p_k}`.  In tournament OCF it is a partition
function or path-homology signature.  In determinant residuals it is the
intersection language of allowed multipliers.  In Collatz it is the rapidity or
linear-form image in which `2^E - 3^k` is visible.

That is what "fundamental" means here: a thing is fundamental when it is the
minimal object that still remembers the obstruction.

## Category + Number Theory

The category-theory side says:

```text
coimage = the minimal quotient that still decides p_0.
Yoneda  = the quotient is canonical because every natural probe recovers it.
```

The number-theory side says which probes matter at the LRC floor.  At
`delta=1/n`, set `C=2n-1`.  The antipodal pairs `{a,-a} mod C`, the unit action
of `(Z/CZ)^x`, and the gcd strata for composite `C` are the resonance observers.
If a unit shell is missed, multiplication by `a^{-1}` exposes a witness at
`a^{-1}/C`, so `p_0` cannot be zero.  If all unit shells are covered, the
remaining question is whether additive resonances align the higher overlaps so
that

```text
p_0 = sum_j (-1)^j S_j
```

cancels exactly.

So `2n-1` is not just a convenient denominator.  It is the number-theoretic
coimage of the floor witness problem.  Prime `C` makes every shell visible to
unit clocks; composite `C` creates nonunit lanes where sporadic cancellation can
hide.  For `n=14`, this means `C=27=3^3`, and the useful object is `Res_27(V)`:
unit-shell coverage plus gcd lanes plus additive-chain relations plus the depth
polynomial cancellation profile.

## Repo Atlas

LRC depth:
The anti-Poisson locus is `p_0=0`.  Generic rows sit near the free baseline with
positive lonely mass.  AP and sporadic additive chains empty the zero-depth
cell.  THM-406 says the mechanism is all-orders inclusion-exclusion, so the
classifier needs the full coimage, not just pair covariance.

Vitali wall:
Measure methods and finite Bonferroni truncations see only initial moments.
The residual can live beyond every fixed moment order.  Anti-Poisson is the
constructive arithmetic side of that wall: classify the exact correlations
that force the ground cell to vanish.

Helly entropy:
Small Helly witnesses are cheap anti-Poisson certificates.  Full-Helly rows are
isostatic: no small subfamily explains the collapse.  The Helly rank is the
certificate entropy of the anti-Poisson event.

Two-block determinant:
The free model is that component languages intersect independently.  The
anti-Poisson event is an empty intersection caused by determinant alignment,
often witnessed by singleton or pair Helly certificates.  If a live strong SCC
survives in the proof-obligation tournament, that is the next genuine branch.

Collatz:
The drift model is the free side.  A putative cycle is anti-Poisson: an exact
tail where powers satisfy the correlated linear form `2^E - 3^k`.  The trivial
cycle has the boundary identity `3+1=4`; nontrivial cycles would require a
coimage-level cancellation against the drift.

Tournament OCF / partition functions:
The hard-core partition function is the coimage.  Random conflict behavior is
the free side.  Forbidden values, non-real-root strata, and path-homology
defects are anti-Poisson-style cancellations in the coefficient/state space.

Additive bases:
High representation entropy is the free side.  Zeckendorf-type uniqueness,
carry constraints, and exact-cover residue phenomena are anti-Poisson
constraints: many apparent representations are collapsed to one canonical
support or to none.

## Why Strong Tournaments Matter

The strong tournament subset matters when the tournament is built from the
right vertices.  A transitive tournament quotient means the proof has a
one-dimensional priority order: peel the first exit, then the next, and so on.
A strong tournament means the quotient has no source/sink simplification.  Its
vertices mutually depend on each other.

That is the tournament analogue of anti-Poisson behavior.  A strong SCC in a
proof-lens tournament says no finite sorted statistic has yet explained the
ground-cell cancellation.  It is the cyclic residue where all-orders structure
can hide.

For LRC, the vertices should be challenged before use.  Runners and raw arcs are
usually too literal.  Better candidates include coverage-depth cells, unit and
nonunit shell gaps, wall-crossing events, residue classes, determinant component
languages, cover arcs, Fourier modes, matroid circuits, and proof obligations.
The chosen quotient should say exactly which predicate it preserves, such as
`p_0>0` versus `p_0=0`, and exactly what it destroys, such as endpoint phase
order.

## What This Means For N=14

For `n=14`, the question is not "is the row AP enough?"  The better question is:

```text
after every cheap projection,
is there still a strong anti-Poisson core?
```

Project through:

```text
depth coimage {p_k},
unit/nonunit C=27 shells,
additive-chain labels,
two-block determinant languages,
Helly certificate rank,
proof-obligation tournament SCCs.
```

If these projections become transitive or small-Helly, the residual is probably
certifiable.  If a strong SCC survives, it should become a named proof object,
not just an annoying computational leftover.

## Next Moves

Add an anti-Poisson signature table to the `p_0` collapse script:

```text
APSig = (p_k, S_j, alternating partial sums, shell floor, Helly rank, tournament SCCs).
```

Compare AP, `(1,3,4,7)`, `(1,3,4,5,9)`, the two `n=8` sporadics, and non-chain
controls by their whole inclusion-exclusion profiles.  The first place where
the profile separates collapsed from non-collapsed rows is the first real
classification invariant.

Then port the same signature to determinant component languages and OCF
partition functions.  If the same anti-Poisson signature keeps reappearing,
the repo gets a shared residual calculus instead of one-off analogies.

## Honest Status

This is a frame, not a proof.  The rigorous anchors are THM-406
(`p_0` is all-orders inclusion-exclusion and `{p_k}` is spectral), HYP-2153
(`p_0` collapse is larger than AP), HYP-2154/HYP-2155 (free baseline and
coimage meaning), and the Helly/logarithm work HYP-2151/HYP-2152.  The new
claim is that anti-Poisson should be promoted to a named coimage-level
signature and used to organize the `n=14` proof residuals.
