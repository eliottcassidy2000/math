# Lonely Runner Feedback Loop: 14, 15, and Failed Disproofs

codex-2026-05-31 S364

The user's rule was useful: attack the 14-runner case, hit a wall, force a new
idea for 15, hit a wall there, then try to build a counterexample and return.
That cycle did not prove a new case, but it corrected the target.

## External Check

As of this session, the public finite frontier remains the reduced cases
`k<=12`, i.e. thirteen total runners.  The next natural case is reduced
`k=13`, fourteen total runners.  The checked anchors were Rosenfeld's
8-runner paper (`arXiv:2509.14111`), Trakulthongchai's 9/10-runner paper
(`arXiv:2511.22427`), Sungkawichai-Trakulthongchai's 11/12/13-runner paper
(`arXiv:2604.23906`), and Jensen's mixed-threshold/Fourier note
(`arXiv:2605.27941`).

The Sungkawichai-Trakulthongchai prime-field lemma is still the most relevant
model: it removes the tight class when `k+1` and the prime grid modulus are odd
primes.  The next frontier has `k+1=14`, so the field move is not directly
available.

## Cycle 1: Fourteen Runners

The S363 micro-staircase idea said: use the full cell arrangement of

```text
floor(14 * {i alpha}),  i=1,...,13,
```

rather than only the coarse samples `alpha=r/14`.

S364 made that exact.  There are `812` cell patterns, hence `11368` pairs
`(s, cell)`.  The first result was a dead end in the best possible way:
every scalar ramp

```text
v_i = m i mod 14
```

blocks every cell.  This includes the unit ramps, not just the nonunit ones.
So the naive lemma "every residue vector has a micro-staircase witness" is
false for a principled reason: scalar ramps are exactly the Dirichlet equality
spine in residue form.

## Cycle 2: Fifteen Runners

Pushing the obstruction to `n=15` gave the same answer.  The arrangement has
`960` patterns and `14400` candidate pairs, and every scalar ramp
`v_i=m i mod 15` blocks all of them.

This made the correction clear.  Unit scalar ramps are just reindexed initial
segments.  Nonunit scalar ramps are the composite quotient cases.  They should
be routed into divisibility/descent arguments, not fought by the generic
micro-staircase lemma.

## Cycle 3: Disproof Construction Attempts

The third route tried to turn the gate condition into actual open covers.  For
the 14-runner case, deterministic families such as

```text
(1,2,...,12,14)
(1,2,...,11,14,28)
(2,4,6,7,8,10,12,14,21,28,35,42,49)
```

all had positive gaps or boundary witnesses.

For the 15-runner case, the analogous `15`-gated and `3x5`-factor ladders also
failed to become open covers.  Random gated searches up to speed `90` for
`n=14` and `105` for `n=15` found only positive-gap cases.  The best sampled
14-gated gap ratio was about `0.043269`; the best sampled 15-gated gap ratio
was about `0.038177`.

That is not a proof, but it is useful negative evidence: simply adding the
forced `n`-divisible speed does not seem to create the endpoint-protection
core a counterexample would need.

## The New Proof Target

HYP-1817 should be read with HYP-1818 attached:

1. Prove and excise the scalar-ramp line for all `n`.
2. Send nonunit ramps through quotient/divisibility descent.
3. Search only non-scalar residue vectors for genuine micro-staircase
   obstructions.
4. Certify the missed cells of the best near-blockers.

The current best non-scalar blockers from S364 are:

```text
n=14: (8,2,10,4,12,13,0,8,2,10,4,12,6), covers 11312/11368.
n=15: (0,0,0,0,0,0,0,0,0,0,0,0,0,10), covers 14280/14400.
```

Those missed cells are the new live edge.  If they have a short structural
description, the 14-runner proof target gets much smaller.  If they do not,
the right next tool is probably Jensen-style mixed thresholds: let different
coordinates have slightly different safe margins inside the missed cells and
look for a Fourier or endpoint-pressure certificate.

## Files

- `04-computation/lonely_runner_feedback_loop_s364.py`
- `05-knowledge/results/lonely_runner_feedback_loop_s364.out`
- `05-knowledge/hypotheses/HYP-1818-lrc-scalar-ramp-excision.md`
