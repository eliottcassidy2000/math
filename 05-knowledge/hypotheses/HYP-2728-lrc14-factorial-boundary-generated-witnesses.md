---
id: HYP-2728
title: LRC14 factorial boundary packets need generated-word witnesses
status: OPEN; exact frontier scout plus Lean finite identity check
source: codex-2026-06-21
depends_on:
  - HYP-2727
  - HYP-2726
  - HYP-2725
  - HYP-2724
  - HYP-2723
  - HYP-2722
  - HYP-2721
  - HYP-2720
  - THM-561
related:
  - THM-558
  - THM-534
  - HYP-2719
  - HYP-2698
  - HYP-2702
  - OPEN-Q-108
---

# HYP-2728: Factorial Boundary Packets Need Generated-Word Witnesses

## Claim

The factorial atom boundary is algebraically closed and easy: if

```text
W_j(q) = sum_{t>=j} binom(t,j) q_t,
```

then

```text
q_0 = sum_j (-1)^j W_j.
```

The finite-difference packet

```text
B_j(t) = (-1)^(j-t) binom(j,t),  t<=j
```

satisfies `W_i(B_j)=delta_ij` and `q_0(B_j)=(-1)^j`.

The proof problem is therefore not the boundary algebra.  It is the
generated-word filter: which formal atom-cone directions can arise from
miss-zeta product words, relation-code packets, and true LRC row contexts?

## Exact Evidence

Script:
`04-computation/lrc14_factorial_boundary_operator_codex_20260621.py`.

Stored output:
`05-knowledge/results/lrc14_factorial_boundary_operator_codex_20260621.out`.

The exact packet table verifies every `B_j` over the seven missed-count atoms.
The abstract cheap LP directions from HYP-2721 show the key obstruction:

```text
r=2,4,5 have W1=W2=0 and U4=0 in the abstract atom cone.
```

Thus low signed leakage does not by itself forbid the dangerous abstract
directions.

The generated frontier tells a different story.  Recomputing the HYP-2722/S71
test bank gives `318` generated-context tests and `0` `q0` failures.  The
robust witnesses on normalized generated moves are:

```text
|W1|+|W2| min = 144154/63487
U4          min = 2187/2005
tail45      min = 182/2005
```

while signed `W1`, `W2`, `W1+W2`, and `B2=-W1+W2` all have nonpositive rows.
So the generated frontier is not controlled by a simple signed `W1/W2/B2`
rule.  The live witnesses are absolute low leakage plus Bonferroni4/tail
readouts.

The nearest generated rows to the first cheap directions remain separated:

```text
r=1 distance 918112/308305
r=2 distance 9682/2479
r=3 distance 528924/61661
```

This is evidence that generated-word compatibility is a real structural
constraint, not a numerical accident.

## Follow-Up: Tail45 Strip Certificate

Follow-up script:
`04-computation/lrc14_generated_separator_certificates_codex_20260621.py`.

Stored output:
`05-knowledge/results/lrc14_generated_separator_certificates_codex_20260621.out`.

The normalized generated frontier lives in an exact `tail45=q5+5q6` strip:

```text
min generated tail45 = 182/2005
max generated tail45 = 10910/21539
```

The HYP-2721 cheap abstract directions have normalized tail45 values

```text
r=1:  1
r=2: -1
r=3:  3/2
r=4: -1
r=5: -1
```

Therefore a single covector separates all five cheap directions from the
generated frontier.  Directions `r=2,4,5` are below the generated tail floor;
directions `r=1,3` are above the generated tail ceiling.  This is sharper than
the earlier witness hierarchy: generated-word exclusion can be attacked as a
bounded `tail45` strip lemma before invoking relation-code/Delsarte packet
classification.

The same run confirms why signed low moments are not the certificate:

```text
W1, W2, and B2 fail on every cheap direction as one-sided named separators.
```

The sign-profile ledger has only five generated profiles for
`(W1,W2,B2,tail45)`, and all have positive tail45.

## Lean Formalization

Lean source:
`04-computation/lean/TournamentH7/TournamentH7/LRCFactorialAtom.lean`.

The module is deliberately self-contained: it defines a local seven-coordinate
binomial table and checks the finite packet identities by `native_decide`.
Theorems currently audited:

```text
basis_moment_delta
originCoeff_delta
basis_q0_sign
U4_basis
low12_basis
cheapScaled_q0
cheapScaled_tail45
```

Direct `lean` checking succeeds in this shell.  `lake build` is still blocked
by an external mathlib clone failure, so the stored lake output is a
source-handoff rather than a full project build.

## Connection To Incoming HYP-2726/HYP-2727

HYP-2726 says the cover bound is a Delsarte/Krawtchouk LP.  HYP-2727 says
relation-code signals are useful but must be applied after generated
miss-zeta word compatibility.  HYP-2728 supplies the finite boundary algebra
for the origin atom and explains why the order matters:

```text
formal B_j boundary -> generated-word exclusion -> relation-code/Delsarte
classification -> factorial odd-L1 tail envelope -> q0 evaluation.
```

The cheap atom cone contains directions invisible to `W1+W2` and `U4`; the
generated frontier does not.  This is the current sharp target.

Follow-up target after the separator scout:

```text
prove every generated miss-zeta frontier move satisfies
    0 < 182/2005 <= tail45 <= 10910/21539 < 1,
then classify the leftover row-level packets by Delsarte/relation-code data.
```

## Tournament Analysis

Vertices:

```text
|W1|+|W2|, U4, tail45, B2, W2, W1+W2, W1
```

Pairwise observable: fewer nonpositive rows on the generated frontier, then
larger exact minimum.  The tournament is transitive:

```text
|W1|+|W2| > U4 > tail45 > B2 > W2 > W1+W2 > W1.
```

Assumption challenge: the tournament vertices are not runners, arcs, or row
shapes.  They are proof-obligation witnesses on the normalized generated atom
frontier.  This quotient preserves which scalar witnesses can certify
generated compatibility and destroys row geometry, relation height, and
sector-label addresses.  Those destroyed data must re-enter through HYP-2727
and the Delsarte/relation-code classification rather than being assumed away.
