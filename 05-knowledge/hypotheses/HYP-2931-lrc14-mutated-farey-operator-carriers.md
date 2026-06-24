---
id: HYP-2931
title: LRC14 mutated Farey operator carriers
status: PROOF-INTERFACE / side-channel hierarchy; not a proof of LRC14
source: codex-2026-06-23-S130
related:
  - HYP-2930
  - HYP-2928
  - HYP-2926
  - HYP-2925
  - HYP-2924
  - HYP-2920
  - HYP-2917
  - HYP-2908
  - HYP-2899
  - HYP-2900
  - THM-572
---

# HYP-2931: Mutated Farey operator carriers

This hypothesis tests the prompt's four mutated Farey operators as possible
LRC14 proof carriers.  For a reduced fraction `p/q`, compare:

```text
sum      = p + q
product  = p*q
denpow   = q^p
numpow   = p^q
```

against the ordinary binding denominator `q`.

The exact LRC14 identity remains:

```text
M(S) = p/q,       e = 14p - q,       M(S)-1/14 = e/(14q).
```

Therefore `q` is the theorem-level binding scale.  The mutated operators are
only safe if they preserve the proof ordering induced by `e/(14q)`, or if they
carry a labelled side channel that explains why an attempted quotient lost
information.

## Computation

The script
`04-computation/lrc14_mutated_farey_tournament_codex_s130.py` stores output at
`05-knowledge/results/lrc14_mutated_farey_tournament_codex_s130.out`.

It audits three objects.

1. **Unit-excess child chain.**
   The fractions above `1/14` with Farey excess `e=1` are

   ```text
   p/q = p/(14p-1).
   ```

   Hence the binding denominator follows the additive lane

   ```text
   q -> q+14.
   ```

   This is the Farey version of the repo's `n+2` recursion: the numerator step
   is `+1`, the denominator step is the apex constant `+14`, and the gap is
   exactly `1/(14q)`.
2. **Farey-order locality in `[1/14,1/13]`.**
   For Farey order `120`, the local interval is too narrow to strongly
   separate the operators: `q`, `p+q`, `p*q`, and `q^p` have the same
   inversion count in the sampled interval, while `p^q` differs by one
   inversion.  This says local Farey adjacency alone does not identify the
   right proof coordinate.
3. **LRC14 row-bank proxy order.**
   On `749` AP/GW/petal/single-replacement rows, the true key is
   `M(S)-1/14`.  The proxy key is `e/W(p,q)`, where `W` is one payload.
   The ordinary denominator `q` and additive payload `p+q` agree with the true
   risk order on every non-tied comparison in this bank:

   ```text
   q        agree=258329 flip=0 tie=21797
   sum      agree=258329 flip=0 tie=21797
   product  agree=213778 flip=43903 tie=22445
   denpow   agree=170851 flip=87478 tie=21797
   numpow   agree=170890 flip=87439 tie=21797
   ```

   Product and power payloads create large inversions, for example between
   the Farey near-miss `12->36` and replacement rows such as `1->59`.

## Tournament Analysis

The explicit tournament-analysis choice is:

```text
vertices: q, sum, product, denpow, numpow
pairwise observables:
  - row-bank pairwise order agreement with M(S)-1/14
  - inverse Farey-order inversion count
switch/gauge: better score beats worse score
tie Hamiltonian path: q -> sum -> product -> denpow -> numpow
```

The risk, locality, and majority tournaments are transitive:

```text
score_hist = ((0,1),(1,1),(2,1),(3,1),(4,1))
c3 = 0
scc = five singletons
hp = 1
majority order = q > sum > product > numpow > denpow
```

This is useful because it orders the payloads as proof carriers rather than
leaving them as analogies.

## Readout

The additive payload `p+q` is the least destructive mutation.  In the audited
row bank it preserves the same non-tied risk order as `q`, and on the
unit-excess chain it remains linear.  It is therefore a good ledger for the
`n+2` / Stern-Brocot recursion.

The product payload `p*q` is the useful multiplicative side channel.  It keeps
both coordinates visible and echoes the repo's product-Mobius / `Div(D) x B_r`
language, but it is not a valid replacement for `q`: it already reverses many
global row-bank comparisons.

The power payloads are magnitude amplifiers.  They are poor local proof
denominators, but good stress tests for fixed-rule tournaments and scalar
quotients.  If a proposed invariant is stable under these power distortions,
then it is probably too magnitude-blind to prove LRC14.

## Proof target

HYP-2930's target is unchanged:

```text
Every non-AP/GW q=14 survivor either has nonunit excess e>1,
or enters a unit-excess Farey child class outside the tight optimum classes.
```

HYP-2931 sharpens how to use mutated coordinates around that target:

```text
q        = theorem-level binding scale
p+q      = additive/Farey recursion ledger
p*q      = multiplicative/coimage side-channel ledger
q^p,p^q  = magnitude-leak detectors for false quotient proofs
```

The missing theorem is still a rigidity theorem for the remaining reduced
LRC14 atom.  The new practical rule is: do not scalarize the fraction address
too early.  Keep `q` as the binding scale, use `p+q` and `p*q` as labelled
side channels, and use power payloads to expose invariants that accidentally
forget magnitude.
