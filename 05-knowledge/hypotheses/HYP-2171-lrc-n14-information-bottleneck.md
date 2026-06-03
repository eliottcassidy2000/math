---
id: HYP-2171
status: SYNTHESIS supported by S612 information-bottleneck audit; sufficiency theorem open
source: user-2026-06-03; codex-2026-06-03-S612
related:
  - HYP-2170
  - HYP-2168
  - HYP-2167
  - HYP-2166
  - HYP-2165
  - HYP-2164
  - HYP-2161
  - HYP-2157
  - HYP-2156
  - HYP-2155
  - HYP-2151
  - THM-401
  - THM-406
  - THM-407
---

# HYP-2171: recent n=14 progress is an information-bottleneck theorem schema

## Claim

The recent LRC n=14 improvements should be understood as an information
bottleneck, not merely as a smaller enumeration.

The proof task is to find the smallest statistic of an integer speed row that
is still sufficient for the floor predicate. HYP-2166 supplies the low-rate
base statistic, while HYP-2165 and HYP-2167 identify two side-information
channels that cannot be discarded:

```text
raw integer row
  -> Res_27 quotient / proof atom
  -> owner route and carry cocycle
  -> floor predicate.
```

The desired theorem shape is a conditional-independence statement:

```text
Floor(V) is independent of the raw integer presentation
given (Res_27 proof atom, owner route, carry cocycle, Cprime window).
```

This is not yet proved. It is the abstract structure suggested by S610/S611
and made explicit by S612.

## S612 Evidence

`04-computation/lrc_n14_information_bottleneck_s612.py` separates address
entropy from relevance entropy.

Address entropy collapses through the recent quotient tower:

```text
raw 13-subsets of nonzero residues:  10,400,600 rows, 23.3102 bits
unit-shell rows:                        340,928 rows, 18.3791 bits
D/U/N survivors:                         27,733 rows, 14.7593 bits
proof-obligation types:                     148 rows,  7.2095 bits
S610 proof atoms:                            11 rows,  3.4594 bits
primitive floor atoms:                        2 rows,  1.0000 bits
```

But proof relevance is not proportional to address entropy:

```text
pinch status on D/U/N survivors: H = 0.0016 bits
owner route through slack <=81:  H = 0.9981 bits
scalar carry shadow route:       H = 0.4138 bits
local carry route, weight <=2:   H = 0.0000 bits
```

The tiny `0.0016` bits for pinch status do not make the floor predicate
unimportant. They say floor rows are rare. Rare events can carry the whole
proof.

The carry side channel is the key example. In the S611 scalar AP/`V*` audit:

```text
scalar probes: 36
route counts: {'floor shadow': 3, 'strict shadow': 33}
carry indicator counts: {'zero carry': 3, 'nonzero carry': 33}
I(carry_nonzero; shadow_route) = 0.4138 bits
```

Within this scalar probe family, `carry_nonzero` perfectly determines whether
the least-positive shadow is floor or strict. So a quotient that forgets carry
can be small and still not conservative.

## Interpretation

HYP-2164 is a low-rate least-positive section. It classifies the base
`Res_27` quotient and proves that the primitive floor shadows are exactly AP
and `V*`.

HYP-2166 is the bottleneck layer. It compresses the proof surface to an
`11`-atom ledger by composing the clock exit, shell fold, pinch quotient, and
owner reattachment.

HYP-2167 identifies the side-channel leak. For `v=r+27k`, the identity

```text
27 == -1 (mod 14), so v == r-k (mod 14)
```

means the carry vector is not decorative metadata. It is the descent datum
that the n-clock sees.

Thus the right abstraction is not "discard as much information as possible."
It is:

```text
discard exactly the information whose conditional mutual information
with the floor predicate is zero, and retain the side information whose
conditional mutual information is nonzero.
```

Owner labels and carry cocycles currently appear to be the essential side
information.

## Anti-Poisson Reading

The anti-Poisson frame says the independent/free baseline predicts a positive
ground cell, while structured arithmetic cancels it at the floor. S612 adds a
rate-distortion warning to that frame.

The positive free term can be huge in address space and tiny in relevance
space. Conversely, a rare bit such as the carry indicator can have low Shannon
entropy but high proof value. So the proof should not optimize only for
compression. It should optimize for sufficient forgetting.

In coimage language, S610 found a candidate minimal sufficient quotient.
S611 showed that the least-positive section alone is not sufficient. S612
therefore proposes the corrected coimage:

```text
corrected statistic =
  (Res_27 proof atom, owner route, carry cocycle, Cprime window).
```

Incoming HYP-2168, the C2b multi-angle assault, fits the same information
picture. It reframes the multiple-of-n branch as `p_0>0`, a positive-measure
lonely interval rather than a measure-zero floor shadow, and isolates
`2n-1` tick data plus small-q clock data as the relevant probes. In S612
language, that says the Cprime window is not an afterthought; it is another
side channel needed for a sufficient statistic.

The renumbered HYP-2170 tight-census enumeration fits as the complementary
finite-evidence side: AP and `V*` are the tight family seen in the bounded
census, while the bottleneck reading explains why that finite fact still needs
owner/carry/Cprime side data before becoming a global proof.

## Tournament Analysis

S612 runs Tournament Analysis on information channels rather than on runners.

```text
vertices: information channels/proof summaries
observable: retained address bits, predicate entropy, side-information flag, role
switch: lower proof burden -> higher proof burden, label tie path
```

The fingerprint is transitive:

```text
vertices: 9
SCC count: 9
largest SCC: 1
directed 3-cycles: 0
score histogram: one vertex at each score 0..8
edge flips between proof-burden and entropy gauges: 6/36
```

The proof-burden Hamiltonian path is:

```text
scalar carry bit
S610 proof atoms
owner route
proof-obligation types
64 fixed classes
THM-407 gcd strata
unit-shell rows
D/U/N survivors
raw Res_27 subsets
```

The edge flips are the signal: entropy-only ranking and proof-burden ranking
do not agree.

## Assumption Challenge

Candidate tournament vertices considered:

```text
runners,
gaps,
fixed circle sections,
section boundaries,
wall-crossing events,
residues,
cover arcs,
Fourier modes,
matroid circuits,
proof obligations,
information channels.
```

This hypothesis chooses information channels because the user asked what the
recent progress fundamentally means. The preserved predicate is sufficiency
for the n=14 floor/certificate route. The destroyed information is raw row
identity, arbitrary unit representative data, and unneeded address labels.

The challenged assumption is that fewer cases automatically means a stronger
proof. S612 says fewer cases only help when the quotient is sufficient.

## Proof Target

A plausible theorem target is:

```text
Let V be a primitive n=14 row after the HYP-2163 no-multiple split.
Let Q(V) be its Res_27 proof atom, owner route, carry cocycle, and Cprime
window. If two rows have the same Q(V), then they have the same floor
classification or discharge to the same certificate class.
```

Equivalently, the remaining proof should show that every raw integer lift
outside AP/`V*` either changes the carry/owner statistic into a known strict
certificate, or falls into the Cprime CRT contradiction.

## Honest Status

What S612 proves by computation:

```text
1. the recent quotient tower is a large address-entropy compression;
2. relevance entropy is concentrated in rare floor/owner/carry predicates;
3. the scalar carry indicator has perfect mutual information with the
   least-positive shadow route inside the S611 scalar probes;
4. proof-burden and entropy-only gauges disagree on 6 of 36 channel pairs.
```

What remains open:

```text
prove the corrected statistic is sufficient for the global n=14 lift theorem.
```

## See

`04-computation/lrc_n14_information_bottleneck_s612.py`,
`05-knowledge/results/lrc_n14_information_bottleneck_s612.out`,
`07-reflections/lrc-n14-information-bottleneck-s612.md`,
`07-reflections/lrc-c2b-multi-angle-assault-s599.md`,
`05-knowledge/hypotheses/HYP-2167-lrc-n14-carry-fiber-conservativity.md`,
`05-knowledge/hypotheses/HYP-2166-lrc-n14-res27-quotient-tower-conservativity.md`,
`05-knowledge/hypotheses/HYP-2161-coimage-yoneda-2nm1-resonance-cancellation.md`.
