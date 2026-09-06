# Independent audit of the missing-interval lemma and three-clock exclusion

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [producer's theorem](continuing6_20260906_lrc_near_resonance.md) passes
independent proof review and an exact implementation that imports no
producer code. The general connected-complement necessary set is reduced
from8201 clocks, maximum14904, to8198 clocks, maximum14886. General LRC14
and the unchanged obligations on the retained clocks remain OPEN.

## Analytic referee

For q=14R+r with1<=r<=13, intersecting the native strict inequalities
||x||<1/14 and ||qx||<1/14 gives exactly2R+1 complete teeth centered at
k/q, -R<=k<=R. Their length is1/(7q). The exclusion of14|q is necessary
because its two extreme intersections are half teeth.

On a translated grid with n=7q-s>0, a tooth has projected length
1-s/(7q)<1. It consequently contributes at most one point. Its absent
phases form a closed interval of length delta=s/(7q) centered at
1/2-sk/q modulo1. This is the exact complement of an open presence arc;
the closed endpoints cannot be discarded.

In unwrapped center order, adjacent gaps are separated by s/q=7delta.
The total span from the left endpoint of the first gap to the right
endpoint of the last is

    2Rs/q+delta = s(14R+1)/(7q).

When this is strictly less than1, all closed gaps are disjoint, including
across the cyclic join. At most one tooth is absent. An absence center
attains2R present teeth, while a phase outside the gaps attains2R+1.
The R=0 case has a single proper gap and gives the same statement.

For q=14R+1,s=7,R>=1, the unwrapped span is exactly1. Only the two
extreme closed gaps touch, at phase0 modulo1. That phase loses two teeth,
giving2R-1. Other gaps remain separated and their total measure is less
than1, so the maximum is2R+1. This independently verifies the stated
sharp endpoint failure, including the genuine atlas pair(1,43).

At q355,R25,s1,2,3 the strict inequality holds, and the exact normalized
minimum is50. A physical pair with sheet gcd6 therefore contributes
overlap credit300, including arbitrary coprime multipliers and translates.
This exceeds even the unconditional full-word maxima113,12,38.

The physical transfer is correctly typed. Under a hypothetical failure,
the inherited joint-shadow constraints provide the complete sheet word;
the connected strict graph provides some edge with sheet gcd<=6. The
entire edge/sheet universe is checked, so that existential edge has a
valid certificate. The proper six-component phase input remains CITED.
Pair overlap bounds the excess dangerous multiplicity pointwise, hence
positive credit minus cost supplies an actual weak-safe lift. A quotient
certificate is not asserted to characterize unsafe physical rows.

## Independent exact referee

The audit regenerates the whole5855-pair atlas by enumerating total sums
and independently factoring them. Native spatial endpoint walls and
midpoint modular predicates reconstruct each danger intersection. This
does not use the producer's pairwise interval-intersection algorithm or
its closed-form length list.

A subset-mask gcd recurrence tests all126 proper positional subsets of
every sorted seven-state multiset. It regenerates all177099 words in the
three unpruned domains and exactly10335 valid words, their ceiling costs,
the complete conditional maxima and the first attaining owners. Repeated
states preserve positional multiplicity. Every valid word and every
candidate is compared with the frozen certificate, not just aggregate
counts or maxima. Native component lengths independently recover all81970
sheet/edge candidates and the single global-cost survivor per clock.

Eight complete translated-grid profiles are regenerated from the actual
modular predicate at every wall and every open cell: three small strict
controls, two equality controls, and the three physical-clock consumers.
At every phase the result is also compared with the number of closed
absence gaps avoided. Full physical grids independently check the
multiplicity and coprime-unit map. Finally the previous necessary-set
array is read in full and exactly the three proved clocks are deleted.

The independent engine passes88593 always-active gates. Its raw LF output
must match normally and with optimization; the continuing6 manifest
records the exact source and output bytes. Reproduction:

    python -B 04-computation/continuing6_20260906_lrc_near_resonance_audit.py
    python -B -O 04-computation/continuing6_20260906_lrc_near_resonance_audit.py

No theorem-ID reservation, Lean proof, additional clock scan, strict-safe
phase, or unrestricted LRC14 conclusion is inferred.
