---
id: THM-976
title: Scale-twelve Hamming-six owner-orthogonality obstruction
status: PROVED FINITE-EXACT — the complete 2,413,458,432-context primitive proper AP-centred common-scale-twelve Hamming-six bank reduces to 64 all-order-twelve projective sign transversals; literal replay of all 262,144 remaining unit words gives zero survivors and pairwise-disjoint owner obligations, independently reconstructed in C++ and Python
source: codex-2026-07-17-S66 scale-twelve exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-963, THM-969, THM-970, THM-974, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twelve_hamming_six_owner_orthogonality_codex_c12.cpp
  - 05-knowledge/results/lrc13_scale_twelve_hamming_six_owner_orthogonality_codex_c12.out
  - 04-computation/lrc13_scale_twelve_hamming_six_referee_codex_c12.py
  - 05-knowledge/results/lrc13_scale_twelve_hamming_six_referee_codex_c12.out
---

# THM-976 — scale twelve is killed by complete owner orthogonality

The primitive proper AP-centred common-scale-twelve Hamming-six sheet bank is
empty.

For `c=12`, the effective state alphabet has orders
`1,2,3,4,6,12` and twelve `(D,e)` states.  Exact leave-one-out lcm enumeration
gives 26,961 divisor words and 2,611,968 state words per support, hence

```text
924*2,611,968 = 2,413,458,432
```

labelled raw contexts.  Scalar capacity leaves 36,830 contexts across 85
order-multiplicity patterns.  Owner-local feasibility then collapses the bank
to exactly 64 contexts on 64 supports, all with six order-twelve coordinates;
the supports are precisely the sign transversals of `F_13^*/{+-1}`.

The remaining fibre has `64*4^6=262,144` global unit words.  Literal replay
finds zero survivors.  Uniformly on every support,

```text
unit words satisfying 0 owners: 3,808,
unit words satisfying 1 owner :   288,
all higher counts             :     0,
|O_o|=48,       O_o intersect O_o'=empty for o!=o'.
```

Thus the owner-obligation intersection graph is empty: scale twelve exhibits
complete owner orthogonality.

## Completeness and independent replay

The C++ certificate independently constructs all twelve effective `(D,e)`
states by literal CRT search, enumerates the 26,961 actual ordered divisor
words, and checks all `924*26,961=24,911,964` labelled order contexts.  Mask
cardinalities are proved unit-independent before scalar capacity is used, so
the weighted state count is an exact expansion of this grammar rather than a
sample.  Owner-local reachability is exact set union on twelve-bit masks.  The
last 262,144 words are then visited literally, with no orbit quotient.
Optimized, unoptimized, and AddressSanitizer/UndefinedBehaviorSanitizer builds
reproduce the frozen output byte-for-byte; strict warnings are clean.

The standard-library Python referee is algorithmically separate.  It derives
CRT bases algebraically, enumerates actual divisor tuples, checks the
provider/owner ratio reduction, uses set-valued owner reachability, and packs
only the unit-independent scalar sums into six carry-free byte lanes.  It
again visits all 262,144 terminal words.  Normal and `python -O` runs reproduce
its frozen output byte-for-byte.

Frozen SHA-256 values:

```text
C++ source     24cf43b8cdc9e8f05f887927c4cb033544ae71e6440b1c2b734b29d419f24652
C++ output     91bf92e3e2e08d75b10b9c9a9f4db57dfbfcd506840881d62216b0886275cbda
Python source  1fe7bf9f1c1ae47e0cb1781ffb3a5fa75360221974992bc7801e7243d976f068
Python output  e76dd23fb2a4d341d292ac9efef2dc9bf0e65cf86ef4b285d5312c40161b5276
```

## Tournament and vertex audit

Use owner obligations as vertices and nonempty pair intersection as the
binary observable.  Every one of the fifteen pairs is a tie.  Gauging by the
projective sign classes and completing ties along `1<2<3<4<5<6` gives the
transitive tournament: score sequence `0,1,2,3,4,5`, no directed triangles,
six singleton SCCs, and one Hamiltonian path.  The full gauge edge-flip
histogram ranges from zero through fifteen, showing that the completion is
coordinate telemetry rather than an invariant obstruction.  Owner
obligations, or dually their exact unit-word sets, are faithful.  Runner,
provider, residue, individual-sheet, and cover-arc vertices each discard the
simultaneous obligation disjointness.

This closes only the primitive proper AP-centred common-scale-twelve
Hamming-six face under THM-860's hereditary-lcm and common-sheet conclusions.
It does not close H5 ramification, non-AP-centred/deep-sheet packets, or global
sporadic emptiness.  Scale thirteen is excluded by primitivity in THM-860, so
the next legal common-scale face is `c=14`.  ∎
