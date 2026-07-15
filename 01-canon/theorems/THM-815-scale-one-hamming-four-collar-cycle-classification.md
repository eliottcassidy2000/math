---
id: THM-815
title: Scale-one Hamming-four collar cycles force a consecutive-doubling label packet
status: RESERVED by codex-2026-07-15-S10; symbolic proof complete, exact replay and full write-up in progress
source: codex-2026-07-15-S10 Hamming-four continuation
depends_on: [THM-806]
related: [THM-770, THM-810, HYP-6820]
---

# THM-815 — Hamming-four collar-cycle classification

For a proper scale-one four-replacement packet

```text
([12] minus R) union {u_r:r in R},
u_r=r (mod 13),                 |R|=4,
```

hypothetical tightness with every `u_r>24` forces the missing-label set to be

```text
R=a{1,2,4,8} in F_13^*.
```

Indeed the THM-806 left-collar argument gives positive cross-handoff indegree
at every owner. Two- and three-cycles are already impossible, so there is a
four-cycle. If `lambda_i` is its speed-ratio word, its exact half-open bands
have integer centres `a_i>=2`, with one `a_i=2`, every `a_i<=17`,
`product a_i=1 (mod 13)`, and

```text
product(a_i-1)<=16<product(a_i+1).
```

The modular/product inequalities force the cyclic band word `(2,2,2,5)`.
Thus three edges double the residue label and the fourth closes by
`5=2^(-3)`, proving the claim. The method-limit row with missing labels
`{1,2,4,8}` and replacements `(79,54,30,34)` realizes the four handoffs but
is loose with exact maximum `3/19`; the collar incidence is a reduction, not
the LRC predicate.

The completed write-up will include the elementary finite-band lemma, the
height-one corollary off the consecutive-doubling packets, analytic finite-box
bounds, and the replayable tournament audit.
