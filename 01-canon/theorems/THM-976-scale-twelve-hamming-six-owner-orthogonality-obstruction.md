---
id: THM-976
title: Scale-twelve Hamming-six owner-orthogonality obstruction
status: CLAIMED — exact exploratory census reduces the 2,413,458,432-context primitive proper AP-centred common-scale-twelve Hamming-six bank to 64 all-order-twelve projective sign transversals; all 262,144 remaining unit words fail and every distinct owner pair is disjoint. Frozen independent certificates are in progress.
source: codex-2026-07-17-S66 scale-twelve frontier probe
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-963, THM-969, THM-970, THM-974, HYP-6820]
---

# THM-976 — scale twelve is killed by complete owner orthogonality

This number is reserved for the scale-twelve closure now under certificate
hardening.

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

The remaining fibre has `64*4^6=262,144` global unit words.  The exploratory
literal replay finds zero survivors.  Uniformly on every support,

```text
unit words satisfying 0 owners: 3,808,
unit words satisfying 1 owner :   288,
all higher counts             :     0,
|O_o|=48,       O_o intersect O_o'=empty for o!=o'.
```

Thus the owner-obligation intersection graph is empty: scale twelve exhibits
complete owner orthogonality.  Primary and independent frozen certificates,
sanitizer replay, hashes, and the exact tournament/vertex audit remain to be
written before promotion to `PROVED FINITE-EXACT`.
