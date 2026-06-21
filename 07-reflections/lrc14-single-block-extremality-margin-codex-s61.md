# LRC14 single-block extremality and the margin ledger

**codex-2026-06-20-S61.**  The useful move this session was to stop treating
HYP-2694 as one statement.  It is really three statements:

1. the decorrelated coherent-block partition function is maximized by one block;
2. the one-block value is below `cap_k`;
3. the finite carrier error is smaller than the remaining margin.

The first two are now exact in the coherent-block quotient.  The anchored speed
`0` is kept as its own cluster, and the `m=k-1` nonzero speeds are partitioned
into far coherent consecutive blocks.  Exact shared-`x` Fraction integration over
all integer partitions of `m=7..11` gives a clean answer: `[m]` is always the
unique top partition.  The closest split is always `[m-1,1]`.

```text
k=8:  D_7=283/1470,     cap-D_7=1111/5880,      split gap=1111/10290
k=9:  D_8=629/2058,     cap-D_8=111019/588588,  split gap=374/5145
k=10: D_9=16969/41160,  cap-D_9=102803/535080,  split gap=6561/96040
k=11: D_10=30551/61740, cap-D_10=184957/802620, split gap=42661/864360
k=12: D_11=71111/123480, cap-D_11=34729/123480, split gap=9047/172872
```

The third statement is also closed for actual shifted single blocks, though
with a conservative cutoff.  For

```text
E_M={0} union {M,M+1,...,M+m-1},
```

write `H(x,phi)` for the event that the shifted block covers all six inner
sectors.  On each `1/M` cell, the carrier phase `Mx` runs exactly once around
the circle.  Freezing `x` and comparing to the cell average leaves only the
`x`-variation of `H`.  For fixed carrier phase, the offset `d` crosses at most
`7d` sector walls, hence

```text
|p0(E_M)-D_m| <= 7*binom(m,2)/M.
```

This proves the single-block branch below cap once
`M > 779,1040,1312,1367,1369` for `k=8..12`.  The exact small samples are much
safer than the bound: at `M=19`, the largest observed error magnitude in the
S61 table is about `0.0125`, while the smallest cap margin is about `0.1886`.

## What changed

The phrase "single-block extremality" is no longer just inspiration.  In the
coherent-block quotient it is a finite partition-function theorem, THM-557.
The proof vertices are not runners or arcs; they are integer partitions of the
nonzero block mass.  The quotient preserves the HYP-2694 cover predicate in the
decorrelated model, but it destroys arbitrary within-cluster shape information.
That is the challenged assumption: consecutive coherent blocks are now proven
best among coherent-block partitions, but arbitrary bounded cluster shapes still
need a compression lemma before the theorem can be applied globally.

## Remaining target

The exact split gaps are now proof currency.  The multi-block branch should not
try to prove a tiny absolute discrepancy.  It only needs a joint carrier
decorrelation bound whose error is below

```text
cap margin + split gap
```

after any nontrivial split.  The single-block diagonal-freeze proof suggests the
right shape: freeze slow coordinates cellwise, let carrier phases sweep, and
pay only the total variation of the block-cover indicator.  The hard extension
is simultaneous carriers; this is exactly the HYP-2684/HYP-2694 joint
Erdos-Turan/Koksma obligation.

No LRC14 proof is claimed.  The branch is now sharper: coherent-block
extremality and single-block large-`M` error are closed; arbitrary-shape
compression plus joint multi-carrier error remain.
