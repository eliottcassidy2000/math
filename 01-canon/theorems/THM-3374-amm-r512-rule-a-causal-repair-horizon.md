---
id: THM-3374
title: "AMM 12592: the R=512 rule-A death has causal repair horizon at least 58"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED
source: kps-s181
depends_on:
  - THM-3329
related:
  - THM-3332
  - THM-3373
companion: 04-computation/amm12592_r512_ruleA_causal_repair_horizon_kps_s181.py
output: 05-knowledge/results/amm12592_r512_ruleA_causal_repair_horizon_kps_s181.out
script_sha256: e2c1a18fe35d74cd445deef549f6d3eda8f87b4a472d340b809f29e3f1935a1d
output_sha256: 8faf5284d478ce1bef1fe23a04c06084ad1221c63ac259855e33207faaf75569
hash_basis: LF-normalized bytes
audit: independent coefficient-capacity, recurrence-indexing, and 57/58 threshold reconstruction
---

# THM-3374 -- the rule-A death cannot be repaired locally

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## Statement

At the exact golden-floor epoch `R=512,D0=0`, THM-3329's rule A survives rows
`0,...,106` and reaches row `107` with error-coordinate constant

```text
-N,
N=199181334599768561751691151979867147451295075845970943582846950031770442839710071820.
                                                               (1)
```

Let another admissible prefix have the same exact degree word and survive row
`107`.  Then it must differ from rule A in at least one row

```text
                              i <= 49.                 (2)
```

Equivalently, agreement through row `49` is impossible.  In particular, no
change confined to the five immediately preceding correction rows
`102,...,106` can repair the death, even if those five changes are allowed to
range over all Lucas-box admissible blocks.  Hence THM-3373's width-five
`R=8` operation law cannot be transplanted as a last-minute row-107 patch.

This does not exclude width-five moves applied earlier or repeatedly.

## Coefficient-capacity lemma

Write a degree-`d` admissible block as

```text
P_Delta(x)=sum_(k=0)^d delta_k x^(d-k)(1-x)^k,
|delta_k|<=C(d,k).                                      (3)
```

For any two such blocks and every `0<=s<=d`,

```text
|[x^s](P_Delta-P_Delta')| <= 2^(s+1) C(d,s).            (4)
```

Indeed, put `h=d-k`.  Only `0<=h<=s` contributes to `[x^s]`; the contribution
from the `h`-term is at most

```text
2 C(d,h) C(d-h,s-h).
```

Vandermonde's identity gives

```text
sum_(h=0)^s 2 C(d,h) C(d-h,s-h)
 =2 C(d,s) sum_(h=0)^s C(s,h)
 =2^(s+1) C(d,s),                                      (5)
```

which proves `(4)`.  Parity can only reduce this coarse capacity.

## Causal unrolling

The exact residual recurrence is

```text
sigma_i=(sigma_(i-1)-P_(Delta_i))/x.                   (6)
```

Suppose two surviving prefixes first may differ among the `L` rows preceding
row `107`.  Repeatedly applying `(6)` shows that the difference of their row-
`107` input constants is bounded by

```text
C_L=sum_(s=1)^L 2^(s+1) C(d_(107-s),s).                (7)
```

The exact floor word and `(7)` give

```text
C_5  = 3478916701028                                    (42 bits),
C_57 = 51445361523493837506256624150514290414698933203924597988799840246920133025721134820
                                                        (275 bits),
C_58 = 418636702487565989501705263678183068399437336707362038805929622565164565084110761700
                                                        (278 bits).             (8)
```

Thus `C_57<N<=C_58`.  A prefix agreeing with rule A through row `49` can
differ only in rows `50,...,106`, so `L=57`; its fatal input remains too far
from the only survivable error constants `0,2`.  This proves `(2)`.  The
five-row conclusion follows from `C_5<N`.

The inequality `N<=C_58` says only that this triangle-inequality obstruction
becomes inconclusive when row `49` is allowed to move.  It does not construct
an admissible repair or prove that horizon `58` is sufficient.

## Exact replay and scope

The standard-library companion independently reconstructs the golden floor
word by Fibonacci--Lucas comparisons, replays the error-coordinate dynamics,
verifies all `107` surviving blocks against the Lucas box, checks the dyadic
parity-fire count is zero, and proves `(5)` directly for every term used in
`(7)`.  Normal and optimized runs are byte-identical.  Its LF-normalized
source hash is

```text
e2c1a18fe35d74cd445deef549f6d3eda8f87b4a472d340b809f29e3f1935a1d.
```

The corresponding stored-output hash is

```text
8faf5284d478ce1bef1fe23a04c06084ad1221c63ac259855e33207faaf75569.
```

An independent audit rederived `(4)` from the reversed Bernstein basis,
unrolled `(6)` with the row-107 indexing from scratch, recomputed `C_5`,
`C_57`, and `C_58`, and confirmed the row-49/50 boundary with no off-by-one.

This theorem is a necessary divergence bound for alternatives to one fixed
rule-A prefix.  It neither proves infeasibility of the `R=512,D0=0` polytope,
nor supplies a new point, nor transfers the twelve `R=8` atoms, nor changes an
AMM deadline constant.  The live positive target is now an early compiler
acting by row `49`, not another end-of-life clamp patch.
