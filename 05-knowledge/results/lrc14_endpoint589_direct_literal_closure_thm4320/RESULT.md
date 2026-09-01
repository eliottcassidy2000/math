# Endpoint-589 direct literal closure

Status: **FINITE-EXACT RELATIVE TO THM-4318**. This packet closes the fixed
endpoint-589 layer of the current thirty-label Haar proof graph. It gives no
physical entry map and does not prove LRC(14).

## Exact carrier exit

The unchanged 3,925-mask endpoint-590 exchanged carrier was replayed on all

```text
28 * C(30,9) = 400,600,200
```

endpoint-589 pair/body cases. Its only failures are

```text
(50,589): 20,025 bodies,
(96,589):     11 bodies.
```

All other bodies have a protected-joint or active-nonjoint carrier witness.
The paired O2/O3 raw, pair, and failure ledgers agree byte-for-byte.

## Direct literal repair

For one hostile pair `(q,589)`, partition the circle by every wall of the
thirty-label pool together with `q` and `589`. On a pair-safe open cell, let
`F` be the pool labels that fail there, and let `w(F)` be the total cell
width. For a failed nine-body `B`, the literal safe-set mass in grid units is

```text
M_q(B) = sum_{F intersect B = empty} w(F).
```

The audits independently compute the smaller quantity

```text
L_q(B) = sum_{F intersect B = empty, |F| <= 9} w(F) <= M_q(B).
```

Both a self-contained C++ wall reconstruction and a clean-room Python wall
reconstruction agree body-for-body. A third raw-cell/class-aggregation audit
also reconstructs the omitted higher-rank mass and checks the displayed
inequality exactly. The extrema are

| row | bodies | grid | low classes | minimum `63 L_q(B)-4 grid` | body |
|---|---:|---:|---:|---:|---:|
| `(50,589)` | 20,025 | 2,827,379,709,554,400 | 2,383 | 14,566,818,763,788,984 | `013c6401` |
| `(96,589)` | 11 | 1,130,951,883,821,760 | 2,352 | 7,172,391,058,639,758 | `0d0c6401` |

Every lower bound is strict. Since endpoint-only wall points have measure
zero, every failed body therefore has literal Haar mass greater than `4/63`.
Together with the carrier witnesses on the complement of the failure ledger,
this closes all twenty-eight endpoint-589 rows without adding or deleting a
mask.

## Fixed-fifty bridge

Using the `C/O` split of THM-4234 and THM-4242, the 20,025 `(50,589)` carrier
failures have exact petal-count distribution

```text
k:      0     1     2     3     4    5
count: 870  4519  7803  5502  1293   38.
```

The first 18,694 bodies (`k<=3`) already lie in inherited fixed-fifty chart
layers; the literal audit supplies all 1,331 remaining carrier failures in
the four- and five-petal layers, whose minimum truncated surplus is
15,777,555,364,138,176 at `20744601`. Consequently THM-4242's complete
fixed-fifty pool ray `r>=590` extends by the single missing row `r=589`.
This is an exact labelled-body bridge, not a claim that petal count preserves
margin, carrier multiplicity, or response structure.

## Labelled failure structure

The q=50 carrier-failure hypergraph has full support, empty common
intersection, and thirty distinct vertex degrees. Therefore every
label-permuting automorphism fixes all thirty vertices: the labelled
automorphism group is trivial, so there is no nontrivial orbit quotient hiding
the 20,025 bodies. Labels 95 and 193 are nevertheless strong hubs: 17,454
bodies contain both, 1,217 contain only 95, 1,347 contain only 193, and only
seven contain neither. Those seven all contain labels 80 and 168. As a hostile
control, 95 and 190 both have gcd 19 with 589 but have failure degrees 18,671
and 3,465, respectively; gcd data alone does not recover this structure.

This census is a diagnostic sidecar, not a dependency of the literal closure.
It explains why forcing another symmetry-quotiented response exchange is a
poor first representation at this endpoint.

## Typed successor

Consuming the twenty-eight closed rows gives

| ledger | count | ordered FNV64 | SHA-256 |
|---|---:|---:|---|
| `post589_typed/typed_union2141.csv` | 2,141 | `c84bb7f7eaa0f230` | `09f6d765afd78c03b60e7c4d047cae7df883ef8b2782f4256ee0e7c7be38ee09` |
| `post589_typed/final_residual20506.csv` | 20,506 | `3cd0863a93c7602e` | `c782f56439163ac6e9f4b5f230cde1397c048e7ae3308bcc08d449aeb4ede2d9` |
| `post589_typed/residual_top588.csv` | 66 | `18cf9a572cf9a5be` | `a490630678d4ba088b98e24f2ca3098fe2b9407e197ace1a4696df823cfeda69` |

The new maximum residual endpoint is 588 with 66 rows. Normal and optimized
typed derivations are byte-identical.

## Scope

The result is finite-exact for the frozen pool, carrier, typed universe, and
endpoint-589 layer. It demonstrates that the rank-eight/rank-nine response
carrier can fail as a representation even when the target literal mass is
comfortably positive. It does not prove a reusable arbitrary-pair theorem,
a physical allocation action, a terminating descent, or LRC(14).
