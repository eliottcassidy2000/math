---
id: THM-3320
title: "Projected-k3 z216 fourth ruler prefix and affine multicover closure"
status: >
  PROVED + VERIFIED-EXACT in the declared projected atlas.  The next
  deterministic intrinsic-ruler prefix after THM-3313 consists of row 141 in
  the complete gcd24/L76440 family and rows 133,219,359 in the complete
  gcd24/L30576 family.  Their exact necessary upper screens are empty:
  230=172 crude +58 status +0 residual.  The integral zero-constant modular
  template reports 55/1/2 status circuits of support 1/2/3, but the full
  common-table dual has 55/3 circuits of support 1/2: two apparent
  support-three packets admit exact half-affine support-two certificates.
  Exact empty- and proper-tail tables prove the stated minimalities, while an
  inherited row-138 packet remains genuinely support three.  The projected
  ledger moves 373157->373153, the z216 wall 353->349, complete families
  31->29, and the cap remains 216.  Physical entry, arbitrary k<=1, the rung,
  and LRC(14) remain OPEN.
audit: >
  The exact companion pins THM-3313 and THM-3308's compiler, reconstructs the
  480-row atlas and every earlier closure, reranks complete live families,
  verifies the strict queue boundary, and independently checks every inherited
  Farkas vector.  It derives affine Boolean majorants by rational vertex
  enumeration, constructs exact feasible empty-tail tables for all 32 distinct
  marginal systems, singleton-tail tables for all 18 proper supports of the
  three new support-two circuits, and all 21 proper supports of the row-138
  support-three control.  Row 64 retains a feasible residual.  Normal and
  optimized outputs byte-match; the source has zero assertion nodes and zero
  floating literals.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3308-threshold-chain-modular-multicovers-and-three-layer-status-circuit
  - THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure
related:
  - THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity
  - THM-3317-unique-protector-cells-and-weighted-scalar-fragility
script: 04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py
output: 05-knowledge/results/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.out
script_sha256: b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213
output_sha256: d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38
semantic_sha256: c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9
hash_basis: LF-normalized bytes
---

# THM-3320 -- projected-k3 z216 fourth ruler prefix and affine multicover closure

**PROVED + VERIFIED-EXACT in the declared projected atlas.**

This theorem both closes the fourth deterministic work prefix and separates
the arity of THM-3308's integral modular template from the arity of the full
common-table dual.

## 1. Complete-family prefix

Reconstruct the state after THM-3313: `353` live wall rows in `31` complete
families.  Rank a family by the inherited intrinsic invoice

```text
cost(g,L)=sum_(rows E in the family) L(E)r(E).              (1)
```

The queue begins

| rank | family | rows | atlas indices | cost |
|---:|---|---:|---|---:|
| 1 | `gcd24/L76440` | 1 | `141` | 2,293,200 |
| 2 | `gcd24/L30576` | 3 | `133,219,359` | 2,629,536 |
| 3 | `gcd72/L7056` | 19 | `21,49,62,...,465` | 3,400,992 |

Select ranks one and two, stopping after the next nonsingleton family.  The
boundary before rank three is strict and all three `L=30576` rows are kept.
Their labelled bodies and component counts are

```text
141: {1,4,6,10,13,14},       r=30;
133: {1,4,6,8,13,14},        r=26;
219: {1,6,8,12,13,14},       r=30;
359: {2,6,8,12,13,14},       r=30.                         (2)
```

Intrinsic cost is only a deterministic work order.

## 2. Empty exact screens

The inherited exact ray/common-status screen gives

| family | states | crude | status | residual |
|---|---:|---:|---:|---:|
| `gcd24/L76440` | 10 | 8 | 2 | 0 |
| `gcd24/L30576` | 220 | 164 | 56 | 0 |
| **total** | **230** | **172** | **58** | **0** |

Rowwise,

```text
141:  10=  8+ 2+0;        133: 10= 0+10+0;
219: 189=147+42+0;        359: 21=17+ 4+0.                 (3)
```

Every status dual is checked over exact arithmetic.  Hence each selected
projected wall row has an empty necessary upper screen and is excluded.

## 3. Integral-template arity is not full-dual arity

For a common status table with mass `q`, marginals `m_i`, subset capacity
`c(P)`, and tail demands `H_t`, THM-3308 uses pointwise majorants

```text
sum_(t in T) 1[c(P)>=t] <= sum_(i in P) w_i               (4)
```

with nonnegative integral `w`.  The `58` new status kills have template
supports

```text
55 of size one,       1 of size two,       2 of size three. (5)
```

The full equality-row dual also permits a rational affine constant:

```text
sum_(t in T) 1[c(P)>=t] <= a+sum_(i in P)w_i.              (6)
```

Summing `(6)` over one common table gives

```text
sum_(t in T)H_t <= aq+sum_i w_i m_i.                       (7)
```

Exact rational vertex enumeration of all sixteen Boolean patterns, followed
by primal feasibility checks on every smaller support, changes `(5)` to

```text
55 of size one,       3 of size two,       0 of size three. (8)
```

Thus the integral-template support is not an intrinsic full-dual circuit
complexity.

## 4. The three support-two circuits

Two packets occur at row `133`, with `q=2548`, histogram

```text
((1,400),(2,860),(3,326),(4,52),(5,754),(6,156)),          (9)
```

and marginal packets

```text
(84,588,728,1274),   m=(1092,1092,728,364),
(84,728,1092,1274),  m=(1092,728,1092,364).               (10)
```

For both, the exact Boolean inequality is

```text
2*(1[c(P)>=1]+1[c(P)>=5])
 <=1+1[0 in P]+1[1 in P]+1[2 in P]+3*1[3 in P].           (11)
```

Equivalently, `a=1/2`, `w=(1/2,1/2,1/2,3/2)`, `T=(1,5)`.
The demand is `3458`, the affine marginal cost is `3276`, and the strict gap
is `182`.  Their zero-constant integral presentations need three layers, so
the affine constant is the missing coordinate.

The third packet is

```text
row 133, divisors=(84,728,1274,2184),
q=2548, m=(1092,728,364,2184),                             (12)
```

with

```text
1[c(P)>=2]+1[c(P)>=5]
 <=1[0 in P]+1[1 in P]+2*1[2 in P].                       (13)
```

Its demand/cost are `3058/2548`, giving gap `510`.

For all three circuits, exact feasible tables exist for every singleton tail
at thresholds `1,...,6`, so support two is minimal.  Exact empty-tail tables
exist for all `32` distinct marginal systems among the `58` status addresses;
there is no hidden support-zero infeasibility.

The controls are sharp in both directions.  Row `64` retains eight residuals,
including an exact feasible common table.  The row-`138` packet at divisors
`(18,49,196,882)` still has the genuine support-three certificate
`T=(1,4,5)`, `w=(1,3,1,1)`, gap `11`; all six singleton and fifteen pair
systems are exactly feasible.

## 5. Ledger and scope

The four exclusions give

```text
projected ledger:       373157 -> 373153,
z1=216 wall rows:           353 -> 349,
complete live families:      31 -> 29,
projected k=3 cap:          216 unchanged.                 (14)
```

The next family is the nineteen-row `gcd72/L7056` fibre.  Capacity-signature
reuse is therefore the next native operation.

Everything here remains in the necessary projected `k=3,z1=216` atlas.  No
physical speed entry, endpoint origin, owner, phase, current, arbitrary
`k<=1`, rung, or LRC(14) conclusion is restored.

QED.
