# LRC(14) `z1=216`: fourth ruler prefix and affine multicover closure

**Status: PROVED + VERIFIED-EXACT DIRECT PROGRESS, PROJECTED-ATLAS SCOPE
ONLY; canonical statement
[THM-3320](../01-canon/theorems/THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure.md).**
After THM-3313, the next deterministic intrinsic-cost prefix is the
`gcd24/L76440` singleton followed by the complete three-row
`gcd24/L30576` family.  All four labelled projected `k=3,z1=216` necessary
upper screens are empty.  The implied ledger change is

```text
373157 -> 373153,
z1=216 wall rows: 353 -> 349,
complete ruler families: 31 -> 29,
projected cap: 216 (unchanged).                              (1)
```

The computation also finds a sharp boundary of the threshold compiler: two
packets which need three layers in the zero-constant nonnegative-integral
modular template already have support-two certificates in the full
common-table dual.  The missing coordinate is an affine constant together
with half-integral weights.

Companions:

- [`lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py`](../04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py);
- its [frozen exact output](../05-knowledge/results/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.out).

The result remains inside a necessary projected quotient and restores no
physical speed entry, endpoint origin, owner, phase, current, arbitrary
`k<=1`, rung, or LRC(14) conclusion.

## Inheritance pass

**Closest proved mechanism.**  [THM-3313 -- projected-k3 z216 third
ruler-cost prefix multicover closure](../01-canon/theorems/THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md)
leaves exactly `353` wall rows in `31` complete ruler families.  THM-3308
supplies the deterministic threshold-layer compiler and the exact
proper-support feasibility machinery.

**Canonical hostile.**  Row `64` has

```text
115 states = 50 crude + 57 status + 8 residual.              (2)
```

Its first residual admits an explicit exact common table, so neither the old
modular compiler nor the affine extension may kill every packet by
convention.

**Corrected near miss.**  MISTAKE-331 and MISTAKE-333 forbid persisting a
solver-selected Farkas basis or its normalization.  Every inherited dual is
checked exactly, but the new frozen certificate is derived instead from the
deterministic capacity table by rational affine-majorant enumeration.

**Least-used relevant sidecar.**  THM-3308 already separated modular-template
support from full-dual support by constructing exact feasible tables for all
smaller tail systems.  Applying that validity gate to the new packets exposes
the half-affine repair below.

## Concept board

| object | representation | preserved predicate | loss / cheapest test |
|---|---|---|---|
| projected wall row | labelled body, ruler and high gate | membership in the maintained necessary atlas | physical entry, owners and current are absent; run the complete screen |
| ruler family | complete fibre of `(gcd(216,L),L)` | every labelled family member remains present | family cost is not a safety statistic; screen all members |
| modular circuit | zero-constant integral coordinate majorant of nested tails | exact contradiction for one common table | may overstate full-dual arity; test every smaller support by exact primal feasibility |
| affine circuit | rational constant plus coordinate weights | full equality-row dual of the common-table LP | still loses physical LRC coordinates; retain the quotient-loss ledger |
| residual hostile | exact feasible status table | prevents false universal elimination | says nothing about later terminal probes |

## Canonical queue reconstruction

The scout pins THM-3313 and its source/output, reconstructs the original
`480`-row atlas, removes every older disjoint closure, and independently
regenerates the first three cost prefixes.  Only then does it rank the `353`
live rows in `31` complete families.  The head is

| rank | family | rows | atlas indices | invoice `sum Lr` |
|---:|---|---:|---|---:|
| 1 | `gcd24/L76440` | 1 | `141` | `2,293,200` |
| 2 | `gcd24/L30576` | 3 | `133,219,359` | `2,629,536` |
| 3 | `gcd72/L7056` | 19 | `21,49,62,...,465` | `3,400,992` |

The selected prefix contains ranks `1--2`; its boundary before the nineteen-row
family is strict.  Its labelled bodies and component counts are

```text
141: {1,4,6,10,13,14},       r=30;
133: {1,4,6,8,13,14},        r=26;
219: {1,6,8,12,13,14},       r=30;
359: {2,6,8,12,13,14},       r=30.                           (3)
```

All three `L=30576` rows are retained.  No representative-row inference is
used.

## Exact screen closure

| family | states | crude kills | exact status kills | residual |
|---|---:|---:|---:|---:|
| `gcd24/L76440` | 10 | 8 | 2 | 0 |
| `gcd24/L30576` | 220 | 164 | 56 | 0 |
| **total** | **230** | **172** | **58** | **0** |

Rowwise,

```text
141:  10 =   8 +  2 + 0;
133:  10 =   0 + 10 + 0;
219: 189 = 147 + 42 + 0;
359:  21 =  17 +  4 + 0.                                  (4)
```

All `58` status duals pass exact independent checks.  There are no residual
passports to classify further.

## Integral modular arity versus full-dual arity

For status mass `x_P`, total `q`, marginals `m_i`, capacity `c(P)`, and tail
demand `H_t`, THM-3308 searches for

```text
sum_(t in T) 1[c(P)>=t] <= sum_(i in P) w_i,                (5)
```

with nonnegative integral `w`.  On the `58` new status kills, this gives

```text
55 support-one,
 1 support-two,
 2 support-three                                        (integral template). (6)
```

The full equality-row dual permits an affine majorant

```text
sum_(t in T) 1[c(P)>=t] <= a + sum_(i in P) w_i.            (7)
```

Summing `(7)` against one common table gives

```text
sum_(t in T) H_t <= a*q + sum_i w_i m_i.                    (8)
```

Rational vertex enumeration of the sixteen Boolean constraints, followed by
exact feasible-table checks on every smaller support, repairs `(6)` to

```text
55 support-one,
 3 support-two,
 0 support-three                                         (full dual).        (9)
```

Thus integral modular arity is an intrinsic complexity statistic for that
template, but it is not automatically the circuit arity of the full
common-table LP.

## The two half-affine circuits

Both new half-affine packets occur at row `133`, have `q=2548`, and share

```text
histogram=((1,400),(2,860),(3,326),(4,52),(5,754),(6,156)). (10)
```

Their divisor/marginal packets are

```text
(84,588,728,1274),   m=(1092,1092,728,364);
(84,728,1092,1274),  m=(1092,728,1092,364).                 (11)
```

For each capacity table, the pointwise Boolean inequality is

```text
2*(1[c(P)>=1] + 1[c(P)>=5])
 <= 1 + 1[0 in P] + 1[1 in P] + 1[2 in P] + 3*1[3 in P].   (12)
```

Equivalently, `(7)` has

```text
a=1/2,       w=(1/2,1/2,1/2,3/2),       T=(1,5).           (13)
```

The demand is

```text
H_1+H_5 = 2548+910 = 3458,
```

while the affine marginal cost is `3276`.  The unscaled strict gap is `182`
for each packet.  The tight Boolean patterns are `(1,2,4,7,8)`.

The zero-constant integral compiler instead needs `T=(1,3,5)` and reports
gap `14`; that certificate remains valid, but it is not support-minimal in the
larger affine dual.

## The integral support-two circuit

The third multi-layer address is

```text
row 133,
divisors=(84,728,1274,2184),
q=2548,
m=(1092,728,364,2184),                                     (14)
```

with the same histogram `(10)`.  Its zero-constant pointwise certificate is

```text
1[c(P)>=2] + 1[c(P)>=5]
 <= 1[0 in P] + 1[1 in P] + 2*1[2 in P].                   (15)
```

Here

```text
H_2+H_5 = 2148+910 = 3058,
m_0+m_1+2m_2 = 2548,                                       (16)
```

so the gap is `510`.

For each of the three circuits, all six singleton tail systems at thresholds
`1,...,6` have explicit exact feasible common tables.  These `18` positive
controls prove support two is minimal in the full dual, not merely in the
affine-majorant search.  Separately, all `32` distinct marginal systems among
the `58` status addresses have exact feasible empty-tail tables.  Hence the
reported support-one circuits are not hiding an infeasible support-zero
marginal system.

## Hostile controls and boundary of the repair

Two controls protect opposite directions:

1. Row `64` still gives `(2)`.  All `57` status kills compile at one layer,
   while its first residual retains an exact feasible common table and no
   modular contradiction.
2. THM-3308's row-`138`, divisors `(18,49,196,882)` packet still needs the
   genuine three-layer certificate `T=(1,4,5)`, `w=(1,3,1,1)`, gap `11`.
   Exact feasible tables are reconstructed for all `6` singletons and all
   `15` pairs.  Hence allowing affine constants does not collapse every
   support-three circuit.

The script additionally verifies family completeness, strict queue boundary,
body/ruler/high fields, cofactor identities, marginals, capacity sets, all raw
duals, every pointwise affine inequality, and every proper-support table.  It
contains no `assert` truth gate and no floating-point literal.  Normal and
optimized transcripts are byte-identical.

## Direction, loss ledger, and stopping boundary

The only row-level implication is

```text
empty exact necessary upper screen
    => selected projected wall row is empty.                (17)
```

No converse from screen feasibility is used.  Intrinsic cost remains a work
order, not a safety invariant.  The source-to-target contract is

```text
source: labelled projected wall row
map: row -> ray states -> common status table -> affine tail circuit
kept: inherited necessary filters and exact infeasibility
lost: physical entry, endpoint, owner, phase, current and chronology
needed sidecar: a lawful physical lift retaining those coordinates
cheapest decisive test: complete exact screen plus hostile feasible table.
```

After `(1)`, `349` wall rows in `29` complete families remain.  The next queue
head is

```text
gcd72/L7056 nineteen-row family,  cost 3,400,992;
gcd72/L144144 singleton,          cost 3,747,744;
gcd24/L64680 two-row family,      cost 3,880,800.             (18)
```

Because the next prefix begins with a nineteen-row family, the natural next
operation is capacity-signature reuse and family-level batching rather than
another unstructured row loop.

## Reproduction and hashes

```text
python3 04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py
python3 -O 04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py
```

```text
script    b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213
output    d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38
semantic  c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9
```
