# LRC(14) `z1=216`: third intrinsic-ruler prefix and multicover closure

**Status:** structural synthesis for
[THM-3313](../01-canon/theorems/THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md),
with projected-atlas scope only.  The next deterministic intrinsic-cost prefix after the two independently audited
prefixes consists of the `gcd72/L72072` singleton and the complete three-row
`gcd24/L25872` family.  All four labelled projected `k=3`, `z1=216` wall rows
have empty exact necessary upper screens.  The implied maintained-ledger
change is

```text
373161 -> 373157,
z1=216 wall rows: 357 -> 353,
complete ruler families: 33 -> 31,
projected cap: 216 (unchanged).                              (1)
```

No maintained ledger or navigation file is changed by this scout.  The exact
companions are

- [`lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py`](../04-computation/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py);
- its [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.out).

This remains a necessary projected quotient.  It does not classify physical
covers and does not restore speed entry, endpoint origin, owners, phase,
current, arbitrary `k<=1`, the rung, or LRC(14).

## Inheritance pass

**Closest proved mechanism.**  THM-3139 supplies the exact ray/status upper
screen, while THM-3281 establishes complete `(gcd(216,L),L)` fibres as lawful
units of projected work.  The [independent audit of the first two
prefixes](lrc14-z216-intrinsic-ruler-cost-prefix-independent-audit-codex-20260803.md)
reconstructs their `23` closures without importing either original scout.

**Canonical hostile.**  Prior row `64` has screen census

```text
115 states = 50 crude + 57 status + 8 residual.              (2)
```

The residuals prevent any convention under which every four-bit status packet
is declared contradictory.

**Corrected near miss.**  MISTAKE-331 and MISTAKE-333 forbid treating a
solver-chosen Farkas basis or its scale as canonical.  The new computation
checks all inherited duals exactly but persists only the deterministic status
instance and the lexicographically selected modular multicover.

**Least-used relevant sidecar.**  The [threshold-layer multicover
compiler](lrc14-z216-threshold-layer-multicover-and-support-minimal-status-circuits-codex-20260803.md)
turns several load-tail constraints on one common status table into a small
Boolean circuit.  This prefix is the first new direct ledger test carried out
with that presentation from the outset.

## Exact queue reconstruction

The computation reconstructs the older disjoint closures, derives the first
and second cost prefixes again, removes their `23` rows, and only then ranks
the `357` live wall rows in `33` complete ruler families.  It does not trust a
prose address.  The exact queue head is

| rank | family | rows | labelled atlas indices | intrinsic invoice `sum Lr` |
|---:|---|---:|---|---:|
| 1 | `gcd72/L72072` | 1 | `157` | `1,729,728` |
| 2 | `gcd24/L25872` | 3 | `53,293,357` | `2,121,504` |
| 3 | `gcd24/L76440` | 1 | `141` | `2,293,200` |

The selected prefix ends strictly before rank `3`; no cost tie crosses its
boundary.  Its four labelled bodies and component counts are

```text
157: {1,4,9,11,12,13},       r=24;
 53: {1,2,6,8,11,14},        r=26;
293: {2,4,6,8,11,14},        r=26;
357: {2,6,8,11,12,14},       r=30.                           (3)
```

The `L=25872` fibre is not sampled: all three rows are taken.  Its bodies are
the common five-set `{2,6,8,11,14}` plus respectively one of `{1,4,12}`.

## Exact screen result

| family | rows | ray states | crude kills | exact status kills | residual |
|---|---:|---:|---:|---:|---:|
| `gcd72/L72072` | 1 | 111 | 93 | 18 | 0 |
| `gcd24/L25872` | 3 | 312 | 196 | 116 | 0 |
| **total** | **4** | **423** | **289** | **134** | **0** |

The rowwise census is

```text
157: 111 = 93 + 18 + 0;
 53:  27 = 17 + 10 + 0;
293:  14 =  9 +  5 + 0;
357: 271 = 170 + 101 + 0.                                  (4)
```

Every one of the `134` returned common-status Farkas certificates is checked
over exact arithmetic.  No residual passport remains to send to a terminal
probe.

## All status kills have modular multicover proofs

For a deterministic four-bit status packet, let `x_P` be the common mass at
pattern `P`, with total `q`, marginals `m_i`, pattern capacity `c(P)`, and
required load tail `H_t`.  For selected thresholds `T`, the compiler searches
for nonnegative integral coordinate weights satisfying the pointwise Boolean
inequality

```text
sum_(t in T) 1[c(P)>=t] <= sum_(i in P) w_i.                 (5)
```

Summing `(5)` against the **one common table** gives

```text
sum_(t in T) H_t <= sum_i w_i m_i.                           (6)
```

A strict reverse inequality is an exact contradiction.  All `134` new status
kills admit this deterministic presentation:

```text
133 one-layer certificates,
  1 two-layer certificate,
  0 uncompiled status kills.                                (7)
```

The older one-threshold taxonomy refines to

| older branch | one layer | two layers |
|---|---:|---:|
| coordinate union | 128 | 0 |
| zero-marginal-reduced union | 2 | 0 |
| two-fan | 3 | 0 |
| weighted core | 0 | 1 |

Thus the sole old “weighted core” is exactly a cross-threshold compatibility
circuit rather than an unstructured LP remainder.

## The new two-layer circuit

The unique multi-layer address is

```text
atlas row 157,
divisors = (99,616,1001,8008),
q = 9009,
m = (1365,1287,1287,1287),
histogram = ((0,99),(1,1338),(2,2382),(3,1992),(4,2874),(5,324)).  (8)
```

Its capacity vector in pattern order `0,...,15` is

```text
(2,8,3,8,8,8,8,8,3,8,4,8,8,8,8,8).                       (9)
```

At thresholds `3` and `4`, the primitive pointwise majorant is

```text
1[c(P)>=3] + 1[c(P)>=4]
 <= 2*1[0 in P] + 1[1 in P] + 2*1[2 in P] + 1[3 in P].     (10)
```

The two required tails have total

```text
H_3 + H_4 = 5190 + 3198 = 8388,
```

whereas the marginal side of `(10)` is

```text
2*1365 + 1287 + 2*1287 + 1287 = 7878.                      (11)
```

The contradiction gap is `510`.  The tight Boolean patterns are
`(0,1,2,4,8,10)`.  All three singleton tail systems at the nontrivial
thresholds `3,4,5` have independently constructed exact feasible common
tables.  Hence support two is minimum in the full common-table dual cone, not
merely minimum among integral modular covers.

## Hostile controls and validity gate

The scout replays both sides of the compiler boundary:

1. All `57` status kills of hostile row `64` have one-layer covers, but its
   first residual retains the previously frozen exact common table and has no
   modular-cover contradiction.
2. The earlier row-`138`, divisors `(18,49,196,882)` packet returns its genuine
   three-layer certificate `T=(1,4,5)`, `w=(1,3,1,1)`, gap `11`.  Thus the
   search is not hard-wired to stop at arity two.
3. Family membership, strict cost boundary, row labels, ruler/high fields,
   marginals, capacity sets, cofactor identities, every raw dual, and every
   pointwise modular inequality are checked exactly.
4. The script has no `assert` truth gates and no floating-point literals.
   Normal and optimized transcripts are byte-identical.

The logical direction remains only

```text
empty exact necessary upper screen
    => selected projected wall row is empty.                (12)
```

Screen feasibility is not used as a converse.  Intrinsic cost is a work
order, not a safety invariant.

## Stopping boundary

After the four exclusions implied by `(1)`, the projected `z1=216` wall has
`353` rows in `31` complete ruler families.  The next exact queue head is

```text
gcd24/L76440 singleton,                  invoice 2,293,200;
gcd24/L30576 three-row family,           invoice 2,629,536;
gcd72/L7056 nineteen-row family,          invoice 3,400,992. (13)
```

The first two entries form the next deterministic prefix through a
non-singleton family.  The third is the next substantially larger fibre and
is a natural boundary for capacity-signature reuse rather than raw row-by-row
screening.

## Reproduction and hashes

```text
python3 04-computation/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py
python3 -O 04-computation/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py
```

```text
script    34710b907f3274057d6cb60a0fd7bd72ae2c3879a90b830a3e376d30e5198646
output    c664acf4e02b40835dfdc76ff135d17332e3d0e57b489e61b19c0f2f7b5998f5
semantic  4d8547f53b287d82bb725cd6b3ed78c390ba03a37a3dcb3ab76ac98389e045a3
```
