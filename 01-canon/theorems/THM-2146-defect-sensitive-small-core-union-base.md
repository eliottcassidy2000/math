---
id: THM-2146
title: "Defect-sensitive exact small-core union base at the 3/41 threshold"
status: >
  FINITE-EXACT, with a proved all-height consequence. Over all 4096 subsets
  C of {1,...,13} of size at most six, two independent exact interval
  algorithms agree on the minimum 3/41-safe mass. The minimizers are unique
  from sizes two through six and are listed below. If a thirteen-speed row
  has defect d>=7 and seven noncore speeds are selected as a body, the
  remaining six-speed mask has safe measure at least B_d. For d=7,...,11,
  B_d improves the separately charged 5/41 base by exact factors
  153101/69300, 1649/900, 19/12, 427/330, and 11/10.
  This does not bound the aggregate covariance with the seven-speed body and
  therefore does not prove near-tight rigidity or LRC(14).
source: codex-2026-07-24-LRC-exact-union-base
depends_on: []
related:
  - THM-735
  - HYP-9024
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
script: 04-computation/lrc14_defect_sensitive_union_base_thm2146.py
output: 05-knowledge/results/lrc14_defect_sensitive_union_base_thm2146.out
---

# THM-2146 -- the exact small-core union base

Put

```text
h=3/41,
D_v={t in R/Z:||vt||<h},
G_h(C)=(R/Z) \ union_(v in C) D_v.                    (1)
```

Changing strict to closed danger bands changes no Haar measure.  We first
compute the exact minimum of `measure(G_h(C))` at every core cardinality.

## 1. Exact finite ledger

For every `0<=j<=6`,

| `j` | inputs | `min measure(G_h(C))` | equality sets |
|---:|---:|---:|---|
| 0 | 1 | `1` | `C=empty` |
| 1 | 13 | `35/41` | every singleton |
| 2 | 78 | `59/82` | `{1,12}` only |
| 3 | 286 | `1615/2706` | `{1,11,12}` only |
| 4 | 715 | `239/492` | `{1,7,8,9}` only |
| 5 | 1287 | `2729/7380` | `{1,5,7,8,9}` only |
| 6 | 1716 | `153101/568260` | `{1,5,7,8,9,11}` only |

The row counts sum to

```text
sum_(j=0)^6 binom(13,j)=4096.                         (2)
```

This is a complete finite statement, not a sample.

### Exact proof and independent path

For fixed `C`, every boundary of (1) is one of the rational points

```text
(k-h)/v or (k+h)/v mod 1,             v in C.         (3)
```

The first implementation sorts (3), takes the exact rational midpoint of
each complementary cell, evaluates every inequality in (1) with
`Fraction`, and sums precisely the safe cell lengths.  The second
implementation independently clips every danger interval to `[0,1]`, merges
the resulting rational intervals, and subtracts their exact total length
from one.  It also compares the longest safe component and the number of
positive-length components; isolated safe endpoints are measure-zero and
are not counted.  (No danger boundaries merely touch in this core universe,
since that would require `41|3(v+/-w)` with `v,w<=13`.)  Both paths agree on
all 4096 inputs; all comparisons and sums are rational.

Reproduce with

```text
python3 04-computation/lrc14_defect_sensitive_union_base_thm2146.py
```

The singleton row is the positive control:
`measure(G_h({v}))=1-2h=35/41` for every `v`.  The equality classifiers in
the last column are the exhaustive hostile boundary output, not chosen
representatives.

## 2. The all-height defect consequence

Let `V` be any set of thirteen distinct positive speeds and define its AP
defect by

```text
d=|V \ {1,...,13}|.                                  (4)
```

Assume `d>=7`.  Select any seven noncore speeds as a body `E`, and let
`F=V\E`.  Then `|F|=6`, and it splits as

```text
C=F intersection {1,...,13},       |C|=13-d,
R=F\C,                              |R|=d-7.          (5)
```

Every extra danger band has measure `2h=6/41`.  The ordinary union bound,
now applied only after retaining the exact core union, gives

```text
measure(G_h(F))
 >= measure(G_h(C))-(d-7)6/41
 >= B_d.                                               (6)
```

The exact values are

| `d` | `B_d` | factor over `5/41` |
|---:|---:|---:|
| 7 | `153101/568260` | `153101/69300` |
| 8 | `1649/7380` | `1649/900` |
| 9 | `95/492` | `19/12` |
| 10 | `427/2706` | `427/330` |
| 11 | `11/82` | `11/10` |
| 12 | `5/41` | `1` |
| 13 | `5/41` | `1` |

No height of a noncore speed occurs in (6); the consequence is uniform.

## 3. The sharpened remaining obligation

Write

```text
m_E=measure(G_h(E)),       m_F=measure(G_h(F)),
kappa=integral (1_(G_h(E))-m_E)(1_(G_h(F))-m_F).      (7)
```

Then the full lonely measure is exactly

```text
measure(G_h(V))=m_E m_F+kappa
              >=m_E B_d-|kappa|.                     (8)
```

Thus the flipped-peel route should no longer pay six danger bands
separately.  Its precise residual is the **one aggregate covariance** in
(7), with a defect-sensitive independent base `m_E B_d`.  At defects
`7,...,11`, the base multipliers above quantify the gain before any Fourier
or relation estimate is used.

Equation (8) also marks the theorem's boundary.  Safe masses alone do not
control the sign or size of `kappa`; an adversarial phase alignment could
make it negative.  THM-2146 neither proves `|kappa|<m_EB_d` nor classifies
the low-frequency relations that govern it.  It is a rigorous strengthening
of the live analytic invoice, not an emptying of the defect-`>=7` branch.

## 4. Perspective ledger

The source object is the six individual residual danger combs.  The target
is their exact union mask `1_(G_h(F))`; the map is Boolean multiplication of
the six safe indicators.  It preserves the predicate needed in (8) and the
full overlap cancellation that separate absolute discrepancies destroy.  It
forgets which runner owns a boundary.  The required sidecar for the next
step is therefore the labelled Fourier/endpoint decomposition of this union
mask, not merely its scalar measure.

There is no intrinsic pair orientation here: forcing a tournament would
discard the higher overlap data.  The faithful finite carrier is the subset
lattice of the six residual bands together with the aggregate covariance.
The independently reserved, still-unproved THM-2145 two-block spectral
crossing proposal targets precisely the missing labelled `6+7` relation
sidecar; it is a compatible next question, not a dependency of this theorem.
The hostile controls are the unique equality cores in the table and the
unrestricted extra speeds in (5).  QED.
