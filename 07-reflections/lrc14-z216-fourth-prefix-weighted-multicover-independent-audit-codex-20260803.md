# LRC(14) `z1=216`: fourth-prefix weighted-multicover independent audit

**Status: INDEPENDENTLY VERIFIED-EXACT FOR THE PREFIX, SCREEN, LEDGER, THREE
NONTRIVIAL SUPPORT-TWO MINIMA, AND CONTROLS IN
[THM-3320](../01-canon/theorems/THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure.md),
PROJECTED-ATLAS SCOPE ONLY.**  This computation was derived concurrently from
THM-3308/3313 and does not import THM-3320.  It exactly reproduces THM-3320's
four selected rows, empty screens, ledger decrement, all three nontrivial
full-dual support-two minima, and hostile controls.  Its additional value is a zero-constant
threshold-**multiset** presentation of the two half-affine circuits.

THM-3320 is the canonical proof source and is stronger: it enumerates the full
affine Boolean-majorant vertices and constructs empty-tail feasible tables for
all `32` distinct status instances.  This audit found no defect in its
statement or proof packet and does not warrant a second theorem.

Here “independent” means independent of THM-3320's source and frozen output.
Both computations deliberately share the older canonical ray/status engine
and the THM-3308/3313 atlas lineage; this is not a second implementation of
that inherited engine.

The independent companions are

- [`lrc14_j7_k3_z216_fourth_intrinsic_ruler_prefix_weighted_multicover_probe_20260803.py`](../04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_prefix_weighted_multicover_probe_20260803.py);
- its [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_fourth_intrinsic_ruler_prefix_weighted_multicover_probe_20260803.out).

## Inheritance and typed operation

The closest proved mechanisms are THM-3308's exact common-table multicover and
THM-3313's complete-family cost queue.  Row `64`, with eight residual masks,
is the canonical nonempty-screen hostile.  The least-used sidecar is the
distinction between the **support** of a tail dual and the **multiplicity** of
its selected layers.

The operation is

```text
labelled projected row
 -> exact ray/common-status screen
 -> nested tail indicators
 -> nonnegative integral threshold multiset
 -> zero-constant modular coordinate majorant.
```

It preserves infeasibility of the maintained necessary relaxation.  It loses
physical speed entry, endpoint origin, owner, phase, current, chronology,
arbitrary `k<=1`, and the rung.

## Exact independent replay

Starting from THM-3313's `353` live wall rows in `31` complete
`(gcd(216,L),L)` families, the independently reconstructed queue begins

| rank | family | rows | indices | cost |
|---:|---|---:|---|---:|
| 1 | `gcd24/L76440` | 1 | `141` | 2,293,200 |
| 2 | `gcd24/L30576` | 3 | `133,219,359` | 2,629,536 |
| 3 | `gcd72/L7056` | 19 | `21,49,...,465` | 3,400,992 |

The complete prefix through the first nonsingleton family is therefore rows
`141,133,219,359`.  Its exact screen is

| row | states | crude | status | residual |
|---:|---:|---:|---:|---:|
| 141 | 10 | 8 | 2 | 0 |
| 133 | 10 | 0 | 10 | 0 |
| 219 | 189 | 147 | 42 | 0 |
| 359 | 21 | 17 | 4 | 0 |
| **total** | **230** | **172** | **58** | **0** |

The selected-packet and canonical-screen hashes are respectively
`b8c63cd3...e27d` and `ea3f2619...387ce`, exactly matching THM-3320.  The
one-way consequence is

```text
projected ledger:       373157 -> 373153,
z1=216 wall rows:           353 -> 349,
complete live families:      31 -> 29,
projected cap:               216 unchanged.
```

No converse from screen feasibility is used, and intrinsic cost is only a
work order.

## Set support versus multiset layer mass

For tail demands `H_t`, status capacities `c(P)`, and marginals `m_i`, the
audit searches the bounded integral language

```text
sum_t alpha_t 1[c(P)>=t] <= sum_(i in P) w_i,               (1)
alpha_t in Z_>=0,  w_i in Z_>=0,  sum_t alpha_t <= 3.
```

Every pointwise inequality is checked on all sixteen patterns.  The `58`
status states have

```text
tail support:       55 of size 1,   3 of size 2;
total layer mass:   55 of mass 1,   1 of mass 2,   2 of mass 3.
```

All three support-two states occur at row `133`.  One is THM-3320's integral
binary circuit

```text
1[c(P)>=2] + 1[c(P)>=5]
 <= 1[0 in P] + 1[1 in P] + 2*1[2 in P],                  (2)
```

with gap `510`.  The other two use the threshold multiset

```text
1[c(P)>=1] + 2*1[c(P)>=5]
 <= 1[0 in P] + 1[1 in P] + 1[2 in P] + 3*1[3 in P],      (3)
```

with gap `364`.  Their divisor packets are

```text
(84,588,728,1274),
(84,728,1092,1274).
```

For each of the three addresses, all six singleton-tail systems have explicit
exact feasible common tables.  Thus their minimum tail support is genuinely
two in the full rational dual, while the first two need total integral layer
mass three in the zero-constant multiset presentation.

## Gauge relation to THM-3320's half-affine form

THM-3320 presents the first two circuits as

```text
2*(1[c(P)>=1] + 1[c(P)>=5])
 <= 1 + 1[0 in P] + 1[1 in P] + 1[2 in P] + 3*1[3 in P].  (4)
```

Equations `(3)` and `(4)` are different pointwise representatives of the same
dual obstruction on the feasible common-table face.  For these two capacity
tables,

```text
1[c(P)>=1] = 1[P is nonempty],     H_1=q.                  (5)
```

Hence every feasible common table has zero mass at the empty pattern.  The
slack in `(4)` is the slack in `(3)` plus
`1-1[c(P)>=1]`, which is supported only at that empty pattern.  After summing
against a feasible table the extra term vanishes.  Both integer-scaled
certificates therefore have gap `364`; THM-3320 divides `(4)` by two and
records the half-affine gap `182`.

This is the useful conceptual addition: the affine constant and a repeated
saturated tail are **gauge-equivalent only after the tail-one constraint has
forced the empty status to vanish**.  They are not the same pointwise Boolean
function on the ambient cube.  Tail support, layer multiplicity, and affine
constant are three distinct coordinates of the certificate passport.

## Controls and failure boundary

- Row `64` replays as `115=50+57+8`; its first residual has an explicit exact
  common table and no weighted multicover through total mass three.
- The inherited row-`138` packet at divisors `(18,49,196,882)` remains genuinely
  support three: `(1,4,5)` with weights `(1,3,1,1)` and gap `11`.  All six
  singleton and fifteen pair systems are exactly feasible.

Thus repeated thresholds do not collapse every three-tail circuit.  The
bounded multiset compiler is sufficient for these `58` states, not a universal
four-bit theorem.

After the four closures, `349` wall rows in `29` families remain.  The next
complete family is the nineteen-row `gcd72/L7056` fibre.  The natural projected
operation is capacity-signature batching; physical progress still requires a
lawful lift restoring endpoint/current rather than another scalar passport.

## Reproduction and hashes

```text
python 04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_prefix_weighted_multicover_probe_20260803.py --processes 4
python -O 04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_prefix_weighted_multicover_probe_20260803.py --processes 4
```

Normal and optimized outputs are byte-identical.  The script contains no
`assert` truth gate and no floating-point literal.  LF-normalized hashes are

```text
script    f69e7e6991caedd572f6438034f3866074bacba22851e547085190464a1c98f0
output    5ed39fd402f76c2d9017af80f21c8d51cb38e617224336eaf762be8c51de5a53
semantic  666458d6122a66a7615a2855a1b0dcf18fb778e1451b44c064bab04726d81b05
```
