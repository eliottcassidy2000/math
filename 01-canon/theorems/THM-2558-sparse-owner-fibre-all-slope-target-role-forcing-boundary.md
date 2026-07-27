---
id: THM-2558
title: "All-slope head visibility and the 202-necklace target-role boundary"
status: >
  PROVED + VERIFIED-EXACT.  For the THM-2531 lexicographic boundary
  selector on a nonconstant thirteen-root mask e, let H(e) be the set of
  selected empty heads obtained from all twelve nonzero slopes.  Exactly
  5,564 of 8,190 labelled masks, or 428 of 630 rotation necklaces, satisfy
  H(e)=e^c.  Every mask of weight 1, 2, or 8--12 does; the 2,626 labelled
  failures are exactly 202 middle-layer necklaces of weights 3--7.  On a
  live depth-one owner fibre the selector mask is the six-role safe mask
  A_0, not the one/two-root mask seen only after Perron rebasing.  Whenever
  its all-empty-visible stratum has positive mass, some slope therefore
  places the selected head on the unique target-active guard/unit failure
  k_a.
  A three-root control hides a singleton target-active failure from every
  slope, so the remaining 202 necklaces cannot be removed by all-slope
  selection or capacity alone.  The forced hit is at the old predecessor
  horizon; no genuinely later target-root map, ancestry/future-digit
  intertwiner, scalar-row exclusion, or LRC(14) conclusion follows.
source: codex-2026-07-27-sparse-owner-head-forcing
depends_on:
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
related:
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2554-translation-quotient-root-displacement-and-endpoint-swap-parity
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
script: 04-computation/lrc14_sparse_owner_head_visibility_thm2558.py
output: 05-knowledge/results/lrc14_sparse_owner_head_visibility_thm2558.out
script_sha256: ae2f57ea26c42c5388e3db897ae89b15449ff36aec649f17e889a0b94f7ed318
output_sha256: ff08c60fb2d0ec446549601ed94981879be880f0428adfae3b7d456a66d4aa74
hash_basis: working-tree bytes (LF)
---

# THM-2558 -- all-slope head visibility has a 202-necklace blind layer

**PROVED + VERIFIED-EXACT.**

The inheritance pass exposes a useful near miss.

```text
closest positive mechanism:  THM-2537's selected occupied-to-empty wall;
closest categorical role:    THM-2461's unique target-active guard/unit k_a;
tempting sparse input:        THM-2305's one/two-root owner packet;
missing coordinate:          that sparsity is in the Perron-rebased owner
                             root, not the selector's predecessor root.
                                                                    (1)
```

The corrected statement is still strong.  All twelve word-tournament
selectors see every empty root on 428 of the 630 prime rotation necklaces.
On those strata the target-active guard/unit failure cannot hide.  The remaining
202 necklaces form an exact finite local obstruction, rather than an
unstructured all-mask problem.

## 1. The two root charts must not be conflated

Use THM-2296/2305 notation.  On a live first-depth-one row,

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

E_j=A_0 intersection D_(c_j)
          minus union_(h!=j)D_(c_h),                         (2)

c_h=13^(lambda_h)u_h,                 lambda_h>=1.
```

The THM-2531/2537 selected wall is read on the direct predecessor fibre

```text
x=iota_r(z)=(z+r)/13,                    r in F_13.          (3)
```

Every blocker status in (2) is constant in `r`, because

```text
c_h iota_r(z)
 =13^(lambda_h-1)u_h z+13^(lambda_h-1)u_h r
 =13^(lambda_h-1)u_h z                       mod 1.          (4)
```

The terminal word at every positive clock is also root-constant on (3).
Consequently, on a fibre where the selected owner-word event is nonempty,
the varying Boolean root mask is exactly

```text
e_r(z)=1_(A_0)(iota_r(z)).                                (5)
```

It is the simultaneous safe mask of the guard and five unit roles.  It need
not have one or two roots.

The one/two-root statement lives in a different chart.  THM-2305 first puts

```text
G_j=P^(lambda_j)1_(E_j),
```

whose support lies in the unit comb `D_(u_j)`, and only then takes the last
inverse branch of `G_j`.  Its label is the owner-ancestry root `r_owner` of
THM-2461, not the direct predecessor root `r` in (3).  Replacing (5) by that
sparse mask without retaining the Perron ancestry would identify two typed
roots which current canon explicitly separates.

## 2. The all-slope selected-head image

Let `R=F_13`, let `e:R->{0,1}` be nonconstant, and fix `tau!=0`.  Read from
each root `a` the cyclic binary word

```text
W_(tau,a)(e)=(e_(a+j tau))_(j=0)^12.                       (6)
```

THM-2531 proves that these thirteen words are distinct.  Let `alpha_tau(e)`
be the root with lexicographically maximal word, let `q_tau(e)` be the length
of its initial run of ones, and put

```text
t_tau(e)=alpha_tau(e)+q_tau(e)tau.                          (7)
```

Then `e_(t_tau(e))=0`.  Define the all-slope head image

```text
H(e)={t_tau(e):tau in F_13^*} subset e^c.                   (8)
```

Call `e` **all-empty-visible** when

```text
H(e)=e^c.                                                   (9)
```

The property is invariant under root translation, and every nonconstant
mask has a free translation orbit of size thirteen.

Two boundary laws are especially transparent.

- For a singleton, the twelve slopes select the twelve empty roots once
  each.
- For a pair, affine covariance reduces to `e={0,1}`.  Its selected heads
  are `2,3,...,12`; the midpoint `7=1/2 mod 13` occurs twice and every other
  empty root once.

Thus every mask of weight one or two is all-empty-visible.

## 3. Exact census: 428 visible and 202 blind necklaces

Exhausting the finite universe of `2^13-2=8,190` nonconstant masks gives the
following complete head-diversity census.  An entry `d:n` means that `n`
labelled masks of that weight have `|H(e)|=d`; the last column counts (9).

| `|e|` | head-diversity histogram | all-empty-visible |
|---:|---:|---:|
| 1 | `12:13` | 13 |
| 2 | `11:78` | 78 |
| 3 | `6:52, 8:156, 10:78` | 78 |
| 4 | `6:52, 7:390, 8:156, 9:117` | 117 |
| 5 | `6:468, 7:468, 8:351` | 351 |
| 6 | `5:78, 6:650, 7:988` | 988 |
| 7 | `5:156, 6:1560` | 1,560 |
| 8 | `5:1287` | 1,287 |
| 9 | `4:715` | 715 |
| 10 | `3:286` | 286 |
| 11 | `2:78` | 78 |
| 12 | `1:13` | 13 |

Therefore

```text
all-empty-visible labelled masks:       5,564;
blind labelled masks:                   2,626;

all-empty-visible rotation necklaces:     428;
blind rotation necklaces:                 202.              (10)
```

Every mask of weight `1,2,8,9,10,11,12` satisfies (9).  Failure occurs only
in the middle weights `3,...,7`.  The exact referee evaluates the definition
(6)--(8) for every mask and every slope; (10) is a finite theorem, not a
sampled extrapolation.

## 4. What a visible stratum forces on a live owner fibre

Choose the THM-2350/2461 pivot with `q_*=k_b`.  Among the five remaining
guard/unit first-failure roles, `k_a` is the unique target-active role.  Its
coefficient is a thirteen-unit.  On every direct fibre (3), its failure set

```text
T_a(z)={r:iota_r(z) in D_(k_a)}                             (11)
```

has one or two roots when `k_a` is ordinary and three or four when it is the
guard; in either case it is nonempty.  Since (5) is safe for all six
guard/unit roles,

```text
T_a(z) subset e(z)^c.                                      (12)
```

If `e(z)` is all-empty-visible, equations (9) and (12) give a slope `tau`
with

```text
t_tau(e(z)) in T_a(z).                                     (13)
```

Thus the selected empty head at that slope carries the old target-active
unit failure `k_a`.  This conclusion is pointwise and needs no Fourier
nonvanishing argument.

There is also a positive-measure form.  Let `g(z)` be any common nonnegative
root-invariant owner/word weight used by THM-2537, and let `B_good` be the
set of active bases on which (9) holds.  Since every selected cut score
`Psi_tau(e)` is a positive integer, (13) gives

```text
sum_(tau!=0) integral_(B_good)
 g(z)Psi_tau(e(z))
 1_(t_tau(e(z)) in T_a(z)) dz
 >=integral_(B_good)g(z)dz.                                (14)
```

If the right side is positive, at least one fixed slope has a positive
old-target-active selected-head packet.  This closes the purely local
“which one of the five failure roles?” question on the visible stratum.

The four neutral roles give a compatible capacity check.  An ordinary unit
failure occupies at most two predecessor roots; the guard failure occupies
at most four.  Hence their union has size at most

```text
4+2+2+2=10.                                                (15)
```

The singleton and pair laws, with respectively twelve and eleven distinct
heads, also force a head outside that neutral union.  Equation (15) does not
extend the conclusion to arbitrary masks because their head diversity can
drop to five or six.

## 5. The middle layer is a sharp obstruction

Take

```text
e={0,1,4}.                                                  (16)
```

Direct evaluation of the twelve word tournaments gives

```text
H(e)={2,7,8,9,11,12},
e^c minus H(e)={3,5,6,10}.                                 (17)
```

A singleton target-active ordinary-unit failure supported at root `3` is disjoint
from `e` and obeys the sharp one/two-root unit capacity, yet every selected
head misses it.  This is a finite support-level hostile: it is not asserted
to be the root fibre of one of the 165 scalar rows.  It proves that the
all-slope selector, mask support, and unit-capacity axioms alone cannot
remove the 202 blind necklaces.

The failed shortcut is now localized exactly:

```text
Perron-rebased owner mask is sparse
  -/-> direct selected-wall mask is sparse;

all twelve slopes are present
  -/-> every empty root is selected on 202 necklaces;

old target-active predecessor hit
  -/-> genuinely later target-root field.                    (18)
```

## 6. Remaining LRC object

The local role-forcing problem has split into two finite/semantic pieces.

1. Prove positive selected mass on one of the 428 visible necklaces, or use
   the physical six-comb incidence constraints to eliminate the 202 blind
   necklaces.
2. On the resulting positive `k_a` packet, place that role genuinely later
   while retaining its ancestry sheet and carry, and prove that its semantic
   later root is the THM-2549 old-sheet root (or construct the future-scale
   intertwiner directly).

Only after the second step does THM-2549 equation (15e) make the Hall graph
diagonal.  The present theorem supplies an old predecessor role on the first
step's visible stratum; it does not define THM-2545's later root map `b`, does
not remove cemetery mass, and does not establish THM-2537 equation (56) in
its later-field sense.  No scalar row is removed.  The LRC(14) ledger remains
`165`.

## Exact referee

Run

```bash
python3 04-computation/lrc14_sparse_owner_head_visibility_thm2558.py
python3 -O 04-computation/lrc14_sparse_owner_head_visibility_thm2558.py
```

Both executions byte-match the stored output.  The referee evaluates all
`98,280` mask--slope selectors, reconstructs every histogram in the table,
checks translation-orbit counts, proves the singleton/pair multiplicities,
exhibits (16)--(17), and checks exact rational unit/guard capacity controls.

**QED.**
