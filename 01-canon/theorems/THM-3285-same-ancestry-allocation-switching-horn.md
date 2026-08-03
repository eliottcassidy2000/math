---
id: THM-3285
title: "Same-ancestry allocation-switching horn"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE AUDIT.
  On the fixed rail-eight, clock-one THM-2806 common-label bank, exactly
  63 labels give a three-address target carrier whose allocation type is
  right-cofiber/common/right-cofiber.  Every positive vertex is one exact
  whole cylinder with the same mass and carry-six coefficient, and all three
  translations preserve the literal THM-2791 rail ancestry sheet.  This
  statement is not promoted into the proved dependency graph until a fresh
  hostile audit independently reconstructs the carrier typing, label census,
  and ancestry chamber.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
script: 04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py
output: 05-knowledge/results/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out
script_sha256: c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee
output_sha256: e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246
hash_basis: LF-normalized bytes
---

# THM-3285 -- same-ancestry allocation-switching horn

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE
AUDIT.**

The exact companion has completed its primary normal/optimized/stored audit
and contains two exact carrier constructions.  The candidate below remains
outside the proved dependency graph until a fresh implementation has checked
the one-sided carrier types, the complete label universe, and the literal
ancestry claim.

## 1. Candidate universe and common-carrier typing

Use the THM-2782/2806 fixed rail-eight geometry.  Fix

```text
physical clock =1,
sigma in S={0,1,2,3,8,9,10,11,12},
tau   in T={3,4,5,6,7,8,9,10,11},                      (1)
```

together with both relative-present safeties, target root one, and the
delayed carry-six coefficient.  This is the complete `9*9=81`-label
clock-one slice of THM-2806's genuine common bank, not the distinct
`tau=12` target-only simplex.

Before any narrow address restriction, write in the source chart

```text
A     = source one-sided carrier,
B^-   = pullback of the target one-sided carrier,
M     = A intersect B^-,
R     = B^- minus A.                                   (2)
```

Equation `(2)` is a literal half-open support decomposition.  For weighted
typing, let `M_B` be the `B^-`-weighted restriction to `supp(A)` and let
`R_B` be the `B^-`-weighted restriction to its complement.  Thus

```text
B^-=M_B disjoint-union R_B,                             (3)
```

and `M_B,R_B` are the target-carrier copies of the common face and right
cofiber.  Push both copies to the target coordinate before cutting by an
address cylinder.  This order is load-bearing: it compares two disjoint
physical pieces in one chart, not two separately integrated scalars.

Use the three THM-2807 addresses

```text
n_0=3454614,
n_+=n_0+13=3454627,
n_a=n_0+689364=4143978,                                (4)
```

whose central indices are `(0,1,53028)`.  Denote the resulting target-chart
cylinder restrictions by `M_n`, `R_n`, and `B_n`.

## 2. Candidate statement: the exact 63-label horn

Put

```text
S_horn={0,3,8,9,10,11,12}.                             (5)
```

For every one of the

```text
1 clock * 7 sigma labels * 9 tau labels =63 labels      (6)
```

in `S_horn x T`, the complete allocation pattern is

```text
                 n_0             n_+             n_a
M                 0                C                0
R                 C                0                C
B=M disjoint R    C                C                C,  (7)
```

where each `C` is one whole open weighted cylinder.  In particular the
allocation word is

```text
R -> M -> R,                                           (8)
```

even though the undecomposed pulled target carrier has constant word
`B=1 -> 1 -> 1`.

Every positive `C` in `(7)` has the same exact rail weight, mass, and actual
delayed carry-six/root-one coefficient:

```text
weight(C)=27581135604,
mass(C)=60781651775958960/371293,
coefficient(C)=790161473087466480.                     (9)
```

No value in `(9)` is assigned after integration; the address restriction is
performed on the physical weighted carrier first.

## 3. Exact whole-cylinder translations

The target-coordinate shifts associated with the three address gaps are

```text
epsilon_1 =7(n_+-n_0)/13^6 =7/371293,
epsilon_2 =7(n_a-n_+)/13^6 =28553/28561,
epsilon_d =7(n_a-n_0)/13^6 =371196/371293.             (10)
```

For every label in `(6)`, exact half-open interval translation gives

```text
T_(epsilon_1) R_(n_0)=M_(n_+),
T_(epsilon_2) M_(n_+)=R_(n_a),
T_(epsilon_d) R_(n_0)=R_(n_a),                         (11)
```

and `epsilon_1+epsilon_2=epsilon_d mod 1`.  Equations `(9)--(11)` are
therefore whole-carrier identities, stronger than equality of the three
integrals.

The mechanism is a double crossing of the one-sided allocation boundary.
The pulled target carrier `B^-` survives at all three addresses, whereas the
source carrier `A` meets it only at the middle cylinder.  Its disjoint
decomposition consequently switches right cofiber, common overlap, right
cofiber.  Scalarizing `B` forgets precisely this switching coordinate.

## 4. Literal same-ancestry audit

THM-2791 types the rail-eight weight before Perron marginalization as the
Boolean Cartesian sheet

```text
U(a,b) x V(e'),
(a,b) in Z/(13^5) x Z/169,       e' in Z/(13^5).       (12)
```

The candidate companion rebuilds every raw `Q`, `E`, and rotated-`E`
contributor wall.  All three cylinders in `(4)` lie strictly inside the one
common chamber

```text
[140890500190440,144190879112280).                     (13)
```

It independently enumerates the complete contributor set at the new middle
vertex.  The two literal factor counts and their typed digest are

```text
|U|=966606,                 |V|=28534,
|U|*|V|=27581135604,
digest(U,V)=15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd.
                                                               (14)
```

The digest in `(14)` equals the independently audited THM-2791 endpoint
digest.  More fundamentally, `(13)` contains no contributor wall, so the
identity map on `(a,b,e')` identifies the complete label sets at all three
vertices.  The explicit Boolean copy

```text
(a,b,e')=(59162,26,56658)                              (15)
```

is active throughout, and the collision/source-root/arrival-root/deep labels
remain `(5,6,12)` at each cylinder.  Thus “same ancestry” means the literal
outer THM-2471 rail sheet is fixed; it is not an inference from equal weights.

## 5. Complete hostile boundary inside the 81-label universe

The remaining `18` labels in `(1)` are exact hostiles, not omitted cases:

```text
sigma=1, every tau in T:       M=000, R=001;
sigma=2, every tau in T:       M=000, R=101.            (16)
```

Together `(7)` and `(16)` exhaust all `81` labels.  Hence `(5)` is the sharp
sigma boundary for this fixed clock-one, rail-eight, three-address universe.
No relabeling or orbit normalization was used to hide `sigma=1,2`.

There is also an external hostile already in canon.  The `tau=12` label is
not in `T`; THM-2806 proves that its selected THM-2807 simplex is target-only
before address restriction.  The present candidate does not revive a common
atom there and does not contradict that no-go.

## 6. Exact proof candidate and reproducibility

The companion uses two calculation paths.

1. It begins with the independently audited full weighted one-sided source
   and target carriers, inserts each full semantic section, pulls the target
   carrier back, and forms `(2)--(3)` by direct set operations.
2. Independently in code structure, it starts from the canonical
   `carrier_objects` support decomposition and applies the inherited
   relative-safety, target push, root-one, and address filters in their
   prescribed order.

The paths agree exactly on all

```text
81 labels * 3 addresses * 2 pieces (M,R)=486 cells.     (17)
```

The script then verifies all `63*3=189` translations in `(11)`, rebuilds the
ancestry walls, enumerates the complete middle contributor set, and checks
the explicit path `(15)`.  Its stable full-universe pattern digest is

```text
9fd5090fb89493d07115ee753a0eadeba82d70d0e52663aa11d60c121f4f9d4f. (18)
```

Run

```bash
python3 04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py
python3 -O 04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out.
```

The LF-normalized source and output hashes are

```text
c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee
e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246.
```

The source has zero Python `assert` nodes and zero floating-point literals.
These checks make the primary candidate reproducible; they do not replace the
pending fresh hostile theorem audit.

## 7. Connection contract and sharp scope

```text
source:
  the 81 clock-one label cells of the fixed rail-eight one-sided bank on the
  three THM-2807 address cylinders;

target:
  a word in the disjoint allocation alphabet {M,R,empty};

map:
  B^- |-> (M=A intersect B^-, R=B^- minus A), then target push and exact
  address restriction;

preserved:
  whole weighted cylinder, exact mass and coefficient, literal rail ancestry,
  clock, sigma, tau, relative safeties, target root, carry, and open endpoints;

lost / not constructed:
  endpoint origin, endpoint current, the full bare/source/target/both K4,
  factorwise translation covariance, and a global packet action;

sidecar:
  the contributor-wall chamber and typed path-set digest;

cheapest next test:
  attach one explicit THM-2772 endpoint-origin fibre to one horn label and
  compute all four literal allocation states at all three vertices.        (19)
```

This candidate supplies a genuine common-carrier middle vertex between the
two THM-2791 `tau=3` endpoint cylinders on one literal ancestry sheet.  It
does **not** supply endpoint allocation, a determinant/root-Cech current,
Fourier noncancellation, row exclusion, or `LRC(14)`.

The procedural synthesis and correction-lineage discussion are recorded in
[`lrc14-same-ancestry-allocation-switching-horn-codex-20260803.md`](../../07-reflections/lrc14-same-ancestry-allocation-switching-horn-codex-20260803.md).
