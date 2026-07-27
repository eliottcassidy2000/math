---
id: THM-2504
title: "Endpoint tournament no-go and thirteen-root chart-holonomy boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the abstract
  untyped affine seven-cube, every
  unordered pair is swapped by a cube translation, so no invariant
  antisymmetric orientation exists. If bit polarity is retained but the
  seven role names are forgotten, the sharp minimum number of ties is
  1,652: every equal-Hamming-weight pair is swapped by a coordinate
  permutation, while weight order orients every unequal-weight pair.
  Hamming parity orients 4,096 pairs and ties 4,032, but on the 128-vertex
  tile graph it is the exact cut delta(parity), with support 4,011 and zero
  cycle part. The THM-2502 lexicographic tournament is complete only after
  role order and polarity are supplied and has zero reversal cochain in
  that gauge. THM-2457's physical co-support matrix cannot repair this:
  its sharp uniform common-chart packet model ties all 8,128 opposite-entry comparisons.
  On thirteen root labels the cyclic half-set tournament has nonzero
  triangle holonomy, but its unique tile-graph cut/cycle support split is
  38/35 and inversion sends the directed tournament to its global converse.
  Holding the numeric path fixed gives a 30/39 split, while transporting the
  path gives the isomorphic 38/35 split. Thus the tournament consumes an
  orientation representative; quotienting by converse loses semantic
  direction. Quadratic residues and guard autocorrelation are
  negation-symmetric and cannot orient. The intrinsic endpoint object is a
  typed directed weighted co-support matrix plus chart, word, and current
  sidecars, not a tournament. No scalar row is removed.
source: mac-mini-2026-07-27-endpoint-holonomy-no-go
depends_on:
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2458-clean-root-guard-danger-thirteen-chart-uniform-offset-hostile
  - THM-2467-bicycle-spaces-of-the-star-flip-split
  - THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary
related:
  - THM-1405-star-quotient-is-the-cycle-space-transverse-to-isomorphism
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2503-sharp-duplicate-probe-active-dipole-sweep
script: 04-computation/lrc14_endpoint_tournament_holonomy_no_go_thm2504.py
output: 05-knowledge/results/lrc14_endpoint_tournament_holonomy_no_go_thm2504.out
script_sha256: a13e2994f1ab5b99c8a73c99fc5cc3d4439e38807569796f51012e1cf077cd66
output_sha256: e0cdfb479da71680a070c45914add0f59b68334f5474cf693787e06e90860a6d
hash_basis: working-tree bytes (LF)
---

# THM-2504 -- endpoint tournament no-go and root-chart holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2502 found a transitive tournament on the `128` complete endpoint
masks. THM-2457 then identified the actual physical relation: a directed
weighted co-support matrix whose two axes have different semantic types.
This theorem determines whether an intrinsic tournament or a tournament
holonomy survives between those descriptions.

It does not. There is a precise sidecar ladder:

| retained structure | strongest immediate relation | exact loss |
|---|---|---|
| abstract untyped affine cube | none | every pair is symmetry-swapped |
| bit polarity, no role/type names | Hamming-weight preorder | `1,652` forced ties |
| Hamming parity | even/odd partial orientation | exact cut, zero cycle part |
| ordered polarized roles | THM-2502 lex tournament | transitive gauge choice |
| typed physical packets | THM-2457 matrix `M` | loops, missing and bidirected pairs |
| oriented `C_13` chart | cyclic half-set tournament | directed only up to global converse |

Thus tournament language becomes proof-bearing only after its pairwise
observable and every orientation/tie sidecar are declared. The surviving
intrinsic object is not an uncoloured tournament.

## 1. Three endpoint sidecar levels

Put

```text
Omega=F_2^7,             |Omega|=128.                         (1)
```

We distinguish three exact quotients.

1. At level `A`, only the abstract untyped affine Hamming-cube structure
   is retained.
2. At level `P`, safe/danger polarity is fixed, but all seven physical
   role names and types (guard/unit/blocker) are forgotten; the allowed
   relabelings are `S_7`.
3. At level `L`, polarity and the full ordered list of roles are retained.

These levels are an audit of information loss, not a claim that XOR
translations or arbitrary role permutations are symmetries of the physically
embedded LRC masks. The safe/danger gates have different measures and semantic
types. That physical typing is exactly sidecar data omitted at levels `A/P`.

For any distinct `x,y in Omega`, translation by

```text
d=x+y=x xor y                                             (2)
```

interchanges them:

```text
x+d=y,                 y+d=x.                            (3)
```

An `A`-invariant orientation would therefore satisfy `x -> y` if and
only if `y -> x`, contradicting antisymmetry. Hence:

> **Affine-cube no-go.** No tournament on `Omega` is invariant under all
> translations of the unlabeled affine seven-cube.

At level `P`, suppose `|x|=|y|`. Pair the coordinates in `x\y` with
those in `y\x`, transpose every paired coordinate, and fix the rest.
This coordinate involution interchanges `x` and `y`. Every such pair is
therefore a forced tie for an `S_7`-invariant antisymmetric partial
orientation. Their number is

```text
sum_(k=0)^7 binom(binom(7,k),2)=1,652.                  (4)
```

This is sharp: orient every unequal-weight pair from smaller to larger
weight. Thus `1,652` is the exact minimum tie count at level `P`, not
merely a lower bound.

THM-2445's labels are still coarser. Its six aggregate first-failure
layers have sizes

```text
4,64,32,16,8,4,                                          (5)
```

and hence leave

```text
sum binom(layer size,2)=2,672                            (6)
```

within-label pairs. Retaining the four blocker states refines these to
the `24` cells with per-state sizes

```text
1,16,8,4,2,1,                                            (7)
```

and leaves exactly

```text
4 sum binom(cell size,2)=620                             (8)
```

within-cell pairs. Equations (6)--(8) quantify why first-failure data are
useful routing data but not a complete pairwise observable.

## 2. Hamming parity is a cut, not a holonomy

There are `64` even- and `64` odd-weight masks. Orienting only
even/odd pairs therefore orients

```text
64*64=4,096                                              (9)
```

pairs and ties

```text
2 binom(64,2)=4,032.                                    (10)
```

Let `pi(x)=|x| mod 2`. On any graph the crossing cochain is literally

```text
c_(xy)=pi(x)+pi(y)=delta pi(xy).                        (11)
```

It is an exact cut. To compare with THM-1405/2467, use the tile graph

```text
H_128=K_128 minus (0--1--...--127).                     (12)
```

The full cube has `4,096` parity-crossing pairs. Incrementing `k` to
`k+1` flips `1+nu_2(k+1)` bits, so the removed successor edge crosses
parity exactly when `nu_2(k+1)` is even. There are

```text
64+16+4+1=85                                            (13)
```

such edges. Hence `c` has tile support `4,096-85=4,011`.

Exact row reduction gives

```text
|E(H_128)|=8,001,       rank_F2 Lap(H_128)=127.          (14)
```

Equivalently, the bicycle space is zero, as also follows from THM-2467
because `128=8 mod 12`. The cut--cycle split is unique, and (11) has

```text
cut support=4,011,      cycle support=0.                 (15)
```

Parity supplies no Wilson-loop or cycle coordinate.

## 3. The complete mask tournament is a chosen total-order gauge

At level `L`, orient a pair at its first differing role from bit `0` to
bit `1`. This is THM-2502's lexicographic tournament. With the same role
order used as binary significance, it is exactly numeric order

```text
0 -> 1 -> ... -> 127.                                   (16)
```

Define the reversal cochain relative to that reference order to be one
on an edge whose orientation is reversed. It is identically zero, so its
tile-graph cycle component and every triangle holonomy are zero.

This does not make the ordered roles useless. THM-2502's edge colours
retain the first-difference/carry depth. It says that after edge colour,
role order, polarity, and the unused tail are erased, the remaining
uncoloured tournament contains no nontransitive residue that could replace
a physical relation current.

## 4. Directed co-support does not canonically complete to a tournament

THM-2457 supplies, after a common oriented root chart and two typed
semantic packets `A,F` are given, the nonnegative matrix

```text
M_(omega,nu)=integral a_omega(y) f_nu(y) dy.             (17)
```

Its direction means `A`-atom to `F`-atom. Loops have physical meaning:
`M_(omega,omega)>0` is exact-atom service. Off-diagonal pairs may be
absent, one-way, or bidirected. This is more information than an
unweighted tournament in some directions and less in others.

The sharp common-chart packet model explicitly constructed in THM-2457 has

```text
M_(omega,nu)=1/16,384             for every ordered pair. (18)
```

Consequently every one of the `8,128` comparisons

```text
M_(omega,nu) ? M_(nu,omega)                              (19)
```

is tied. This is a sharp model inside THM-2457's packet axioms, not a full
scalar LRC counterexample row. Thus antisymmetrizing `M` cannot universally produce even one
tournament edge. A lexicographic or arbitrary completion of (19) is extra
data, not a consequence of co-support. In particular, a sufficient
tournament certificate could never be treated as equivalent to the
directed weighted service object.

## 5. Thirteen roots: genuine holonomy spends the chart orientation

Now let the root labels be `C_13`. Choosing the positive half-system

```text
H={1,2,3,4,5,6}                                        (20)
```

defines the cyclic tournament

```text
u -> v       iff       v-u in H.                       (21)
```

Unlike the endpoint lex tournament, (21) has a genuine cycle part. On

```text
H_13=K_13 minus (0--1--...--12),                       (22)
```

its reversal cochain relative to numeric order is one exactly when the
positive numeric difference is in `{7,...,12}`. Therefore

```text
|E(H_13)|=66,            cochain support=21.             (23)
```

The exact Laplacian rank is `12`; THM-2467 gives the same zero bicycle
dimension because `13=1 mod 12`. Its unique decomposition has

```text
cut support=38,           cycle support=35.              (24)
```

For the tile triangle `(0,2,7)`, the three reversal bits are

```text
0,0,1,                                                     (25)
```

so its holonomy is one.

But THM-2457 permits a harmless common chart relabeling `u -> -u`.
It exchanges `H` and `-H` and sends the directed tournament (21) to its
global converse. If the numeric path in (22) is held fixed, the cochain
becomes `c+1`. Its
unique decomposition has

```text
cut support=30,           cycle support=39.              (26)
```

In particular its fixed-path cycle vector differs from (24). This support
change is a coordinate statement, not a gauge invariant: transporting the
path together with `u -> -u` gives an isomorphic tile graph and preserves
the original `38/35` support split. What fails to descend is the **directed
choice** between the tournament and its global converse. One may quotient by
that pair, but the quotient no longer orients the semantic `A -> F` current.
Thus a directed cyclic tournament is available only after an oriented root
chart, or an equivalent half-system representative, is retained.

Two tempting canonical replacements fail for the same exact reason.
The nonzero quadratic residues modulo `13` are

```text
QR={1,3,4,9,10,12}=-QR,                                 (27)
```

because `-1` is a square. The relation `v-u in QR` is symmetric, not a
Paley tournament. For the THM-2458 four-root guard AP with start zero
and step five,

```text
S={0,5,10,2},                                            (28)
```

the translate autocorrelation is

```text
|S intersect (S-t)|, t=0,...,12
 =4,0,1,2,0,3,0,0,3,0,2,1,0.                           (29)
```

It satisfies `A(t)=A(-t)`. Guard-translate overlap therefore cannot orient
root differences. Other starts translate (28), and other AP steps permute
the same negation symmetry. This statement is only about the guard-overlap
observable; it does not exhaust every typed field of THM-2458's full
`(q,a,b,c)` atlas.

## 6. Exact consequence and next proof object

The no-go is scoped:

```text
unlabeled endpoint cube
  -> no invariant tournament;

ordered endpoint labels
  -> a transitive tournament with no cycle residue;

root half-system
  -> nonzero cycle residue and a converse pair, but a directed choice only
     with oriented-chart sidecar;

physical endpoint service
  -> typed directed weighted co-support, not pairwise completeness. (30)
```

Nothing here forbids a tournament argument when an intrinsic binary
observable really exists. It proves that an LRC endpoint tournament must
declare at least:

```text
vertices; A/F semantic type; pairwise observable; root chart/section;
orientation gauge; tie completion; preserved target/current; lost data. (31)
```

The next proof-bearing object should therefore be a **typed current
cochain**: retain the directed weights `M`, the service loops, a common
root-character vector, and the semantic owner/terminal word before taking
any cut/cycle quotient. A successful holonomy must transform covariantly
under common chart reversal and still preserve the THM-2350 target/deep
character. Neither an uncoloured mask tournament nor the antisymmetric
scalar part `M-M^T` meets that specification.

This refines the structural bottleneck after THM-2466. Delayed mixing can
retain drift and service once a common owner-supported root base and word
are supplied; the missing step is still to produce that typed common base
from the physical scalar-cover coupling. No scalar row is excluded, the
open ledger remains `165`, and LRC(14) is not proved.

## 7. Exact companion

Run

```text
python3 04-computation/lrc14_endpoint_tournament_holonomy_no_go_thm2504.py
python3 -O 04-computation/lrc14_endpoint_tournament_holonomy_no_go_thm2504.py
```

The dependency-free exact companion:

- exhausts all `8,128` mask pairs and constructs every translation and
  equal-weight coordinate swap;
- reproduces all parity, first-failure-layer, and `24`-cell tie counts;
- verifies the lexicographic zero cochain and the uniform co-support tie
  hostile;
- computes both tile-graph Laplacian ranks and unique cut--cycle splits;
- checks the nonzero `(0,2,7)` triangle, fixed-path global converse, and
  fully transported inversion control; and
- verifies the quadratic-residue and guard-autocorrelation symmetries.

Normal and optimized runs must byte-match

```text
05-knowledge/results/lrc14_endpoint_tournament_holonomy_no_go_thm2504.out
```

QED.
