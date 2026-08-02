# Rank-four forest holotopy: literal positive lifts die on the first faces

Date: 2026-08-02  
Status: **FINITE-EXACT STOPPING RESULT / REFLECTION**, not a theorem and not a
proved dependency

## Inheritance pass

The anchor is proved
`THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction.md`.
It replaces arbitrary anchored width-three product-Gamma coefficients by two
exact labelled Ewens currents `W_1,W_2`, supported on rank-four set partitions
of eight and nine macro packets.  Its closest near miss is the uniform forest
lift recorded in
`arbitrary-anchored-product-gamma-n6-jucys-murphy-boundary-codex-20260802.md`:
deleting one edge leaves many nonzero three-forest boundaries.

The new question is stronger and choice-free inside the natural positive
class.  Can the weight of each macro partition be redistributed among all
spanning forests in its fibre, without changing its sign, so that the result
is a genuine four-cycle?  The exact answer is no.  Every such coefficient is
forced to zero by one-sign three-faces.

The niche connection is the cycle-weighted Young-subgroup decomposition in
`THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary.md`.
The forest obstruction explains why its projection gaps, rather than a
literal positive Orlik--Solomon chain, are the viable carrier.  The wildcard
is an oriented local-system lift: the obstruction leaves that possibility
open and in fact identifies the intra-fibre sign cancellation it must carry.

## 1. The literal forest lift

Let `E_N` be the lexicographically ordered edges of `K_N`.  Write
`F_4(N)` for the canonically oriented four-edge forests, with boundary

```text
partial[e_1,e_2,e_3,e_4]
 =sum_(j=1)^4 (-1)^(j-1)[e_1,...,omit e_j,...,e_4].    (1)
```

The component partition `pi(F)` of a four-forest has rank four.  Let `W(pi)`
be the exact THM-3110 zeta weight.  A **literal same-sign distribution lift**
is a chain

```text
C=sum_(F:W(pi(F))!=0) sign(W(pi(F))) x_F [F],
x_F>=0,
sum_(F:pi(F)=pi) x_F=|W(pi)|.                          (2)
```

The last condition says that forgetting the tree inside each component
recovers the macro current exactly.  It is deliberately stronger than an
arbitrary signed lift and weaker than the failed uniform lift: the mass may
be redistributed in any way inside every forest fibre.

For a three-forest `G`, the equation `[G] partial C=0` is

```text
sum_(F covers G) sign(W(pi(F))) epsilon(G,F) x_F=0,    (3)
```

where `epsilon(G,F)` is the deletion sign in `(1)`.  If all currently live
coefficients in `(3)` have the same sign, nonnegativity forces every one of
their `x_F` to vanish.  Substituting these exact zeros may make further faces
one-signed.  This is a finite zero-propagation proof, not an LP tolerance or
floating-point infeasibility claim.

## 2. Exact propagation census

The companion reconstructs `W` directly from the live THM-3110 `24/25` banks,
enumerates every four-forest, builds every incident oriented three-face, and
at each round simultaneously fires all currently one-signed faces.

```text
                         I1 / K8                 I2 / K9
nonzero W fibres          480                     1620
positive / negative       285 / 195               720 / 900
all four-forests           18865                   55755
literal variables          1440                    4860
incident three-faces       1860                    3960

round: forcing faces / newly killed variables
1                         1203 / 1368             2106 / 4038
2                          120 /   72              496 /  551
3                                                    228 /  244
4                                                     73 /   27
remaining variables          0                        0.       (4)
```

Consequently `partial C=0` implies `x_F=0` for every forest in `(2)`.  This
contradicts the positive fibre sum `|W(pi)|` on every nonzero macro fibre.
Thus neither bank has a nonzero literal same-sign forest-cycle lift.

This strictly strengthens the old uniform-lift hostile.  The failure is not
caused by choosing equal weights among the spanning trees.  It persists for
every nonnegative redistribution preserving the macro sign.

## 3. A three-face local witness

One fibre makes the mechanism visible without the census.  In `K_8` take

```text
pi=((0,1,2),(3,4),(5,6),(7)),       W_1(pi)=1/30.      (5)
```

Its three forests differ only by the spanning tree on `{0,1,2}`.  Deleting
the edge `(3,4)` produces the faces

```text
{(0,1),(0,2),(5,6)},
{(0,1),(1,2),(5,6)},
{(0,2),(1,2),(5,6)}.                                  (6)
```

Each face has exactly three incident nonzero literal variables and all three
oriented signs are `+`.  Hence all three forests in `(5)` vanish.

For `K_9`, append the singleton `(8)`.  Now

```text
W_2(pi)=-1/60.                                         (7)
```

Each face in `(6)` has six incident nonzero variables.  Five lie over
negative `-1/60` fibres with positive deletion sign; the sixth lies over a
positive `1/30` fibre with negative deletion sign.  All six **oriented**
coefficients are therefore `-`.  Again the entire displayed fibre dies.
The global propagation in `(4)` is an iteration of precisely this local
phenomenon.

## 4. What the obstruction preserves and destroys

The relevant connection can now be typed exactly.

```text
source:       canonically oriented four-forest chain
target:       rank-four macro-partition Ewens current
map:          forget the tree, retain its component partition
preserved:    fibre mass, macro sign, rank, symmetric relabelling
destroyed:    root, edge order beyond its sign, insertion history,
              and cancellation between presentations of one macro flat
needed sidecar: an oriented/rooted local system allowing opposite signs
                inside a single macro fibre.                         (8)
```

The no-go is therefore not “topology cannot help.”  It says that topology
cannot be added *after* symmetrization by merely spreading positive macro
mass among trees.  A successful lift must be created before the row labels
are collected, or must carry an additional orientation/root character whose
opposite sheets cancel on the faces in `(6)`.

This also clarifies why the THM-3112 Young-subgroup carrier is different.
For a set partition `pi`, `I-P_(H_pi)` is already a positive projection gap;
when `pi` is refined by `kappa`, subgroup inclusion orders these gaps.  It
forgets the alternating forest boundary and retains the monotone quotient
that `(4)` does not obstruct.  Grouping the colourings by occupancy type
suggests the exact next object

```text
beta_N(S)=sum_(mu partition N) m_mu(S) Lbar_mu,         (9)
```

where `Lbar_mu` is the uniform average of `I-P_(H_pi)` over set partitions
of block-size type `mu`.  A biregular refinement incidence couples the two
uniform measures and gives

```text
mu refines nu  ==>  Lbar_mu <= Lbar_nu.                (10)
```

Equations `(9)--(10)` are the quotient/star--mesh lane untouched by the
forest no-go.  At this checkpoint they are a proof target being tested, not
a new proved dependency of this reflection.

## 5. The surviving signed-holotopy lane

A separate sparse numerical solve gives useful but explicitly noncanonical
evidence.  If all four-forests are allowed arbitrary signed coefficients and
zero macro fibres may carry cancelling null-pairs, the simultaneous system

```text
partial C=0,                 sum_(pi(F)=pi) C_F=W(pi)  (11)
```

has numerical solutions in both banks: relative residuals about `7e-12` in
`K_8` and `1.5e-10` in `K_9`.  This is not an exact existence theorem and is
not part of the companion transcript.  Its value is diagnostic: together
with `(4)`, it points specifically to **intra-fibre signed null-pairs**, not
to a missing amount of positive mass.

The cheapest decisive continuation is therefore one of:

1. construct a rational rooted/NBC lift of `(11)` by deletion--contraction;
2. prove the occupancy-type refinement flow for the coefficients in `(9)`
   symbolically in every degree; or
3. find a small exact dual witness showing that even a chosen rooted local
   system is too small.

A further uniform positive redistribution or another unrooted forest basis
is ruled out by `(4)`.

## 6. Reproduction and scope

Run

```text
python 04-computation/gmc_product_gamma_rank4_forest_holotopy.py
python -O 04-computation/gmc_product_gamma_rank4_forest_holotopy.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_product_gamma_rank4_forest_holotopy.out.
```

The script uses exact integers and rational numbers only.  It checks the
THM-3110 zeta counts, all forest and face counts, every propagation round,
and the local witness `(5)--(7)`.

LF-normalized SHA-256 at this checkpoint:

```text
script  5b40ba16ee89709bbc62212e576a192a37303940a5d6ed6bee3a8a0c3345f175
output  f50cc7852b694f8066d31a3efea4ab26bcdb91455454f791d367d2d31240a20f
```

This reflection does not prove the all-degree row extremum, arbitrary
anchored product-Gamma goodness, SFC in width three, GMC(2), NC2, LRC(14),
JC(2), or DC(2).  It proves only that the canonical same-sign forest lift is
impossible and locates the missing datum inside oriented fibre cancellation.
