---
id: THM-2780
title: "Marked D3 weight-colouring torsor descent and affine inertia gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The four absolute determinant weight-colourings of the six D3
  root lines are W(D3)-equivariantly the four even-sign states, hence a
  regular V4-set and exactly the abstract generic fibre of the quartic
  Kummer torsor.  A retained A2 chamber gives a chiral tournament, but its
  weighted switching class is self-converse after the chamber gauge is
  forgotten, so no orientation character descends.  For a realized Kummer
  plane the colouring cover is quasi-etale exactly when every normalized
  divisor row is zero.  Nevertheless the marked fibre, chamber, common root
  sum, and even the complete divisor-row collection do not determine the
  torsor: two distinct S3-standard unit Kummer planes on one four-torus give
  an exact unramified twin hostile.  The minimal missing coordinate is the
  embedded Q-equivariant Kummer H1-class, not another finite D3 marking.
source: d3-torsor-descent/marked-weight-colouring-affine-gate-2026-07-28
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
  - THM-2777-marked-d3-six-root-determinant-tournament-and-binary-ternary-edge-spectrum
related:
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2775-modular-s4-to-weyl-d3-generator-frame-and-affine-parity-blindness
script: .scratch/d3_torsor_descent/thm2780_weight_colouring_descent.py
output: .scratch/d3_torsor_descent/thm2780_weight_colouring_descent.out
script_sha256: 456e87aa2411798ac6678d03cb8e69b4e741a944cb1f6d1e80c29d9d8d736385
output_sha256: c04b424af35ac8dccc48e518f851d85c4f370f30c71e8df8dc0459ed52d16279
hash_basis: LF-normalized bytes
---

# THM-2780 -- the four weight colourings are the torsor fibre, not its twist

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2777 leaves two apparently competing facts.  Its marked determinant
tournament is a genuine chiral object, but its orientation needs a chamber.
Its absolute `1,2,3` weights need no root orientations.  The exact survivor
is the latter: the four possible weight colourings are precisely the four
states of the quartic `V4` torsor.

That identification is a generic-fibre theorem, not an affine
trivialization.  A marked state is a point after the torsor has been pulled
back to itself.  The twist of the four-state fibre is a Kummer cohomology
class.  Divisor rows detect its ramification but not its unramified part.

## 1. Four absolute colourings of the six `D3` root lines

Let

```text
L={R(e_i-e_j), R(e_i+e_j): 1<=i<j<=3}                  (1)
```

be the six unoriented `D3` root lines in the oriented lattice `Z^3`.  Put

```text
Omega={delta in {+1,-1}^3: delta_1 delta_2 delta_3=1}.  (2)
```

For `delta in Omega` and two distinct lines `l,m in L`, choose any root
representatives `alpha in l`, `beta in m` and define

```text
kappa_delta({l,m})=|det(alpha,beta,delta)|.             (3)
```

Absolute value makes `(3)` independent of both representative signs.  For
every `delta`, the complete classification is

```text
kappa_delta=3  iff l and m both lie in delta^perp;
kappa_delta=2  iff l and m are orthogonal root lines;
kappa_delta=1  otherwise.                              (4)
```

In particular there are no zero pairs and the spectrum is

```text
1^9 2^3 3^3.                                            (5)
```

### Proof

The even signed-permutation group

```text
W(D3)=V4 semidirect S3                                  (6)
```

acts transitively on `Omega`.  It is therefore enough to take
`delta=(1,1,1)`.  The three perpendicular lines are represented by
`e_1-e_2,e_1-e_3,e_2-e_3`; their three pair determinants have magnitude
three.  The three orthogonal root-line pairs are

```text
{e_i-e_j,e_i+e_j},                                      (7)
```

and have magnitude two.  The remaining nine direct evaluations have
magnitude one.  This proves `(4)--(5)`.

For `g in W(D3)`,

```text
kappa_(g delta)({g l,g m})
 =|det(g)| kappa_delta({l,m})
 =kappa_delta({l,m}).                                   (8)
```

Thus the construction is `W(D3)`-equivariant.  Moreover, the three vertices
of the weight-three graph are exactly `L intersect delta^perp`.  Their
common perpendicular line determines `+/-delta`, and the product-one
condition in `(2)` selects exactly one sign.  Hence

```text
delta |-> kappa_delta                                   (9)
```

is injective.  Its image `C` has four elements.  Coordinatewise
multiplication by the diagonal `V4=Omega` acts simply transitively, so

```text
C is a regular V4-set, W(D3)-equivariantly isomorphic to Omega. (10)
```

This is stronger than the one-colouring statement in THM-2777: the complete
four-colouring orbit recovers the abstract quartic sheet torsor.

## 2. Why the directed tournament does not descend with it

Fix `delta=(1,1,1)`, orient the three roots in `delta^perp` by the chamber
`1<2<3`, orient the other three by positive dot product with `delta`, and
fix ambient volume.  The signs of the determinants in `(3)` give the
THM-2777 tournament.  It has score sequence

```text
(1,2,2,3,3,4)                                           (11)
```

and is not isomorphic to its converse.  Thus it is genuinely chiral while
the chamber and representative choices are retained.

Changing one root representative is vertex switching.  Write
`a_ij=e_i-e_j` and `b_ij=e_i+e_j`.  Switch the vertex `a_12`, then relabel

```text
a_23 <-> a_13,             b_13 <-> b_23.              (12)
```

The result is the global converse, and `(12)` preserves every absolute
weight.  Therefore

```text
retained chamber:       chiral weighted tournament;
forgotten chamber:      self-converse weighted switching class.        (13)
```

Odd `W(D3)` elements globally reverse the retained tournament, so that
tournament realizes the quartic sign character only before the chamber
gauge is removed.  No sign character survives in `(13)`.  The survivor is
exactly the four-state absolute colouring torsor `(10)`.

## 3. Exact field-level map to quartic roots

Let `E` be a characteristic-zero field.  Suppose

```text
tau_1 tau_2 tau_3=c^2,                                  (14)
dim_F2 <[tau_1],[tau_2],[tau_3]>=2                      (15)
```

in `E*/E*2`.  Choose `s_i^2=tau_i` and orient their product by

```text
s_1s_2s_3=c.                                            (16)
```

Then

```text
M=E(s_1,s_2,s_3)
```

is a connected `V4` Kummer extension.  Its Galois group is the even sign
plane, acting on the `s_i`.  Given a common root sum `e_1 in E`, put

```text
r_delta=e_1/4+(delta_1s_1+delta_2s_2+delta_3s_3)/2,
delta in Omega.                                         (17)
```

The four values in `(17)` have sum `e_1`.  In the centered coordinate
`Z=Y-e_1/4`, their polynomial is

```text
product_delta (Z-(delta.s)/2)
 =Z^4-(tau_1+tau_2+tau_3)Z^2/2-cZ
  +((tau_1+tau_2+tau_3)^2
    -4(tau_1tau_2+tau_1tau_3+tau_2tau_3))/16.           (18)
```

Also

```text
r_+++ + r_+-- - e_1/2=s_1,
r_+++ + r_-+- - e_1/2=s_2,
r_+++ + r_--+ - e_1/2=s_3.                             (19)
```

Thus `tau_i` are exactly the three squared complementary-pair-sum
coordinates, and `(18)` is the inverse quartic construction.

For an even sign vector `epsilon`,

```text
epsilon(r_delta)=r_(epsilon delta).                     (20)
```

Combining `(9)` and `(20)` gives the promised equivariant bijection

```text
kappa_delta |-> r_delta.                                (21)
```

It identifies the four-colouring regular set with the generic quartic root
torsor over the full matching field.

This is an actual Galois descent statement.  If `P=Spec(M)` is viewed as a
right `V4`-torsor, twist the constant regular set `C` by `P`.  Equivariance
of `(21)` gives an `E`-isomorphism

```text
P times^(V4) C  ->  the four-root finite E-scheme.       (21a)
```

Because two distinct even states differ in two signs, equality of their
values in `(17)` would make two `tau_i` have the same squareclass and the
third square; that contradicts `(15)`.  Hence the orbit is separable and
regular, each root generates `M/E`, and `(21a)` is the `V4` torsor itself,
not merely a quotient of it.

There is a vital type distinction in `(21)`.  The coefficient vector
`delta` or its colouring is an abstract fibre state.  The evaluated number
`delta.s` uses the actual Kummer square roots and exists upstairs on
`Spec(M)`.  If the symbol `h` means that evaluated number, then `h` already
contains the missing Kummer data and is not a sidecar on `Spec(E)`.

## 4. Divisor inertia acts by translation on the colourings

Let `R` be a normal integral affine complex domain with fraction field `E`,
let `Z` be its normalization in `M`, and let `D` run over the height-one
prime divisors of `R`.  Define

```text
nu_D=(v_D(tau_1),v_D(tau_2),v_D(tau_3)) mod 2.          (22)
```

Equation `(14)` places `nu_D` in the even code

```text
{000,110,101,011}=V4.                                   (23)
```

Over the DVR at `D`, tame Kummer inertia flips precisely the square roots
whose valuations are odd.  Hence the inertia element represented by
`nu_D` acts on both sides of `(21)` by the free translation

```text
delta |-> (-1)^nu_D delta.                              (24)
```

For `nu_D!=000`, this is a double transposition on the four states and has
no fixed point.  It preserves the family of four weight colourings and all
their `1,2,3` incidence, so the finite frame cannot reveal whether that
element has occurred as divisor inertia.

The normalization criterion is exact:

```text
Z -> Spec(R) is unramified in codimension one
iff nu_D=000 for every D.                               (25)
```

Equivalently, the colouring torsor extends quasi-etale exactly under
`(25)`.  A single marked colouring does not descend as a section of a
connected rank-two torsor: because the `V4` action is regular, such a
section would trivialize the generic torsor.  What descends is the
four-element cover, not one globally preferred `h`.

Neither `e_1` nor a chamber changes `(22)`: `e_1` is fixed by the Kummer
group, while a chamber only labels and orients the three coordinates.
THM-2769's row `110` is the sharp ramified hostile to any attempted repair
by these sidecars.

## 5. Even all divisor rows do not determine the torsor

The residue test `(25)` decides ramification of a **given** Kummer plane.  It
does not reconstruct that plane.  This loss is sharp even with full `S3`
equivariance.

Work on the four-dimensional torus

```text
R=C[x_1^+-,x_2^+-,x_3^+-,y_1^+-,y_2^+-,y_3^+-]/
  (x_1x_2x_3-1, y_1y_2y_3-1),                          (26)
```

where `+-` means adjoining both a variable and its inverse.  Let
`Q=S3` permute the three indices simultaneously.  In
`R*/R*2`, define

```text
W_x=<[x_1],[x_2],[x_3]>,
W_y=<[y_1],[y_2],[y_3]>.                               (27)
```

Both planes have rank two because the product in each triple is one.  Their
three nonzero classes are permuted transitively by `S3`, so each is the
standard `F2[S3]` plane.  They are nevertheless distinct.  Indeed

```text
R*/R*2 = <[x_1],[x_2],[y_1],[y_2]> = F2^4,             (28)
```

with `x_3=x_1x_2` and `y_3=y_1y_2` modulo squares, and the two planes in
`(27)` meet only in zero.

Every `x_i,y_i` is a unit of `R`.  Consequently all normalized
height-one-divisor rows of both Kummer planes are identically zero.  With
the same common sum `e_1=0`, the two inverse quartics are

```text
P_x(Z)=Z^4-(sum x_i)Z^2/2-Z
        +((sum x_i)^2-4 sum_(i<j)x_ix_j)/16,

P_y(Z)=Z^4-(sum y_i)Z^2/2-Z
        +((sum y_i)^2-4 sum_(i<j)y_iy_j)/16.            (29)
```

Each coefficient packet is `S3`-invariant, each generic Kummer plane is
connected, and each has the same abstract marked `D3` fibre, chamber
options, common sum, and complete zero residue collection.  But the covers
are not isomorphic over the fixed base, even after changing the `V4` basis:
connected `V4` Kummer torsors up to basis change are classified by their
embedded two-dimensional squareclass planes, and `(27)` are different.

Thus the tuple

```text
(marked D3 state, chamber, e_1, every divisor row)       (30)
```

does not determine the quartic `V4` torsor.  The exact missing coordinate is

```text
the embedded Q-equivariant Kummer class
W -> H^1_et(R_reg,mu_2),                                (31)
```

or equivalently the actual rank-two squareclass plane together with its
unit/class-group realization.  The boundary map records only the residues
of `(31)` and has a large unramified kernel; `(26)--(29)` isolate the unit
branch of that kernel.

## 6. Descent ledger

```text
SOURCE:        D3 root-line frame, its four absolute determinant colourings,
               and a realized rank-two product-oriented Kummer plane.

TARGET:        quartic four-root cover over the full matching normalization.

MAP:           kappa_delta -> e_1/4+(delta.s)/2.

PRESERVED:     regular V4 action; W(D3)=S4 covariance; opposite-edge
               matching; marked orthogonal triangle; weights 1,2,3;
               normalized inertia translation.

DESTROYED:     after chamber quotient, tournament chirality/sign;
               after taking residues, the unramified Kummer H1-class;
               after retaining only a marked fibre, the twisting cocycle.

NEEDED:        actual Q-equivariant squareclass plane W; all divisor rows
               to test quasi-etaleness; then its unit or Cl[2] carrier.

CHEAP HOSTILES:
               row 110 from THM-2769 for ramified descent;
               twin unit planes (27) for unramified nonuniqueness.       (32)
```

The theorem gives an exact interpretation of the four weight colourings and
an exact stopping rule.  It does not force the affine inertia rows to vanish,
identify the Kummer plane of an unknown graph quartic, produce a global
torsor section, or exclude `A4/S4` Keller monodromy.  It proves no `JC(2)`,
`DC(2)`, Graceful Tree, or `LRC(14)` result.

## 7. Exact companion

Run

```bash
python3 .scratch/d3_torsor_descent/thm2780_weight_colouring_descent.py
python3 -O .scratch/d3_torsor_descent/thm2780_weight_colouring_descent.py
```

The no-`assert`, integer/rational exact companion enumerates all four
colourings and all `24` Weyl elements; checks injectivity, covariance,
weight spectra, and the regular `V4` action; exhausts all `6!` tournament
relabelings to prove retained-chamber chirality; verifies the explicit
weighted switching self-converse witness; reconstructs the quartic
symbolically; and checks the two distinct simultaneous-`S3` standard
planes in the four-torus unit squareclass lattice.  Normal and optimized
runs byte-match the stored scratch transcript.

QED (candidate).
