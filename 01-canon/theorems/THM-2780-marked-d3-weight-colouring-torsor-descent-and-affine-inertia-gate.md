---
id: THM-2780
title: "Marked D3 weight-colouring torsor descent and affine inertia gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The four absolute
  determinant-weight colourings of the six D3 root lines are W(D3)-
  equivariantly the four quartic Kummer states: the weight-three triangle
  recovers the mark and V4 acts simply transitively.  A retained chamber
  gives a chiral tournament, but its weighted switching class is self-
  converse after chamber gauge, so directed sign does not descend.  For a
  given Kummer plane, all normalized divisor rows vanish iff the colouring
  cover is quasi-etale.  A nonzero row is the unique mod-two null of its
  determinant-two frame; determinant-three frames are mod-two units.  Even
  the complete zero-row collection, common root sum, chamber, and finite
  frame do not reconstruct the torsor: two distinct simultaneous-S3 unit
  Kummer planes on one four-torus are an exact unramified twin.  The missing
  coordinate is the embedded equivariant Kummer H1-class, not another finite
  mark.  No A4/S4, JC(2), or DC(2) exclusion follows.
source: jc2-v4-affine-gate/marked-weight-descent-2026-07-28
audit: >
  d3-torsor-descent/THM-2780-audit-2026-07-28 independently reconstructed
  all four colourings, Weyl covariance, stabilizers, inertia/Smith gates,
  retained-chamber chirality, the weighted self-converse witness, inverse
  quartic descent, and the distinct unramified S3-standard unit planes;
  replayed both companions under normal and optimized Python against stored
  output with zero assert nodes: ACCEPT.
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2777-marked-d3-six-root-determinant-tournament-and-binary-ternary-edge-spectrum
related:
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2775-modular-s4-to-weyl-d3-generator-frame-and-affine-parity-blindness
script: 04-computation/marked_d3_weight_torsor_descent_thm2780.py
output: 05-knowledge/results/marked_d3_weight_torsor_descent_thm2780.out
script_sha256: 22cbd2e0097e9a48f5695deae2bbd881142cedafcd5048d8001d11bb1be599f1
output_sha256: 2ab9c0ceca0993263f62c3637b58d13f0101af7c1c3d810b3ff8a1565642fe76
secondary_script: 04-computation/marked_d3_weight_torsor_descent_hostile_thm2780.py
secondary_output: 05-knowledge/results/marked_d3_weight_torsor_descent_hostile_thm2780.out
secondary_script_sha256: 20104ad6ddc7bfa83c743ef6d41977965c6b41bd62aeaaddc3b4fce0a24b8dbe
secondary_output_sha256: d11e6703f207d8c7b93fbcaf416b5bdb326fcb0c1573dac3b539c3be9f62a2ec
hash_basis: LF-normalized bytes
---

# THM-2780 -- marked `D3` weight colourings and affine inertia

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2777 obtains a genuine weighted tournament only after marking one of the
four half-Hadamard states.  The mark is not a harmless gauge.  Its
weight-three triangle recovers the marked state exactly, so the four possible
weightings form the original quartic `V4` torsor.  This gives a sharp affine
answer:

```text
an individual marked weighting = a quartic root-state section;
local inertia fixes one         iff the Kummer divisor row is zero.       (1)
```

Thus the marked `D3/A2` structure gives a combinatorial presentation of the
known parity gate.  It does not provide an independent reason for the parity
row to vanish.

## 1. The four determinant-weight colourings

Let

```text
Omega={1,2,3,4},                  E=binom(Omega,2).       (2)
```

The six elements of `E` are the six unoriented `D3` root lines under
THM-2777's half-Hadamard dictionary.  For a mark `h in Omega`, define a
weight on unordered pairs of distinct edges by

```text
w_h(e,f)=2,   if e and f are disjoint;
          3,   if e and f meet and h is not in e union f;
          1,   otherwise.                               (3)
```

Equation `(3)` is exactly THM-2777's determinant magnitude:

```text
w_h=2: the three opposite K4-edge pairs;
w_h=3: the three pairs of edges in the triangle Omega\{h};
w_h=1: the remaining nine pairs.                        (4)
```

In particular every `w_h` has spectrum

```text
1^9 2^3 3^3.                                            (5)
```

The weight-three support is the complete graph on the three edge-vertices

```text
binom(Omega\{h},2).                                     (6)
```

Taking the union of the endpoints in `(6)` gives `Omega\{h}`.  Therefore its
unique omitted point is `h`, and

```text
h |-> w_h                                                (7)
```

is injective.  It is visibly `S4`-equivariant:

```text
w_(g h)(g e,g f)=w_h(e,f).                              (8)
```

Consequently

```text
Stab_S4(w_h)=Stab_S4(h)=S3,
Stab_A4(w_h)=Stab_A4(h)=C3.                             (9)
```

The diagonal `V4` acts simply transitively on `Omega`; by `(7)--(8)` it acts
simply transitively on the four colourings as well.  Marking the weight-three
triangle is therefore literally choosing one point of the `V4` torsor.

The oriented tournament of THM-2777 contains still more choices: an `A2`
chamber and ambient volume orientation.  But it contains `w_h` as its arc
magnitude.  Hence descent of that weighted tournament already implies
descent of `h`; the chamber and orientation cannot repair a missing mark.

### 1.1 The chamber-free weighted switching class loses sign

The directed information has a sharper failure than orbit symmetrization.
Fix `h=(1,1,1)`, the chamber `1<2<3`, and the root names

```text
a_ij=e_i-e_j,                 b_ij=e_i+e_j.            (9a)
```

With these choices retained, the determinant tournament is chiral: an
exhaustion of all `6!` vertex relabellings finds no isomorphism to its
converse.  Now switch the single vertex `a_12`, meaning reverse every arc
incident to it, and then relabel

```text
a_23 <-> a_13,                b_13 <-> b_23.           (9b)
```

The result is the global converse, and `(9b)` preserves all absolute arc
weights.  Therefore

```text
retained chamber:       chiral weighted tournament;
chamber quotient:       self-converse weighted switching class.        (9c)
```

This refines THM-2777's labelled boundary: its converse is not in the same
labelled switching orbit, but it becomes switching-isomorphic after the
allowed chamber relabelling.  Thus directed sign is not an additional
chamber-free torsor coordinate; only the absolute colouring `(7)` survives.

## 2. Exact field-of-definition consequence

Let `L/K` be a separable quartic splitting field with transitive group

```text
G=A4 or S4,
```

and let

```text
E_0=L^V4                                                   (10)
```

be its cubic matching field.  Identify the four quartic root states with
`Omega` and their six difference lines with `E`.  Equivariance `(8)` and
the stabilizers `(9)` give

```text
L^Stab_G(w_h)=L^Stab_G(h).                               (11)
```

Thus one marked colouring has the same degree-four field of definition as
one quartic root.  More sharply,

```text
Stab_V4(w_h)=1,                  E_0(w_h)=L.             (12)
```

So over the cubic matching field an individual colouring contains the
entire rank-two Kummer extension.  It is not a new invariant living on
`E_0`.

An explicit inverse makes the descent map literal.  Suppose in a field `E`

```text
s_i^2=tau_i,            s_1s_2s_3=c,
dim_F2 <[tau_1],[tau_2],[tau_3]>=2.                   (12a)
```

For the four even sign states `delta`, and a common root sum `e_1`, put

```text
r_delta=e_1/4+(delta_1s_1+delta_2s_2+delta_3s_3)/2.   (12b)
```

Their centered polynomial is

```text
prod_delta (Z-(delta.s)/2)
 =Z^4-(sum tau_i)Z^2/2-cZ
  +((sum tau_i)^2-4sum_(i<j)tau_itau_j)/16.           (12c)
```

Conversely, complementary pair sums recover each `s_i`.  Even sign change
sends `r_delta` to `r_(epsilon delta)`, so

```text
w_delta |-> r_delta                                    (12d)
```

is a `V4`-equivariant isomorphism from the twist of the four-colouring set
to the quartic root scheme.  Rank two prevents collisions.  The abstract
label `delta` is base-independent; the evaluated quantity `delta.s` already
uses the actual Kummer square roots upstairs and is not a new finite sidecar.

This is a useful typing test.  A proposed *canonical* marked determinant
colouring on the matching normalization has already supplied a section of
the quartic `V4` torsor.  For a connected nontrivial `A4/S4` extension no
such generic section exists.  Keeping the orbit of all four colourings is
lawful, but that orbit is just the torsor again.

## 3. Divisor parity is the local fixed-colouring test

Now use the affine model of THM-2685.  Let `R` be a normal affine model of
`E_0`, and let its normalization in `L` be the rank-two Kummer cover.  At a
prime divisor `D` of `R`, write the three squared opposite differences as
`tau_1,tau_2,tau_3` and form

```text
nu_D=(v_D(tau_1),v_D(tau_2),v_D(tau_3)) mod 2
     in {000,110,101,011}.                              (13)
```

THM-2685 identifies a nonzero row with the corresponding tame `V4` inertia
element.  In the half-Hadamard states of THM-2775, the three cases act as
nonzero translations.  For example,

```text
110: (h_1 h_4)(h_2 h_3).                               (14)
```

It fixes no state and, by `(7)--(8)`, fixes no marked colouring.  The other
two nonzero words have the same `2+2` cycle type.  Hence

```text
nu_D=000
 iff the inertia fixes every marked colouring
 iff the inertia fixes at least one marked colouring.   (15)
```

Over the strict henselization at `D`, `(15)` says that the four-point
colouring cover splits into four unramified local sections exactly when the
Kummer row is zero.  If the row is nonzero, the colourings occur in two
ramified pairs.  Applying `(15)` at every prime divisor recovers precisely
THM-2685's quasi-etale criterion.

The qualification **local** is essential.  A connected quasi-etale `V4`
torsor can have zero parity row at every divisor and still have no global
section.  The class-group model

```text
Spec C[a,b,c,d]/(d^2-abc)                               (16)
```

from THM-2685 is the sharp control: its standard `V4` cover is quasi-etale
and connected, while its class-group plane is nontrivial.  Therefore
codimension-one extendability of all four local colourings does not choose
one globally.

## 4. The inertia direction is a binary root-line pair

Use THM-2777's oriented root names

```text
a_ij=e_i-e_j,                     b_ij=e_i+e_j.          (17)
```

The sign flip `110=diag(-1,-1,1)` fixes the two root **lines**

```text
{a_12,b_12}                                             (18)
```

and interchanges

```text
a_13 <-> b_13,                    a_23 <-> b_23.         (19)
```

Likewise `101` fixes `{a_13,b_13}` and `011` fixes
`{a_23,b_23}`.  Each fixed pair is one of the three weight-two opposite
`K4`-edge pairs.  Thus the three nonzero parity rows have an exact
binary-edge avatar:

```text
nonzero Kummer inertia
 -> unique pointwise-fixed weight-two root-line pair.    (20)
```

There is an equally exact Smith/parity formulation.  For the marked frame

```text
M_ij=[a_ij; b_ij; h],                                   (20a)
```

all four half-Hadamard states satisfy `h=111 mod 2`, while

```text
a_ij=b_ij=e_i+e_j mod 2.
```

Therefore

```text
ker(M_ij mod 2)=<e_i+e_j>.                              (20b)
```

This is exactly the parity word which flips coordinates `i,j`.  It is the
unique mod-two null direction of THM-2777's determinant-two frame.  By
contrast, every determinant-one frame and every determinant-three `A2`
frame has odd determinant and is invertible over `F_2`.  Thus the ternary
`P(A2)/Q(A2)=Z/3` defect supplies no binary null line and no independent
mod-two vanishing gate.  The only Smith rank drop aligned with Kummer inertia
is the already-known weight-two pair.

But `(20)` uses the **action** of inertia on labelled root lines.  The static
set of three weight-two pairs is present for every mark and every divisor,
so it does not select which row occurs.  The ternary `A2` index-three pairs
do not constrain this two-primary translation through their Smith defect:
their weight-three triangle moves with `h`, and their frames are mod-two
units.

## 5. Orbit symmetrization is exactly the stopping boundary

Forget `h` but retain, for each root-line pair, the multiset of its four
weights.  Equations `(3)--(4)` give the complete answer:

```text
e,f opposite:   {w_h(e,f):h in Omega}={2,2,2,2};
e,f adjacent:   {w_h(e,f):h in Omega}={1,1,1,3}.        (21)
```

There are three pairs of the first type and twelve of the second.  Thus
orbit symmetrization recovers only the ordinary adjacent/opposite
association scheme on `E`.  It forgets the marked state and the placement of
inertia.  In particular the same static data `(21)` occurs with row `000`
and with row `110`; one must retain the local monodromy action or the
normalized divisor row.

THM-2769 is the sharp affine hostile.  Its generic group is full `S4`, its
pair-sum product and six-root discriminant are squares, and the same modular
`W(D3)` frame is present.  At `t=0`, however,

```text
nu_D=110.                                                (22)
```

Equations `(14)` and `(18)--(19)` then pair all four marked colourings while
fixing the binary pair `{a_12,b_12}`.  Thus the complete marked-weight
construction coexists with failure of quasi-etaleness; only a separately
supplied inertia-free extension makes it an affine positive tool.

## 6. Complete divisor rows do not reconstruct the unramified torsor

Criterion `(15)` tests ramification of a **given** Kummer plane; it does not
identify that plane.  This loss persists with full `S3` covariance.  Work on
the four-torus

```text
R=C[x_1^+-,x_2^+-,x_3^+-,y_1^+-,y_2^+-,y_3^+-]/
  (x_1x_2x_3-1, y_1y_2y_3-1),                         (22a)
```

where `+-` means that each variable is inverted, and let `S3` permute the
three indices simultaneously.  In `R*/R*2`, define

```text
W_x=<[x_1],[x_2],[x_3]>,       W_y=<[y_1],[y_2],[y_3]>. (22b)
```

Both are rank-two standard `F2[S3]` planes because each triple has product
one.  Yet

```text
R*/R*2=<[x_1],[x_2],[y_1],[y_2]>=F2^4,
W_x intersect W_y=0,                                  (22c)
```

so they are distinct even after changing the `V4` basis.  All six defining
classes are units.  Consequently every height-one divisor row for both
planes is zero, while the connected `V4` covers are nonisomorphic over the
identity of `R`.  Taking `e_1=0` and `c=1` in `(12c)` gives two explicit
`S3`-invariant quartics with the same common sum, abstract marked frame,
chamber choices, and complete residue collection, but different Kummer
torsors.

Thus

```text
(marked D3 fibre, chamber, e_1, every divisor row)     (22d)
```

does not determine the cover.  The missing coordinate is the embedded
equivariant class

```text
W -> H^1_et(R_reg,mu_2),                               (22e)
```

equivalently the actual squareclass plane and its unit or class-group
realization.  Divisor rows record only the residues of `(22e)` and have a
large unramified kernel.

## 7. Connection contract and exact scope

The proved candidate connection is

```text
source:       one marked D3 determinant-weight colouring;
target:       one point of the quartic V4 root-state torsor;
map:          recover the point omitted by the weight-three triangle;
preserved:    S4/A4 action, stabilizer, and codimension-one inertia;
destroyed
 by orbiting: the marked point and inertia placement;
sidecar:      the embedded equivariant Kummer H1-class, together with its
              normalized divisor rows and unit/class-group realization;
hostiles:     THM-2769's full-S4 row 110 (ramified), and the distinct
              unit planes W_x,W_y in (22b) (unramified). (23)
```

Run

```bash
python 04-computation/marked_d3_weight_torsor_descent_thm2780.py
python -O 04-computation/marked_d3_weight_torsor_descent_thm2780.py
python3 04-computation/marked_d3_weight_torsor_descent_hostile_thm2780.py
python3 -O 04-computation/marked_d3_weight_torsor_descent_hostile_thm2780.py
```

The integer/set-only companion uses explicit exceptions and no truth-bearing
assertions.  It checks all `24` sheet permutations, all four marks, all
fifteen root-line pairs, both `S4/A4` stabilizers, simple transitivity of
`V4`, every even sign word on the half-Hadamard states and six root lines,
all sixty marked determinant frames and their complete mod-two kernels, and
the complete symmetrized profile `(21)`.  The independent companion
additionally exhausts retained-chamber chirality, checks the explicit
weighted self-converse witness, reconstructs `(12c)`, and verifies the twin
unit planes `(22b)--(22c)`.  Both normal/optimized pairs byte-match their
stored transcripts, with zero Python `assert` nodes.

```text
PROVED HERE:              the four marked weightings are the V4 torsor;
                          exact S3/C3 stabilizers and field degrees;
                          explicit inverse quartic descent map;
                          chamber-free weighted self-converse boundary;
                          local fixed-colouring iff zero parity row;
                          fixed binary root-line pair for each nonzero row;
                          determinant-two null-word alignment;
                          determinant-three mod-two-unit stopping rule;
                          complete orbit-symmetrized no-go;
                          THM-2769 ramified hostile;
                          unramified equivariant unit-plane twins.

NOT PROVED:               a global marked section from quasi-etaleness;
                          vanishing of the Kummer boundary rows;
                          exclusion of A4 or S4 Keller monodromy;
                          a new Keller/Jelonek realization;
                          JC(2), DC(2), Graceful Tree, or LRC(14).       (24)
```

QED.
