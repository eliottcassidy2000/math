---
id: THM-2777
title: "Marked D3 six-root determinant tournament and binary/ternary edge spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Marking one
  even-sign body diagonal h in D3 identifies its
  stabilizer with W(A2)=S3.  After choosing an A2 chamber and volume
  orientation, the six positive D3 root lines form a tie-free determinant
  tournament.  Its fifteen edge magnitudes are nine 1s, three 2s, and three
  3s: weight two is exactly the three opposite K4-edge matches, while weight
  three is exactly the three A2 spanning-tree pairs in the triangle opposite
  the marked half-Hadamard state.  The corresponding frame cokernels are
  trivial, Z/2, and Z/3.  Forgetting root orientations leaves only a switching
  class; no canonical unmarked tournament or graceful/Keller/LRC conclusion
  follows.
source: root/marked-d3-six-root-tournament-2026-07-28
audit: >
  thm2777-hostile-audit/2026-07-28 independently reconstructed the Weyl
  orbit, K4 dictionary, all 15 determinants and Smith forms, tournament
  census, chamber and reorientation gauges, normal/-O/stored replay, LF
  hashes, and documentation gate: ACCEPT after the three installed scope
  repairs.
depends_on:
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
related:
  - THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects
  - THM-2775-modular-s4-to-weyl-d3-generator-frame-and-affine-parity-blindness
script: 04-computation/marked_d3_six_root_tournament_thm2777.py
output: 05-knowledge/results/marked_d3_six_root_tournament_thm2777.out
script_sha256: 4b02d65554fa2c941beb83f6926e2c6c9f2818c8f5527657eb2e8e1ef5c06ed3
output_sha256: ca2ab9c0f4b82b1b608513a7962dbda82d881633e26b09ecd524aeb528fdffaa
hash_basis: LF-normalized bytes
---

# THM-2777 -- a marked `D3` frame gives a genuine six-root tournament

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

There is a natural six-object tournament in the rank-three binary/ternary
frame, but it is not gauge-free.  The vertices are the six `D3` root lines.
One must first mark a half-Hadamard state, orient the orthogonal `A2` roots by
a chamber, and orient volume.  With those coordinates retained, every pair
has a nonzero integral determinant and the edge magnitude is exactly its
selected-frame Smith index.

## 1. The marked root frame

Work in the oriented lattice `Z^3` with basis `e_1,e_2,e_3` and mark

```text
h=e_1+e_2+e_3.                                          (1)
```

Choose the standard chamber `1<2<3` of the `A2` subsystem orthogonal to
`h`.  Use the six positive `D3` roots

```text
a_12=e_1-e_2,    b_12=e_1+e_2,
a_13=e_1-e_3,    b_13=e_1+e_3,
a_23=e_2-e_3,    b_23=e_2+e_3.                          (2)
```

The `a_ij` are exactly the three `D3` root lines orthogonal to `h`, one
orientation from each line; they form the positive roots of `A2`.  The
`b_ij` are oriented canonically by `b_ij dot h>0`.

The even-sign Weyl group is

```text
W(D3)=V4 semidirect S3,                                 (3)
```

of order `24`.  Its orbit of `h` is

```text
(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1),              (4)
```

and the diagonal `V4` acts simply transitively on `(4)`.  If a signed
permutation fixes `h`, all three row signs are positive, so

```text
Stab_(W(D3))(h)=S3=W(A2).                               (5)
```

Thus marking `h` is literally choosing one state of the `V4` torsor, and
the surviving ternary coordinate symmetry is its `A2` Weyl stabilizer.

## 2. A tie-free determinant tournament

For distinct roots `alpha,beta` in `(2)`, orient the pair by

```text
alpha -> beta  iff  det(alpha,beta,h)>0,                (6)
```

and give the arc the positive integer weight

```text
w(alpha,beta)=|det(alpha,beta,h)|.                       (7)
```

Direct evaluation of the fifteen pairs gives no zero.  Hence `(6)` is a
genuine tournament, not a relation with discarded ties.  In the vertex
order

```text
(a_12,b_12,a_13,b_13,a_23,b_23),                        (8)
```

its outdegrees are

```text
(3,2,3,4,1,2),          sorted=(1,2,2,3,3,4).           (9)
```

It has `6` directed triangles and `23` directed Hamilton paths.  In the
fixed gauge `(1)--(2)` its automorphism group is trivial.  These last
statistics record the displayed oriented object; they do not classify a
six-vertex tournament and are not invariants of an unoriented root-line set.
In particular, the nonisomorphic converse has the same displayed census.

## 3. The weights are the marked `K4` edge geometry

Use the half-Hadamard states

```text
h_1=( 1, 1, 1), h_2=( 1,-1,-1),
h_3=(-1, 1,-1), h_4=(-1,-1, 1).                         (10)
```

Each `D3` root line is one edge of the complete graph on these four states:

```text
a_12 <-> 23,   b_12 <-> 14,
a_13 <-> 24,   b_13 <-> 13,
a_23 <-> 34,   b_23 <-> 12.                            (11)
```

Indeed each difference `h_i-h_j` is twice the displayed root, up to sign.
The marking `h=h_1` therefore splits the six edges into the incident star
`{b_12,b_13,b_23}` and the opposite triangle
`{a_12,a_13,a_23}`.

The complete determinant classification is

```text
w=2  iff the two K4 edges are opposite
     iff {alpha,beta} is one of
        {a_12,b_12}, {a_13,b_13}, {a_23,b_23};

w=3  iff both edges lie in the triangle opposite h_1
     iff alpha,beta are distinct A2 positive roots;

w=1  in the other nine cases.                           (12)
```

So the weighted spectrum is exactly

```text
1^9 2^3 3^3.                                            (13)
```

This refines the six-edge `S4` action: the ordinary matching quotient sees
the three weight-two opposite pairs, while marking one torsor state exposes
the three weight-three pairs in its opposite triangle.

## 4. Every arc weight is its cyclic frame cokernel

For an arc pair form the row matrix

```text
M_(alpha,beta)=[alpha; beta; h].                         (14)
```

Every matrix `(14)` has a unit entry and a unit `2 x 2` minor.  Its Smith
form is consequently

```text
diag(1,1,w(alpha,beta)).                                 (15)
```

The nine weight-one frames are unimodular.  For an opposite pair
`{a_ij,b_ij}`, if `k` is the remaining coordinate then

```text
h-b_ij=e_k.                                              (16)
```

Thus integer row operations turn `(14)` into the pure binary frame
`[e_i-e_j;e_i+e_j;e_k]`, with quotient `Z/2`.

For any two distinct `A2` roots, their rows generate the root lattice

```text
Q(A2)={z in Z^3:sum z_i=0},                              (17)
```

while adjoining `h` gives

```text
rowspan_Z M={z in Z^3:sum z_i=0 mod 3}.                 (18)
```

Hence the three weight-three frames have quotient `Z/3`, explicitly by
the coordinate-sum map.  Equivalently, after quotienting `Z h`, the ambient
lattice is the `A2` weight lattice and `(17)` is the root lattice, so this
is the familiar `P(A2)/Q(A2)=Z/3` defect.

The weighted tournament therefore records all one-marked-state full frame
cokernels at once:

```text
9 trivial,              3 binary Z/2,              3 ternary Z/3.       (19)
```

## 5. The orientation gauge is load-bearing

The determinant sign in `(6)` needs three choices that the unmarked root
arrangement does not contain by itself:

1. the state `h` from the four-element torsor `(4)`;
2. an `A2` chamber orienting the three roots orthogonal to `h`;
3. an orientation of the ambient volume form.

Permuting the ordered coordinate chamber by an even element of `S3`
relabels the same tournament.  An odd permutation reverses every arc; this
is the volume-orientation change.  More generally, changing the chosen
representative of one root line reverses every incident arc.  With `h` and
the ambient volume orientation fixed, the six labelled unoriented lines
therefore give a `32`-member reorientation orbit.  Reversing volume gives
the distinct labelled converse orbit; the companion checks that it is not
one of those `32` switchings.  The weight function `(12)` survives all of
these changes, but the individual scores and directed cycles do not.

This is the sharp stopping rule.  The construction supplies an intrinsic
binary relation only on the **marked and oriented** frame.  It does not turn
the six roots into a canonical tournament before those sidecars are chosen.

## 6. Exact verification and scope

Run

```bash
python 04-computation/marked_d3_six_root_tournament_thm2777.py
python -O 04-computation/marked_d3_six_root_tournament_thm2777.py
```

The integer-only companion uses explicit exceptions and no truth-bearing
Python assertions.  It enumerates `W(D3)`, its four-state orbit and
stabilizer; checks the half-Hadamard/K4 dictionary; evaluates all fifteen
determinants and Smith minors; classifies the weight-two and weight-three
arcs; computes the tournament census; tests all six chamber permutations;
and enumerates all `32` labelled root-line switchings.  Normal and optimized
runs byte-match the stored transcript.

```text
PROVED HERE:              marked-state W(D3)/W(A2) decomposition;
                          tie-free six-root determinant tournament;
                          exact 1^9 2^3 3^3 edge spectrum;
                          K4 opposite/marked-triangle classification;
                          cyclic Smith quotients 1, Z/2, Z/3;
                          tournament score/cycle/path census;
                          chamber reversal and switching-class boundary.

NOT PROVED:               a chamber-free canonical tournament;
                          a modular free-product action on these six roots;
                          a signed graceful coefficient selector;
                          a Keller or quartic-monodromy exclusion;
                          Graceful Tree, JC(2), DC(2), or LRC(14).       (20)
```

QED.
