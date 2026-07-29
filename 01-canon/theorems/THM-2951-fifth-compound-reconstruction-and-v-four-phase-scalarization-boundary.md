---
id: THM-2951
title: "Fifth-compound reconstruction and V-four phase scalarization boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  On a reduced real length-six algebra with three
  conjugate pairs, the sixth, fifth, and balanced third compounds of
  multiplication carry respectively the nonnegative norm, the
  norm-scaled inverse operator, and the eight half-system weights
  underlying the V_4 trace frame.  The full invertible fifth-compound
  operator reconstructs multiplication and hence loses no phase
  information.  However, there is no nonzero linear map from the fifth
  compound to the balanced third compound equivariant under the full
  signed-pair group (C_2^3) semidirect S_3.  No nonzero pure
  two-form contraction lands in the balanced sector; after balanced
  projection the contraction family has rank 12 and the three
  within-pair area forms form its kernel.  The norm, fifth-compound
  energy, and one selected cofactor can agree while all four V_4
  traces change.  No SFC(4), LRC, knot, modular-group, Jacobian, or
  Dixmier consequence is claimed.
source: root-gmc-compound-torsor-ladder-2026-07-29
audit: Pending independent hostile audit.
depends_on:
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
related:
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
script: 04-computation/gmc_fifth_compound_v4_phase_boundary_thm2951.py
output: 05-knowledge/results/gmc_fifth_compound_v4_phase_boundary_thm2951.out
script_sha256: 22186ea829aa2f5f7f4ac23d138def3bdd05772a9a124d7a2d42d7f90b359573
output_sha256: 7d62ae55f82099ebced31fab669ea8465aae482365ff1debf2872a45c68e186f
hash_basis: LF-normalized bytes
---

# THM-2951 -- fifth-compound reconstruction and V-four phase scalarization boundary

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. The three compound levels

Let

```text
V=E_1 (+) E_2 (+) E_3,             dim_R E_i=2,       (1)
```

be the reduced real length-six algebra carrier from THM-2947.  After
complexification, `E_i` splits into the two primitive lines exchanged
by conjugation.  Let multiplication by a real unit `f` be `M`, with
complex eigenvalues

```text
z_1,zbar_1,z_2,zbar_2,z_3,zbar_3.                    (2)
```

Put

```text
N=det(M)=product_i |z_i|^2>0.                         (3)
```

Then the top compound is the scalar

```text
Lambda^6 M=N.                                         (4)
```

There is a canonical isomorphism

```text
kappa:Lambda^5 V -> det(V) tensor V^*,
kappa(v_1 wedge ... omit v_a ... wedge v_6)
 =(-1)^(a-1)(v_1 wedge ... wedge v_6) tensor v_a^*.  (5)
```

Under `(5)`,

```text
C_5:=Lambda^5 M=N(M^-1)^*.                            (6)
```

Thus its eigenvalues are

```text
N/z_i,                         N/zbar_i.               (7)
```

Finally the pair decomposition canonically defines the balanced
middle sector

```text
B=E_1 wedge E_2 wedge E_3             subset Lambda^3 V,
dim_R B=8.                                                (8)
```

After temporarily choosing one primitive line in each pair, its eight
complex eigenweights are

```text
h_epsilon=product_i z_(i,epsilon_i),
epsilon in F_2^3,       z_(i,0)=z_i, z_(i,1)=zbar_i.  (9)
```

Complementary weights are conjugate.  Their four sums

```text
T_[epsilon]=h_epsilon+h_(epsilon+(1,1,1))             (10)
```

are precisely the four real half-system traces of THM-2950.  Notice
that `B` is eight-dimensional; conjugation pairs its eight spectral
weights but does not turn `B` itself into a four-dimensional quotient.

## 2. The full fifth compound is faithful

Taking determinants in `(6)` gives

```text
det(C_5)=N^5.                                         (11)
```

Over the reals, `N` is the unique real fifth root of `(11)`, and

```text
M^*=N C_5^-1.                                         (12)
```

Consequently the **full operator** `C_5` reconstructs `M`.  Once the
pair decomposition `(1)` is retained, it therefore reconstructs
`Lambda^3 M|_B` and the unordered four-trace frame `(10)`.

The same fact is visible spectrally.  If

```text
c_(i,epsilon)=N/z_(i,epsilon),                        (13)
```

are the paired eigenvalues of `C_5`, then

```text
h_epsilon
 =N^3/product_i c_(i,epsilon_i).                      (14)
```

Thus the obstruction below is not information-theoretic loss from the
full fifth compound.  THM-2947 and THM-2949 use a positivity readout or
one matrix coefficient of it, not the entire invertible operator.

## 3. No signed-pair-equivariant linear descent

Choose temporary bases `e_(i,0),e_(i,1)` of the three pairs and let

```text
W=(C_2^3) semidirect S_3                              (15)
```

act by flipping the two entries in any pair and permuting the three
pairs.  This is the signed-pair automorphism group used in THM-2950.
Both `Lambda^5 V` and `B` are intrinsic `W`-modules, despite the
temporary bases.

Let `H=C_2^3`, and write `chi_i` for the sign character of the `i`th
flip.  With

```text
s_i=e_(i,0)+e_(i,1),       a_i=e_(i,0)-e_(i,1),       (16)
```

one has

```text
Lambda^5 V|_H
 =3 chi_123 (+) chi_12 (+) chi_13 (+) chi_23,         (17)
B|_H
 =direct_sum_(S subset {1,2,3}) chi_S.                (18)
```

Indeed `Lambda^5 V=det(V) tensor V^*`; `det(V)` has
`H`-character `chi_123`, while a permutation of two two-dimensional
pair blocks acts trivially on `det(V)`.

Any `W`-map must preserve the characters in `(17)--(18)`.  On the
`chi_123` block, the multiplicity space in `Lambda^5 V` is the
three-point permutation representation of `S_3`, whereas the
`chi_123` line in `B`, spanned by

```text
a_1 wedge a_2 wedge a_3,                              (19)
```

is the sign representation.  There is no map between them.  On the
orbit of `chi_23`, its stabilizer contains the transposition `(23)`.
That transposition fixes the `chi_23` line of `Lambda^5 V`, but
negates the line

```text
s_1 wedge a_2 wedge a_3                               (20)
```

in `B`.  The same holds cyclically.  Hence

```text
Hom_W(Lambda^5 V,B)=0.                                (21)
```

Equivalently, the exact character inner product is zero:

```text
(1/48) sum_(w in W)
  chi_(Lambda^5 V)(w) chi_B(w)=0.                     (22)
```

Equation `(21)` is the precise naturality obstruction.  It does not
contradict the nonlinear reconstruction `(11)--(14)`.

## 4. The contraction boundary is sharper

For a two-form

```text
beta in Lambda^2 V^*,                                 (23)
```

let `i_beta:Lambda^5 V->Lambda^3 V` be contraction.  A two-form, not a
bivector, contracts multivectors under this convention.

There is no nonzero `beta` for which

```text
i_beta(Lambda^5 V) subset B.                          (24)
```

To see this, use the six basis vectors from Section 3.  If a nonzero
coefficient of `beta` uses the two entries of one pair, contract a
five-vector missing a coordinate in another pair.  The resulting
three-vector omits the first pair and doubles a different pair, so it
is outside `B`.  If the coefficient uses entries from two different
pairs, contract a five-vector missing the mate of the first entry.
The result again has pair multiplicities `(0,1,2)`, up to permutation.
For a fixed input and output the removed pair of coordinates is
unique, so these offending terms cannot cancel.  Every coefficient of
`beta` is therefore zero.

The pair decomposition does supply a canonical projection

```text
pi_B:Lambda^3 V->B.                                   (25)
```

For the projected contraction, decompose

```text
Lambda^5 V
 =direct_sum_i E_i tensor det(E_j) tensor det(E_k),
{i,j,k}={1,2,3}.                                      (26)
```

On the `i`th summand, only the cross-pair component

```text
beta_jk in E_j^* tensor E_k^*                         (27)
```

can survive `pi_B`.  Hence

```text
rank[beta |-> pi_B i_beta]=3*2*2=12,
kernel=direct_sum_i Lambda^2 E_i^*,       dim=3.       (28)
```

In particular, a sum of natural within-pair area forms is in the
kernel, not a solution.  A nonzero projected descent requires a
chosen cross-pair contraction tensor.  An oriented half-system selects
a weight in `(9)`, but by itself does not canonically supply the
two-form `(27)`.

Nor is every linear map `Lambda^5 V->B` a contraction: the full Hom
space has dimension `6*8=48`, while the projected contraction family
has dimension `12`.  Thus contraction is a particular sidecar, not a
classification of arbitrary linear maps.

## 5. Sharp scalarization hostile

Let

```text
M_0=I_6,
M_1=J (+) I_2 (+) I_2,       J=[[0,-1],[1,0]].        (29)
```

Both are real pair-diagonal units with

```text
N=1.                                                  (30)
```

Both fifth compounds are orthogonal, so the squared fifth-compound
energy from THM-2947 is

```text
||Lambda^5 M_0||_F^2
 =||Lambda^5 M_1||_F^2=6.                             (31)
```

A diagonal fifth-cofactor coordinate belonging to either unchanged
pair equals `1` for both.  Nevertheless their eigenvalue triples may
be written

```text
(1,1,1),                         (i,1,1),             (32)
```

and their four trace vectors are

```text
(2,2,2,2),                       (0,0,0,0).            (33)
```

Thus `N`, the positive energy, and one nonzero cofactor can all agree
while the whole `V_4` trace frame changes.  The full operators in
`(29)` are different, exactly as reconstruction predicts.

## 6. Singular and nonreduced boundaries

The reconstruction formula requires `C_5` invertible.  In general,

```text
rank M<=4  ==> Lambda^5 M=0,                           (34)
rank M=5   ==> rank Lambda^5 M=1.                      (35)
```

THM-2947's no-real-support parity excludes rank `5`: on its singular
physical branch the entire fifth compound vanishes.  If a real local
factor is allowed, `diag(0,1,1,1,1,1)` realizes `(35)`, so neither the
binary rank gap nor `(12)` survives.

The exterior-power identities and invertible reconstruction apply to
any six-dimensional real vector space.  The balanced sector, the
eight weights, and the `V_4` frame require the reduced three-pair
decomposition.  They have no asserted extension through collision of
the conjugate-pair idempotents.  No conclusion about SFC(4), LRC,
knots, the modular group, the planar Jacobian conjecture, or the
Dixmier conjecture follows from this representation-theoretic bridge.

## 7. Exact companion

Run

```text
python 04-computation/gmc_fifth_compound_v4_phase_boundary_thm2951.py
python -O 04-computation/gmc_fifth_compound_v4_phase_boundary_thm2951.py
```

The companion uses exact integer and rational arithmetic and explicit
`require` gates.  It checks:

1. the `48` signed-pair elements and the zero character inner product;
2. the `15`-rank unbalanced contraction constraint matrix;
3. projected-contraction rank `12` and kernel dimension `3`;
4. the canonical fifth-compound reconstruction on two controls;
5. all eight instances of the spectral identity `(14)`;
6. the scalarization hostile `(29)--(33)`; and
7. the singular ranks `(34)--(35)`.

Normal and optimized executions reproduce the stored transcript
byte-for-byte.
