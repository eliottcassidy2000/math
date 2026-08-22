---
id: THM-3206
title: "Heterogeneous factorial exterior reflection groupoid and fixed-plane holonomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every off-discriminant normalized factorial exterior block factors as
  iota F pi through one parameter-independent two-plane, with pi iota=I and
  F a projective involution.  Hence arbitrary heterogeneous products retain
  the same kernel, image plane, and complementary chart atlas.  Their
  internal transport is the ordered PGL2 product of the reflections.  A
  two-block product is scalar exactly when the two parameter pairs agree;
  its eigenvalue discriminant is sixteen times the binary resultant of the
  two underlying quadratics.
audit: >
  The assertion-independent integer companion exhausts 300 off-wall
  one-block systems and 32,304 ordered pairs over p=3,5,7,11,13, including
  3,316 affine eigenline defects, 3,616 binary fixed-line identities, every
  pair invariant, a directly constructed Sylvester resultant, and the sharp
  scalar classification.  It also checks 119,256 heterogeneous words through
  lengths 5,4,3 over p=3,5,7 and two exact p=5 hostiles.  Normal and optimized
  replay agree with the stored transcript.  An independent immutable audit
  rederived the split carrier, arbitrary-word atlas, binary-root involution,
  Sylvester resultant, commutator and parabolic walls, including the
  projective-infinity boundary; replayed both modes; and accepted the hashes
  and scope.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on:
  - THM-3191-factorial-block-exterior-clifford-law-and-global-carry-smith-profile
related:
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
script: 04-computation/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.py
output: 05-knowledge/results/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.out
script_sha256: 9eed2306e6fc34735e32bfccfb6c2804dfee2181e4b5908aca70243f77c26217
output_sha256: 72b1ec241cf924b18c634ce9fbd7b3b54889a4b7c4ac0ce94b0ef4ad76a53b75
hash_basis: LF-normalized bytes
---

# THM-3206 -- heterogeneous factorial exterior reflection groupoid and fixed-plane holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3191 closes repeated factorial blocks with fixed residues by the cubic
law `D^3=s^2 Delta D`.  Parameter variation does not destroy the exterior
carrier.  More rigidly, every normalized block passes through the same split
two-plane.  What variation creates is genuine projective holonomy *inside*
that plane.

This separates two questions which the fixed-parameter law conflates:

```text
existence of one global exterior atlas:       always, off Delta=0;
collapse to the fixed D/D^2 parity orbit:     only for a constant block.   (1)
```

## 1. One split carrier for every parameter pair

Let `p` be odd.  For `s in F_p^*` and `v in F_p`, put

```text
Delta(s,v)=1-4sv !=0.                                         (2)
```

In the THM-3191 exterior basis `(w_01,w_02,w_12)`, remove the nonzero
quadratic-character and carry scalar from its normalized block and write

```text
E(s,v)=
[0    2s+1             -1  ]
[0   -(1+2v)            2v ]
[0   -2(s+v+1)         2v+1].                              (3)
```

Define

```text
pi(y_0,y_1,y_2)=(y_1,y_2),
iota(x,y)=(x-y,x,y),

F(s,v)=
[-(1+2v)             2v ]
[-2(s+v+1)          2v+1].                                (4)
```

Then direct multiplication gives the split factorization

```text
pi iota=I_2,                    E(s,v)=iota F(s,v) pi.      (5)
```

Moreover,

```text
tr F=0,                 det F=-Delta,       F^2=Delta I_2. (6)
```

Thus `[F(s,v)]` is an involution in `PGL_2(F_p)`.  The word
*reflection* below means precisely this projective order-two statement; no
positive or metric reflection is being asserted.

### Proof

The last two rows of `(3)` are `F pi`.  The first row is the first row of
`F` minus its second row, so `(5)` follows.  The trace and determinant in
`(6)` are immediate, and Cayley--Hamilton then gives the square.  QED.

### The reflections are the quadratic root involutions

Homogenize the reduced factorial quadratic as

```text
q_(s,v)^h(X,Z)=vX^2-XZ+sZ^2.                              (7)
```

The projective change of coordinates

```text
[X:Z] |-> [Z:Z+X] in P^1                                 (8)
```

identifies its roots with the fixed lines of `[F(s,v)]`.  Indeed,

```text
det([Z,Z+X]^t, F(s,v)[Z,Z+X]^t)=-2q_(s,v)^h(X,Z).         (9)
```

On the affine chart `Z=1`, this refines to the defect identity

```text
F(s,v)[1,1+x]^t
 =q'(x)[1,1+x]^t-2q(x)[0,1]^t,

q(x)=vx^2-x+s,                    q'(x)=-1+2vx.            (10)
```

Thus a root line has eigenvalue `q'(x)` in the affine normalization.  The
two fixed lines collide exactly at `Delta=0`.  When `v=0`, the second root is
the projective point at infinity, and `(9)` still applies without a degree
convention.

## 2. Arbitrary heterogeneous words keep the same flags

Choose any word of off-wall pairs

```text
(s_1,v_1),...,(s_K,v_K)                                  (11)
```

and arbitrary nonzero block scalars `rho_i`.  Put

```text
D_i=rho_i E(s_i,v_i),
Q_K=F(s_K,v_K)...F(s_1,v_1),
rho=rho_K...rho_1.                                        (12)
```

Repeated use of `pi iota=I_2` gives

```text
D_K...D_1=rho iota Q_K pi,                                (13)

det Q_K=(-1)^K product_(i=1)^K Delta(s_i,v_i) !=0.        (14)
```

Consequently every nonempty product in `(13)` has

```text
rank=2,
right kernel=span(w_01),
image=H:={y_0-y_1+y_2=0},
left kernel=span(1,-1,1).                                 (15)
```

The restriction to `H` is the invertible transport `rho Q_K` in the split
coordinates `(y_1,y_2)`.  In particular the two projective charts

```text
y_1!=0,                         y_2!=0                     (16)
```

cover every heterogeneous output.  The parameter-dependent image planes do
not become transverse: they are exactly the same plane.

For an actual list of local factorial blocks, THM-3191 has

```text
p^(-h_i) Lambda^2 B_i ==rho_i E(s_i,v_i)          (mod p), (17)
```

where its carried specialization is `rho_i=-u_i chi_i s_i`.  Therefore, if
`H_tot=sum h_i`, compound multiplicativity gives the formal heterogeneous
concatenation

```text
p^(-H_tot)Lambda^2(B_K...B_1)==rho iota Q_K pi     (mod p). (18)
```

Equation `(18)` is a theorem about concatenating locally typed blocks.  A
single factorial moment system has fixed `s,v`; arbitrary variation in `(11)`
is not silently asserted to arise from one such system.

## 3. Exact two-reflection resultant holonomy

For two pairs abbreviate

```text
F_i=F(s_i,v_i),             Delta_i=1-4s_i v_i,
Q=F_2F_1.                                                   (19)
```

The complete characteristic data are

```text
tau=tr Q=-2(2s_1v_2+2s_2v_1-1),
delta=det Q=Delta_1 Delta_2,

Q^2-tau Q+delta I_2=0.                                    (20)
```

Let

```text
R_12=Res_(X:Z)(q_(s_1,v_1)^h,q_(s_2,v_2)^h)              (21)

 =s_1^2v_2^2-2s_1s_2v_1v_2+s_1v_1-s_1v_2
   +s_2^2v_1^2-s_2v_1+s_2v_2.
```

This is the determinant of the literal `4x4` degree-two Sylvester matrix,
including roots at infinity.  The pair discriminant and additive commutator
are exactly

```text
tau^2-4delta=16R_12,
det(F_2F_1-F_1F_2)=-16R_12.                               (22)
```

Consequently `R_12=0` iff the two binary quadratics have a common
projective root over the algebraic closure, iff `Q` has a repeated
eigenvalue.  If the parameter pairs are distinct this is a non-scalar
parabolic projective holonomy.  If they agree, the return is scalar.

Most sharply,

```text
F_2F_1 is scalar     iff     (s_2,v_2)=(s_1,v_1).         (23)
```

### Proof of the sharp return criterion

If the pairs agree, `(6)` gives `F_2F_1=Delta_1 I`.  Conversely suppose
`F_2F_1=cI`.  Since `F_1^(-1)=Delta_1^(-1)F_1`, one has

```text
F_2=(c/Delta_1)F_1.                                       (24)
```

The two entries in the first row of every normalized `F(s,v)` sum to `-1`.
Thus the scalar in `(24)` is one.  Equality of the upper-right entries gives
`v_2=v_1` because `p` is odd, and equality of the lower-left entries then
gives `s_2=s_1`.  This proves `(23)`.  Expanding the trace, determinant, the
Sylvester determinant, and the commutator proves `(20)--(22)`.  QED.

Thus the constant-parameter case is the scalar-return diagonal of the
heterogeneous parameter square.  Away from that diagonal, even two blocks
already leave THM-3191's scalar even return.  Ordered words of the projective
reflections carry a genuine `PGL_2` holonomy on `P(H)`.

## 4. Two exact boundary controls

Over `F_5`, take

```text
(s_1,v_1)=(1,1),             Delta_1=2,
(s_2,v_2)=(3,1),             Delta_2=4.                   (25)
```

Then

```text
F_1=[2 2;4 3],               F_2=[2 2;0 3],

Q=F_2F_1=[2 0;2 4],          tr Q=1,       det Q=3,
Q^2=[4 0;2 1].                                               (26)
```

Neither `Q` nor `Q^2` is scalar.  The block discriminants and the pair
eigenvalue discriminant are nonzero, so this is not a wall artifact.
Including THM-3191's character scalars merely multiplies `Q` by a unit and
cannot restore the scalar return.

The projective-root wording after `(22)` is load-bearing.  Still over `F_5`,
take `(s_1,v_1)=(1,0)` and `(s_2,v_2)=(2,0)`.  Their affine linear roots are
distinct, but their binary quadratics share `[1:0]`.  Accordingly,

```text
F(2,0)F(1,0)=[1 0;2 1]                                   (27)
```

is non-scalar unipotent and `R_12=0`.  The simultaneous leading-degree drop
does not refute the resultant law; it realizes its root at infinity.

## 5. Holotopy meaning and the remaining selector debt

The exact object supplied by `(13)` is a rank-two local system with a fixed
ambient carrier and reflection-valued transition words:

```text
fixed data:       kernel line, image plane H, conormal, two charts;
moving data:      the ordered PGL_2 word Q_K;
lost simplifier:  the fixed-parameter D/D^2 parity orbit.                (28)
```

This is better news for adaptive transport than a transverse-plane hostile:
there is no Cech gluing problem for the exterior plane itself.  Equations
`(21)--(22)` also identify a genuine resultant gate: common projective roots
of the block quadratics are exactly parabolic pair holonomy.  But this is the
resultant of the *input quadratics in `x`*, not yet the coefficient-degree
PRS of the factorial moment polynomials in `v`.

One still needs a proved map from a PRS-selected Schur-complement coordinate
to a vector or covector on `H`, together with permission to switch between
`(16)`.  THM-3186 shows that one visible exit-time coordinate can cancel
while the complementary transverse pair remains nonzero.  THM-2624 is the
complementary warning that multi-chart tomography need not produce a
physical transport.

The connection contract is therefore

```text
source:       the fixed exterior carrier H with holonomy Q_K;
target:       a selected coefficient-degree PRS pivot;
map needed:   an index-compatible PRS-to-H vector/covector assignment;
preserved:    both exterior charts and all ordered block holonomy;
destroyed:    the fixed two-state parity reduction;
cheapest test: identify one PRS leading row with a nonzero covector on H
               and check its covariance under one heterogeneous block.  (29)
```

No such map is constructed here.  No arbitrary-support `NC(2)`, Gaussian
Moment, Euclidean-depth staircase, or `LRC(14)` conclusion follows.

## 6. Exact evidence

Run

```text
python 04-computation/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.py
python -O 04-computation/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.py
```

and compare LF-normalized bytes with the declared output.  The companion is
assertion-independent and uses integer arithmetic modulo `p` only.  It
exhausts all `300` off-wall one-block systems, `3,316` affine defect
identities `(10)`, `3,616` binary fixed-line identities `(9)`, and all
`32,304` ordered pairs over `p=3,5,7,11,13`.  For every pair it constructs
the `4x4` Sylvester matrix and checks `(20)--(23)`.  It also checks all
`119,256` words through respective maximum lengths `5,4,3` over `p=3,5,7`,
including the full physical character scalar, factorization `(13)`,
determinant, kernel, conormal, and image-plane chart.  Finally it checks every
entry of `(25)--(27)`.  Normal and optimized executions agree with the stored
transcript.

**QED.**
