---
id: THM-3206
title: "Heterogeneous factorial exterior reflection groupoid and fixed-plane holonomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  Every off-discriminant normalized factorial exterior block factors as
  iota F pi through one parameter-independent two-plane, with pi iota=I and
  F a projective involution.  Hence arbitrary heterogeneous products retain
  the same kernel, image plane, and complementary chart atlas.  Their
  internal transport is the ordered PGL2 product of the reflections.  A
  two-block product is scalar exactly when the two parameter pairs agree;
  an explicit p=5 pair therefore disproves heterogeneous two-state closure.
audit: >
  The assertion-independent integer companion exhausts 300 off-wall
  one-block systems and 32,304 ordered pairs over p=3,5,7,11,13, including
  every pair invariant and the sharp scalar classification.  It also checks
  119,256 heterogeneous words through lengths 5,4,3 over p=3,5,7 and the
  exact p=5 hostile.  Normal and optimized replay agree with the stored
  transcript.  Independent immutable audit is pending.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on:
  - THM-3191-factorial-block-exterior-clifford-law-and-global-carry-smith-profile
related:
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
script: 04-computation/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.py
output: 05-knowledge/results/factorial_heterogeneous_exterior_reflection_holonomy_thm3206.out
script_sha256: 547e841c259fab54f7ae352708b50d438468b11e985e0cf2e2fd656aac83002f
output_sha256: ddf9275eed048d9c7619d4989d3d3f2ac0340217971265c6a4074a59bf2e7e90
hash_basis: LF-normalized bytes
---

# THM-3206 -- heterogeneous factorial exterior reflection groupoid and fixed-plane holonomy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT
PENDING.**

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

## 2. Arbitrary heterogeneous words keep the same flags

Choose any word of off-wall pairs

```text
(s_1,v_1),...,(s_K,v_K)                                  (7)
```

and arbitrary nonzero block scalars `rho_i`.  Put

```text
D_i=rho_i E(s_i,v_i),
Q_K=F(s_K,v_K)...F(s_1,v_1),
rho=rho_K...rho_1.                                        (8)
```

Repeated use of `pi iota=I_2` gives

```text
D_K...D_1=rho iota Q_K pi,                                (9)

det Q_K=(-1)^K product_(i=1)^K Delta(s_i,v_i) !=0.        (10)
```

Consequently every nonempty product in `(9)` has

```text
rank=2,
right kernel=span(w_01),
image=H:={y_0-y_1+y_2=0},
left kernel=span(1,-1,1).                                 (11)
```

The restriction to `H` is the invertible transport `rho Q_K` in the split
coordinates `(y_1,y_2)`.  In particular the two projective charts

```text
y_1!=0,                         y_2!=0                     (12)
```

cover every heterogeneous output.  The parameter-dependent image planes do
not become transverse: they are exactly the same plane.

For an actual list of local factorial blocks, THM-3191 has

```text
p^(-h_i) Lambda^2 B_i ==rho_i E(s_i,v_i)          (mod p), (13)
```

where its carried specialization is `rho_i=-u_i chi_i s_i`.  Therefore, if
`H=sum h_i`, compound multiplicativity gives the formal heterogeneous
concatenation

```text
p^(-H)Lambda^2(B_K...B_1)==rho iota Q_K pi         (mod p). (14)
```

Equation `(14)` is a theorem about concatenating locally typed blocks.  A
single factorial moment system has fixed `s,v`; arbitrary variation in `(7)`
is not silently asserted to arise from one such system.

## 3. Exact two-reflection holonomy

For two pairs abbreviate

```text
F_i=F(s_i,v_i),             Delta_i=1-4s_i v_i,
Q=F_2F_1.                                                   (15)
```

The complete characteristic data are

```text
tau=tr Q=-2(2s_1v_2+2s_2v_1-1),
delta=det Q=Delta_1 Delta_2,

Q^2-tau Q+delta I_2=0.                                    (16)
```

The additive commutator records the repeated-eigenvalue wall exactly:

```text
det(F_2F_1-F_1F_2)=4delta-tau^2.                          (17)
```

Most sharply,

```text
F_2F_1 is scalar     iff     (s_2,v_2)=(s_1,v_1).         (18)
```

### Proof of the sharp return criterion

If the pairs agree, `(6)` gives `F_2F_1=Delta_1 I`.  Conversely suppose
`F_2F_1=cI`.  Since `F_1^(-1)=Delta_1^(-1)F_1`, one has

```text
F_2=(c/Delta_1)F_1.                                       (19)
```

The two entries in the first row of every normalized `F(s,v)` sum to `-1`.
Thus the scalar in `(19)` is one.  Equality of the upper-right entries gives
`v_2=v_1` because `p` is odd, and equality of the lower-left entries then
gives `s_2=s_1`.  This proves `(18)`.  Expanding the trace, determinant, and
commutator proves `(16)--(17)`.  QED.

Thus the constant-parameter case is a zero-curvature diagonal of the
heterogeneous parameter square.  Away from that diagonal, even two blocks
already leave THM-3191's scalar even return.  Ordered words of the projective
reflections carry a genuine `PGL_2` holonomy on `P(H)`.

## 4. Smallest explicit loss of the parity orbit

Over `F_5`, take

```text
(s_1,v_1)=(1,1),             Delta_1=2,
(s_2,v_2)=(3,1),             Delta_2=4.                   (20)
```

Then

```text
F_1=[2 2;4 3],               F_2=[2 2;0 3],

Q=F_2F_1=[2 0;2 4],          tr Q=1,       det Q=3,
Q^2=[4 0;2 1].                                               (21)
```

Neither `Q` nor `Q^2` is scalar.  Both discriminants are nonzero, so this is
not a wall artifact.  Including THM-3191's character scalars merely
multiplies `Q` by a unit and cannot restore the scalar return.

## 5. Holotopy meaning and the remaining selector debt

The exact object supplied by `(9)` is a rank-two local system with a fixed
ambient carrier and reflection-valued transition words:

```text
fixed data:       kernel line, image plane H, conormal, two charts;
moving data:      the ordered PGL_2 word Q_K;
lost simplifier:  the fixed-parameter D/D^2 parity orbit.                (22)
```

This is better news for adaptive transport than a transverse-plane hostile:
there is no Cech gluing problem for the exterior plane itself.  But it does
not yet solve the coefficient-degree PRS problem.  One still needs a proved
map from a PRS-selected Schur-complement coordinate to a vector or covector
on `H`, together with permission to switch between `(12)`.  THM-3186 shows
that one visible exit-time coordinate can cancel while the complementary
transverse pair remains nonzero.  THM-2624 is the complementary warning that
multi-chart tomography need not produce a physical transport.

The connection contract is therefore

```text
source:       the fixed exterior carrier H with holonomy Q_K;
target:       a selected coefficient-degree PRS pivot;
map needed:   an index-compatible PRS-to-H vector/covector assignment;
preserved:    both exterior charts and all ordered block holonomy;
destroyed:    the fixed two-state parity reduction;
cheapest test: identify one PRS leading row with a nonzero covector on H
               and check its covariance under one heterogeneous block.  (23)
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
exhausts all `300` off-wall one-block systems and all `32,304` ordered pairs
over `p=3,5,7,11,13`.  For every pair it checks `(16)--(18)`.  It also checks
all `119,256` words through respective maximum lengths `5,4,3` over
`p=3,5,7`, including the full physical character scalar, factorization
`(9)`, determinant, kernel, conormal, and image-plane chart.  Finally it
checks every entry of `(20)--(21)`.  Normal and optimized executions agree
with the stored transcript.

**QED (candidate pending independent immutable audit).**
