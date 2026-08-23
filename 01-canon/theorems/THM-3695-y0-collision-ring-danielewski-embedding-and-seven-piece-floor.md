---
id: THM-3695
title: "Danielewski embedding, line-jet and filtration floors, and a W003 ray-family nonentry for the y=0 collision ring"
status: >
  PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The y=0 collision ring of THM-3686 is a graded subalgebra of the squarefree
  exponent-two Danielewski Poisson ring with Sigma(b)=1-b^2; the ambient
  embedding into the source is Poisson, but the collision ring itself is not
  asserted bracket-closed.  Cited injective-line rigidity and the source jet
  force at least three retained weights in each output.  Repaired THM-3592
  then forces at least seven in total, so the first live cell is 3x4 up to
  exchange.  A cited reduced-height theorem forces one output to have minimal
  target-filtration degree at least 22.  The collision ring's extra b=0
  divisor closes one complete W003 anchor-20 placement family, including all
  odd scales left by THM-3613.  No Darboux pair or JC(2) counterexample is
  constructed.
source: root / planar-Jacobian counterexample long session, 2026-08-22
audit: >
  PASS after terminology repair and two independent valuation derivations.
  Audits checked injectivity, all Poisson signs, grading-support preservation,
  both source-line restrictions, the transverse origin jet, the repaired
  scope of THM-3592, the n=2 retained-zero boundary, and the all-scale W003
  contradiction.  Normal and optimized companion runs byte-match the stored
  754-gate transcript.  The false phrase "Poisson subalgebra" was removed: R
  is only a graded subalgebra of D, and bracket compatibility in D suffices.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3609-three-by-four-size-nine-euler-factor-nonentry
  - THM-3613-three-by-four-size-seven-ray-parity-gate
related:
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3691-y0-collision-ring-two-weight-darboux-no-go
  - THM-3693-y0-collision-ring-two-by-three-weight-darboux-no-go
  - THM-3694-one-pure-x-q-two-by-three-support-closure
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
  - "Guccione--Guccione--Horruitiner--Valqui, arXiv:2204.14178v1, reduced height at least 108."
script: 04-computation/jacobian_y0_danielewski_transfer_support_floor_thm3695.py
output: 05-knowledge/results/jacobian_y0_danielewski_transfer_support_floor_thm3695.out
script_sha256: 28c30bf1733e2270e2339207a593c1856d7e10cd17e9756c0b454ddeb4949bfb
output_sha256: e5a07fb645f73d626b0f5c100fe10c0db39df9b94ca2262f353b0167b3d5edf7
hash_basis: LF-normalized bytes
---

# THM-3695 -- the collision ring inherits stronger Danielewski gates

**PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  The original inheritance theorem is strengthened here by
a direct source-line/jet argument, a target-filtration floor, and the first
collision-ring-specific `3 x 4` ray-family nonentry.  These are necessary
conditions for a counterexample; no Darboux pair is constructed.

All rings are over `C`.  On the source plane with coordinates `(x,z)`, put

```text
b=1-x^2z,             c=x,             e=z(2-x^2z).    (1)
```

Then

```text
c^2e=1-b^2.                                             (2)
```

Thus the assignments in `(1)` define a homomorphism

```text
D=C[b,c,e]/(c^2e-(1-b^2))  -->  C[x,z].                 (3)
```

It is injective.  Indeed, `c^2e+b^2-1` is primitive and linear in `e`, hence
irreducible, so `D` is a domain.  After inverting `c`, `(3)` is the
isomorphism

```text
D[c^-1]=C[b,c,c^-1]  ~=  C[x,x^-1,z],
x=c,                 z=(1-b)/c^2.                      (4)
```

The localization map from the domain `D` is injective, proving the claim.

## 1. The ambient embedding preserves the Poisson object

Use the source bracket `{F,G}=F_xG_z-F_zG_x`.  Direct differentiation gives

```text
{b,c}=c^2,             {c,e}=2b,             {b,e}=-2ce. (5)
```

For `Sigma(b)=1-b^2`, these are exactly

```text
{b,c}=c^2,       {c,e}=-Sigma'(b),       {b,e}=-2ce,   (6)
```

the bracket used in THM-3569/3583/3592.  Moreover `Sigma` is squarefree of
degree two, so every hypothesis of those universal theorems holds.  Give `D`
their grading

```text
wt(b,c,e)=(0,1,-2).                                    (7)
```

The three collision-ring generators of THM-3686 become

```text
A=3e,                 B=2ce,                 C=bc.     (8)
```

Consequently

```text
R=C[A,B,C]=C[3e,2ce,bc]  subset D                      (9)
```

is a graded subalgebra, with `wt(A,B,C)=(-2,-1,1)`.  It is important that
`R` is not asserted to be closed under the bracket; for example `{A,B}` need
not lie in `R`.  What `(5)--(6)` prove is that the injection `D -> C[x,z]` is
Poisson.  Hence for `P,Q in R subset D`, their ambient bracket maps exactly to
their source bracket.  In particular, a source equation `{P,Q}=1` is the same
equation in `D`, by injectivity.

Because the inclusion `R subset D` is injective and grading-preserving, the
nonzero homogeneous pieces
of an element of `R` are exactly the pieces seen when that element is regarded
in `D`; neither quotienting nor cancellation changes the support count.

## 2. A direct every-line and source-jet gate

Let `P,Q in R` and suppose

```text
{P,Q}=1.                                               (10)
```

Then `H=(P,Q)` is a planar Keller map.  It is noninjective because every
element of `R` identifies the fixed collision

```text
(x,z)=(1,0),(-1,2)  -->  (A,B,C)=(0,0,1).              (11)
```

Gwozdziewicz's **CITED** injective-line theorem says that a complex planar
Keller map injective on one affine line is an automorphism.  Thus the
restriction of `H` to every affine source line is noninjective.  It is also
immersed, since the invertible differential of `H` kills no nonzero line
tangent.

An immersed noninjective polynomial curve `(f,g):A1_C->A2_C` has both
coordinate degrees at least two.  Indeed, a degree-one coordinate is
injective.  If one coordinate is constant, immersion makes the derivative
of the other nowhere zero; over `C` that derivative is constant, again giving
an injective curve.

Apply this lemma to

```text
z=0:       (A,B,C)=(0,0,x),
x=0:       (A,B,C)=(6z,0,0).                           (12)
```

After scalar weight-zero pieces are removed, each of `P,Q` therefore has a
pure-`C` weight at least `2` and a pure-`A` weight at most `-4`.  At the
origin,

```text
A=6z+higher terms,       B=terms of degree >=2,
C=x+higher terms.                                      (13)
```

Each gradient row of `(P,Q)` is nonzero.  Its linear jet can only come from
weight `-2` (`A`) or weight `1` (`C`), neither of the two extreme pieces just
found.  Hence

```text
#supp_wt(P)>=3,                    #supp_wt(Q)>=3.      (14)
```

This ring-specific argument, not the finite `2 x 3` rectangle of THM-3693,
supplies the coordinatewise floor.

## 3. The seven-piece floor and first live cell

Regard `P,Q` as elements of `D`.  The grading and every nonzero homogeneous
piece are preserved.  Repaired THM-3592 applies to the squarefree quadratic
`Sigma=1-b^2` and excludes `3 x 3`.  Combining it with `(14)` gives

```text
#supp_wt(P)+#supp_wt(Q)>=7,
first live partitions: (3,4) and (4,3).                (15)
```

The inherited THM-3569 and THM-3583 exclusions remain independent controls,
but `(14)` also removes their ambient `2 x 5` survivor from the collision
ring.  This is the correction recorded in MISTAKE-442.

## 4. A cited target-filtration floor

Because the hypersurface relation is not homogeneous in ordinary target
degree, let `F_N R` consist of elements admitting representatives by target
monomials `A^iB^jC^k` with `i+j+k<=N`.  Their source degrees are

```text
deg_(x,z)(A,B,C)=(4,5,4).                              (16)
```

If `P,Q in F_N R`, every member of their target pencil has source degree at
most `5N`; cancellation can only lower that bound.  The **CITED**
Guccione--Guccione--Horruitiner--Valqui result recorded in
[CORE-PAPERS](../../05-knowledge/reference/CORE-PAPERS.md) gives reduced
pencil height at least `108` for a planar Jacobian counterexample.  The
collision `(11)` therefore implies

```text
P,Q cannot both lie in F_21 R,
max(delta_R(P),delta_R(Q))>=ceil(108/5)=22,             (17)
```

where `delta_R(F)=min{N:F in F_N R}`.  At common cap `22`, raw height is at
most `110<125`, so the same cited classification forces reduced pair
`(72,108)` up to order.  This asserts neither existence nor the displayed
coordinate degrees: a linear pencil cancellation may expose degree `72`.

## 5. The collision-ring coefficient sidecar

A target monomial `A^iB^jC^k` of weight

```text
r=-2i-j+k                                               (18)
```

has `c`-homogeneous coefficient proportional to

```text
(1-b^2)^(i+j)b^k.                                      (19)
```

Consequently every retained weight coefficient in `R` obeys

```text
r>0  => b^r divides the coefficient,
r=0  => b(1-b^2) divides the nonscalar coefficient,
r<0  => (1-b^2)^ceil(-r/2) divides the coefficient.    (20)
```

The roots `b=+1,-1` are the ambient arm divisor.  The divisor `b=0` is the
extra information lost when one enlarges `R` to `D`.  THM-3696 strengthens
`(20)` to the complete coefficient modules and a synchronized three-branch
scalar law; the W003 argument below needs only the displayed coarse divisors.

### 5.1 An all-scale W003 placement nonentry

On the THM-3603/3613 word `W003`, choose scalar anchor `20` with orientation
`(1,-2)`.  At scale `n>=2`,

```text
supp(P)=(1-3n,1-2n,1),
supp(Q)=(-2,n-2,2n-2,3n-2).                            (21)
```

Its ordered fibres are

```text
00; 01+10; 02+11; 03+12+20; 13+21; 22; 23,            (22)
```

with the middle triple scalar.  Write `f_i,g_j` for the coefficient
polynomials and `nu` for order at `b=0`.  The positive-weight floor in `(20)`
makes address `20` the only possible scalar-unit supplier, forcing

```text
nu(f_2)=1,                         nu(g_0)=0.            (23)
```

The singleton UFD rows `00`, `22`, and `23` then force

```text
nu(f_0)=0,        nu(g_2)=2n-2,       nu(g_3)=3n-2.    (24)
```

Put

```text
alpha=nu(f_1)>=0,                  beta=nu(g_1)>=1.    (25)
```

At `n=2`, `g_1` has retained weight zero, so the middle clause of `(20)`
still supplies the bound on `beta`.  The row `02+11` has two nonresonant
initial terms; equality of their orders gives

```text
alpha+beta=2n-2.                                        (26)
```

In `01+10`, either `alpha=0` and the resonant possibility is conservatively
retained, or both initial terms are nonresonant and cancellation forces
`alpha=beta`.  Thus

```text
(alpha,beta)=(0,2n-2)  or  (n-1,n-1).                 (27)
```

Finally, the two terms of `13+21` have exact initial orders

```text
alpha+3n-3,                         beta.               (28)
```

They are respectively `(3n-3,2n-2)` or `(4n-4,n-1)`.  The second is uniquely
lowest, with nonzero multiplier `-n` or `-1`, and cannot cancel.  Therefore

```text
the W003 anchor-20 placement family is empty in R for every n>=2.          (29)
```

The placement is absent at `n=1`.  THM-3613 already rejects every even scale;
`(23)--(28)` newly close every odd scale `n>=3`, and output exchange gives the
transposed family.  This closes one complete placement family, not the whole
word `W003` and not the full `3 x 4` frontier.  The proof is symbolic for all
`n`; finite checks through `n=64` are controls only.

## 6. Exact construction frontier and scope

The ambient THM-3603/3606/3609 atlas has `144` scalar schemes on `27` oriented
words before the THM-3613 parity refinement.  The next subring-specific test
is to impose THM-3696's exact coefficient modules and synchronized scalar law
at `b=+1,-1,0`; passing the resulting three-place valuation system is still
necessary only and does not lift orders to coefficient polynomials.

Any construction inside `R` must clear

```text
support:           3 x 4 or wider, outside (29),
target filtration: at least 22 in one output,
boundary:          THM-3686 normalization and axis-jet laws,
subring:           all three divisor clauses in (20).                  (30)
```

The every-line and height steps are **CITED**.  THM-3592 is used only in its
repaired post-MISTAKE-432 form.  No statement here constructs a Darboux pair
on `D` or `R`, a polynomial counterexample, or a proof of `JC(2)`.

## 7. Exact reproduction

```bash
python3 -B 04-computation/jacobian_y0_danielewski_transfer_support_floor_thm3695.py
python3 -B -O 04-computation/jacobian_y0_danielewski_transfer_support_floor_thm3695.py
```

Both modes byte-match the stored transcript.  The companion checks the
ambient Poisson signs and injection identities, collision, line restrictions,
origin jet, `164` coefficient-module monomials, `126` W003 controls through
scale `64`, support arithmetic, and `(16)--(17)`.  The cited inputs and the
all-scale arguments remain proof-driven, not finite extrapolations.

**QED.**
