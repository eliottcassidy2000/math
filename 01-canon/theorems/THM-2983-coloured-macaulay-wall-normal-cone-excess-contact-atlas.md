---
id: THM-2983
title: "Coloured Macaulay wall normal-cone and excess-contact atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The natural
  row-deficit-coloured width-six Macaulay chart has exact top-n quadratic
  sections, first normal cones, and direction-indexed liftable-kernel
  barcodes. At root six, the projectivized first normal cone is a non-SNC
  K2 death-jump arrangement with a shared line and singular quadric point.
  Arc contact, ambient multiplicity, and ambient divisibility remain
  rigorously separated; no coloured strict-transform or all-arc theorem is
  claimed.
source: codex-gmc-coloured-wall-normal-cones-2026-07-30
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
related:
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
script: 04-computation/gmc_m6_coloured_wall_normal_cone_excess_contact_atlas_thm2983.py
output: 05-knowledge/results/gmc_m6_coloured_wall_normal_cone_excess_contact_atlas_thm2983.out
script_sha256: f68624bf90099eabda35214795f48d9505d39c9a867fd4283aecaf3f8a967e06
output_sha256: 0df63db9a79bbc441bcd084e7094893cf2e583832149b45a2819d3403e7549ed
constructor_dependency_sha256: 5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5
hash_basis: LF-normalized bytes
---

# THM-2983 -- coloured Macaulay wall normal-cone and excess-contact atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The theorem identifies what survives when the one-variable width-six wall
invoice is restored to the natural row-deficit-coloured Macaulay chart.  The
answer is not a coloured copy of the diagonal factor.  It is a stratified
normal-cone and arc-barcode atlas.

## 1. Coloured chart and its top n-jet

Give a deficit of `k` in a quadratic, cubic, or quartic Macaulay row the
weight `q^k,c^k,f^k`, respectively, and let

```text
P(n,q,c,f)                                                (1)
```

be the resulting homogeneous width-six determinant.  It has degree `312`.
The canonically coloured THM-2942 flag `q_200^6 c_300 K` has degree `62` and
divides `(1)`; write the pure resultant quotient as `R`, of degree `250`.

The exact quadratics obtained from `partial_n^(310)P` and
`partial_n^(248)R`, up to their positive common scalars, both have inertia

```text
(1 positive, 3 negative).                                (2)
```

The full-chart and quotient top two diagonal coefficients independently
match the one-variable interpolation.  On the colour-transverse plane
`q+c+f=0`, the quotient colour block has

```text
trace = 244051230785617/5627343750,
det   = -10882890806368883653776584/376988067626953125.  (3)
```

Equation `(2)` is one exact top-jet Lorentzian section, not a proof that the
four-variable polynomial is Lorentzian.  A factor depending only on the
common colour cannot change `(3)`.  The Smith/seam wall is defined only after
the diagonal specialization, so no three-colour strict-transform identity
follows from `(2)` or `(3)`.

## 2. Half-wall and integer first normal cones

At the half-wall point `(-1,2,2,2)`, projectively the diagonal root
`n=-1/2`, the chart has rank `30` and corank `6`; the coloured flag is a
unit.  In tangent coordinates `(z,q,c,f)`, THM-2985's canonical first normal
map has determinant

```text
unit * (f+2z)^6.                                         (4)
```

Thus the half-wall has the expected ambient multiplicity six, but its first
normal cone is one repeated hyperplane rather than six independently
coloured factors.

Let

```text
W_6=B_6 E_6/H_6
```

be THM-2969's diagonal Smith/seam wall product after its selected-flag
overlap is removed.  At the six integer points `(-r,1,1,1)`, the exact
rank/order census is

| `r` | corank | flag order | diagonal `W` order | diagonal chart order |
|---:|---:|---:|---:|---:|
| 1 | 2 | 0 | 6 | 6 |
| 2 | 10 | 1 | 24 | 25 |
| 3 | 15 | 1 | 26 | 27 |
| 4 | 10 | 0 | 24 | 24 |
| 5 | 9 | 0 | 23 | 23 |
| 6 | 15 | 1 | 20 | 21 |

Exact rational ranks computed during the nullspace/tangent and block-Toeplitz
calculations give this table; agreement at three large primes is a separate
rank check, not the proof of the rational upper bound.  The first normal
determinants are:

```text
r=1:  unit * (q+f-2c)^2;                                (5)

r=2:  0, with common right/left kernel dimensions 5/2;
r=3:  0, with common right/left kernel dimensions 8/10;
r=4:  0, with common right/left kernel dimensions 9/8;
r=5:  0, with common right/left kernel dimensions 8/8;   (6)

r=6:  unit * (6c+z)^2
             * (7748852908937866554c
                -2517055682762701080q
                +871966204362527579z)
             * (-6cf-6cq-2cz+12fq+fz+qz)^6.             (7)
```

The degrees in `(7)` are `2+1+12=15`, exactly the corank.  In particular,
root six has a genuine three-component first normal cone.  Roots two through
five instead begin with excess ambient contact certified by their kernel
flags.  The common-kernel certificate is sufficient, not logically
necessary in general, by THM-2985's alternating-pencil hostile.

The apparently unstructured coefficients in `(7)` hide an exact incidence.
Put

```text
A_6=z+6c,             u=c-q,
B=2517055682762701080, C=871966204362527579.             (7a)
```

Since

```text
7748852908937866554=B+6C,
```

the linear and quadratic factors in `(7)` are exactly

```text
L=C A_6+B u,
Q=A_6(f-2c+q)+12(c-q)(c-f).                              (7b)
```

Consequently the three projective components are not SNC: `A_6=L=0`
implies `u=0` and hence `Q=0`.  They share the projective line

```text
z=-6c, q=c,
```

and the quadric `Q=0` has its unique projective singular point on that line
at `(-6,1,1,1)`.

## 3. Arc holotopy: liftable-kernel barcodes

For a line through one of the integer points, expand

```text
A(t)=A_0+tA_1+t^2A_2+... .                               (8)
```

Let `B_j` be the `36j by 36j` lower block-Toeplitz matrix expressing
`A(t)v(t)=0 mod t^j`.  THM-2985 gives, and the companion checks exactly,

```text
dim K_j = 36-rank(B_j)+rank(B_(j-1)).                    (9)
```

The affine diagonal wall crossing is the `n`-only direction `(1,0,0,0)`.
Its liftable dimensions, including the final zero, are

```text
r=1: (2,2,1,1,0),       contact 6;
r=2: (10,9,4,2,0),      contact 25;
r=3: (15,10,2,0),       contact 27;
r=4: (10,9,4,1,0),      contact 24;
r=5: (9,8,5,1,0),       contact 23;
r=6: (15,6,0),          contact 21.                    (10)
```

Their sums exactly reproduce the flag-plus-wall invoice.  Along the explicit
coloured line `(z,q,c,f)=(1,2,3,5)`, the barcodes become

```text
r=1: (2,0),             contact 2;
r=2: (10,5,2,0),        contact 17;
r=3: (15,10,2,0),       contact 27;
r=4: (10,9,4,1,0),      contact 24;
r=5: (9,8,5,1,0),       contact 23;
r=6: (15,0),            contact 15.                    (11)
```

Two further exact directions `(2,3,5,7)` and `(1,-1,2,-3)` reproduce the
root-two through root-six sequences in `(11)`, except that the former lies
on the hyperplane in `(5)` and has root-one barcode `(2,1,0)`, contact
three.  This is the equality boundary of the first cone, not noise in the
rank calculation.

Equations `(10)` and `(11)` are a direction-indexed kernel **holotopy**:
the special-fibre kernel is fixed, but its liftable classes die at different
grades under different deformations.  No claim is made that the three tested
directions classify every arc.

Every displayed sequence terminates at zero.  This also verifies the
nonzero-determinant hypothesis in THM-2985: if `det A(t)` were identically
zero, a nonzero formal kernel vector, divided by its least `t`-power, would
give a nonzero special-fibre class lifting to every depth.  The terminal zero
rules that out.  Hence summing `(9)` is lawful for every contact in
`(10)--(11)`.

### 3.1 The root-six first-death arrangement

For a line direction `v`, the first lifting condition is

```text
A_1(v)v_0 in im(A_0).
```

Thus `K_2(v)` is exactly the kernel of THM-2985's first normal map evaluated
at `v`.  The projectivized first normal cone `(7)` is therefore precisely the
first death-jump arrangement for the root-six barcode, not merely a
determinant factorization.

Exact block-Toeplitz ranks give the following rational points on its strata.
Coordinates are `(z,q,c,f)`, and `D=B+12C`.

| stratum | direction | liftable dimensions | contact |
|---|---|---|---:|
| off the cone | `(1,2,3,5)` | `(15,0)` | `15` |
| `A_6` only | `(-6,2,1,3)` | `(15,2,0)` | `17` |
| `L` only | `(B,C,0,1)` | `(15,1,0)` | `16` |
| `Q` only | `(1,1,1,1)` | `(15,6,0)` | `21` |
| smooth `L cap Q`, off `A_6` | `(BD,CD,0,-BC)` | `(15,6,1,0)` | `22` |

This table also shows why factor exponent and total arc contact must not be
identified.  On the smooth `L cap Q` point the first surviving space has
dimension six, while one class lives for one additional grade.

The non-SNC boundary behaves differently.  At the `A_6 cap Q` point
`(-6,0,1,1)` the exact checked prefix is

```text
(15,6,6,6,6,6),
```

and this is a genuine resultant-null arc, not merely a long finite prefix.
Indeed, along the whole affine line put `s=1+t`, so that

```text
n=-6s, q=1, c=f=s.
```

For `ell=x_0+x_1+x_2` and

```text
H=63299 x_0^2+121066 x_0x_1+114169 x_0x_2
  +57772 x_1^2+108669 x_1x_2+50919 x_2^2,             (11a)
```

the three root-six forms at `s=1` factor exactly as

```text
Q_1=-24H,
C_1=641571840 ell H,
F_1=-93261243646771200 ell^2 H,                         (11b)
```

while on the full arc

```text
C_s=641571840 s^11 ell H,
F_s=-93261243646771200 s^17 ell^2 H.                    (11c)
```

For `s` nonzero, the projective plane conics `Q_s=0` and `H=0` have a
common point over the algebraic closure, and that point also kills `C_s`
and `F_s`.  At `s=0`, the last two forms vanish identically and the same
conclusion is immediate.  Hence the resultant, and therefore this Macaulay
determinant, vanishes identically on the `A_6 cap Q` arc.

There is a second, logically different null mechanism.  On the common triple
line, the generic point `(-6,1,1,2)` and the singular-quadric point
`(-6,1,1,1)` have checked prefixes

```text
(15,9,9,9,9,9),             (15,15,15,15,15,15),       (11d)
```

respectively.  There the coloured THM-2942 curvature flag vanishes
identically on the whole affine arc.  Thus `(11b)--(11c)` certify a
resultant-null common-zero mechanism, whereas the triple line is
flag-null.  In both cases the displayed prefixes are null-arc boundary data,
not finite Smith barcodes.

## 4. Exact separation of the diagonal wall from an ambient factor

At roots one, two, and six the diagonal contacts are respectively

```text
6,25,21,
```

while the explicit coloured line has contacts

```text
2,17,15.                                                 (12)
```

Consequently no ambient divisor through the corresponding point can occur
in `P` with the full diagonal multiplicity: such a power would force at
least that order on every arc through the point.  The missing order in
`(12)` is tangency/lifetime along the diagonal arc, not an ambient power that
may be divided from the coloured determinant.

At roots three through five, the tested transverse contacts happen to equal
the diagonal contacts.  This finite agreement is not promoted to an ambient
factorization or an all-arc classification.  At the half-wall `(4)` gives
the initial sixth power, but an initial-form power alone likewise does not
prove divisibility of the complete germ.

The repo's noncanonical `S_3`-equivariant SNC reflection supplies the
complementary Gram control.  On the iterated blow-up it takes
`D_i=E_i-F_i`; the three classes are primitive and independent even modulo
two, while their Gram matrix is `-2I_3` and hence has a full mod-two zero
mode.  Thus neither `(2)` nor a Gram/Hessian kernel can replace the normal
map or the liftable-kernel tower.  This geometric control is proved on its
reflection surface, but is not a canonical theorem dependency and has no
quartic-resolvent realization attached.

## 5. Scope

This theorem is an exact width-six local atlas.  It proves:

- the natural coloured top-`n` quadratic sections `(2)`;
- the first normal cones `(4)--(7)`;
- the explicit arc barcodes `(10)--(11)`;
- the distinct resultant-null and flag-null boundary mechanisms
  `(11a)--(11d)`; and
- the arc-versus-factor separation `(12)`.

It does **not** prove that the coloured determinant is Lorentzian, construct
a coloured lift of the complete Smith/seam wall, classify every arc, or
extend the strict-ULC range of the separately proposed THM-2982 candidate.

## 6. Exact evidence

The consolidated companion pins the LF hash of the imported PROVED THM-2969
constructor, uses exact `QQ/fmpz` ranks as the load-bearing wall and
block-Toeplitz calculations, and uses the three large-prime ranks only as
independent cross-checks.  Reproduce with

```text
python 04-computation/gmc_m6_coloured_wall_normal_cone_excess_contact_atlas_thm2983.py --output .scratch/thm2983.normal.out
python -O 04-computation/gmc_m6_coloured_wall_normal_cone_excess_contact_atlas_thm2983.py --output .scratch/thm2983.opt.out
```

The normal, optimized, and stored transcripts are byte-identical UTF-8/LF
files of `3,410` bytes.  Their frozen hashes are

```text
script=f68624bf90099eabda35214795f48d9505d39c9a867fd4283aecaf3f8a967e06,
output=0df63db9a79bbc441bcd084e7094893cf2e583832149b45a2819d3403e7549ed,
THM-2969 dependency=5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5.
```

An independent hostile audit rederived the top-jet quotient subtraction,
all exact first-normal factors, the rational and modular rank scopes, the
block-Toeplitz lifetime identity, and the terminal-zero determinant gate. It
also exposed and verified the non-SNC root-six death arrangement, checked the
exact `A_6 cap Q` factorization `(11a)--(11c)` and its algebraic-closure
common-zero argument, separated that resultant-null arc from the triple
flag-null line, and accepted the arc-versus-factor and noncanonical
Gram-control boundaries. Normal,
optimized, and stored replay matched exactly; no load-bearing defect remains
within the stated width-six scope.

**QED.**
