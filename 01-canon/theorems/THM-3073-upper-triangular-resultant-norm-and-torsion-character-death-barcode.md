---
id: THM-3073
title: "Upper-triangular resultant norm and torsion-character death barcode"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A pure
  lower complete intersection followed by an arbitrary upper-triangular
  homogeneous block has resultant S^(upper degree product) times T^(lower
  degree product), with canonical positive sign.  Hence its scalar phase is
  the corresponding power-map quotient of the two input phases.  Applied to
  THM-3063 and THM-3069 this gives exact C7, C13, C91, and Kummer death
  barcodes, and separates genuine pullback holotopy from norm contraction.
  It is an information-loss theorem, not an LRC, JC, or GMC closure.
source: root-gmc-lrc-jc-norm-barcode-2026-08-01
audit: >
  Two independent hostile audits rederived the universal vanishing-set and
  multidegree proof, canonical sign, both physical GMC exponent maps, every
  direct and iterated torsion-death row, the phase-placement and relative
  cancellation controls, and the mixed-lower hostile.  The immutable audit
  additionally checked eighteen small normalized resultants, every proved
  dependency, the pullback/Kummer and LRC scope boundaries, normal and
  optimized replay, stored bytes, declared LF hashes, and documentation.
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2962-planar-suspension-resolvent-and-v4-kummer-descent
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
related:
  - THM-2543-augmentation-norm-relative-phase-local-system-dichotomy
  - THM-2865-gamma-transverse-null-holotopy-and-uniform-fourth-exit
  - THM-2983-coloured-macaulay-wall-normal-cone-excess-contact-atlas
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
script: 04-computation/gmc_upper_triangular_resultant_norm_barcode_thm3073.py
output: 05-knowledge/results/gmc_upper_triangular_resultant_norm_barcode_thm3073.out
script_sha256: 4d5274f372cfd16a0e9adb7238ea340a6d56cde5a2bb40569c23bb23dd9395b0
output_sha256: c81120eb64cc5792ffd2d70384ee12052fd8a5a9b69921f5daf73822c60b9fbc
hash_basis: LF-normalized bytes
---

# THM-3073 -- upper-triangular resultant norm and torsion-character death barcode

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Several maintained arguments add a coordinate and then ask whether an old
phase, owner, root origin, or Kummer class survives.  There are two different
operations hiding under that sentence:

1. a **pullback** carries the old class to the new space with winding one;
2. a **norm/resultant contraction** multiplies the old phase over a finite
   fibre and therefore raises it to a degree.

The second operation can prove nonvanishing or positivity while destroying
exactly the torsion label needed downstream.  The distinction is not
heuristic.  It follows from one universal resultant identity.

## 1. The upper-triangular product theorem

Let `K` be a field, let `a,b>=1`, and let

```text
y=(y_1,...,y_a),       z=(z_1,...,z_b),                    (1)
```

and take positive integers `d_i,e_j`.  For `1<=i<=a`, let `H_i(y)` be a
homogeneous form of degree `d_i` depending **only** on `y`.  For `1<=j<=b`,
let `G_j(y,z)` be an arbitrary homogeneous form of degree `e_j`, and put

```text
A_j(z)=G_j(0,z),
D_y=prod_i d_i,         D_z=prod_j e_j,
S=Res_y(H_1,...,H_a),   T=Res_z(A_1,...,A_b).              (2)
```

Order variables as `(y,z)` and forms as `(H,G)`.  Normalize every homogeneous
resultant by

```text
Res(y_1^d1,...,y_a^da,z_1^e1,...,z_b^eb)=1.               (3)
```

Then

```text
Res_(y,z)(H_1,...,H_a,G_1,...,G_b)=S^D_z T^D_y.            (4)
```

In particular, every mixed `y,z` coefficient in the upper forms `G_j` is
invisible to the scalar resultant once the pure lower block and the normal
restrictions `A_j` are fixed.

### Proof

It is enough to work over an algebraic closure.  Suppose first that `T!=0`
and `S=0`, and choose a common base direction `y_0`.  On the projective space
with coordinates `[lambda:z]`, the equations

```text
G_j(lambda y_0,z)=0                                      (5)
```

are `b` positive-degree hypersurfaces in `P^b`.  Their intersection number is
`D_z`.  At `lambda=0` they become the `A_j`, which have no common projective
zero because `T!=0`.  Hence their intersection is nonempty and lies in the
affine chart `lambda!=0`.  Rescaling gives a common zero `(y_0,z)` of the full
system.  If instead `T=0`, a common normal direction `z_0` gives the full zero
`(0,z_0)` directly.

Conversely, a full common projective zero `(y,z)` has either `y!=0`, in which
case the pure equations `H_i(y)=0` force `S=0`, or `y=0`, in which case
`z!=0` and `A_j(z)=0` force `T=0`.  Thus the zero divisor of the specialized
full resultant is exactly the union of the irreducible divisors `S=0` and
`T=0`.  It follows in the universal coefficient ring that

```text
Res_(y,z)(H,G)=c S^alpha T^beta.                            (6)
```

The full resultant has degree

```text
(D_y/d_i)D_z
```

in the coefficients of `H_i`, whereas `S` has degree `D_y/d_i`.  Hence
`alpha=D_z`.  The degree in the coefficients of `G_j` is
`D_y(D_z/e_j)`, whereas `T` has degree `D_z/e_j`; hence `beta=D_y`.
Finally the coordinate-power system `(3)` gives `c=+1`.  This proves `(4)`
over the universal integer coefficient ring and therefore over every field.

The purity of the lower block is load-bearing.  If both equations mix both
blocks, then

```text
F_1=y(y+cz),            F_2=z(z+dy)                         (7)
```

have nominal axis restrictions `y^2,z^2`, each with resultant one, but

```text
Res_(y,z)(F_1,F_2)=1-cd.                                   (8)
```

Thus `(4)` is upper-triangular, not a statement about an arbitrary block
decomposition.

## 2. The phase torus and its exact kernel

Write the two nonzero input phases as `(u,v)`.  Formula `(4)` induces the
monomial map

```text
Phi_(D_z,D_y): G_m^2 -> G_m,
(u,v) |-> u^D_z v^D_y.                                    (9)
```

On fundamental groups its winding row is `(D_z,D_y)`.  On `mu_q`, in additive
exponent coordinates, it is

```text
(r,s) |-> D_z r + D_y s mod q.                            (10)
```

Put `g=gcd(q,D_z,D_y)`.  Then

```text
|image Phi|=q/g,             |kernel Phi|=qg.              (11)
```

A phase of exact order `q` placed only in the lower block leaves with exact
order

```text
q/gcd(q,D_z),                                             (12)
```

and the analogous normal phase leaves with order `q/gcd(q,D_y)`.  It dies
precisely when the relevant exponent is divisible by its **actual** order.
Calling an element a `C91` label is not enough: `zeta_91^13` has order seven,
whereas `zeta_91^7` has order thirteen.

Separate survival also does not imply scalar survival.  In the two-normal
`k=4` map below the exponents are `(12,2)`.  If
`u=v=zeta_7`, then neither input is killed separately, but

```text
u^12 v^2=zeta_7^14=1.                                    (13)
```

The scalar carrier has forgotten a nontrivial relative phase.

## 3. The two physical GMC suspension maps

THM-3069 has lower degrees `2,...,k-1` and one normal degree `k`.  Thus

```text
D_y=(k-1)!,             D_z=k,
(u,v) |-> u^k v^((k-1)!).                                 (14)
```

Its physical normal coefficient is `v=(-1)^k`; since `(k-1)!` is even for
`k>=3`, only `u^k` remains.  This is the phase-level reason that the remote
four-slot result has sign `sign(S^4)=+1`.

THM-3063 has lower degrees `2,...,k-2` and two normal degrees `k-1,k`.  Hence

```text
D_y=(k-2)!,             D_z=k(k-1),
(u,v) |-> u^(k(k-1)) v^((k-2)!).                           (15)
```

The base exponent is always even, and the physical binary normal resultant
has the positive standard orientation.  Again, the same contraction which
proves eventual positivity has deliberately removed base orientation.

For a phase placed in exactly one block, the first stage `k>=3` at which the
whole character of the displayed order dies is:

| phase location | `C7` | `C13` | primitive `C91` |
|---|---:|---:|---:|
| THM-3069 base, exponent `k` | 7 | 13 | 91 |
| THM-3069 normal, exponent `(k-1)!` | 8 | 14 | 14 |
| THM-3063 base, exponent `k(k-1)` | 7 | 13 | 14 |
| THM-3063 normal, exponent `(k-2)!` | 9 | 15 | 15 |

Placement therefore matters.  At `k=8`, for example, a `C7` phase in the
THM-3063 base is dead because `56=0 mod 7`, while a normal `C7` phase survives
because `6!=6 mod 7`.  In THM-3069 the roles reverse between `k=7` and `k=8`:
the base dies and normal survives at seven, while the base survives and normal
dies at eight.  An unlabelled positive scalar cannot reconstruct which input
carried the phase.

## 4. The iterated one-normal death barcode

Iterate the THM-3069 suspension through one fixed lexicographically ordered
tower, starting with a width-`m` carrier and ending at width `k`.  Repeated use
of `(14)` gives

```text
R_m |-> R_m^(k!/m!),
U_j |-> U_j^(k!/j),             m<j<=k.                    (16)
```

Thus a phase of exact order `q` on `R_m` has order

```text
q/gcd(q,k!/m!)                                             (17)
```

at width `k`.  Define its death duration by

```text
delta_q(m)=min{d>=1: q divides (m+1)(m+2)...(m+d)}.         (18)
```

For a prime `p`, with residue zero interpreted as `p`, one has

```text
delta_p(m)=p-(m mod p).                                    (19)
```

Since seven and thirteen are coprime,

```text
delta_91(m)=max(delta_7(m),delta_13(m)).                   (20)
```

Consequently:

| starting width | `C7` duration | `C13` duration | `C91` duration |
|---|---:|---:|---:|
| `m=2` | 5 | 11 | 11 |
| `m=3` | 4 | 10 | 10 |

The corresponding terminal widths are seven and thirteen.  Once a phase is
trivial in this **same composable tower**, later power maps keep it trivial.
There is no such monotonicity across unrelated fixed-`k` families: the direct
exponent `k` can vanish modulo `q` and then become a unit again at `k+1`.

Equation `(16)` is also only lexicographic.  If lower supports move on the
same scale, their sidecar can approach zero as fast as the discarded terms.
THM-3069's ordered limit and this barcode then supply no simultaneous chamber.

## 5. Pullback bars do not die

THM-2962 is the opposite operation.  Its direct `A1` suspension is a literal
base change, and pullback gives a natural isomorphism

```text
H^1_et(U,mu_2) ~= H^1_et(U x A1,mu_2).                     (21)
```

The `V4` torsor and all three nonzero Kummer characters therefore persist
with winding one.  Applying `(17)` to them would contradict `(21)` because
no resultant norm has been taken.

The contrast is already visible at order two.  In THM-3063 the base exponent
`k(k-1)` is always even, and the normal exponent `(k-2)!` is even for `k>=4`;
the scalar phase map is then trivial on `mu_2^2`.  In THM-3069 both exponents
are even whenever `k>=4` is even.  Such a scalar norm cannot replace the
Kummer plane preserved by THM-2962.

THM-2542 has a related but internal LRC power law: pulling its chart holonomy
to an `n`-fold clock cover multiplies the order-thirteen class by `n`, and the
minimal trivializing cover has degree thirteen.  Formula `(12)` describes the
same elementary power-map arithmetic.  There is, however, **no proved map**
identifying that clock-cover degree with a GMC factorial-slot exponent.
Likewise, THM-2543's cyclic norm is a positive product of seven phase means,
not the torsion-character norm `(9)`.

The valid cross-frontier conclusion is an information audit:

> Any proposed bridge which stores an LRC `C91` root-clock label or a JC
> Kummer label only in the scalar resultant after `(4)` must first quotient
> that label by `(10)`.  Nonvanishing and positivity of the scalar do not
> restore anything in the kernel.

A surviving bridge therefore needs a pre-norm sidecar: the ordered pair
`(S,T)`, a determinant/root line, a pointed branch or root origin, the Kummer
character plane, or an independently typed semantic cospan.

## 6. Norm holotopy and the zero-sidecar boundary

On phase circles, `(9)` is a holotopy contraction whose winding map is the row
`(D_z,D_y)`.  Across the one-normal slot tower, `(17)--(20)` form a literal
torsion persistence barcode: a bar ends when its winding becomes divisible by
the character order.  This phase barcode is complementary to THM-2865 and
THM-3069's magnitude holotopies.  Those theorems deform a positive carrier to
the physical resultant; the present theorem records what their scalar norm
forgets during that deformation.

If `S=0`, formula `(4)` correctly says that the leading scalar carrier
vanishes.  It does not identify the next physical factor.  The honest object
is then the direction-dependent kernel-to-cokernel normal map and liftable
kernel barcode of THM-2983/2985.  A factor seen on one arc need not extend to
the whole normal cone.

## 7. Exact evidence and scope

The exact companion verifies:

1. all `7,056` degree-block pairs with one to three equations per block and
   degrees one through four, comprising `38,304` multidegree equalities;
2. `324` exact phase image/kernel cells for the selected nine-modulus bank
   through `91` and both GMC exponent families;
3. every direct first-death stage in the table;
4. `76` iterated tower cells and the prime/`C91` duration formulas;
5. the `k=8` placement hostile and `k=4` two-input cancellation hostile; and
6. the exact symbolic Sylvester determinant `1-cd` in `(8)`.

Run

```text
python 04-computation/gmc_upper_triangular_resultant_norm_barcode_thm3073.py
python -O 04-computation/gmc_upper_triangular_resultant_norm_barcode_thm3073.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem proves an upper-triangular resultant identity and its exact
phase quotient.  It does not construct a transport from GMC phases to LRC or
JC labels, a semantic vertical edge, a root origin, a Kummer owner, a
simultaneous moving-scale cone, arbitrary-radial GMC(2), NC2, LRC(14), JC(2),
or DC(2).  Its purpose is to prevent a positive norm from being mistaken for
the sidecar it has erased.

**QED.**
