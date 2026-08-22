---
id: THM-3188
title: "Quadratic-character pre-reset holonomy and exterior-flag rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Before every factorial prime reset, the homogeneous Gauss--Manin plane has
  the exact Legendre-scaled unipotent holonomy
  -chi(Delta)[[1,0],[-1,1]].  Its determinant wedge is fixed, so the complete
  block retains the same missing-wedge and conormal flags and the same
  squared-scale return.  On the discriminant wall the homogeneous plane
  collapses instead.
audit: >
  The pure-integer companion proves the closed continuant formula through
  length twelve, exhausts 948 unit/wall parameter triples over seven odd
  primes, and checks 300 exact full-block compound flags.  Normal and
  optimized replay agree with the stored transcript.  Two independent
  immutable audits rederived the continuant, all boundary cases, compound
  flags, squared return, wall collapse, and scope.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on:
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
related:
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
script: 04-computation/quadratic_character_pre_reset_holonomy_thm3188.py
output: 05-knowledge/results/quadratic_character_pre_reset_holonomy_thm3188.out
script_sha256: 47356320287682419a6d5a6e0b65c9522f15e0b738b45e875d5cc93b36a216bd
output_sha256: 6347861d5abb398c1f76b92db787768bae32443e5284ceafe8aff58d240191e9
hash_basis: LF-normalized bytes
---

# THM-3188 -- quadratic-character pre-reset holonomy and exterior-flag rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3185 identifies the final transfer in each prime block as a rank-one
reset.  The invertible tail immediately before that reset is not arbitrary:
its homogeneous part is one fixed Jordan shear, scaled only by the quadratic
character of the discriminant.  This is a finite-field Gauss--Manin
holonomy, and it rigidifies the oriented exterior flag through the whole
block.

## 1. The pre-reset transfer

Use THM-3185's exact state and transfer

```text
z_n=(M_n,X_n,d^n)^t,                    z_(n+1)=T_n z_n,

T_n=
[-(n+1)                  2v(n+1)                    d       ]
[-(n+1)(2n+3+2d)        (n+1)(1+2v(2n+3))         d(2n+3)]
[0                       0                           d       ]. (1)
```

Let `p` be odd, let `d==s (mod p)` with `s!=0`, and put

```text
Delta=1-4sv in F_p,             chi=Delta^((p-1)/2) in {+1,-1}. (2)
```

Assume first that `Delta!=0`.  Define the pre-reset and complete-block maps

```text
A=T_(p-2)...T_0,                       B=T_(p-1)A.             (3)
```

Then there are residues `beta,gamma` such that

```text
Abar=
[-chi       0          beta ]
[ chi      -chi        gamma]
[  0        0            1  ].                               (4)
```

In particular the homogeneous plane `H={D=0}` has exact monodromy

```text
Abar|H=-chi [1 0; -1 1].                                    (5)
```

The unit `s` has disappeared from `(5)` except through `Delta`.  Its separate
block multiplier remains the `s` in THM-3185's rank-one formula for `Bbar`.

## 2. Continuant proof of the holonomy

Set the boundary coordinate to zero and put

```text
u_n=M_n/n!,                  x=2v,
u_0=M_0,                     u_1=-M_0+2vX_0.                 (6)
```

For `1<=n<p-1`, the scalar recurrence becomes

```text
u_(n+1)=(2n+1)x u_n+Delta u_(n-1).                          (7)
```

Let `K_m` be the continuant with diagonal
`3x,5x,...,(2m+1)x` and dimer weight `Delta`.  Direct induction gives the
closed formula

```text
K_m=sum_(0<=j<=floor(m/2))
    binom(m-j,j) (2m-2j+1)!!/(2j+1)!!
    x^(m-2j) Delta^j.                                      (8)
```

For `m=p-3`, every coefficient except the last contains the factor `p` in
the numerator double factorial and none in its denominator.  For `m=p-2`,
every coefficient contains that factor.  Hence

```text
K_(p-3)==Delta^((p-3)/2),                  K_(p-2)==0.       (9)
```

For `p>=5`, the shifted continuant on `a_2,...,a_(p-3)`, where
`a_j=(2j+1)x`, also vanishes.  Indeed reflection

```text
j |-> p-1-j                                                    (10)
```

negates every monomer weight.  The interval has odd size.  Nonfixed
matchings cancel in reflected pairs, while a fixed matching must leave the
central vertex as a monomer, whose weight is `px=0`.  Thus

```text
u_(p-2)==Delta^((p-3)/2)u_1,                               (11)
```

and the coefficient of `u_1` in `u_(p-1)` is zero.

It follows that the homogeneous `2x2` block of `Abar` has form
`[alpha,0;lambda,mu]`.  Its determinant is

```text
prod_(n=0)^(p-2) (-(n+1)^2 Delta)
  =(-Delta)^(p-1)((p-1)!)^2=1.                             (12)
```

Wilson also gives `(p-2)!=1`.  For `v!=0`, THM-3183's exact gauge identity at
`n=p-1` reduces to

```text
2vX_(p-1)=(1-2v)M_(p-1)-Delta M_(p-2).                    (13)
```

Apply `(11)` first to initial data `(M_0,X_0)=(0,1)`.  Equation `(13)` gives
`mu=-chi`; `(12)` then gives `alpha=-chi`.  Applying `(11)` to `(1,0)` gives
`lambda=chi`.  This proves `(4)--(5)`.

When `v=0`, the top blocks in `(1)` are triangular.  The two diagonal
products are `-1`, and the alternating sum of `2n+3+2s` is `-1` before the
Wilson factor, giving the lower-left entry `1`; this is again `(5)` because
`Delta=chi=1`.  For `p=3`, direct multiplication of `T_1T_0` gives the same
formula.  Thus every odd prime is covered.

## 3. Exterior holonomy and the canonical flag

Use the ordered wedge basis

```text
w_01=M wedge X,             w_02=M wedge D,
w_12=X wedge D.                                             (14)
```

Taking the compound matrix of `(4)` gives

```text
Lambda^2 Abar=
[1   -chi(beta+gamma)    chi beta]
[0        -chi              0   ]
[0         chi             -chi ].                         (15)
```

Therefore

```text
(Lambda^2 Abar)w_01=w_01.                                  (16)
```

On the quotient by `w_01`, the same Legendre-scaled shear reappears.  More
invariantly, `Lambda^2 H` is fixed exactly, while

```text
Abar|H+chi I                                                 (17)
```

is a nonzero rank-one nilpotent whose image equals its kernel.  The scalar
`-chi`, the nontrivial Jordan class, and that invariant line survive a change
of basis in `H`; the displayed unit shear coefficient belongs to the named
factorial frame.  The extension coordinates `beta,gamma` depend on the chosen
boundary-coordinate splitting and are not canonical.

## 4. Whole-block flag rigidity

Let

```text
C_0=p^(-1)Lambda^2T_(p-1) mod p.                            (18)
```

THM-3185 proves that `C_0` has rank two, right kernel `w_01`, left kernel
`(1,-1,1)`, and squared-scale return

```text
p^(-2)Lambda^2T_(p-1)(w_01)==-Delta w_01.                  (19)
```

Compound multiplicativity and `(15)--(16)` give

```text
p^(-1)Lambda^2B mod p=C_0(Lambda^2Abar).                   (20)
```

Consequently the complete block has the same right missing wedge `w_01`, the
same left conormal `(1,-1,1)`, and rank two.  Moreover `Lambda^2A` preserves
the line `w_01` integrally with multiplier equal to the determinant of its
homogeneous block, which is `1` modulo `p`.  Hence the complete block also has
the unchanged return

```text
p^(-2)Lambda^2B(w_01)==-Delta w_01.                        (21)
```

Thus the invertible pre-reset tail may shear the two visible exterior
coordinates, but it cannot rotate or rephase the named missing-wedge line or
the named left conormal.  A chosen visible coordinate or path sum can still
cancel across the invertible tail; only loss of these canonical exterior flags
is forced to occur downstream.

## 5. Sharp discriminant wall

If `Delta=0` in the residue field, `(7)` becomes first order and its central
factor `(2n+1)x` vanishes.  Since `s!=0`, the wall equation forces `v!=0`.
The recurrence first kills `M_(p-2),M_(p-1)`, and the `Delta=0` specialization
of `(13)` then kills `X_(p-1)`.  Both homogeneous coordinates are therefore
zero by the end of the pre-reset tail:

```text
Abar|H=0,                      Lambda^2Abar=0.               (22)
```

So `(4)--(21)` do not extend across the wall by putting `chi=0`.  The
quadratic-character shear and the exterior determinant flag die together.
This is the exact failure boundary.

## 6. Scope and next bridge

The canonical object here is the determinant line of the homogeneous plane
inside the named augmented factorial state.  An arbitrary rational gauge
which mixes the boundary coordinate with that plane can destroy the displayed
splitting; Smith data alone still cannot recover a projected moment row.

The theorem supplies a rigid tail flag to the arbitrary-length exterior
continuants targeted by THM-3186.  It does not prevent different visible
paths from cancelling after projection, identify the PRS walls `H,J,K`, or
prove the empirical `floor(s/2)` depth law.  The appearance of a quadratic
character is formally analogous to parity/central-involution quotients in the
LRC stream, but no object map or implication is asserted.  No `NC(2)`,
`GMC(2)`, arbitrary-support, or `LRC(14)` conclusion follows.

## 7. Exact evidence

Run

```text
python 04-computation/quadratic_character_pre_reset_holonomy_thm3188.py
python -O 04-computation/quadratic_character_pre_reset_holonomy_thm3188.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer arithmetic only.  It checks `(8)` coefficientwise through length
twelve, exhausts every `s!=0` and every `v` over
`p=3,5,7,11,13,17,19`, separates 880 unit cases from 68 discriminant-wall
cases, and reconstructs 300 exact whole-block compound flags and squared
returns over the first five primes.  There is no floating point, random
sampling, imported executable, or assertion-sensitive test.

**QED.**
