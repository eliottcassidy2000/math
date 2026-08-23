---
id: THM-3883
title: "Marked-root carrier descent reduces locally to sign matching and one cube"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For an arbitrary reduced irreducible affine
  plane curve with a regular normalization root z^2=1+(2/3)AC, a fixed-sign
  polynomial carrier exists exactly when: z has the chosen sign above A=0;
  all nonzero root residues agree on each singular fibre; and at every
  zero-root singularity z^3 belongs to the completed local branch ring.
  Matching roots at A=0 descend by an all-order binomial cancellation, and
  common nonzero roots descend by Hensel uniqueness.  This subsumes the
  node/A2 criterion of THM-3880 and isolates one square/cube sidecar at every
  higher singularity.
source: jc_sparse_direct_search / post-THM-3880 arbitrary-singularity endpoint, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit rederived the
  fixed-sign rational carrier, checked Hensel uniqueness in the possibly
  multibranch completed local ring, proved the all-order `A^2` binomial
  cancellation, and verified that the zero-root seam is equivalent exactly
  to cube descent.  It also checked finite-normalization/completion
  globalization, all smooth/vertical/unit boundaries, and both plane-branch
  controls.  The assertion-free exact companion checks
  the cusp identity, the completed binomial quotient through order ten with
  every A-order nonnegative, the simple-root factorization, the zero-root
  carrier reduction, and the sharp k[[t^2,t^3]] / k[[t^2,t^5]] cube controls.
  Normal and optimized runs byte-match the frozen 43-gate transcript.  The
  all-order proof is the formal binomial series and Hensel uniqueness, not
  the bounded replay.
depends_on:
  - THM-3880-marked-root-opposite-sign-normalization-fibre-nondescent
related:
  - THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate
script: 04-computation/jc2_marked_root_complete_local_square_cube_thm3883.py
output: 05-knowledge/results/jc2_marked_root_complete_local_square_cube_thm3883.out
script_sha256: 10cc1364cc90c21e6848ca9b9914c06619b6bfa19c22767a757251dce2d3c630
output_sha256: d0e9558e6242147d678ba860627aff96f45fd6d8de50d27ca0238ba3f3cc2ba3
semantic_sha256: f091ccb803702bc6b3bb6acec3ab700f705bf7f8017c36692093c7cbc181ceb0
hash_basis: raw LF bytes
---

# THM-3883 -- the complete-local sign/square/cube criterion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b,                       (1)
P=1+(2/3)AC.                                                   (2)
```

Let `Gamma subset A2_(A,C)` be an arbitrary irreducible reduced affine
curve, not the vertical line `A=0`, and let

```text
nu: X -> Gamma                                                  (3)
```

be its normalization.  Assume there is a regular marked root

```text
z in O(X),                         z^2=P.                       (4)
```

Fix `epsilon in {+1,-1}`.  For a singular point `x in Gamma`, write

```text
R_x=O_hat_(Gamma,x),
S_x=product_(p in nu^-1(x)) O_hat_(X,p),                       (5)
```

with the natural injection `R_x -> S_x`.  Then there is a polynomial
`b_epsilon in k[A,C]` satisfying

```text
Delta_(b_epsilon)|_Gamma=0,
1+AC+A^2b_epsilon=epsilon z^3 on X                             (6)
```

if and only if all three conditions hold:

```text
(A0) for every p in X with A(p)=0,               z(p)=epsilon;

(NZ) for every singular x with A(x)!=0 and P(x)!=0, all residues z(p),
     p in nu^-1(x), are equal;

(Z3) for every singular x with P(x)=0,             z^3 in R_x. (7)
```

In `(Z3)`, the function `z^3` is viewed in the normalized completed ring
`S_x`.  Since `z^2=P` already belongs to `R_x`, this is exactly one missing
square/cube descent condition, not an unspecified higher-jet family.

The theorem makes no smoothness, nodal, cuspidal, seminormality, or degree
assumption on `Gamma`.  It does assume that the normalization square root
`z` exists.  It classifies the marked component itself, not the projective
normalizations of the other components of `V(Delta_b)` and not a planar
Jacobian mate.

## 1. Global sign and the forced rational carrier

The universal identity is

```text
A^2 Delta_b=27[P^3-(1+AC+A^2b)^2].                            (8)
```

Indeed, on the normalization the identity `(8)` factors as
`(u-z^3)(u+z^3)=0`, where `u=1+AC+A^2b`.  Since `X` is integral, any carrier
has one global sign.  On the dense open set `A!=0`, a sign-`epsilon` carrier
must therefore have the unique value

```text
B_epsilon=(epsilon z^3-1-AC)/A^2 in k(X).                      (9)
```

Once `(9)` is regular on `X`, its defining identity and `(8)` imply
`Delta_(B_epsilon)=0`, because `A` is a nonzero function in the domain.
Thus the proof has two exact tasks: regularity on `X`, then descent from
`O(X)` to `O(Gamma)`.

## 2. The A=0 fibre: matching sign gives every completed jet

At a normalization point above `A=0`, one has `P=1` and `z=+1` or `-1`.
THM-3880 proves the exact finite pole criterion

```text
B_epsilon regular above A=0             iff z=epsilon there.  (10)
```

For the arbitrary-singularity theorem we need more: the matching sign must
place `B_epsilon` in the completed **branch ring**, not merely in every
normal branch.

Fix a point `x in Gamma` with `A(x)=0` and assume `(A0)`.  In the complete
local ring `R_x`, put

```text
T=(2/3)AC in m_x.                                              (11)
```

There is a unique Hensel root congruent to one,

```text
w=(1+T)^(1/2) in R_x.                                         (12)
```

In every factor of `S_x`, the function `epsilon z` is also a root of
`W^2-(1+T)` with residue one.  The derivative is two, a unit; uniqueness
gives

```text
w=epsilon z in S_x.                                           (13)
```

Now use the all-order binomial expansion:

```text
w^3-1-AC
 =(1+T)^(3/2)-1-(3/2)T
 =sum_(j>=2) binom(3/2,j) T^j.                                (14)
```

Substitute `(11)` and divide by `A^2`:

```text
B_epsilon
 =sum_(j>=2) binom(3/2,j)(2/3)^j A^(j-2) C^j in R_x.          (15)
```

Every exponent `j-2` is nonnegative, so `(15)` proves all conductor jets at
once, for an arbitrary singularity and arbitrary orders of `A,C`.  The first
term is `C^2/6`, agreeing with the direct `A=0` specialization of `Delta`.

At the wrong sign, the numerator of `(9)` is the unit `-2`, so there is a
pole of order `2 ord_p(A)`.  This proves both necessity and sufficiency of
`(A0)` and closes the finite regularity gate.

## 3. Nonzero marked roots descend by Hensel uniqueness

Let `x` be singular with `A(x)!=0` and `P(x)!=0`.  Condition `(NZ)` gives a
common nonzero residue

```text
z(p)=lambda in k^*                         for all p|x.        (16)
```

The polynomial `W^2-P in R_x[W]` has the simple residue root `lambda`.
Hensel's lemma supplies a unique

```text
w in R_x,                         w^2=P, w mod m_x=lambda.     (17)
```

Its image in each factor of `S_x` equals `z`: both are roots with the same
residue, and

```text
(w-z)(w+z)=0                                                   (18)
```

with `w+z` a unit.  Hence `z in R_x`.  Since `A` is a unit at `x`, formula
`(9)` now gives

```text
B_epsilon in R_x.                                             (19)
```

Necessity of `(NZ)` is equally exact.  If a carrier descends, then
`epsilon z^3=1+AC+A^2b` has one residue on the singular fibre.  The possible
nonzero values of `z` are `lambda,-lambda`, whose cubes are opposite, so all
root residues must agree.

Points with `A(x)=0` were already treated by `(A0)` and Section 2, so this
covers every nonzero-root singularity.

## 4. At a zero marked root, the cube is the whole obstruction

Suppose `P(x)=0`.  Then

```text
A(x)C(x)=-3/2,                                                (20)
```

so both `A` and `C` are units in `R_x`, and every normalization branch has
residue `z=0`.  Formula `(9)` becomes

```text
B_epsilon=A^-2[epsilon z^3-(1+AC)].                           (21)
```

The second term lies in `R_x` and `A^-2` is a unit.  Therefore

```text
B_epsilon in R_x                         iff z^3 in R_x.       (22)
```

This proves sufficiency of `(Z3)`.  It also proves necessity directly from
`(6)`: if `b` descends, then

```text
epsilon z^3=1+AC+A^2b in R_x.                                 (23)
```

Although `z^2=P` always descends, `z^3` need not.  Equation `(22)` is the
precise higher-contact sidecar which the pointwise sign packet cannot see.

## 5. Completion globalizes the three local tests

Put `R=O(Gamma)` and `S=O(X)`.  Normalization is finite, so `S/R` is a
finite `R`-module supported at the singular locus.  Condition `(A0)` first
makes `(9)` a regular element of `S`.  Sections 2--4 prove that its image
belongs to `R_x` in every completed singular stalk.

Completion is faithfully flat.  The class of `B_epsilon` in `S/R` therefore
vanishes after localization and completion at every maximal ideal in its
support, hence vanishes globally.  Thus

```text
B_epsilon in R.                                               (24)
```

Because `Gamma` is a closed affine plane curve, every element of `R` is the
restriction of a polynomial in `k[A,C]`.  Choose a lift `b_epsilon`; equations
`(8)--(9)` give `(6)`.  This proves sufficiency of `(7)`.  Necessity was
proved in Sections 2--4, completing the iff.

## 6. Sharp complete-local controls

The zero-root cube condition is neither automatic nor a disguised demand
that `z` itself descend.

For the positive control, take

```text
R_+=k[[t^2,t^3]],
A=1+t^3,
C=(3/2)(t^2-1)/(1+t^3),
z=t.                                                          (25)
```

Then `P=t^2=z^2` and `z^3=t^3 in R_+`, so `(Z3)` holds even though `z` is
not in the cusp ring.  These are genuinely plane-curve coordinates: writing
`v=A-1=t^3` and `c=(2/3)(C+3/2)`, one recovers

```text
t^2=(1+v)c-v.                                                 (26)
```

For the hostile control, take

```text
R_-=k[[t^2,t^5]],
A=1+t^5,
C=(3/2)(t^2-1)/(1+t^5),
z=t.                                                          (27)
```

Again `P=t^2`, but the semigroup `<2,5>` has gaps `1,3`, so

```text
z^2=t^2 in R_-,                         z^3=t^3 notin R_-.     (28)
```

The corresponding forced carrier has a nonzero `t^3` coefficient and does
not descend.  All exponents at least four lie in `<2,5>`, so `(28)` locates
the first and only relevant low cube gap rather than relying on a bounded
coefficient search.  Here too the displayed coordinates generate the full
completed plane branch: with `v=A-1=t^5` and the same definition of `c`,
formula `(26)` again recovers `t^2`.

THM-3864 supplies a global multi-cusp realization of the positive
square/cube phenomenon: its seminormal defect elements need not lie in the
branch ring even though both their squares and cubes do.  THM-3883 explains
why the cube, rather than the root itself, is the exact carrier observable.

## 7. Replay and scope

Run

```text
python3 04-computation/jc2_marked_root_complete_local_square_cube_thm3883.py
python3 -O 04-computation/jc2_marked_root_complete_local_square_cube_thm3883.py
```

and compare both streams byte-for-byte with the frozen output.  The
companion checks the binomial quotient through depth ten and the complete
low semigroup packet.  Equations `(11)--(15)` are the all-order proof, and
Hensel uniqueness plus the finite-module completion argument give the
unbounded singularity scope.

This theorem completes affine descent for every reduced irreducible marked
curve once a regular normalization square root exists.  It does not decide
whether `P` admits such a root on a proposed curve, does not classify the
other discriminant components forced by the chosen polynomial lift, and
does not produce or exclude a planar Keller map by itself.
