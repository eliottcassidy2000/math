---
id: THM-2681
title: "THM-1310 S3 normalization and quartic A4/S4 resolvent exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the target
  chart g!=0 (written c!=0 below), the full Galois normalization of the
  actual THM-1310 cubic root field is the smooth factorial threefold
  Spec C[c^+-1,r,s]=G_m x A^2, with S3 acting by
  rho(r,s)=(zeta r,zeta^-1 s) and tau(r,s)=(s,r).  The explicit root
  x=6c/(r^2+rs+s^2) has stabilizer <tau> and therefore has this entire
  degree-six field as its Galois closure.  The chart has procyclic etale
  fundamental group, unit-squareclass rank one, and trivial class group,
  so it admits no connected etale V4 torsor.  THM-2655 therefore excludes
  this actual cubic field as the classical resolvent-root field of a
  dimension-three quartic S4 Keller map under a polynomial target
  automorphism or explicit base-ring identification.  The A4 case is
  excluded even earlier: its canonical V4-fixed matching field is cyclic
  Galois of degree three, whereas the actual THM-1310 cubic root field is
  nonnormal with S3 closure.  This field-type obstruction is over the
  original target field and need not survive a discriminant quadratic base
  change.  Arbitrary transport and matches of only the discriminant, cusp,
  Jelonek, or odd-valuation shadows remain open; no general A4/S4 exclusion
  or Jacobian conjecture follows.
source: root-2026-07-28-thm1310-resolvent-normalization
depends_on:
  - THM-1310-jacobian-counterexample-fiber-geometry-S3-resolvent-jelonek
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2570-jacobian-cusp-submersion-and-unramified-fold-law
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
script: 04-computation/jacobian_thm1310_resolvent_galois_open_v4_no_go.py
output: 05-knowledge/results/jacobian_thm1310_resolvent_galois_open_v4_no_go.out
script_sha256: b1a4ad05e86bf449fb8f6ef54dc22c6d7661ad4521140972b052d325848b6b4d
output_sha256: d16de1daf7b82c42a936c86fdfe309f617aa7e0a8c52ab0e47b9db6f3b91fc4a
secondary_script: 04-computation/jacobian_thm1310_a4_resolvent_field_type_referee.py
secondary_output: 05-knowledge/results/jacobian_thm1310_a4_resolvent_field_type_referee.out
secondary_script_sha256: 55c72c3832c6cd7ff62208f29568b71a5b4672a2369e8541b4739923467d96bc
secondary_output_sha256: f41000269e672f4ac66a5294801345494912e6d31859f155d5a2bd62a2ccd933
hash_basis: LF-normalized bytes
---

# THM-2681 -- the THM-1310 resolvent chart cannot carry the quartic kernel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2655 turns a live quartic `S4` Keller map into a topological demand on
the **full** degree-six Galois normalization of its classical cubic
resolvent: that normalization must carry a connected etale `V4` torsor on
its regular locus.  The discriminant of the THM-1310 cubic had long suggested
that it might furnish this resolvent layer.  Here the actual cubic root field,
not merely its discriminant, is normalized explicitly.  Its simplest dense
chart is too topologically small to meet the demand.

Throughout, `c` denotes the target coordinate called `g` in THM-1310.

## 1. Recoordinate the target chart

Let

```text
K=C(a,b,c),
p=4-3bc,
Q=27ac^2-9bc+8,                                         (1)
```

and let

```text
L=27a^2c^2-18abc+16a+b^3c-b^2.                         (2)
```

THM-1310's generic fibre coordinate has minimal polynomial

```text
N(X)=L X^3+pX-2c.                                       (3)
```

On `c!=0`, the changes of variables

```text
b=(4-p)/(3c),
a=(Q-3p+4)/(27c^2)                                     (4)
```

are polynomial in `c^+-1,p,Q`.  Thus the target coordinate ring is exactly

```text
A=C[a,b,c,c^-1]=C[c^+-1,p,Q].                           (5)
```

The THM-1310 cusp identity becomes

```text
Q^2-p^3=27c^2L,                                        (6)
```

and `Disc_X(N)=-4Q^2L`.

## 2. The full `S3` normalization

Choose constants `eta,zeta in C` with

```text
eta^2=-3,             zeta=(-1+eta)/2.                  (7)
```

On the discriminant cover put `w^2=-L` and define the Cardano factors

```text
U=Q+3 eta c w,          V=Q-3 eta c w.                  (8)
```

Equation (6) gives

```text
UV=p^3.                                                   (9)
```

Normalize by adjoining `r,s` with

```text
r^3=U,       s^3=V,       rs=p.                         (10)
```

Then

```text
B=C[c^+-1,r,s]                                          (11)
```

carries the faithful action

```text
rho(r,s)=(zeta r,zeta^-1 s),
tau(r,s)=(s,r),
rho^3=tau^2=1,          tau rho tau=rho^-1.              (12)
```

This is not just a convenient degree-six extension.  A monomial fixed by
`<rho>` has exponent difference divisible by three, so

```text
B^<rho>=C[c^+-1,r^3,s^3,rs]
        =C[c^+-1,U,V,p]/(UV-p^3).                       (13)
```

The transposition swaps `U,V` and fixes `p`; its invariant ring in (13) is
generated by `p` and `U+V=2Q`.  Therefore

```text
B^S3=C[c^+-1,p,Q]=A.                                    (14)
```

Since `B` is a normal finite `A`-algebra and its fraction field has the six
distinct transformations (12), it is the normalization of (5) in this
`S3` extension.

## 3. It is the *actual* THM-1310 splitting field

Inside `Frac(B)`, equations (8)--(10) give

```text
w=(r^3-s^3)/(6 eta c),
x=(r-s)/(eta w)=6c/(r^2+rs+s^2).                        (15)
```

Direct substitution proves

```text
Lx^3+px-2c=0.                                           (16)
```

Moreover `tau(x)=x`, whereas `rho(x)!=x` as a rational function.  Hence

```text
Stab_S3(x)=<tau>,       [K(x):K]=3.                      (17)
```

The core of a transposition subgroup in `S3` is trivial.  Consequently the
Galois closure of the degree-three field `K(x)` is all of `Frac(B)`.  This
is the load-bearing upgrade over a shared-discriminant calculation: (11) is
the full Galois normalization of the specific root field in THM-1310.

## 4. The topology is too small for `V4`

The normalization chart is

```text
Spec(B)=G_m x A^2.                                      (18)
```

It is smooth and factorial, with

```text
B^*=C^* c^Z,             Cl(B)=0,                       (19)
H^1_et(Spec(B),mu_2)=B^*/B^(*)2=F_2.                    (20)
```

Equivalently,

```text
pi_1^et(Spec(B))=Z-hat,                                 (21)
```

which has no quotient isomorphic to the noncyclic group `V4`.  Thus (18)
admits no connected etale `V4` torsor.

Now suppose a degree-four polynomial Keller map `A^3_C -> A^3_C` has
monodromy `S4` and its classical cubic-resolvent root field is identified,
over the target coordinate ring, with the THM-1310 field above.  THM-2655
produces a connected quasi-etale `V4` cover of the full `S3` Galois-resolvent
normalization, etale over its regular locus.  Restrict it to `c!=0`.
Equations (11)--(18) say that this entire chart is regular.  The restricted
total space remains nonempty and integral--it is a localization of the
integral full normalization--and hence remains connected.  It would
therefore be a connected etale `V4` torsor over (18), contradicting
(20), or equivalently (21).

Therefore:

```text
the actual THM-1310 cubic root field cannot be the classical
resolvent-root field of a dimension-three quartic S4 Keller map.           (22)
```

## 5. The `A4` resolvent is the wrong field type

Let `M/K` be the Galois closure of a separable transitive quartic with
generic group `A4`.  Its action on the three perfect matchings of four roots
has kernel `V4` and regular image

```text
A4/V4=C3.                                                   (23)
```

Consequently the canonical matching-resolvent field

```text
E=M^V4                                                     (24)
```

is cyclic Galois of degree three.  Equivalently, a root of the canonical
cubic resolvent already generates its own cubic-resolvent splitting field
`E`.  By contrast,
Section 3 proves that the THM-1310 field `K(x)` has transposition stabilizer
inside its `S3` closure.  It is a nonnormal cubic whose Galois closure has
degree six.  Hence

```text
K(x) is not K-isomorphic to an A4 canonical cubic-resolvent field.         (25)
```

There is no reducible-resolvent loophole under the stated generic `A4`
hypothesis: the quotient `C3` acts transitively on the three matchings.
Specializations where the group drops are outside `(25)`.  Nor is `(25)`
invariant under arbitrary base extension.  After adjoining the discriminant
quadratic field `K(sqrt(-L))`, the `S3` closure in Section 2 becomes cyclic
of degree three.  Thus `(25)` is an exact target-field/base-ring obstruction,
not a claim that the two cubics remain distinct after every extension.

## 6. Sharp local Kummer hostile at `a!=0`

The global obstruction in Section 4 cannot be checked only at the generic
root torus.  A useful presentation of the same full `S3` normalization is

```text
B'=C[c^+-1,r_1,r_2,r_3]/(r_1+r_2+r_3-1),                  (26)
```

where the `r_i` are the ordered roots of

```text
f(T)=T^3-T^2+(bc/4)T-ac^2/4.                              (27)
```

Indeed, with `h=6T-2`, direct expansion gives

```text
216 f((h+2)/6)=h^3-3ph-2Q,                                (28)
```

the Cardano cubic normalized in Section 2.  Hence uniqueness of
normalization identifies `(26)` with `(11)`.

On `D(a)` inside this chart,

```text
r_1 r_2 r_3=ac^2/4                                       (29)
```

is a unit, so all three `r_i` are units.  The two squareclasses

```text
[r_1/r_3], [r_2/r_3]                                     (30)
```

are independent and span the standard even-weight plane in `F_2^3`.
Adjoining their square roots therefore gives a connected finite-etale
`S3`-equivariant `V4` Kummer torsor on `D(a)`.  This is the sharp generic
hostile to any proof which inspects only the all-roots-nonzero chart.

It does not extend as a quasi-etale carrier over the full smooth chart.
The divisors

```text
D_i={r_i=0}                                               (31)
```

are irreducible, since `B'/(r_i)` is a Laurent polynomial ring in `c` and
one remaining root.  Moreover

```text
div(r_1/r_3)=D_1-D_3,       div(r_2/r_3)=D_2-D_3.         (32)
```

Thus the normalization of `B'` in the local `V4` field is ramified in
codimension one along `D_1 union D_2 union D_3`; it is neither etale nor
quasi-etale over `Spec(B')`.  Restoring precisely the projection-completion
divisors above `a=0` kills the tempting standard plane.  Globally `(19)`--
`(20)` leave only the one-dimensional `S3`-trivial squareclass line generated
by `[c]`.  This is a
local-carrier/global-divisor boundary, not an `A4` hostile and not a Keller
example.

## 7. Exact scope and surviving shadows

Statement (22) is stable under a polynomial target-coordinate automorphism,
and under any explicit base-ring isomorphism which identifies the relevant
normalizations and the chart (5).  It is **not** a statement about arbitrary
birational transport of function fields: such a transport need not preserve
the affine normalization or the open topology used in Section 4.

In THM-2465 language, this closes the intended exact `N1`--`N5` realization
in which `N5` identifies the actual degree-six THM-1310 Galois-resolvent
layer.  The weaker `N2`--`N4` shadows--shared discriminant squareclass,
embedded Jelonek equation, and odd divisorial valuation--remain viable.
Neither the identity (6), the cusp shape, the factor `-L`, nor matching
Jelonek data alone identifies the field or supplies the topology in (18).

The theorem does not exclude arbitrary quartic `A4` or `S4` monodromy, an
unrelated `C3`/`S3` resolvent normalization, or quartics whose resolvents only
resemble THM-1310 through those shadows.  It proves no general degree-four
exclusion, `G1`, `JC(2)`, `JC(3)`, general Jacobian conjecture, or Dixmier
conjecture.

## 8. Reproduction and audit

Run

```bash
python3 04-computation/jacobian_thm1310_resolvent_galois_open_v4_no_go.py
python3 -O 04-computation/jacobian_thm1310_resolvent_galois_open_v4_no_go.py
```

Both executions must byte-match

```text
05-knowledge/results/jacobian_thm1310_resolvent_galois_open_v4_no_go.out.
```

The companion performs `38` explicit hard checks: the coordinate inverse,
cusp and discriminant identities, Cardano norm, cyclic invariant residues,
normalization coordinates, actual-root formula and cubic equation, the six
`S3` transformations and stabilizer test, and the Kummer-rank obstruction.
The invariant-ring, integral-normalization, stabilizer-core, class-group,
fundamental-group, localization-connectedness, and THM-2655 steps were also
independently rederived by two hostile audits.  Normal and optimized
executions byte-match the stored transcript and the hashes above.

The independent `A4` field-type and local-divisor referee is

```bash
python3 04-computation/jacobian_thm1310_a4_resolvent_field_type_referee.py
python3 -O 04-computation/jacobian_thm1310_a4_resolvent_field_type_referee.py
```

Both runs byte-match
`05-knowledge/results/jacobian_thm1310_a4_resolvent_field_type_referee.out`.
It enumerates the exact `S4/A4/V4` matching actions and stabilizers, checks
the nonnormal `S3` root-field fingerprint, and verifies the standard
two-dimensional divisor-parity plane in `(30)`--`(32)`.  The field and
codimension-one extension arguments were independently hostile-audited.

QED.
