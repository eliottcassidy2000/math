---
id: THM-3855
title: "Formal inverse-discriminant lift and algebraization gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At the THM-3808
  rational four-ray packet, the four coefficient
  gradients of the binary-cubic discriminant form a basis of all homogeneous
  cubics.  More strongly, the two base-coordinate gradients form a cubic
  complete intersection with socle degree four.  Hence every discriminant
  deformation beginning in total degree five is formally right-equivalent
  to the four-line packet by a tangent-identity base automorphism.  The
  one-place THM-3853 deformation therefore has a connected normal
  nonmonogenic formal S3 cubic completion, explicitly the base-changed
  completion of THM-3808.  Polynomial algebraization, constant units on a
  global etale open, and a Keller atlas remain open.
source: root / inverse binary-cubic discriminant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).  The
  audit independently rederived the coefficient
  gradient matrix, whose determinant is 640000, identifies its ideal with
  `(A,C)^3`, and gives a homogeneous right-inverse recursion.  It also checks
  the base-gradient resultant, complete-intersection degree-five determinant,
  explicit tangent-identity quadratic correction, and a second homogeneous
  right-equivalence recursion.  The companion replays both lifts through
  total degree twelve, with nonzero degree-thirteen residuals, and contains
  no Python asserts.  The independent audit also checked the
  degree safety of every nonlinear recursion term, the one-place
  parametrization, all height-one index-zero normality cases, henselian
  connectedness, the nonsquare/S3 step, and the index nonunit gate.  Normal
  and optimized runs byte-match the frozen 99-gate transcript and both
  recorded hashes.  A second independent 173-gate audit
  (thm3855_hostile_audit) builds
  every homogeneous coefficient multiplication map, solves independent RREF
  pivot minors, lifts both the one-place target and a dense all-monomial
  target, and rechecks normalization, branch, index, normality, connectedness
  and S3 gates.  Together they isolate lambda=0, target order four,
  characteristics 2 and 5, and finite termination as the genuine boundaries.
  THM-3855 proves no arbitrary finite-degree coefficient obstruction;
  THM-3853 supplies only the quadratic-depth exclusion.
depends_on:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3853-quadratic-depth-inverse-discriminant-one-place-gluing-obstruction
script: 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
output: 05-knowledge/results/jc2_formal_inverse_discriminant_lift_thm3855.out
script_sha256: 79be066d6e7ddd459ce3d7d3e9d81c041cbf0804561013691a73aed27c546874
output_sha256: bcbff3d657957c1e2e242621388618ccb2d9734308a1535e6a47f07b9803fc99
semantic_sha256: d8ba779da7618b29472ca255bb717a3a2c833559223d41e0404360aa47eb2b6e
independent_audit_script: 04-computation/jc2_formal_inverse_discriminant_lift_independent_audit_thm3855.py
independent_audit_output: 05-knowledge/results/jc2_formal_inverse_discriminant_lift_independent_audit_thm3855.out
independent_audit_script_sha256: ff37c33f6d1417184399e669adc3e25f90561d75b1af4f4288a1badade2dc216
independent_audit_output_sha256: 1d5b355e7576b8376a276747b29c1bb2074f6ac5690e812e6c5d1ed80c7ebe21
independent_audit_semantic_sha256: 586228efe890015d8d3e0d99eebf9996db0539033da109538a8ef95f9b374dc1
hash_basis: raw LF bytes
---

# THM-3855 -- the one-place inverse discriminant has no formal obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of
characteristic zero.  Put

```text
R=k[A,C],               m=(A,C),               Rhat=k[[A,C]].   (1)
```

For a binary cubic coefficient row `f=(a,b,c,d)`, write

```text
Disc(f)=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.                    (2)
```

At the THM-3808 row

```text
f_1=(A,C,7A,-3A),                                                (3)
```

one has

```text
Delta_0=Disc(f_1)
       =A(C+5A)(4C+19A)(3C-17A).                                (4)
```

The inverse-discriminant map at `(3)` is formally surjective with a
three-degree shift.  More strongly, it is surjective through changes of the
two base coordinates: for every `Phi in m^5 Rhat`, there is a
tangent-identity formal automorphism

```text
phi=(P,Q),             P=A mod m^2,              Q=C mod m^2,  (4a)

Delta_0(P,Q)=Delta_0(A,C)+Phi.                                 (4b)
```

Consequently, taking

```text
f_phi=(P,Q,7P,-3P)                                             (4c)
```

already gives the following coefficient lift:

> For every `Phi in m^5 Rhat`, there is an
> `eta=(eta_a,eta_b,eta_c,eta_d) in (m^2 Rhat)^4` such that
>
> ```text
> Disc(f_1+eta)=Delta_0+Phi.                                    (5)
> ```

In particular, for every `lambda in k*`, the irreducible one-place target

```text
delta_lambda=Delta_0+lambda C^5                                 (6)
```

has an exact formal binary-cubic lift with fixed linear part.
The Delone--Faddeev algebra of the specific row `(4c)` over `Rhat` is a
connected normal rank-three domain, is globally nonmonogenic over `Rhat`,
and has generic Galois closure `S3`.  It is the completed THM-3808 algebra
after the base automorphism `phi`, not a new completed singularity type.

For completeness, the one-place assertion in `(6)` is visible without using
the inverse-discriminant calculation.  On `C!=0`, put `t=A/C` and

```text
F(t)=t(1+5t)(4+19t)(3-17t).                                    (6a)
```

Then `(6)` is `C^4(F(t)+lambda C)=0` on this chart, and the curve has the
polynomial normalization

```text
C=-F(t)/lambda,                    A=-tF(t)/lambda.             (6b)
```

The ratio `A/C=t` is a rational inverse.  The four distinct roots of `F`
map to the origin, while every finite `t` remains affine and the polynomial
parametrization has only `t=infinity` on its projective completion.  Thus
`delta_lambda` is irreducible with normalization `A1` and one place at
infinity.  Its degree-four tangent cone is exactly the four distinct lines
in `(4)`.

This does **not** produce a polynomial cubic over `k[A,C]`.  THM-3853's
independently audited saturated computation excludes a lift consisting only
of `(3)` plus homogeneous quadratic corrections.  Finite termination from
degree three onward, the unit group of a global deleted-ramification open,
and any planar Keller map remain open.

## 1. The base Jacobian is a cubic complete intersection

Differentiate `(4)` in the base coordinates.  One obtains

```text
Delta_A=-2(3230A^3+567A^2C-49AC^2-6C^3),
Delta_C=-2A(189A^2-49AC-18C^2).                               (6c)
```

These two cubics are coprime.  The exact resultant and the first square
homogeneous multiplication determinant are

```text
Res_A(Delta_A,Delta_C)=36864000000 C^9,

det[(R_2)^2 --(Delta_A,Delta_C)--> R_5]=-36864000000.          (6d)
```

Thus `(Delta_A,Delta_C)` is a homogeneous complete intersection of type
`(3,3)`.  Its quotient has Hilbert series

```text
(1-z^3)^2/(1-z)^2=(1+z+z^2)^2,                                (6e)
```

so its socle degree is four and

```text
m^5 subset (Delta_A,Delta_C).                                  (6f)
```

Multiplying `(6f)` by arbitrary homogeneous forms gives a surjection

```text
(R_(n-3))^2 -> R_n,
(alpha,beta) |-> Delta_A alpha+Delta_C beta                    (6g)
```

for every `n>=5`.  For the target `C^5`, the unique first quadratic pair is

```text
alpha_2=-(1722409941/16000000)A^2
        +(3586046429/72000000)AC+(1/12)C^2,

beta_2= (26492305283/14400000)A^2
        -(2460621541/48000000)AC
        -(1211682143/72000000)C^2,                             (6h)

Delta_A alpha_2+Delta_C beta_2=C^5.
```

The large constants are not structural; `(6d)`--`(6f)` are.  They say that
the ordinary four-line singularity is four-determined under formal right
equivalence.

## 2. Homogeneous recursion gives formal right-equivalence

Write `Phi=sum_(n>=5)Phi_n`.  Suppose a tangent-identity pair `(P,Q)` has
been chosen so that `(4b)` holds through degree `n-1`, and let `E_n` be the
remaining homogeneous degree-`n` error.  By `(6g)`, choose

```text
(alpha_(n-3),beta_(n-3)) in (R_(n-3))^2,
Delta_A alpha_(n-3)+Delta_C beta_(n-3)=-E_n.                   (6i)
```

Replace `(P,Q)` by `(P+alpha_(n-3),Q+beta_(n-3))`.  At degree `n`, the
change is exactly the left side of `(6i)`.  Replacing the derivative at the
identity by the derivative at the current tangent-identity pair costs one
additional degree, and every quadratic Taylor term costs at least one more
degree when `n>=5`.  Thus no lower equation is disturbed.  Completeness of
`Rhat` gives `(4a)`--`(4b)`, and its identity linear part makes `phi` an
automorphism.

This sharpens the geometric typing.  The global polynomial curve `(6)` is
irreducible and has affine normalization `A1`, but its completed germ at the
origin is an ordinary four-branch point.  The formal automorphism identifies
that germ with the completed union of the four lines `(4)`; it does not
identify the global curves.

Nor can this particular formal solution terminate as a polynomial base
map.  If `P,Q in k[A,C]` have the identity linear part, then

```text
Delta_0(P,Q)=P(Q+5P)(4Q+19P)(3Q-17P)                           (6j)
```

is a product of four nonconstant polynomials with distinct linear parts.
It cannot equal the irreducible `delta_lambda`.  Hence any eventual
polynomial coefficient lift must leave the base-change orbit and become a
genuine four-coefficient deformation.

## 3. The coefficient gradients span every cubic

Differentiate `(2)` with respect to `(a,b,c,d)` and evaluate at `(3)`.  In
the ordered monomial basis

```text
(A^3,A^2C,AC^2,C^3),                                            (7)
```

the four resulting homogeneous cubics are the columns of

```text
M=[[-1858, -378, -588,  162],
   [ -378,   98,  -54,  126],
   [    0,   36,   14,    0],
   [    0,    0,    0,   -4]].                                  (8)
```

Directly,

```text
det(M)=640000 !=0.                                               (9)
```

Thus the coefficient gradients form a basis of `R_3`, and their homogeneous
ideal is exactly

```text
(Disc_a(f_1),Disc_b(f_1),Disc_c(f_1),Disc_d(f_1))=m^3.          (10)
```

This is stronger than a single favorable tangent direction: multiplication
by arbitrary forms gives a surjection

```text
(R_(n-3))^4 -> R_n,       u |-> sum_i Disc_i(f_1)u_i            (11)
```

for every `n>=3`.

## 4. An independent coefficient recursion proves the formal lift

Write `Phi=sum_(n>=5) Phi_n` by total degree.  Suppose corrections have been
chosen so that `(5)` is true through degree `n-1`, and let `E_n in R_n` be
the remaining degree-`n` error.  By `(11)`, choose a homogeneous coefficient
correction

```text
u_(n-3) in (R_(n-3))^4,
sum_i Disc_i(f_1)(u_(n-3))_i=-E_n.                              (12)
```

The degree-`n` change in the discriminant is exactly the left side of
`(12)`.  Indeed, all previous corrections have degree at least two, so
replacing the derivative at `f_1` by the derivative at the current row costs
at least one extra total degree; terms quadratic in `u_(n-3)` do likewise.
The new row therefore solves `(5)` through degree `n` without disturbing
lower degrees.

Starting at `n=5` and using completeness of `Rhat` gives a convergent
`m`-adic coefficient row `f_1+eta` satisfying `(5)` exactly.  No analytic or
finite-degree convergence is asserted.

For the target `(6)`, the inverse of `(8)` gives the particularly concrete
first correction

```text
eta_2=lambda C^2(
  91449/40000,
 151263/40000,
-194481/20000,
          -1/4).                                                (13)
```

It creates exactly `lambda C^5` at degree five.  The deterministic companion
uses a fixed monomial right inverse of `(8)` to continue through degree
twelve.  The degree-thirteen error of that finite truncation is nonzero, as
it should be: only the infinite recursion has been proved exact.

## 5. The formal cubic passes the index, normality, and S3 gates

For the coordinate lift `(4c)`, these properties can first be read by
inheritance.  Let `S_0` be the completion at the vertex of the normal
third-Veronese cubic algebra in THM-3808.  Base change along the automorphism
`phi` gives

```text
S_phi=S_0 completed_tensor_(Rhat,phi) Rhat.                     (13a)
```

Because `phi` is an automorphism, `(13a)` is the same completed algebra with
its base coordinates relabelled.  In particular it is a normal local domain,
and its generic cubic and index packet are inherited.  The intrinsic argument
below independently checks the same conclusion and applies to any formal
coefficient lift with discriminant `(6)` and coefficients in `m`.

The binary index form of the Delone--Faddeev algebra is

```text
I(x,y)=-(a x^3+b x^2y+c xy^2+d y^3).                            (14)
```

Every coefficient of `f_1+eta` belongs to `m`, because `(3)` is linear and
`eta` starts in `m^2`.  Hence

```text
I(Rhat^2) subset m,                                             (15)
```

so the index form represents no unit and the formal cubic algebra is not
monogenic.

For `(6)`, the tangent cone `(4)` is a product of four distinct lines.
Consequently `delta_lambda` is reduced in `Rhat`, with four smooth formal
branches of multiplicity one.  The Delone--Faddeev algebra is finite free,
hence `S2`, and is generically etale.  Away from `(6)` it is etale in
codimension one; at each branch the discriminant valuation is one, so the
discriminant-index formula forces height-one index zero in the integral
closure.  It is therefore normal.

Modulo `m`, all four coefficients vanish and the special fibre is

```text
k[omega,theta]/(omega,theta)^2,                                 (16)
```

which is local.  Since `Rhat` is henselian, the finite algebra has no
nontrivial idempotent.  A connected normal finite algebra is a domain.
Finally, `(6)` is not a square in `Frac(Rhat)` because each of its four
height-one factors has odd valuation.  The connected generic cubic therefore
has nonsquare discriminant and Galois closure `S3`.

## 6. Exact boundary and next construction step

The formal packet now passes four gates simultaneously:

```text
one-place global discriminant target,
connected normal cubic order,
nonmonogenic index form,
generic S3 monodromy.                                           (17)
```

The missing implication is algebraization.  A finite polynomial row would
have to realize `(6)` exactly while retaining `(15)`.  THM-3853 excludes
maximum coefficient degree two; degree-three and interacting higher-degree
rows are the first live search.  Even after such a row is found, its maximal
etale open must still have only scalar units and admit the required plane
atlas.  Thus `(5)` is a positive local existence theorem and a sharp
algebraization gate, not a Jacobian counterexample.

Reproduction:

```bash
python3 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
python3 -O 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
python3 04-computation/jc2_formal_inverse_discriminant_lift_independent_audit_thm3855.py
python3 -O 04-computation/jc2_formal_inverse_discriminant_lift_independent_audit_thm3855.py
```

The primary modes byte-match the strengthened frozen 99-gate transcript.
The independent 173-gate RREF recursion lifts the target through degree twelve and
a dense target through degree nine, while independently recovering every
formal-order gate above.  Its normal and optimized runs byte-match its
frozen output.  The separate hostile audit of the right-equivalence extension
also passes.  **QED.**
