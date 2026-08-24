---
id: THM-3923
title: "Binary-cubic four-ray tangent packet cannot algebraize to a cubic Keller atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a reduced
  irreducible plane branch, r distinct tangent-cone factors give at least r
  normalization addresses over the marked point. A degree-three Keller
  completion permits at most three such addresses. Consequently an
  irreducible binary-cubic discriminant with four distinct tangent rays at
  one affine point cannot occur in a planar cubic Keller completion. The
  THM-3808 fixed linear packet always has the four rays
  A(C+5A)(4C+19A)(3C-17A), so every irreducible polynomial algebraization of
  the THM-3855 formal coefficient germ is excluded. In particular the two
  THM-3853 one-place inverse-discriminant targets and every tangent-identity
  formal deformation retain the fatal four-address packet. Formal
  inverse-discriminant surjectivity remains true; it is now proved to point
  into the wrong local branch grammar for a cubic Keller atlas.
source: jc_zero_debt_lift / post-THM-3920 tangent-cone address lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (incoming_thm3923_3924_audit/root,
  2026-08-23).  The audit independently reconstructed the completed-germ
  Hensel/normalization bridge, checked the reduced/irreducible and
  per-component scope, derived the degree-three address contradiction, and
  verified its application to both THM-3853 targets and all THM-3855 fixed-
  linear-packet algebraizations.  The assertion-free companion verifies the
  binary-cubic discriminant, all six pairwise ray determinants, invariance
  under arbitrary quadratic coefficient jets and tangent-identity base jets,
  and both four-address parametrizations.  LF-normalized normal and optimized
  streams match the frozen LF output in all 28 gates; script, output, and
  semantic hashes agree.  No repair was required.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3853-quadratic-depth-inverse-discriminant-one-place-gluing-obstruction
  - THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
related:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3865-one-place-inverse-discriminant-resolvent-class-group
  - THM-3887-binary-cubic-common-zero-quintic-one-tangent-obstruction
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
script: 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
output: 05-knowledge/results/jc2_binary_cubic_four_ray_tangent_address_thm3923.out
script_sha256: 3232643254947a0970ec60355d95056b023be4943dff09ca0f04c01e1624c878
output_sha256: 9856a88626944ab9821a889554b85f29d38ccc0d153995f179beb016880e5fbb
semantic_sha256: 264aef4f00430ba6db9215f5d26e9322f8d06813a7b643e692570eb92359c9e8
hash_basis: raw LF bytes
---

# THM-3923 -- four formal rays are one ray too many

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero.

The theorem isolates the missing local invoice in the inverse-discriminant
program:

```text
formal solvability of the discriminant equation
       does not change the normalization-address packet;
cubic Keller completion
       can carry at most three addresses on one branch component.       (1)
```

## 1. Tangent directions inject into normalization addresses

Let `Gamma=V(delta) subset A2_(A,C)` be a reduced irreducible plane curve,
and let `o=(0,0)`. Write

```text
delta=delta_m+terms of order >m,             delta_m!=0,              (2)
```

and suppose the homogeneous tangent cone `delta_m` has `r` distinct linear
factors over `k`. Then the normalization `Gamma^nu` has at least `r`
distinct points above `o`.

Indeed, factor the reduced completed germ in `k[[A,C]]`. The initial form of
an irreducible plane branch is a power of one linear form: if it contained
two coprime homogeneous factors, Hensel factorization would split the germ.
Thus one formal branch can account for only one tangent direction. Distinct
tangent factors therefore require distinct formal branches, and the branches
of a reduced excellent plane curve germ are exactly the points of its
normalization over `o`.

When `m=r=4` and all four factors are distinct, the conclusion is sharper.
Hensel lifting splits `delta` into four factors having pairwise distinct
linear terms. Each factor is smooth, so `o` is an ordinary four-branch point
and has exactly four normalization addresses.

Both adjectives in the statement matter. If the curve is nonreduced, a
factor of the displayed equation need not be a branch of its reduced
support. If the global branch is reducible, different local rays may belong
to different global components, and the per-component cubic address bound
below need not see all of them at once.

## 2. The cubic address obstruction

Suppose `Gamma` is an irreducible component of the branch divisor of the
normal finite-flat completion of a hypothetical degree-three Keller map

```text
f:A2_source -> A2_(A,C).                                      (3)
```

THM-3801 proves that the prime over `Gamma` omitted from the etale source
has decomposition type `(2,1)` and residue degree one. Hence its
normalization has the same function field as `Gamma^nu`, and every
normalization address of `Gamma` over an affine target point is an address
of that one irreducible boundary divisor. THM-3920's finite-flat boundary
cap gives

```text
#(Gamma^nu over o) <=3.                                      (4)
```

Combining Sections 1 and 2 proves the general local gate

```text
reduced irreducible cubic branch in a Keller completion
       ==> at most three distinct tangent directions at every affine point.
                                                                    (5)
```

In particular, an ordinary four-branch affine point is impossible. This is
a statement about the normalization of one global branch component, not
merely about the number of factors in a tangent polynomial.

## 3. The fixed THM-3808 linear packet always carries four rays

Consider a binary cubic index form

```text
Phi=aX^3+bX^2Y+cXY^2+dY^3                                (6)
```

whose coefficient row at `o` begins with the THM-3808 linear packet

```text
(a,b,c,d)=(A,C,7A,-3A) mod (A,C)^2.                       (7)
```

The discriminant is homogeneous of coefficient degree four:

```text
Disc(Phi)=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.             (8)
```

Therefore its degree-four base jet depends only on `(7)` and is

```text
Delta_0=A(C+5A)(4C+19A)(3C-17A).                          (9)
```

The four displayed lines are pairwise distinct. Arbitrary coefficient terms
in `(A,C)^2` can change only degree five and higher in the discriminant.
Consequently:

> If a polynomial binary-cubic order has linear packet `(7)` and its
> discriminant is reduced and irreducible, then its branch has four
> normalization addresses over `o`; its finite completion cannot contain a
> same-field affine-plane Keller source.

Normality or existence of the polynomial order is not inferred from the
discriminant equation. The implication says that even a successful normal
polynomial algebraization inside this fixed linear packet would fail the
independent boundary-address gate.

## 4. The THM-3853 one-place targets are excluded

For `L=C` or `L=A+C`, THM-3853 considers

```text
delta_(L,lambda)=Delta_0+lambda L^5,             lambda!=0.   (10)
```

It proves that `(10)` is irreducible, its affine normalization is `A1`, and
its projective normalization has one place at infinity. In complementary
coordinates `(M,L)`, write

```text
Delta_0(M,L)=L^4 D_L(M/L).                                  (11)
```

The two address polynomials are

```text
D_C(t)    =-t(5t+1)(17t-3)(19t+4),
D_(A+C)(t)=-t(4t+1)(15t+4)(20t-3).                          (12)
```

Each has four distinct roots. The normalization is explicitly

```text
L=-D_L(t)/lambda,                    M=tL,                  (13)
```

so all four roots map to the affine origin. Thus both targets violate `(4)`.
No normal finite-flat cubic algebra with either discriminant can be the
finite completion of a degree-three Keller map. This conclusion is
independent of THM-3853's bounded coefficient-lift no-go and THM-3865's
quadratic-resolvent class-group obstruction.

More generally, every reduced irreducible polynomial

```text
Delta_0+Phi,                         Phi in (A,C)^5,         (14)
```

has the same four formal branches at `o` and is excluded by `(5)`. It need
not have one place at infinity; that global property is irrelevant to the
local address contradiction.

## 5. Why THM-3855 formal surjectivity cannot escape

THM-3855 proves two compatible formal statements:

1. every target perturbation `(14)` is obtained from `Delta_0` by a
   tangent-identity formal base automorphism; and
2. every formal binary-cubic coefficient row with linear part `(7)` is,
   after such a base automorphism and an identity `SL2` binary-variable
   gauge, the completed THM-3808 row.

If

```text
P=A mod (A,C)^2,                    Q=C mod (A,C)^2,         (15)
```

then

```text
Delta_0(P,Q)=P(Q+5P)(4Q+19P)(3Q-17P).                      (16)
```

The four factors in `(16)` are smooth formal branches with the four distinct
linear terms in `(9)`. A tangent-identity automorphism preserves their
branch count and tangent packet. Hence the strength of the formal lifting
theorem is now also its obstruction: all its lifts remain in the ordinary
four-branch completed orbit.

The formal result itself is not retracted. It proves that there is no local
coefficient-equation obstruction. What `(5)` proves is that **every
polynomial algebraization with this completed linear packet lies outside the
cubic Keller geometry**. Higher coefficient depth cannot pay a debt already
visible in the degree-four tangent cone.

## 6. Exact boundary and the next positive grammar

The obstruction is sharp at the address count. Three distinct rays meet the
cubic cap and are not excluded. If the nonzero quartic tangent cone of a
common-zero binary cubic has at most three distinct factors, its possible
root-multiplicity partitions are

```text
3+1,                 2+2,                 2+1+1,            (17)
```

or the single-ray partition `4`. THM-3887 proves that the last case forces a
rank-one linear coefficient pencil and cannot give a reduced irreducible
unibranch discriminant of total degree at most five. Thus the least-used
positive common-zero designs are precisely the three collision partitions
in `(17)`, together with the higher-order case in which the quartic jet
vanishes.

The theorem does not exclude:

- a reducible discriminant whose four rays lie on different global branch
  components;
- the equality case of three addresses;
- a degree at least four field extension;
- a unit-ideal nonmonogenic packet such as THM-3907, whose coefficients have
  no common zero; or
- the general planar Jacobian conjecture.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
python3 -O 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
