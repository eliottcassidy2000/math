---
id: THM-3812
title: "Nodal arm-coefficient second-normal profiles cannot be Darboux pairs"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  On the
  cubic pseudo-plane, no Darboux pair with the normalized nodal arm jet can
  have all of its higher canonical coefficients confined to
  r g(e)+z^2 f(e) in each output.  Monic reduction in z makes the pure-z
  bracket coefficient a Bezout law for f,k, while the pure-r^3 coefficient
  forces their Wronskian to vanish.  Proportionality then makes a
  nonconstant quadratic-linear factor a polynomial unit, a contradiction.
  This is an exact canonical-profile theorem, not a claim about every
  element of the arm ideal I^2.  The first omitted mixed coefficient is rz.
source: jc_zero_debt_lift / nodal second-normal direct-slice lane, 2026-08-23
audit: >
  PROOF CANDIDATE.  The exact companion has 26 active gates checking the
  nodal first-normal Bezout row, Poisson/Casimir signs, monic canonical
  reduction, the pure-z and pure-r^3 coefficients, the Wronskian
  logarithmic derivative, all degree-at-most-eight hostile controls, and
  the precise arm-ideal sidecar.  Normal and optimized runs byte-match the
  frozen transcript.  Independent hostile audit remains required.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
  - THM-3810-affine-r-profile-constant-z2-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_nodal_second_normal_profile_thm3812.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_nodal_second_normal_profile_thm3812.out
script_sha256: 4cd8f6570a777c647b244a7ac1987b8e787defd99bf5816ec5c64529e9f8369d
output_sha256: 5d2b4df1a06a578570b7a02c28a3ffbf2739acb85d01fb414c1f001b27f48d2e
semantic_sha256: fa9ebcecfa1bdfc74b02412ae611a38810b261949c27c010ae742aded1f2e1d2
hash_basis: raw LF bytes
---

# THM-3812 -- the nodal pair cannot stop at arm-coefficient r/z2 profiles

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.**  Let `k`
be a field of characteristic zero, fix `c in k*`, and put

```text
B=k[r,z,e]/(r^2e-z^3+c^3r),                         (1)
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re. (2)
```

For arbitrary polynomials `f,g,h,kappa in k[e]`, define

```text
A=e^2-z/(3c^3)+r g(e)+z^2 f(e),                     (3)
C=e^3-e-ez/(2c^3)+r h(e)+z^2 kappa(e).              (4)
```

Then

```text
{A,C}!=1.                                            (5)
```

Thus a Darboux lift of the normalized nodal arm profile cannot stop with
arm-dependent coefficients in only the canonical `r` and `z^2` slots.

## 1. Exact scope and canonical normal form

Let `I=(r,z)` be the arm ideal.  The restrictions and first-normal
coefficients of `(3)--(4)` are

```text
(A|L,C|L)=(e^2,e^3-e),
(a_1,b_1)=(-1/(3c^3),-e/(2c^3)).                    (6)
```

They satisfy the exact conormal law

```text
3c^3[a_1(e)(e^3-e)'- (e^2)'b_1(e)]=1.              (7)
```

The added terms do not alter `(6)`: `z^2 in I^2`, while the surface
relation gives successively

```text
c^3r=z^3-r^2e,
r in I^2,
r^2e in I^3,                 hence r in I^3.         (8)
```

The theorem's wording is deliberately narrower than “all corrections in
`I^2`.”  It quantifies exactly the profiles in `(3)--(4)`, whose `r` and
`z^2` coefficients are pulled back from the arm `k[e]`; coefficients
depending on `r` or `z` and higher canonical monomials are not included.

There is no ambiguity or quotient loss in the coefficient calculation.
Indeed

```text
B=k[r,e][z]/(z^3-r^2e-c^3r)                         (9)
```

is free over `k[r,e]` with basis `(1,z,z^2)`, because its defining
polynomial is monic in `z`.  Every bracket therefore has a unique normal
form of `z`-degree at most two.  A bracket is the scalar one exactly when
all nonscalar coefficients in this normal form vanish.

## 2. Two necessary buckets

Compute `{A,C}-1` from `(2)` and reduce every `z^3` using `(9)`.  The
coefficient of the pure monomial `z` is

```text
[r^0z^1]({A,C}-1)
 = [36c^6e^2f-24c^6e kappa-12c^6f+1]/(2c^3).        (10)
```

Consequently `(5)` could fail only if

```text
(3e^2-1)f-2e kappa=-1/(12c^6).                      (11)
```

Independently, the coefficient of the pure monomial `r^3` is

```text
[r^3z^0]({A,C}-1)=12e^2(f kappa'-kappa f').          (12)
```

The profiles `g,h` do not occur in either `(10)` or `(12)`.  Since `k[e]`
is a domain, vanishing of `(12)` forces

```text
f kappa'-kappa f'=0.                                 (13)
```

Equations `(11),(13)` are necessary polynomial identities, obtained from
two distinct coefficients in the unique normal form.  No saturation,
division by `e`, or truncation of the other bracket equations has occurred.

## 3. Bezout and Wronskian are incompatible

If `kappa=0`, equation `(11)` says that the nonconstant polynomial
`3e^2-1` times `f` is a nonzero scalar, impossible in `k[e]`.

Suppose `kappa!=0`.  In the rational function field `k(e)`, equation `(13)`
is

```text
(f/kappa)'=0.                                        (14)
```

Characteristic zero makes the constant field of `k(e)` equal to `k`, so
there is a `lambda in k` with

```text
f=lambda kappa.                                      (15)
```

Substitution in `(11)` gives

```text
kappa(e)[lambda(3e^2-1)-2e]=-1/(12c^6).             (16)
```

A product in `k[e]` is a nonzero scalar only when both factors are units.
But the coefficient of `e` in the second factor of `(16)` is always `-2`;
it is nonconstant for every `lambda`.  This contradiction proves `(5)`.

## 4. The first positive design boundary

Equation `(8)` also types the next canonical layer accurately:

```text
r in I^3,                         rz in I^4.          (17)
```

Accordingly the theorem does **not** call `rz` an `I^3` term.  It says that
at least one output of a live nodal lift must leave the exact profile
`(3)--(4)`.  In the canonical `r,z` monomial grammar, the cheapest omitted
mixed coefficient is

```text
rz p(e),                                               (18)
```

in one or both outputs; a construction may instead jump to a still higher
coefficient.  This is a design gate, not an existence claim for `(18)`.
It is orthogonal to THM-3787's Euler-support census: the present failure is
the collision between a second-normal Bezout equation and a higher normal
Wronskian, not merely a count of active weights.

The exact companion named in the metadata checks `(7)--(13)`, the monic
normal reduction, all Poisson signs, and the proportional-profile factor.
As hostile controls only, it independently computes unit Groebner ideals
for all pairs `f,kappa` of degree at most eight after the good specialization
`c=2`; the all-degree proof is `(14)--(16)`, not that finite scan.  Normal
and optimized runs byte-match the frozen 26-gate transcript.  **QED.**
