---
id: THM-3809
title: "Smallest mixed R/Z2 nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  c=1 cubic pseudo-plane, every carrier A=e^2-z/3+b r+h z^2, for arbitrary
  b,h, has a critical point.  In the genuinely mixed cell b*h!=0, the slope
  u=r/z compresses criticality to two equations whose residual resultant is
  a degree-twelve polynomial R with R(0)=-9 and nonzero value on the only
  other forbidden boundary 1-6hz=0.  Specialized residuals close both
  parameter axes and the origin without generic-specialization arguments.
  Hence the smallest simultaneous r/z^2 repair cannot admit a Darboux mate.
source: jc-cohn-boundary / mixed-cell slope-resultant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-23).  The
  Hamiltonian packet, Casimir syzygy, slope compression, compact five-term
  resultant, all four parameter strata, forbidden-boundary values, fixed
  projective u-degrees, and denominator-free source/critical reconstruction
  were independently rederived.  The 34 active gates pass normally and
  under Python -O, both runs byte-match the frozen transcript, raw script
  and output hashes match, and the agent-facing documentation check passes.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_smallest_mixed_carrier_thm3809.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_smallest_mixed_carrier_thm3809.out
script_sha256: d9821d9c9f4d551723d9fd19fbc7c9d570477f39a69a0de649e3aded407a0a29
output_sha256: 79249d6c80b96fd9de4c8d5db63e999d454f1f1f8023b3447c9d4f250e5e4444
semantic_sha256: d0b8cf13c5c1e68f0b1a76e19c36b444983ad9b95b389ec23a924c79a5c631f7
hash_basis: raw LF bytes
---

# THM-3809 -- every smallest mixed r/z2 carrier remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero.  On the `c=1`
member of the THM-3785 cubic pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),                         (1)
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.      (2)
```

For arbitrary `b,h in k`, the regular function

```text
A=e^2-z/3+b r+h z^2                                  (3)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

This is a genuinely mixed statement.  THM-3795 closes the `b=0`
`r`-independent cell and THM-3799 closes the `h=0` monomial-`r` cell, but
criticality is nonlinear and those two theorems do not control their sum.
The proof below treats all of `(3)` uniformly and recomputes every
degree-drop axis.

## 1. Slope compression

Put

```text
L=1-6hz,                         K=1+2re.             (4)
```

Directly from `(2)`, the three Hamiltonian components are

```text
C_r={A,r}=Lr^2-18ez^2,
C_z={A,z}=3br^2-6eK,
C_e={A,e}=9bz^2-LK.                                  (5)
```

They satisfy the surface Casimir syzygy

```text
K C_r-3z^2 C_z+r^2 C_e=0.                            (6)
```

A critical point cannot have `z=0`: then `C_r=r^2=0`, hence `r=0`, but
`C_e=-1`.  On `z!=0`, introduce the slope

```text
u=r/z.                                                (7)
```

The equation `C_r=0` gives `e=Lu^2/18`.  After this substitution, surface
membership and `C_e=0` become exactly

```text
p1=18z^2-u^4Lz-18u=0,
p2=81bz^2-u^3L^2z-9L=0.                              (8)
```

Thus no division by a parameter or by `u` enters the compressed equations.

## 2. A five-term residual resultant

Exact elimination of `u` from `(8)` gives

```text
Res_u(p1,p2)=729z^3(6hz-1)^3 R(z),                   (9)
```

where the residual has the compact five-term form

```text
R(z)=59049b^4z^8+13122b^3z^6L-162bz^2L^3
     -8z^7L^5-9L^4.                                  (10)
```

The apparent simplicity of `(10)` is the reason the slope is the useful
coordinate: all mixed expansion terms are stored in powers of the single
first-normal factor `L`.

Suppose first that `bh!=0`.  Then

```text
deg_z R=12,                 LC_z(R)=62208h^5,
R(0)=-9,                    R(1/(6h))=9b^4/(256h^8). (11)
```

Hence `R` is nonconstant and every one of its roots avoids both forbidden
values `z=0` and `L=0`.

The parameter axes must not be inferred by specializing `(11)`, because its
degree drops there.  Recompute them instead:

```text
b=0, h!=0:   R=L^4 S,       S=48hz^8-8z^7-9,
             deg S=8,       S(0)=S(1/(6h))=-9;

h=0, b!=0:   deg R=8,       LC(R)=59049b^4,  R(0)=-9;

b=h=0:       R=-8z^7-9.                                  (12)
```

In each row algebraic closedness supplies a root of `R` (of `S` in the
first row) away from `zL=0`.  Therefore `(9)` vanishes at some

```text
z=z_0,                    z_0 L(z_0)!=0.              (13)
```

This is an existence proof over every allowed parameter pair, not a generic
parameter argument.

## 3. Every safe resultant root is an actual critical point

At `(13)`, the leading `u`-coefficients in `(8)` are

```text
LC_u(p1)=-Lz,                    LC_u(p2)=-L^2z.       (14)
```

Both are nonzero, so the quartic and cubic retain their degrees.  The
vanishing affine resultant therefore gives a finite common root `u_0`.
Moreover `p1(0)=18z_0^2!=0`, so `u_0!=0`.  Define

```text
r_0=u_0z_0,                    e_0=L(z_0)u_0^2/18.    (15)
```

No denominator occurs in this reconstruction.  Before setting `p1=p2=0`,
direct substitution gives the polynomial identities

```text
r_0^2e_0-z_0^3+r_0=-z_0p1/18,
C_r(r_0,z_0,e_0)=0,
C_e(r_0,z_0,e_0)=p2/9,
C_z(r_0,z_0,e_0)=u_0^2p2/27.                         (16)
```

Thus `(r_0,z_0,e_0)` lies on `Y` and all three Hamiltonian components
vanish.  Since the Poisson structure of `(2)` is symplectic on the smooth
surface, this is an actual critical point of `A`.  At such a point
`{A,B}=0` for every regular `B`, so `{A,B}=1` is impossible.

The theorem closes the first simultaneous `r`/higher-normal cell while
leaving a sharply defined next boundary: carriers

```text
A=e^2-z/3+r g(e)+z^2 h(e)                            (17)
```

with nonconstant `g` or `h`, higher `z`-normal powers, different arm
profiles, and rational mates with poles are not covered.  The exact
companion named in the metadata checks `(5)--(16)`, including independently
specialized resultants on both axes and at the origin.  Normal and optimized
executions byte-match the frozen 34-gate transcript.  **QED.**
