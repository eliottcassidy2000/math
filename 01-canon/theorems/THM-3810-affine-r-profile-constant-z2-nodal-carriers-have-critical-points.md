---
id: THM-3810
title: "Affine R-profile constant Z2 nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  c=1 cubic pseudo-plane, every carrier
  A=e^2-z/3+r(ae+b)+h z^2, for arbitrary a,b,h, has a critical point.
  The slope u=r/z again compresses criticality to a quartic and cubic.  If
  their residual resultant had support only on 1-6hz=0, its first three
  coefficients would force degree four, b=0, and a=0, after which a nonzero
  32z^7(1-6hz)^5 term remains.  The h=0 seam is recomputed directly.
  Hence no affine r-profile plus constant z^2 repair has a Darboux mate.
source: jc-cohn-boundary / affine mixed-cell slope-resultant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-23).  The
  Hamiltonian packet, Casimir syzygy, affine slope compression, compact
  residual resultant, pure-power boundary implication and all four decisive
  coefficient comparisons, separately recomputed h=0 seam, fixed
  quartic/cubic u-degrees, and denominator-free critical reconstruction were
  independently rederived.  The 26 active gates pass normally and under
  Python -O, both runs byte-match the frozen transcript, raw hashes match,
  and the agent-facing documentation check passes.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3803-affine-linear-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3809-smallest-mixed-r-z2-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_affine_mixed_carrier_thm3810.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_affine_mixed_carrier_thm3810.out
script_sha256: b955f43a83b035e2cbc9f43df8b3b52025844f251bb29ed2b4d778a862c1df2e
output_sha256: 1294a32311a7a69536a4275f33f1bf584dfd5f1132c0d80c1a8fc831311de7ef
semantic_sha256: 8914826aa78e2fa2ce140322d68d4fe752ebcbc3a66c0184aa34e4271bcdd03e
hash_basis: raw LF bytes
---

# THM-3810 -- every affine-profile mixed carrier remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero.  Retain the smooth
symplectic `c=1` cubic pseudo-plane

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),                         (1)
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.      (2)
```

For arbitrary `a,b,h in k`, the regular function

```text
A=e^2-z/3+r(ae+b)+h z^2                              (3)
```

has a critical point on `Y`.  Consequently it has no regular Darboux mate.

THM-3803 closes `(3)` when `h=0`, and THM-3809 closes it when `a=0`.
Neither boundary theorem controls the genuinely three-parameter interior.
Here the slope method survives the `ae` term because the first Hamiltonian
remains linear in `e`.

## 1. The affine slope equations

Put

```text
L=1-6hz,                  K=1+2re,                  g=ae+b.           (4)
```

Direct calculation gives

```text
C_r={A,r}=Lr^2-9z^2(2e+ar),
C_z={A,z}=3gr^2-3K(2e+ar),
C_e={A,e}=9gz^2-LK.                                  (5)
```

As always, the surface Casimir supplies

```text
K C_r-3z^2 C_z+r^2 C_e=0.                            (6)
```

A critical point cannot lie on `z=0`: there `(5)` first gives `r=0`,
whereas the last equation reads `-1=0`.  On `z!=0`, set

```text
u=r/z,                     M=Lu^2-9auz.              (7)
```

The equation `C_r=0` is exactly `e=M/18`.  After this substitution, surface
membership and `C_e=0`, with the latter multiplied by `18`, are

```text
p1=18z^2-u^2zM-18u=0,
p2=162bz^2+9aMz^2-18L-2LuzM=0.                      (8)
```

This remains a quartic/cubic pair in `u`; no parameter has been inverted.

## 2. The residual and its forbidden-root contradiction

Exact elimination gives

```text
Res_u(p1,p2)=-2916z^3(6hz-1)^3 T(z),                 (9)
```

where

```text
T=32L^5z^7+36L^4
  +324L^3a^2z^5-3888L^3abz^8+648L^3bz^2
  +729L^2a^4z^10-2916L^2a^2bz^7
  -13122La^3bz^6+104976La^2b^2z^9-52488Lb^3z^6
  -59049a^5bz^11+118098a^3b^2z^8-236196b^4z^8.      (10)
```

The residual has

```text
T(0)=36.                                              (11)
```

Suppose first that `h!=0`.  Its linear coefficient is

```text
[z]T=-864h,                                          (12)
```

so `T` is nonconstant.  If all its roots lay on the forbidden boundary
`L=0`, factorization over the algebraically closed field and `(11)` would
force

```text
T=36L^d,                       d=deg_z T>0.           (13)
```

The linear coefficient of the right side is `-216dh`; comparison with
`(12)` forces `d=4`.  Now compare successively with `36L^4`:

```text
[z^2](T-36L^4)=648b,
[z^5](T-36L^4)|_(b=0)=324a^2.                       (14)
```

Thus `(13)` would force `b=0` and then `a=0`.  But on that terminal stratum

```text
T-36L^4=32L^5z^7!=0,                                (15)
```

a contradiction.  Therefore `T` has a root away from `L=0`; by `(11)` that
root also avoids `z=0`.

The degree changes when `h=0`, so recompute this seam rather than
specializing the preceding degree argument.  Here `L=1`, `T(0)=36`, and

```text
b!=0:       [z^2]T=648b!=0,
b=0:        [z^7]T=32!=0.                            (16)
```

Thus `T` is again nonconstant and has a nonzero root.  Combining both rows,
there is always a residual root

```text
z=z_0,                         z_0L(z_0)!=0.          (17)
```

## 3. Denominator-free reconstruction

At `(17)` the leading coefficients of the two polynomials in `(8)` are

```text
LC_u(p1)=-Lz,                    LC_u(p2)=-2L^2z.     (18)
```

They are nonzero, so the quartic and cubic retain their actual degrees and
the vanishing resultant yields a finite common root `u_0`.  Also
`p1(0)=18z_0^2`, hence `u_0!=0`.  Define

```text
r_0=u_0z_0,                       e_0=M(u_0,z_0)/18.  (19)
```

Before imposing `p1=p2=0`, exact polynomial identities give

```text
r_0^2e_0-z_0^3+r_0=-z_0p1/18,
C_r(r_0,z_0,e_0)=0,
C_e(r_0,z_0,e_0)=p2/18,
C_z(r_0,z_0,e_0)=u_0^2p2/54.                         (20)
```

Thus `(19)` is a point of `Y` at which every Hamiltonian component of `A`
vanishes.  Since `(2)` is symplectic, it is an actual critical point.  No
regular `B` can satisfy `{A,B}=1` there.

This theorem closes all constant `z^2` corrections paired with an affine
`r`-profile.  The first mixed profile not covered is therefore

```text
A=e^2-z/3+r g(e)+z^2h(e)                             (21)
```

with `deg g>=2` or nonconstant `h(e)`.  In that range `C_r=0` is no longer
uniformly denominator-free: for quadratic `g`, it remains linear in `e`
but its leading coefficient is the moving factor `1+a_2r`, while genuine
nonlinearity begins at `deg g>=3`.  The slope still organizes the equations,
but a new sidecar is needed to control that denominator divisor and the
extra `e`-branches.  The exact companion named in the
metadata checks `(5)--(20)`, including the specialized `h=0` resultant and
all reconstruction identities.  Normal and optimized executions byte-match
the frozen 26-gate transcript.  **QED.**
