---
id: THM-3792
title: "Pure first-normal nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  THM-3785 cubic pseudo-plane, every pure first-normal lift
  A_f=e^2+z f(e) of the quadratic coordinate of the minimal nodal arm
  immersion has a critical point whenever f(0) is nonzero.  In particular
  every polynomial Bezout completion required by THM-3790 fails at the
  carrier level, before solving for a mate.  The critical locus is produced
  uniformly by the nonconstant polynomial
  c^6(f-3ef')^3+864e^7.  Corrections in the square of the arm ideal remain
  open, so no general Darboux obstruction or JC(2) conclusion is claimed.
source: jc_zero_debt_lift / nodal carrier critical-locus elimination, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The three Hamiltonian
  components, restriction to K=0, compatibility in both directions,
  nonvanishing of H_f at a critical-polynomial root, actual point recovery,
  all degree seams, and the exact scope modulo I^2 were rederived.  The
  proof derives all three Hamiltonian
  equations, isolates the K=c^3+2re sheet, proves the critical polynomial is
  nonconstant in every degree (including the impossible f-3ef'=0 branch),
  and reconstructs r and z exactly at every root.  The deterministic
  companion checks those identities over Q(c), 257 degree controls, the
  forced Bezout constant, and the THM-3790 canonical seven-point boundary.
  Normal and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
script: 04-computation/jc2_cubic_pseudoplane_pure_nodal_carrier_thm3792.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_pure_nodal_carrier_thm3792.out
script_sha256: 1b5f077050e5be032d213981f7559ecb36e243bf6df5b1f0a612eba11ddc235e
output_sha256: ec3b6d4dba447e3bf25543a1f7bb42a0a01cd18bace0542defe6b22211c55cb3
semantic_sha256: 10a84bf8ce701c2806bd60a64cf978ea335ad966ff5e525f59963443d6105d3e
hash_basis: raw LF bytes
---

# THM-3792 -- a first-normal nodal correction always has a critical point

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero, fix `c in k*`, and
retain the smooth symplectic surface

```text
Y=Spec k[r,z,e]/(r^2e-z^3+c^3r)                       (1)
```

with Poisson packet

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re.  (2)
```

For every polynomial `f in k[e]` with `f(0)!=0`, put

```text
A_f=e^2+z f(e).                                       (3)
```

Then `A_f:Y->A1` has a critical point.  More explicitly, define

```text
H_f(e)=f(e)-3e f'(e),
G_f(e)=c^6 H_f(e)^3+864e^7.                           (4)
```

For every root `eta` of `G_f`, the values

```text
e_0=eta,
r_0=-c^3/(2eta),
z_0=6eta^2/H_f(eta)                                  (5)
```

are defined and give a critical point of `A_f` on `Y`.

The qualifier *pure first-normal* is exact.  If `I=(r,z)` is the arm ideal,
THM-3790 gives

```text
B/I^2=k[e,z]/(z^2).                                   (6)
```

Thus `(3)` chooses the arm value `e^2` and an arbitrary polynomial first
normal coefficient, but sets every correction in `I^2` to zero.  The theorem
does not treat carriers

```text
e^2+zf(e)+R,                         R in I^2.         (7)
```

Those second-normal corrections are the new live coordinate.

## 1. Why every nodal Bezout first jet is included

The minimal immersed noninjective arm profile of THM-3790 is

```text
(a_0(e),b_0(e))=(e^2,e^3-e).                         (8)
```

If `f` and `g` are the two first-normal coefficients of a hypothetical
Darboux lift, the exact conormal law is

```text
3c^3[f(e)(3e^2-1)-2e g(e)]=1.                        (9)
```

Setting `e=0` forces

```text
f(0)=-1/(3c^3)!=0.                                   (10)
```

Conversely, every polynomial `f` with the value in `(10)` admits a polynomial
`g`: the numerator `3c^3f(3e^2-1)-1` vanishes at zero and can be divided by
`6c^3e`.  Equivalently all solutions are obtained from the THM-3790
particular row by

```text
f=-1/(3c^3)+2e h(e),
g=-e/(2c^3)+(3e^2-1)h(e),             h in k[e].      (11)
```

Hence the theorem closes every all-polynomial choice of the first-normal
Bezout sidecar for the nodal coordinate, not merely the constant choice used
in THM-3790.

## 2. Exact Hamiltonian reduction

Directly from `(2)`, the three Hamiltonian components are

```text
{A_f,r}=-3f r^2-9(2e+zf')z^2,
{A_f,z}=-3(2e+zf')(c^3+2re),
{A_f,e}= 3f(c^3+2re).                                (12)
```

It is enough to look on the sheet

```text
K=c^3+2re=0.                                         (13)
```

Because `c!=0`, this sheet has `r,e!=0`.  Combining `(13)` with the surface
equation `(1)` gives

```text
r=-c^3/(2e),                    z^3=-c^6/(4e).        (14)
```

The last two equations in `(12)` now vanish.  Dividing the first by `-3`
and using `(14)` reduces the remaining equation to

```text
c^6 H_f(e)+24e^3z^2=0.                               (15)
```

For `e!=0`, equations `(14),(15)` are compatible if and only if

```text
c^6 H_f(e)^3+864e^7=0.                               (16)
```

Indeed, cubing the value of `z^2` in `(15)` and comparing it with the square
of `z^3` in `(14)` gives `(16)`.  Conversely, if `(16)` holds, then
`H_f(e)!=0` and the value `z=6e^2/H_f(e)` simultaneously satisfies both
`(14)` and `(15)`.  Thus `(4)--(5)` are not a resultant-only necessary
condition: they reconstruct an actual point where all three expressions in
`(12)` vanish.  Since the Poisson structure on the smooth surface is
symplectic, this is exactly a critical point of `A_f`.

## 3. The critical polynomial always has a nonzero root

First,

```text
G_f(0)=c^6f(0)^3!=0.                                 (17)
```

So every root of `G_f` is nonzero, as required in `(5)`.  It remains to show
that `G_f` is nonconstant.

If `m=deg f=0`, then `H_f=f` is a nonzero constant and the term `864e^7`
makes `G_f` degree seven.  If `m>=1`, then

```text
deg H_f=m,
lc(H_f)=(1-3m)lc(f)!=0.                              (18)
```

The two terms in `G_f` have degrees `3m` and `7`.  These degrees never agree
for an integer `m`, so their leading terms cannot cancel.  In particular the
formal exceptional equation `H_f=0` has no nonzero polynomial solution:
coefficientwise it would require `(1-3j)[e^j]f=0` for every integer `j>=0`.
Therefore `G_f` is nonconstant in every case.  Algebraic closedness supplies
a root `eta`, and `(17)` makes it nonzero.  Section 2 then supplies the
critical point `(5)`.

## 4. Sharp boundary and remaining repair coordinate

For the constant completion in THM-3790,

```text
f=-1/(3c^3),
G_f(e)=864e^7-1/(27c^3).                              (19)
```

Under `z=6e^2/f=-18c^3e^2`, equation `(19)` is exactly

```text
8z^7+9c^15=0,                                        (20)
```

recovering all seven critical points of that theorem.  The new content is
that no polynomial change of the Bezout parameter `h(e)` in `(11)` removes
the critical locus.

Thus the next nodal construction cannot merely improve the first-normal
completion.  It must use a genuine class in `I^2`, change the arm immersion,
or leave the minimal nodal bidegree.  This is a carrier obstruction, not a
proof that arbitrary nonlinear carriers or Darboux pairs on `Y` do not exist.

The deterministic companion named in the metadata verifies `(9)--(20)`, the
degree split through `m=256`, and the canonical seven-point specialization
with active assertion-free gates.  **QED.**
