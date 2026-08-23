---
id: THM-3800
title: "Sharp torus-escaping nodal carrier has fourteen critical points"
status: >
  PROVED + VERIFIED-EXACT + CORRECTION OF RESERVED DIRECTION + INDEPENDENTLY
  HOSTILE-AUDITED.  The THM-3795 sharp carrier
  A=e^2-z/(3c^3)+(4/c^6)z^3e^3 does escape the K=c^3+2re=0 torus, but it is
  not smooth: for every c!=0 it has exactly fourteen reduced critical
  points on two explicit K!=0 branches.  Therefore it has no regular
  Darboux mate, and the originally reserved four-weight mate-support
  enumeration is superseded before its first step.  No claim about rational
  mates with poles is made.
source: root reserved mate-support lane / jc_zero_debt_lift carrier validity audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The K=0 exclusion,
  K!=0 derivative reduction, seven-root law, both explicit r-branches,
  nonzero K values, set-theoretic exhaustiveness, general-c triangular
  Groebner basis, squarefreeness, discriminant, transverse minor, and exact
  reduced length fourteen were independently rederived.  The c=1 replay
  agrees; normal and optimized runs byte-match the frozen transcript, all
  hashes match, and documentation checks pass.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
related:
  - THM-3796-cubic-pseudoplane-first-two-by-four-support-peel
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_sharp_nodal_carrier_critical_points_thm3800.py
output: 05-knowledge/results/jc2_sharp_nodal_carrier_critical_points_thm3800.out
script_sha256: 401c48db675a151fd0438c03906dfab05385186e57b24632071471ca8aaa9216
output_sha256: 9c3537b698a7d97fd6fef2e9109e661bbe44f75bdd48e22acdd6b9526e041c5e
semantic_sha256: b9b81fd0e16042bf6bab0b44fa00ca91d5591032af857bffe97385a8360ea679
hash_basis: raw LF bytes
---

# THM-3800 -- the sharp torus escape is still a critical carrier

**PROVED + VERIFIED-EXACT + CORRECTION OF RESERVED DIRECTION + INDEPENDENTLY
HOSTILE-AUDITED.**  Work over an algebraically closed field `k` of
characteristic zero, fix `c in k*`, and retain the smooth symplectic surface

```text
Y=Spec k[r,z,e]/(r^2e-z^3+c^3r)                       (1)
```

with Poisson packet

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re.  (2)
```

The sharp boundary left by THM-3795 is

```text
A_sharp=e^2-z/(3c^3)+(4/c^6)z^3e^3.                  (3)
```

Its critical scheme on `Y` consists of exactly fourteen reduced points:
choose any of the seven roots `eta` of

```text
1296c^3 eta^7=1,                                      (4)
```

and either sign in

```text
e=eta,
z=-6c^3eta^2,
r=c^3(-1 +/- 1/sqrt(3))/(2eta).                       (5)
```

Consequently `A_sharp` has no regular Darboux mate.

## 1. Correction lineage: the missed validity gate

THM-3800 was first reserved for a four-weight enumeration of possible mates
of `(3)`.  The motivating observation was correct but insufficient:
`A_sharp` removes the THM-3795 critical points on the torus

```text
K=c^3+2re=0.                                          (6)
```

The reserved direction implicitly treated escape from that one witness as
evidence that the carrier might be smooth.  Testing carrier criticality is
logically prior to enumerating mates; the explicit points `(4)--(5)` refute
that premise.  The reserved stub had no proved force, so no proved theorem
is retracted.  Its proposed support enumeration is simply superseded by the
stronger carrier-level obstruction here.

The distinction between regular and rational mates is essential.  A
critical point rules out `{A_sharp,C}=1` for regular `C`, because every
Hamiltonian bracket vanishes there.  It does not by itself rule out a
rational `C` having poles at the critical points; no such claim is made.

## 2. Exhaust the escaped torus first

Since `(3)` is independent of `r`, write

```text
A_z=-1/(3c^3)+(12/c^6)z^2e^3,
A_e=2e+(12/c^6)z^3e^2.                               (7)
```

The Hamiltonian components are

```text
{A,r}=-3r^2A_z-9z^2A_e,
{A,z}=-3K A_e,
{A,e}= 3K A_z.                                       (8)
```

On `K=0`, the torus parametrization from THM-3795 is

```text
e=-c^6/(4z^3),               r=2z^3/c^3,             (9)
```

and `(3)` restricts exactly to

```text
A_sharp|_(K=0)=-z/(3c^3).                             (10)
```

Thus its torus derivative is the nonzero constant `-1/(3c^3)`.  Equivalently
the remaining component in `(8)` is

```text
{A,r}|_(K=0)=r^2/c^3!=0.                              (11)
```

There are no critical points on the sheet it was designed to escape.

## 3. The two off-torus branches

At every critical point one therefore has `K!=0`.  The last two equations
in `(8)` force

```text
A_z=A_e=0.                                            (12)
```

The first equation in `(7)` excludes `e=0`; solving `(12)` gives

```text
z^2e^3=c^3/36,                 z^3e=-c^6/6.           (13)
```

Their quotient and substitution yield exactly

```text
z=-6c^3e^2,                    1296c^3e^7=1.          (14)
```

For each of the seven distinct roots `eta` in `(4)`, the surface equation
is a quadratic in `r`.  Its two roots are precisely `(5)`.  At them

```text
K=c^3+2re=+/-c^3/sqrt(3),                             (15)
```

so neither branch falls back into the excluded case.  Conversely `(7)`
vanishes at every point `(4)--(5)`, hence all three expressions in `(8)`
vanish.  This proves both existence and set-theoretic exhaustiveness.

## 4. Exact count and reducedness

Over the rational function field `Q(c)`, the lexicographic Groebner basis
of the surface equation and the three Hamiltonian components is

```text
r^2+1296c^6e^6r+216c^9e^5,
z+6c^3e^2,
e^7-1/(1296c^3).                                     (16)
```

All leading coefficients remain units for `c!=0`.  The last polynomial is
squarefree and has seven nonzero roots.  Modulo it, the discriminant of the
first polynomial is

```text
432c^9e^5!=0.                                         (17)
```

Thus the two `r`-roots over every `e`-root are distinct.  The `z`-equation
is linear.  Equivalently the determinant of the `(z,e)` derivatives of
`(A_z,A_e)` reduces to `-1008e^5/c^3`, also nonzero at all seven roots.
Hence the critical scheme is reduced and has length `7*2=14`.

As an independent hostile specialization, at `c=1` the basis `(16)` becomes

```text
r^2+1296e^6r+216e^5,
z+6e^2,
1296e^7-1,                                            (18)
```

where the last row has merely been cleared of its unit denominator.  This
was the cheapest computation that exposed the failed mate-support premise.

The deterministic companion named in the metadata checks `(7)--(18)`, the
two explicit branches modulo `sqrt(3)^2=3`, both general-`c` and `c=1`
Groebner bases, the discriminant, transverse minor, and correction scope
with 36 active assertion-free gates.  Normal and optimized executions
byte-match the frozen transcript.  **QED.**
