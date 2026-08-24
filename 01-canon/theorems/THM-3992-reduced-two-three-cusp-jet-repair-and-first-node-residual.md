---
id: THM-3992
title: "Reduced 2:3 cusp jets force a centered nodal boundary and mixed-generator residual"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In THM-3989's
  first reduced cusp-pole cell, determinant-one target shears and translations
  erase the tau^-4, tau^-3, and tau^-2 integration constants. The tau^-1 row
  leaves a constant residual I in a specified gauge slice, while the tau^0
  row and THM-3989's scalar moment give the exact repair equation
  2*h^2*W-I*q0=s. Polynomiality and normality force h=gamma*s. The deleted
  source line is then forced to the cubic A=X^2+a, C=X^3+(3a/2)X. The cusp
  a=0 contradicts the Keller condition, so a!=0 and the image has a node with
  two normalization addresses. Those addresses center X at gamma*x, force
  b in s^2*k[s], and yield the exact mixed-generator normal form
  C^2-A^3+(3a^2/4)A+a^3/4=gamma*u+(3a/(2gamma))*p+R(p,y)
  with R in (p^2,y). This sharply narrows but does not close the 2:3 cell.
source: root + frontier_transfer_scout / Hopf repair-quotient transfer, 2026-08-24
audit: >
  PASS (root, 2026-08-24). The proof was independently checked from the full
  Laurent convolution law. The companion script separately expands the
  Laurent Jacobian, verifies every row and factorization, and passes the cusp
  and node controls. Normal and optimized executions are byte-identical.
depends_on:
  - THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
script: 04-computation/jc2_reduced_23_cusp_jet_repair_thm3992.py
output: 05-knowledge/results/jc2_reduced_23_cusp_jet_repair_thm3992.out
script_sha256: 985b4833162ac0080b305fba5c2847a0736997940914033c3f2cda11cfd906dc
output_sha256: 87c92173c5fb9893b2bb1ee68165938f24ae375705d4677271793a853c1c6cf5
hash_basis: raw LF bytes
---

# THM-3992 -- reduced 2:3 cusp jets force a centered nodal boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero in the logarithmic
coordinates of THM-3989:

```text
s=xt,                         tau=t,
J_(x,t)(A,C)=tau*(A_s C_tau-A_tau C_s).                 (1)
```

Assume that a hypothetical Darboux pair has the first reduced pole depths
`(2,3)`. After a determinant-one constant target scaling, its leading rows
may be written

```text
A=h^2 tau^-2+b tau^-1+A0+A1 tau+A2 tau^2+A3 tau^3+...,
C=h^3 tau^-3+d tau^-2+e tau^-1+C0+C1 tau+C2 tau^2+...,

h,b,Ai,d,e,Ci in k[s],         h!=0,        s|h.         (2)
```

The divisibility `s|h` is inherited from THM-3989. The theorem proves that
the bracket equations force the reverse divisibility and hence

```text
h=gamma*s                 for some gamma in k*.          (3)
```

Before proving `(3)`, we identify the exact repair quotient through the
constant row.

## 1. The first three rows are gauge constants

For coefficient dictionaries `A=sum a_i tau^i`, `C=sum c_j tau^j`, put

```text
E_k=sum_(i+j=k)(j a_i' c_j-i a_i c_j').                 (4)
```

THM-3989 gives `E_k=delta_(k,0)`. At row `tau^-4`,

```text
E_-4=h^4*(2d/h^2-3b/h)'.                                (5)
```

Thus `lambda=2d/h^2-3b/h` is constant. The determinant-one target shear

```text
C -> C-(lambda/2)A                                      (6)
```

does not change the depths or leading coefficients and permits the relabeling

```text
d=(3/2)h b.                                              (7)
```

Put `B=b/h`, `E=e/h`. The next row is

```text
E_-3=h^3*(2E-(3/4)B^2-3A0)'.                            (8)
```

Its bracketed expression is another constant `mu`. The translation
`A->A+mu/3`, followed by relabeling `A0`, makes

```text
e=h*((3/8)B^2+(3/2)A0).                                 (9)
```

Define rational approximate-root coefficients in `k(s)` by

```text
q0=B/2,
q1=(A0-q0^2)/(2h),
q2=(A1-2q0 q1)/(2h),
C0*=q0^3+6h q0 q1+3h^2 q2.                             (10)
```

If `R=C0-C0*`, direct expansion gives

```text
E_-2=2h^2 R'.                                           (11)
```

Hence `R` is constant, and `C->C-R` kills it without changing any negative
row. Equations `(5)`, `(8)`, and `(11)` show that all apparent obstructions
through `tau^-2` lie in the image of allowed polynomial target repairs.

## 2. The first residual and the scalar repair equation

Let `q3 in k(s)` be arbitrary and put

```text
U=A2-(q1^2+2q0 q2+2h q3),
C1*=3q0^2 q1+3h q1^2+6h q0 q2+3h^2 q3,
V=C1-C1*,
I=3h^2 U-2h V.                                          (12)
```

The `q3` terms cancel from `I`, and the next row is exactly

```text
E_-1=-h I'.                                             (13)
```

Thus `I in k`. Choose

```text
q3=(A2-q1^2-2q0 q2)/(2h),                               (14)
```

so that `U=0` and `V=-I/(2h)`. Next put

```text
q4=(A3-2q0 q3-2q1 q2)/(2h),
C2*=3h^2q4+6hq0q3+6hq1q2+3q0^2q2+3q0q1^2,
W=C2-C2*.                                               (15)
```

Then the constant bracket row and the scalar moment of THM-3989 are

```text
E_0=(2h^2W-Iq0)',
M=sum_i i a_i c_-i=Iq0-2h^2W.                           (16)
```

Since `M=-s`, the integration constant is not free:

```text
2h^2W-Iq0=s.                                            (17)
```

Consequently `(17)` also gives `E_0=1`, as required. The first three
constants were gauge; `I` is the first surviving constant in this displayed
slice, and it is coupled to the next coefficient by the nonconstant moment.

## 3. Polynomiality forces the leading root to be linear

Put

```text
r=A0-q0^2=2h q1.                                        (18)
```

Eliminating `q1,q2,q3` from `(12)` gives the exact identity

```text
I=3h^2A2-2hC1+(3/4)r^2+(3/2)bA1.                        (19)
```

All terms except `r^2` visibly belong to `k[s]`, so `r^2 in k[s]`. Since
`r in k(s)` and the UFD `k[s]` is integrally closed, `r in k[s]`. Equation
`q0^2=A0-r` and the same argument give

```text
q0 in k[s],                    b=2h q0,                 (20)
```

so in particular `h|b`. This divisibility was not assumed.

The expression in `(15)` also simplifies to

```text
C2*=(3/2)*(hA3+q0A2+q1A1-q0q1^2).                      (21)
```

Multiplying `(17)` by four, substituting `(19)--(21)`, and using
`b=2hq0` gives

```text
8h^2 C2=4s+hK,                                          (22)

K=12h^2A3+24q0hA2+6rA1+12q0^2A1-8q0C1 in k[s].         (23)
```

The left side and `hK` are divisible by `h`; hence `h|4s`, and `4` is a unit.
Therefore `h|s`. Together with THM-3989's `s|h`, this proves `(3)`.

In particular, the leading approximate root is no longer merely a Laurent
fraction:

```text
h tau^-1=gamma*s/tau=gamma*x in B_2.                    (24)
```

This is a sharp collapse of the first reduced cell. It is not yet a depth
reduction: subtracting a function of the source variable `x` from `A` or `C`
is not a polynomial target automorphism and need not preserve the bracket.

## 4. The target cusp polynomial has one exact mixed residual

Put

```text
D=C^2-A^3+IA.                                           (25)
```

The coefficients already used through `(17)` determine every negative
Laurent row of `D`. Direct expansion gives

```text
[tau^k]D=0                         for k<=-2,
[tau^-1]D=h*(2h^2W-Iq0)=h*s=gamma*s^2.                 (26)
```

Since `u=s^2/tau`, the element `D-gamma*u` lies both in `B_2` and in
`k[s,tau]`. THM-3989's exact depth-zero intersection therefore supplies a
unique polynomial `F in k[p,y]` such that

```text
C^2-A^3+IA=gamma*u+F(p,y).                              (27)
```

This is a ring identity, not only a leading asymptotic. It is the direct
payoff of taking the Laurent defect modulo all earlier target repairs.

## 5. Every hypothetical pair has a centered nodal deleted-line image

Let

```text
qbar=q0(0),                 a=A0(0)-qbar^2,
X=gamma*x+qbar.                                           (28)
```

Evaluate the normalized Laurent expansions on the source line `t=0`.
Positive powers of `tau=t` vanish; using `(7)`, `(9)`, and `(10)` gives

```text
C0=q0^3+(3/2)q0*r+(3/2)h*A1.                            (28a)
```

Together with `(19)` and `h(0)=b(0)=0`, this yields

```text
A(x,0)=X^2+a,
C(x,0)=X^3+(3a/2)X,
I=3a^2/4.                                                (29)
```

Consequently

```text
C(x,0)^2=(A(x,0)-a)(A(x,0)+a/2)^2,
F(0,0)=-a^3/4.                                          (30)
```

If `a=0`, take the source point `t=0,X=0`. Equation `(29)` gives
`A_x=C_x=0` there, so `J_(x,t)(A,C)=0`, contradicting the Keller equation.
Therefore

```text
a!=0,                         I!=0.                     (31)
```

The deleted-line image is forced to be the node with target singular point
`(A,C)=(-a/2,0)` and two source addresses

```text
X^2=-3a/2.                                               (32)
```

Thus the one-place cusp branch is closed inside the reduced `(2,3)` cell.
The raw scalar `I` depends on the residual fifth-root normalization, but its
vanishing and the one-address-versus-two-address type do not.

The two node addresses also remove the translation `qbar`. Put

```text
G=C^2-A^3+IA+a^3/4.                                     (33)
```

Equation `(27)` and `(30)` show that `G=tQ` for some `Q in k[x,t]`. Since

```text
u_t(x,0)=x^2,          p_t(x,0)=1,          y_t(x,0)=0,
```

one has

```text
Q(x,0)=G_t(x,0)=gamma*x^2+F_p(0,0).                     (34)
```

At both addresses `(32)`, the target gradient of `G` vanishes, so the
pullback differential vanishes and both addresses are roots of `(34)`. The
even quadratic `(34)` has root sum zero, whereas the sum of the two roots in
`(32)` is `-2qbar/gamma`. Hence

```text
qbar=0,                    F_p(0,0)=3a/(2gamma).         (35)
```

It follows that `q0 in s*k[s]` and

```text
b=2gamma*s*q0 in s^2*k[s].                              (36)
```

Writing the constant and `p`-linear terms of `F` explicitly, `(27)` becomes
the final normal form

```text
C^2-A^3+(3a^2/4)A+a^3/4
 =gamma*u+(3a/(2gamma))*p+R(p,y),

R in (p^2,y) subset k[p,y].                             (37)
```

The boundary Bezout identity gives one further exact jet. At `(x,t)=(0,0)`,
equation `(29)` with `qbar=0` has `A_x=0`, `C_x=3a gamma/2`, and
`A_t=A1(0)`. Thus

```text
A1(0)=-2/(3gamma*a).                                    (38)
```

Finally, `(34)--(35)` give

```text
Q(x,0)=gamma*(x^2+3a/(2gamma^2)).                        (39)
```

Its two roots are simple. Hence the pullback of the target node has exactly
two transverse clutch points between the line `t=0` and the companion branch
`Q=0`. Whether those two local branches belong to one global component or two
is the new discrete conductor-incidence bottleneck.

The family `A=x^2+a`, `C=x^3+(3a/2)x` remains a useful hostile control for
`(29)--(30)`: it has the same cusp/node geometry and `I=3a^2/4`, but bracket
zero and failure of `(17)`. Thus the boundary normal form alone is not a
Keller construction.

## 6. The one-place cusp cotangent sidecar

For the cusp ring and its normalization

```text
R=k[P,Q]/(Q^2-P^3)  subset  S=k[z],
P=z^2,                         Q=z^3,                    (40)
```

the conductor and cotangent sequence are

```text
(R:S)=(P,Q)=z^2S,

0 -> (R/(P^2,Q))*theta -> Omega_(R/k) -> S dz
  -> k dz -> 0,

theta=3Q dP-2P dQ,              Ann_R(theta)=(P^2,Q).   (41)
```

Indeed, `dP` and `dQ` map to `2z dz` and `3z^2 dz`, so the image is `zS dz`
and the cokernel is `k dz`. The relation in `Omega_R` is
`2Q dQ-3P^2 dP=0`, and `theta` maps to zero. If
`u dP+v dQ` maps to zero, then `2u+3zv=0`. The condition `zv in R` forces
`v in (P,Q)`; splitting `v=P r0+Qq` reduces the kernel element to a multiple
of `theta` plus the defining relation. Finally, `r0 theta=0` exactly when
`r0/z in R`, equivalently `r0 in (P^2,Q)`.

This is the one-branch counterpart to THM-3955's node. The cusp has a genuine
conductor-torsion differential and a one-dimensional normalization cokernel;
the node instead has two normalization addresses. “One projective point” does
not choose between these local modules.

## 7. Scope and reproduction

The theorem does **not** prove:

1. that the scalar representative `I` is fixed by every target automorphism;
2. that `q1,q2,...` belong to `B_2` rather than only `k(s)`;
3. that the source term `gamma*x` can be removed by an allowed target repair;
4. that `R in (p^2,y)` vanishes or is forced into a previously closed cell;
5. that the two companion branches in `(39)` lie on one global component;
6. that the reduced `(2,3)` cell is empty, or that `JC(2)` holds.

The residual fifth-root normalization acts by
`h->zeta*h`, `A->zeta^2A`, `C->zeta^3C`, with `zeta^5=1`; the scalar `I`
has weight four. The robust invariant is not its chosen value but the forced
node type and its two labelled normalization addresses.

The next exact object is now the coordinated residual `R in (p^2,y)` together
with the two transverse node clutches `(39)`, not an amorphous `(2,3)` cell.
The cheapest decisive test is its **complete** oriented node-address graph and
integral class group.  The roots of `(39)` enumerate only intersections with
the known line `t=0`, not every source address over the target node.  A
two-owner split is not itself a boundary forest obstruction: THM-3951 concerns
completion-boundary primes, whereas these pullback curves lie inside the
source.  One must first determine whether further node addresses occur and
whether the target node lies in the nonproper-value locus.  See MISTAKE-469.
THM-3996 proves the resulting exact alternative: distinct companion owners
force an additional address or a nonproper target value.  More strongly, the
Keller degree-two exclusion says that the complete fibre has an additional
address whenever the node is proper; a two-edge cycle can only be one
connected subpacket.

Reproduce from the repository root:

```text
python3 04-computation/jc2_reduced_23_cusp_jet_repair_thm3992.py
python3 -O 04-computation/jc2_reduced_23_cusp_jet_repair_thm3992.py
```

Both executions must byte-match
[the frozen output](../../05-knowledge/results/jc2_reduced_23_cusp_jet_repair_thm3992.out)
and end with `ALL EXACT CHECKS PASSED`.
