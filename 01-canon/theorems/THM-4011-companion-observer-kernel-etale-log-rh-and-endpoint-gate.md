---
id: THM-4011
title: "Companion observer kernel, etale log-Riemann--Hurwitz, and endpoint gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the THM-3999
  height-two completion put
  Delta=p^3-y^2=t*p^2. For every normalized companion G=tQ and every
  T in p*k[p,y], multiplication by 1+T*Delta stays in the same mixed-residual
  normal-form class while preserving the boundary endpoint, boundary order,
  total strict class, two visible clutches, and arbitrarily many prescribed
  source-normal rows. Thus those observers alone cannot prove Q irreducible.
  The prime factors H_M=1+p^M*Delta have class zero, avoid both D and L, and
  have projective genus floor((M+3)/2), with two punctures for odd M and three
  for even M. This insertion is ambient and need not preserve Darboux
  realizability. If H_M actually divides the companion of an etale pair,
  log Riemann--Hurwitz excludes every odd M; even M forces covering degree
  at least M+3 and two distinct finite normalization-side escape values.
  Finally the extended boundary D maps into the nonproperness locus, so every
  THM-3999 endpoint is a nonproper value on the forced nodal cubic. A boundary
  endpoint over the target node is necessarily repeated and makes the strict
  companion singular there; repeated endpoints alone do not identify the
  node.
source: root + companion_geometry / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (companion_geometry primary exact audit + next_hasse_row independent
  SymPy-free hostile audit, 2026-08-24). The B2 repair
  identities, factor monoid, exact changed-row orders, endpoint/clutch/class
  blindness, prime-curve smoothness, branch and puncture parity, odd/even
  log-Riemann--Hurwitz invoices, live b=d=0 control, and node/smooth-tangency
  hostiles were reconstructed. The valuative endpoint lemma and actual-etale
  lift are proved geometrically below. The independent audit also caught and
  repaired an initially boundary-ramified converse hostile: the final
  repeated-endpoint control is an etale identity map with a smooth tangent
  pullback. Normal and optimized certificate streams are byte-identical.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3985-cusp-plane-rational-time-residue-and-height-two-mixed-submersion
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-3999-companion-divisor-boundary-endpoint-and-class-ledger
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
related:
  - THM-3994-double-resultant-collision-separates-two-address-and-length-two-seams
  - THM-4005-reduced-two-three-live-seam-invariant-support-atlas
script: 04-computation/jc2_companion_observer_kernel_thm4011.py
output: 05-knowledge/results/jc2_companion_observer_kernel_thm4011.out
script_sha256: cdc634f98ba0dc4aa6e705d1a56192bd553e71f2cfa4c173f5690a03f029048b
output_sha256: 44f3b587a56fbb8793fa8dc1e77d2167a42884403c32b0ea14896ff4e98af10a
independent_audit_script: 04-computation/jc2_companion_observer_kernel_thm4011_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_companion_observer_kernel_thm4011_independent_audit.out
independent_audit_script_sha256: 6d0c176f90daad2a0ea9ad86b3082bd196af3a47ecefb162eb052205e2adfa52
independent_audit_output_sha256: 41e4ab052f371280c95f040efe6f57e43a77205f748a3187c8d0d6acce512427
independent_audit_semantic_sha256: 66cd205bdb868b546f6e187d8fbda46377b8f070860fc8d28d628f6a03a2be49
hash_basis: raw LF bytes
---

# THM-4011 -- the companion observer has an exact class-zero kernel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. On the height-two
completion of THM-3973 put

```text
z=1+x^2t,                 p=zt,                 y=xtp,
u=x^2t=z-1,               Delta=p^3-y^2=tp^2,
B_2=k[x,z,p,y],           X_2=Spec(B_2),
U=A2_(x,t),               D=X_2 minus U,         L=closure(V(t)).       (1)
```

The affine-modification identities of THM-3985 include

```text
x=yp/Delta,        u=y^2/Delta,        t=Delta/p^2,
u*Delta=y^2.                                                     (2)
```

Begin with any algebraic normal-form datum of the type isolated in THM-3992
and THM-3999:

```text
G=gamma*u+alpha*p+R(p,y)=tQ,
alpha=3a/(2gamma),       a*gamma!=0,
R in (p^2,y).                                                     (3)
```

When `(3)` comes from a Darboux pair, `G` is the pullback of the fixed nodal
target polynomial

```text
P(A,C)=C^2-A^3+(3a^2/4)A+a^3/4.                         (4)
```

Section 1 first works in the larger space of algebraic data `(3)`. It does
**not** assert that its multiplication operation sends Darboux pairs to
Darboux pairs. Sections 2 and 3 then state what the actual etale-map
hypothesis recovers.

## 1. An arbitrarily deep class-zero factor is invisible to the companion observer

### 1.1 The exact factor-insertion monoid

For any

```text
T in p*k[p,y]
```

define

```text
H_T=1+T*Delta,
R_T=R+T*(gamma*y^2+(alpha*p+R)*Delta).                  (5)
```

Using `(2)` once gives the ring identity

```text
(gamma*u+alpha*p+R)H_T
 =gamma*u+alpha*p+R_T.                                  (6)
```

The right side is again in normal form, and

```text
R_T in (p^2,y),                R_T(0,y)=R(0,y).          (7)
```

Indeed, the only apparent `u`-term created by multiplication is
`gamma*T*u*Delta=gamma*T*y^2`. The assumption `p|T` makes every added term
vanish under `p=0`. Dividing `(6)` by `t` proves

```text
Q_T=Q*H_T.                                               (8)
```

These factors form a monoid because

```text
(1+T_1 Delta)(1+T_2 Delta)
 =1+(T_1+T_2+T_1T_2 Delta)Delta,                        (9)
```

and the new coefficient still lies in `p*k[p,y]`.

The two boundary colors both miss this operation. On `D`, one has `p=0`,
while on `L` one has `p=y=Delta=0`. Hence

```text
H_T|D=H_T|L=1.                                          (10)
```

Consequently multiplication by `H_T` preserves, scheme-theoretically,

```text
ord_D(Q),
Spec k[y]/(gamma-R(0,y)),
Q(x,0)=gamma*x^2+alpha,
the two labelled clutch roots on L.                     (11)
```

It also preserves THM-3999's total class. Every nonconstant `H_T` is a
regular nonunit of `B_2`, has no boundary component by `(10)`, and satisfies

```text
[div(H_T)]=0 in Cl(X_2).
```

Thus `(8)` merely appends a class-zero strict divisor to the total class
`-2[D]`.

### 1.2 No finite source-normal packet detects the factors `H_M`

For an integer `M>=1`, specialize to

```text
H_M=1+p^M Delta.                                        (12)
```

At fixed `x` near `t=0`,

```text
ord_t(p)=1,        ord_t(Delta)=3,       ord_t(G)=1.
```

The first changed rows are exact:

```text
G_M-G=t^(M+4)(gamma*x^2+alpha)+O(t^(M+5)),
Q_M-Q=t^(M+3)(gamma*x^2+alpha)+O(t^(M+4)).               (13)
```

Since `gamma*x^2+alpha` is nonzero, these are equalities, not lower bounds.
Choosing `M` arbitrarily large preserves any prescribed finite collection of
source-normal or Hasse rows. The residual fifth-root action also fixes
`p,y,Delta,H_M`, so this blindness survives the THM-3997 invariant quotient.

Combining `(10)--(13)`, the observer

```text
Q |-> (ord_D Q, endpoint scheme, total strict class,
       clutch polynomial, any fixed finite normal-row packet)             (14)
```

cannot determine irreducibility, connectedness, or component ownership on
the ambient normal-form space.

### 1.3 The inserted divisor is one prime curve of explicitly growing genus

Let

```text
C_M=V(H_M) subset X_2.                                  (15)
```

Modulo `H_M`, both `p` and `Delta` are units because `p^M Delta=-1`.
Equations `(1)--(2)` therefore identify the whole curve, not merely a dense
open, with

```text
C_M=Spec k[p,p^-1,y]/(y^2-p^3-p^-M).                    (16)
```

The rational function

```text
p^3+p^-M=p^-M(p^(M+3)+1)                               (17)
```

is not a square in `k(p)`: every root of `p^(M+3)+1` has valuation one.
Thus `H_M` is prime in `B_2`. The affine curve is smooth. If its `y`-partial
vanished, then `y=0`; equations `(16)` and its `p`-partial would respectively
give `p^(M+3)=-1` and `(M+3)p^(M+2)=0`, an impossibility.

The double cover of `P1_p` has branch places

```text
M+3 simple roots of p^(M+3)+1,
p=0 exactly when M is odd,
p=infinity for every M.                                 (18)
```

Riemann--Hurwitz gives the smooth projective genus

```text
g_M = (M+2)/2   if M is even,
g_M = (M+3)/2   if M is odd
    = floor((M+3)/2).                                   (19)
```

The affine curve `(16)` removes only the places over `p=0` and `p=infinity`.
For even `M`, the even pole at zero gives two places and the odd pole at
infinity gives one, so there are three punctures. For odd `M`, both poles are
odd and give one place each, so there are two punctures:

```text
r_M=3 if M is even,             r_M=2 if M is odd.       (20)
```

In particular, `C_M` is a smooth prime class-zero curve, is disjoint from
both `D` and `L`, and has genus at least two. The completion class and the
two clutch points are blind even to an infinite family of pairwise distinct
prime curves of growing genus.

### 1.4 Sharp ambient controls, including the current live coefficients

At `R=0`, THM-3999 proves that

```text
Q_0=gamma*x^2+alpha*z
```

is irreducible. For every `M`, equations `(5)--(8)` give a reducible datum

```text
R_M=p^M*(gamma*y^2+alpha*p*Delta),
Q_M=Q_0 H_M                                                (21)
```

with the same endpoint, class, clutch polynomial, and any chosen finite row
packet once `M` is large. The two factors are distinct: on `C_M`,
`t=-p^(-M-2)` and

```text
G_0=-gamma*(p^(M+3)+1)+alpha*p
```

is not the zero function. Thus `(21)` is a reduced factorization hostile,
not merely an inserted multiplicity.

There is also a control retaining every residual coefficient currently fixed
on the live minimal-support stratum. Set

```text
a=1,       gamma=-1/2,       alpha=-3,       b=d=0,
R_*= (8/3)p^2-(1376/135)p^3+(5696/105)p^4.                (22)
```

The first two coefficients in `(22)` are THM-3997's forced unnormalized
`p^2,p^3` rows. THM-4007 further forces, at `A5=1`,

```text
[p^4](R/gamma)+(6/7)[y^2](R/gamma)=-11392/105.
```

The choice `[y^2]R=0` in `(22)` therefore fixes the raw coefficient
`[p^4]R=5696/105`. With `M=1`, formula `(5)` gives

```text
R_host=R_*+p*(-y^2/2+(-3p+R_*)*(p^3-y^2)),
Q_host=Q_* *(1+p(p^3-y^2)).                              (23)
```

It remains boundary-disjoint and leaves

```text
([p^4],[p^2y],[y^2])R_host=(5696/105,0,0),               (24)
```

or, after division by `gamma`, `(-11392/105,0,0)`. Thus it satisfies the
THM-4007 residual comparison. The factor insertion first changes `Q` at
`t^4` and `G` at `t^5`; in particular it leaves the `G`-row at `t^4` forced
by the third source-normal `A,C` row. It does not supply those `A,C` rows.

THM-4008 now excludes the **entire** pure-`p` base datum `(22)` from an
actual Keller pair. Equations `(21)--(24)` are ambient observer hostiles,
not constructions or formal lifts. This firewall is load-bearing.

## 2. Actual etaleness recovers an odd-factor no-go and an even-factor invoice

Now suppose that `(A,C)` is an actual hypothetical Darboux pair in the
normalized reduced `(2,3)` cell, and assume

```text
H_M divides Q in B_2.                                    (25)
```

This is a new hypothesis on an actual pair, not the ambient operation of
Section 1. Since `H_M` is prime and disjoint from `L`, its curve `C_M` is an
irreducible component of the pullback of the nodal cubic `N_a=V(P)`. The
ambient map

```text
F=(A,C):A2 -> A2
```

is etale and quasi-finite, so it cannot contract `C_M`. The dominant map
`C_M->N_a` lifts uniquely to the normalization

```text
phi_M:C_M -> Ntilde_a=A1_X.                              (26)
```

This lift is unramified at every affine point. Away from the node this is
base change of the etale surface map. At a point over the node, the completed
etale local map identifies each selected source branch with one normalized
target branch. Equivalently, `(26)` lies in the etale base change
`A2 times_(A2) Ntilde_a -> Ntilde_a`.

Let `Cbar_M` be the smooth projective completion and extend `(26)` to a map

```text
phibar_M:Cbar_M -> P1_X
```

of degree `e`. All ramification lies among the `r_M` punctures in `(20)`.
Riemann--Hurwitz reads

```text
sum_(P puncture)(e_P-1)=2g_M-2+2e.                      (27)
```

### 2.1 Odd `M` is impossible

For odd `M`, there are only two punctures. Their total contribution is at
most `2e-2`, while the right side of `(27)` is

```text
2e+2g_M-2 > 2e-2
```

because `g_M>=2`. Therefore

```text
M odd  ==>  H_M does not divide an actual companion Q.  (28)
```

This is an actual Keller-cell exclusion, unlike Section 1's observer hostile.

### 2.2 Even `M` forces high degree and two finite exits

Let `M` be even. There are three punctures. At least one maps to target
infinity. If two or three did, then the fibre-degree identity over infinity
and the remaining maximum ramification would give respectively at most
`2e-3` or `e-3`, both smaller than `(27)`. Hence exactly one puncture maps to
infinity.

The other two punctures map to finite values of `A1_X`. They cannot map to
the same value: their ramification indices would sum to at most `e`, and the
total contribution in `(27)` would again be at most `2e-3`. Thus they give
two distinct finite **normalization-side** values. Finally the three-point
bound in `(27)` gives

```text
2g_M-2+2e <= 3(e-1),
e >= 2g_M+1=M+3.                                        (29)
```

Equality holds exactly when all three punctures are totally ramified of
index `e=M+3` over their three distinct normalization values.

Each finite puncture is a source escape with a finite target limit, so its
image in `N_a` lies in the nonproperness locus `S_F`. The two normalization
values can fold to one target point only if they are the two preimages of the
node; otherwise they give two distinct target values. The finite-envelope
degree ledger also gives

```text
d=[k(x,t):k(A,C)] >= e >= M+3.                           (30)
```

Thus even ghost factors are not excluded, but they must pay a growing field-
degree invoice and two labelled finite escapes. No factor `H_M` is asserted
to occur.

## 3. The endpoint does certify boundary-visible nonproperness

For an actual pair, both `A,C` lie in `B_2`, so the source map extends to

```text
Phi=(A,C):X_2 -> A2.                                    (31)
```

Let `d` be any closed point of `D`. Choose a DVR arc in the smooth surface
`X_2` transverse to `D` at `d`; its generic point lies in `U=A2`. If
`Phi(d)` had a neighborhood over which `F` were proper, the valuative
criterion would extend that generic source arc to a point of `U`. Composing
with the open immersion `U->X_2` would give a second extension of the same
generic arc to the separated surface `X_2`, one centered in `U` and one at
`d`. This is impossible. Therefore

```text
Phi(D) subset S_F.                                       (32)
```

Write the boundary restriction as

```text
phi_D(y)=(A_D(y),C_D(y)).                                (33)
```

Near `D`, THM-3999's strict companion equation is

```text
q=Q/x^2=G/u,                    u|D=-1,
q|D=E(y)=gamma-R(0,y).                                  (34)
```

Thus

```text
P(phi_D(y))=R(0,y)-gamma=-E(y),                         (35)
```

and the endpoint scheme is the literal fibre product

```text
D times_(A2) N_a = Spec k[y]/(E(y)).                    (36)
```

Combining `(32)` and `(36)` proves the new geometric content:

```text
every strict-companion endpoint maps into N_a intersect S_F.              (37)
```

If `R(0,y)` is nonzero, then `(35)` makes `phi_D` nonconstant. Its image
closure is an irreducible curve contained in `S_F`, hence a component of the
nonproperness curve, and `(36)` is its parametrized intersection with `N_a`.
In particular, the live condition `b!=0` names a boundary-generated Jelonek
component meeting the forced nodal cubic.

### 3.1 A node endpoint is repeated and singular

Let `y_0` be an endpoint and suppose

```text
phi_D(y_0)=o=(-a/2,0),                                  (38)
```

the unique node of `N_a`. Since `dP_o=0`, the chain rule gives

```text
dG_(0,y_0)=dP_o compose dPhi=0.
```

At an endpoint `G=0`, while `u` is a unit, so `(34)` yields

```text
dq_(0,y_0)=0,                 E'(y_0)=0.                 (39)
```

Hence a `D`-visible target node forces both

```text
y_0 is a repeated root of E,
the strict companion V(q) is singular at (0,y_0).       (40)
```

Consequently, if `E` is squarefree, every endpoint is a smooth transverse
point of the strict companion mapping to a smooth nonproper value of `N_a`;
none pays THM-3996's target-node alternative. If `R(0,y)=0`, equivalently
`E=gamma` or `E/gamma=1`, the strict companion is boundary-disjoint and
`N_a` is disjoint from `Phi(D)`. Any occurrence of the node in `S_F` must
then be witnessed by a valuation with no center on this finite boundary `D`,
such as an escape to infinity inside `X_2`.

The scalar converse "repeated `E` implies target node" is false even when
the map remains etale across the displayed boundary. In local coordinates
take

```text
(x,y) |-> (r,s)=(x,y),                  Jac=1,           (41)
D=V(x).
```

The smooth target curve `r-s^2=0` pulls back to the smooth strict curve
`q=x-y^2=0`, with `dq_(0,0)=dx`, but its restriction to `D` is `-y^2`.
Thus an etale map can produce a repeated endpoint purely by tangency to the
boundary, while both the target curve and its pullback remain smooth. Scalar
endpoint multiplicity therefore does not distinguish a target node from
smooth tangency, just as THM-3994's double resultant does not distinguish two
addresses from one length-two center.

The missing sidecar is now exact:

```text
endpoint polynomial E
  preserves: parametrized intersection multiplicities with N_a;
  destroys: target labels, node versus smooth tangency, component owners;
  restore with: phi_D=(A_D,C_D) and the normalized-address graph;
  cheapest node test: gcd(E,A_D+a/2,C_D), then dq at its roots.             (42)
```

## 4. Scope and next boundary

This theorem proves neither `JC(2)` nor emptiness of the reduced `(2,3)`
cell. More precisely:

1. Section 1's factor insertion preserves algebraic normal-form data, not
   bracket equations or Darboux realizability.
2. The observer kernel disproves an inference from endpoint/class/clutch and
   finite rows to irreducibility; it does not prove an actual companion is
   reducible.
3. Section 2 excludes only the named odd factors `1+p^M(p^3-y^2)` and gives
   a necessary invoice for their even analogues. It is not a factorization
   classification.
4. Endpoint existence proves a nonproper value on `N_a`, not that the value
   is its node, that all node addresses are known, or that `S_F` has no other
   component.
5. A repeated endpoint is necessary for a `D`-visible node and insufficient
   for identifying it.

The cheapest continuation is to compute the boundary map `(33)`, rather than
only its scalar composite `(35)`, and to test actual companion factors by
their genus, punctures, covering degree, and two normalization-side exit
values. **QED.**

## Reproduction

```bash
python3 04-computation/jc2_companion_observer_kernel_thm4011.py
python3 -O 04-computation/jc2_companion_observer_kernel_thm4011.py
python3 04-computation/jc2_companion_observer_kernel_thm4011_independent_audit.py
python3 -O 04-computation/jc2_companion_observer_kernel_thm4011_independent_audit.py
sha256sum 04-computation/jc2_companion_observer_kernel_thm4011.py \
  05-knowledge/results/jc2_companion_observer_kernel_thm4011.out \
  04-computation/jc2_companion_observer_kernel_thm4011_independent_audit.py \
  05-knowledge/results/jc2_companion_observer_kernel_thm4011_independent_audit.out
python3 agents/check_docs.py
```
