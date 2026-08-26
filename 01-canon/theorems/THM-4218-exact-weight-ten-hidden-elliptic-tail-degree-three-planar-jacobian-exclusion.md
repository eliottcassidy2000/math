---
id: THM-4218
title: "Exact-weight-ten hidden elliptic-tail degree-three planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4007/4012/4045/4103/4120/4122/4147
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the inherited b=d=0 reduced
  (2,3) seam, the complete exact-M=10 locus with nonzero p^5, p^2y^2, y^3,
  and separated top roots contains no nonautomorphic planar Keller pair. A
  multiplicity-one j=0 tail has two attachments differing by exact
  three-torsion, so specialization forces cover degree divisible by three;
  the only carrier responses have degrees 31 and 25. JC(2), DC(2), entry,
  other cells, and the coefficient walls remain OPEN.
source: codex-planar-jacobian-weight-ten-session-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
related:
  - THM-4217-complete-mixed-off-antidiagonal-delta-zero-planar-jacobian-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the audited use in
  THM-4045. The external theorem supplies the face/edge model and rational
  toric chains; it proves none of the in-repo face arithmetic, attachment,
  carrier-response, or planar-Jacobian claims.
script: 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218.py
output: 05-knowledge/results/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218.out
script_sha256: 438c8b33c16b74fb5fd6997169ad22256a288a99353b82dc5ea86faa1276804e
output_sha256: 06afd580fbe0999c2c050b98609450c353f47ea49d2e666d6f45b1b7d24fe05d
independent_audit_script: 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218_independent_audit.out
independent_audit_script_sha256: 90aa3481a25111d2cd2b1f6c5cf0377142901b45f6fc584c87636435d40b769f
independent_audit_output_sha256: 46c4b3dd85ca7a2b209fa54492b22387d174446fcb95a2d37c4d260975f69a5a
---

# THM-4218 -- exact-weight-ten hidden elliptic-tail degree-three exclusion

**PROVED RELATIVE TO THM-3992/3997/4007/4012/4045/4103/4120/4122/4147
+ VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Let

```text
upsilon=[p^5]H,       xi=[p^2y^2]H,       zeta=[y^3]H. (1)
```

> **Theorem.** The complete exact-weight-ten coefficient locus
>
> ```text
> upsilon*xi*zeta*(upsilon+xi) != 0                    (2)
> ```
>
> contains no nonautomorphic planar Keller pair.

Every lower coefficient allowed by the inherited normal form is arbitrary.
The theorem begins inside the named seam; it does not prove entry or `JC(2)`.

The inheritance pass is deliberately component-sensitive:

- the closest proved mechanism is THM-4045's complete lower-Newton regular
  model, not a highest-face proxy;
- the canonical hostile is the side elliptic tail: the new genus-two main
  component has no elliptic quotient, but that alone does not close the map;
- the corrected near miss is to call equal attachment images impossible.
  Here they are possible through an Eisenstein degree-three isogeny;
- the least-used sidecar is the *difference* of the two labelled attachments,
  retained together with the prime-carrier response degree.

The decisive obstruction is therefore not `Hom=0` by itself. The tail
attachment forces divisibility by three, while both complete responses have
degree one modulo three.

## 2. Complete exact-weight-ten source

Use the rational source coordinates

```text
s=XT,             p=T+s^2,             y=sp,
P=T+X^2T^2,       Y=XTP,               t=p-s^2=T.       (3)
```

THM-3992/3997 gives

```text
G=-s^2/(2t)+H(p,y),       [y]H=[py]H=0,                (4)
```

after the maintained gauge. Enumerating `0<2i+3j<=10` and deleting exactly
the two forbidden rows leaves

```text
p, p^2, p^3, y^2, p^2y, p^4, py^2, p^3y, y^3,
p^5, p^2y^2.                                             (5)
```

Thus the complete polynomial is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta y^3
  +upsilon p^5+xi p^2y^2,                              (6)

K=2848/45-(7/6)Delta.                                  (7)
```

There is no ellipsis in `(6)`. In source coordinates the last two terms are
`upsilon p^5+xi s^2p^4`. Condition `(2)` makes the maximum residual weight
exactly ten and retains the weight-nine cubic owner.

For a generic pencil value put `Q=q^-1` and clear the source denominator:

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (8)
```

On the torus of `(8)`, `p=s^2` would imply `F_Q=-Qs^2/2!=0`. Hence `t` is
nonzero there and

```text
t=p-s^2,                    X=s/t                       (9)
```

recovers the original source function field. This is the actual Keller
source curve, not an auxiliary Newton curve.

## 3. Universal two-face lower model

A monomial `p^i y^j` in `H` contributes the two height-one endpoints

```text
(j+2,i+j,1),                 (j,i+j+1,1).              (10)
```

Put `(r,k)` for the first two coordinates. The only lower height planes are

```text
M: nu_M=(r+2k-2)/10,
T: nu_T=(r+k-2)/6.                                    (11)
```

For either endpoint in `(10)`, its gap above `M` is
`1-(2i+3j)/10`. Thus every term below weight ten is strictly above `M`, and
the two top terms lie on it. On `T`, equality occurs only for the first
endpoints of `zeta y^3` and `xi p^2y^2`; all other allowed endpoints lie
strictly above it.

The possible coincident coefficient of `(2,5,1)` is an interior point of an
`M` edge. Even if it cancels, the three outer vertices remain uniquely owned
by `upsilon,xi,zeta`. Consequently all `2^8` lower-support subsets have the
same two-face hull under `upsilon*xi*zeta!=0`. The two lower faces and global
polygon are

```text
M: (0,1),(2,0),(4,4),(0,6),
T: (2,0),(5,3),(4,4),

P_global=(0,1),(2,0),(5,3),(4,4),(0,6).               (12)
```

Their Pick ledgers are

| polygon | `2Area` | boundary | interior |
|---|---:|---:|---:|
| `M` | 30 | 10 | 11 |
| `T` | 6 | 6 | 1 |
| global | 36 | 12 | 13 |

The outer primitive faces, with monomials and units suppressed, are

```text
AB: rational linear,
BC: q-1/2=K W^2+zeta W^3,
CD: zeta s+xi p,
DE: (p-s^2)(upsilon p+xi s^2).                        (13)
```

The vertical edge `EA` lies on the affine source divisor `s=0` and is not an
infinity puncture. THM-4103 gives the complete candidate packet

```text
(9,9,6,2,2,2,1),       n=31,       defect=24=2*13-2. (14)
```

The two index-nine roots are distinct exactly when `upsilon+xi!=0`. The
index-six root is simple when `xi*zeta!=0`. The index-two face is one prime
separable cubic closed point over `C(q)` because the polynomial map
`W |-> 1/2+KW^2+zeta W^3` has degree three. No lower coefficient can alter
these conclusions.

The generic outer restrictions are linear on `AB,CD`; the two separated
linear roots on `DE`; the cubic `H_3(W)-q` on `BC`; and a polynomial
`H_5(W)-q` of degree five and leading coefficient `upsilon` on `EA`.
In characteristic zero, `H_d(W)-q` is coprime to its derivative over `C(q)`:
a common derivative root is constant and cannot make the transcendental `q`
constant. Hence the toric generic completion is smooth. By `(9)` and
uniqueness of the smooth projective model, it is the actual source
completion.

## 4. Multiplicity-one regular special fibre

Take the single base change

```text
Q=sigma^30.                                             (15)
```

The integral lower heights are

```text
M: 3r+6k-6,                  T: 5r+5k-10.              (16)
```

Their primitive three-dimensional normals have last coordinate one. Thus
Dokchitser's face/edge model, in the exact form already audited in THM-4045,
gives both face components multiplicity one. The six special edge schemes
are, up to a nonzero monomial and scalar,

```text
AB: 1-Z,                 BC: 1-zeta Z^3,
CD: zeta+xi Z,           DE: (1-Z)(upsilon+xi Z),
EA: 1-upsilon Z^5,       BD: 1-xi Z^2.                (17)
```

Every scheme in `(17)` is reduced under `(2)`, and no root meets a corner.
The exact adjacent slopes at the internal edge `BD` are `-5>-6`, so no
intermediate toric chain occurs; its two reduced roots give two direct,
labelled attachment paths.

The `M`-face polynomial factors exactly as

```text
g_M=(S^2-P)(1-upsilon P^5-xi S^2P^4)=R*C.             (18)
```

The rational branch `R` and the other branch `C` meet at

```text
P=S^2,                    S^10=(upsilon+xi)^-1.        (19)
```

There are ten such points. Their gradient determinant is

```text
-10S^9(upsilon+xi) != 0,                               (20)
```

so all intersections are transverse and there are no others at the
compactified boundary.

In the main chart

```text
s=sigma^-3 S,              p=sigma^-6 P,
H_M=sigma^30 H(sigma^-6P,sigma^-9SP),                 (21)

U=S^2-P,                   V=(1-H_M)/S^2.
```

The exact scaled equation is

```text
UV=sigma^30/2.                                           (22)
```

Thus each point in `(19)` is an `A_29` smoothing. Regular resolution inserts
only multiplicity-one rational curves along each of the ten paths.

For `C`, put `Y=SP^2`; its normalization is

```text
xi Y^2=1-upsilon P^5.                                  (23)
```

This is a smooth genus-two curve. The `T` face normalizes under
`T_0=SP,W=SP^2` to

```text
xi W^2=1-zeta T_0^3,                                   (24)
```

a smooth `j=0` elliptic curve `E`. The two derivative equations for the
untransformed face would force simultaneously
`3zeta S+2xi P=3zeta S+4xi P=0`, impossible on its torus, so no singularity
is hidden before normalization.

The two roots of the internal edge `BD` attach `E` to `C` at two distinct
points. The normalized special fibre therefore has three positive/core
vertices `R,C,E`, ten `R--C` paths, and two `C--E` paths. Its dual graph has

```text
b_1=12-3+1=10,
g_special=0+2+1+10=13.                                (25)
```

All other edge, toric, and resolution components are rational. This matches
the global Pick genus in `(12)`, proves connectedness, and proves that `(14)`
is the complete boundary packet rather than a candidate subset.

## 5. The genus-two and attachment sidecars

Let `rho` be a primitive fifth root. On `(23)`, the automorphism

```text
P |-> rho P
```

acts on the holomorphic differentials `dP/Y,PdP/Y` with characters
`rho,rho^2`. Hence `Q(zeta_5)` embeds in `End^0(J(C))`.

The Jacobian `J(C)` is simple. If it were isogenous to two nonisogenous
elliptic curves, a degree-four field could not embed into either elliptic
endomorphism algebra. If it were isogenous to `E_1^2` with `E_1` non-CM,
`M_2(Q)` has maximal commutative degree two. If `E_1` has imaginary-quadratic
CM field `L`, then `L` is central in `M_2(L)`. The only quadratic subfield of
`Q(zeta_5)` is the real field `Q(sqrt(5))`, so the compositum would be an
eight-dimensional commutative `Q`-algebra inside `M_2(L)`, whose commutative
maximum has dimension four. Every alternative is impossible. Therefore

```text
Hom(J(C),E_0)=0                                        (26)
```

for every elliptic curve `E_0`, and every morphism from `C` to an elliptic
curve is constant.

After scaling `(24)` to

```text
E_0: y^2=x^3+1,
```

the two labelled attachments are

```text
P_+=(0,1),                    P_-=(0,-1).              (27)
```

The tangent at `P_+` is horizontal, so `2P_+=P_-`; their difference is a
nonzero point of exact order three. If a nonconstant morphism
`f:E_0->E_0` satisfies `f(P_+)=f(P_-)`, its homomorphic part kills that
three-torsion point. It is an isogeny with kernel order, hence degree,
divisible by three. This is sharp: the Eisenstein isogeny `1-zeta_3` has
kernel `{O,P_+,P_-}`.

## 6. Good target reduction and degree conservation

Since `q=sigma^-30`, scale the inherited target coordinates by

```text
A=sigma^-10 X,                  C=sigma^-15 Y.          (28)
```

The exact target model becomes

```text
Y^2=X^3+1-(3a^2/4)sigma^20X-(a^3/4)sigma^30,           (29)
```

with smooth special fibre `E_0:y^2=x^3+1`.

A hypothetical Keller pair gives a finite nonconstant generic morphism from
the source completion of Section 3 to this target. Resolve its rational
extension over the regular source model. Every new exceptional component is
rational. The branch `R`, every toric/resolution chain, and every exceptional
component maps constantly to `E_0` by Riemann--Hurwitz. Equation `(26)` makes
`C` constant as well. Constants agree across their connected subcurve.

Both attachments of the multiplicity-one tail `E` meet that same constant
subcurve, so the specialized tail map satisfies

```text
phi_E(P_+)=phi_E(P_-).                                  (30)
```

Let `L` be a relative target line bundle of degree one. Proper flat degree
conservation and the multiplicity-one ledger give

```text
deg(phi_generic)
 =sum_i m_i deg(phi_i^*L)
 =deg(phi_E).                                           (31)
```

The right side is zero or divisible by three by Section 5. The generic map
is nonconstant, so

```text
3 | deg(phi_generic).                                  (32)
```

## 7. The two carrier responses

The cubic in `(13)` is a prime separable degree-three residue extension.
THM-4120 gives `E_q(C(q))={O}`, and THM-4122/4147's source-independent
finite-carrier lemma leaves exactly two responses.

In the full response, all places in `(14)` map to the target origin, so

```text
deg(phi_generic)=9+9+6+2+2+2+1=31.                    (33)
```

In the finite response, the carrier maps birationally to one horizontal
degree-three target point. Its three geometric conjugates, each of index two,
are removed from the origin fibre, leaving

```text
deg(phi_generic)=31-3*2=25.                            (34)
```

Both `31` and `25` are congruent to one modulo three, contradicting `(32)`.
Thus neither response exists and the theorem follows.

As an independent positive control on the response method, the exact rational
coefficient choice

```text
Delta=1, K=5591/90, Phi=2, Theta=5, eta=7, zeta=11,
upsilon=13, xi=17                                      (35)
```

has residual source-critical degree `25`; restoring four universal Morse
points gives `L=29`. The inherited full and finite permutation caps are
`4<24` and `23<24`. This proves the critical-open subchamber is nonempty, but
the control is not needed for the all-lower-coefficient divisibility proof.

## 8. Sharp boundary and next fronts

Every unit in `(2)` performs a visible proof step:

```text
upsilon=0:          loses the degree-five outer owner;
xi=0:               loses the tail vertex and attachment pair;
zeta=0:             makes the side face rational;
upsilon+xi=0:       merges the two index-nine roots and makes (19) nontransverse.
                                                                    (36)
```

These walls require new polygons or strict transforms. The result proves
nothing about them, another reduced cell, seam entry, exact `M>=11`, `JC(2)`,
or `DC(2)`.

The typed connection is

```text
source:       complete exact-M=10 lower Newton model;
target:       cover degree modulo three;
map:          Q=sigma^30 regular specialization of the Keller fibre map;
preserved:    component multiplicity, degree, and labelled attachments;
destroyed:    lower coefficients and rational resolution-chain lengths;
sidecars:     genus-two simplicity, attachment difference, carrier response;
decisive test: 31,25 are nonzero modulo 3;
hostile:      any wall in (36).                         (37)
```

## 9. Verification

The primary certificate expands the complete weight-ten support, enumerates
the two lower faces, checks Pick and packet ledgers, all special and generic
edge schemes, face smoothness, ten transverse nodes, the exact `A_29` chart,
component genera, the order-five and order-three sidecars, target scaling,
and both response degrees. The independent audit uses a separate support and
Euclidean path and rechecks the regular-model and degree-conservation inputs.

Replay with

```bash
python3 -B 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218.py
python3 -B -O 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218.py
PYTHONHASHSEED=4218 python3 -B 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218.py

python3 -B 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218_independent_audit.py
python3 -B -O 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218_independent_audit.py
PYTHONHASHSEED=4218 python3 -B 04-computation/jc23_weight10_hidden_elliptic_tail_degree3_exclusion_thm4218_independent_audit.py
```

Compare stdout with the two frozen outputs. **QED.**
