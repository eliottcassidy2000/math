---
id: THM-3237
title: "Degree-nine Jacobian infinity wall and square-root escape"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In both THM-3212 cubic accessory algebras, the nonlinear clutch
  B_t=1+t*x^9 is the first monomial degree B=1+t*x^d that changes the
  infinity degree of the THM-3225 critical resultant.  Its saturated
  residual has generic degree 55.  At the exact wall t=2/Gamma its first
  two infinity coefficients vanish, the nonzero kappa_infinity carry leaves
  degree 53, and two roots escape with a reciprocal square-root law.  Exact
  good reductions certify squarefree, boundary-disjoint generic and wall
  controls.  This produces critical-point obstructions, not a Keller mate or
  a conclusion about JC(2).
source: root/creative-reframes/2026-08-03
audit: >
  The self-contained assertion-independent companion derives the universal
  y-resultant, checks the degree ledger d=0 through 9, divides every t
  coefficient by the passport boundary in both cubic accessory algebras,
  verifies the two leading cancellations and kappa carry, computes root-free
  wall and carry norms, and checks the reciprocal local law.  Good reductions
  (4111,113,85) and (3211,101,64) certify the t=2 and t=2/Gamma controls as
  squarefree and disjoint from the owner boundary.  Normal, optimized, and
  stored transcripts agree byte-for-byte.  An independent hostile audit
  rederived the gradient/resultant, complete degree-one-through-nine ledger,
  wall coefficients, carry norms, degree drop, reciprocal law, and
  good-reduction implications through a separate symbolic path.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape
related:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3172-shear-invariant-differential-owner-filtration-and-transverse-recurrence
  - THM-3229-hasse-pluecker-simple-root-contact-gcd-flag-and-degree-termination
script: 04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py
output: 05-knowledge/results/jc_heptic_degree_nine_infinity_wall_thm3237.out
script_sha256: 92518f258afeca233e90790fa2f713fcfd375c295271ff85dc4c5c66c0057d81
output_sha256: 38599c5d2a9b527098274d0dfccd427b5f091528671df168ca3a7ffb31c3b9cb
hash_basis: LF-normalized bytes
---

# THM-3237 -- degree-nine Jacobian infinity wall and square-root escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem isolates the first nonlinear `B`-degree in the monomial
family `B=1+t*x^d` that changes the infinity face of THM-3225's critical
resultant.  It does not remove the critical obstruction: the wall loses two
roots at infinity but leaves a squarefree residual of degree `53`.

## 1. Universe and theorem statement

Let `K_i` be either cubic accessory field of THM-3212 and let `K_i -> K_0`
be any embedding into an algebraically closed characteristic-zero field.
Retain

```text
Gamma=-7A_0,
V=4SDT^2/Gamma^2,                A=A_src=2SET/Gamma,
g=ST,                            2VA'-AV'=2V.             (1)
```

The embeddings cover the current one `(4,1,1,1)` and three `(3,2,1,1)`
unmarked heptic covers.  Consider

```text
B_t=1+t x^9,
P_t(x,z)=(V(x)z^2+B_t(x)z)^2+A(x)z+x.                    (2)
```

Let

```text
boundary_(4111)=S^3T^8x^9,
boundary_(3211)=S^3T^8x^6(x-1)^3.                        (3)
```

Both boundaries are monic of degree `44`.  The exact companion proves,
coefficient by coefficient in `t`, that `(3)` divides the critical
resultant for `(2)`.  If `H_t` denotes the quotient after the inherited
`V^3` saturation, then

```text
deg_x H_0=52,
deg_x H_t=55       for t!=0 and Gamma*t!=2,
deg_x H_(2/Gamma)=53.                                    (4)
```

The two-unit loss in `(4)` is an infinity escape with a square-root local
law.  Exact good reductions prove that `t=2` and `t=2/Gamma` are respectively
squarefree degree-`55` and degree-`53` critical controls, disjoint from `g`.

Assertions `(1)--(4)` and the local law below are **PROVED**.  The printed
modular rows are **FINITE-EXACT** controls; their good-reduction consequences
are used only for the two named characteristic-zero specializations.

## 2. Gradient reduction and universal resultant

On `V!=0` put

```text
y=Vz,                    L=y^2+B_t y,
J_t=2VB_t'-B_tV'.                                          (5)
```

Exactly as in THM-3225, multiplying the two gradient equations by powers of
`V` and subtracting `(V'y/2)` times the first gives

```text
R_1=2L(2y+B_t)+VA,
R_2=V^3+V^2y+J_t yL.                                      (6)
```

These generate the localized gradient ideal.  Direct symbolic elimination
gives

```text
Res_y(R_1,R_2)=V^3K_t,                                    (7)
```

where, writing `B=B_t` and `J=J_t`,

```text
K_t=
 -A^3J^3+12A^2J^2V^2-4AB^3J^2V+4AB^2JV^2
 +24ABJV^3-48AJV^4-16AV^4-8B^4JV^2+8B^3JV^3
 +32B^2V^4-96BV^5+64V^6.                                 (8)
```

Exact division of every parameter coefficient gives residual-degree lists,
indexed by `t^0,...,t^5`,

```text
B=1+t x^8:  (52,50,48,46,50,48),
B=1+t x^9:  (52,52,52,52,55,55),                         (9)
```

in both accessory fields.

## 3. Why degree nine is minimal in this family

For the monomial family `B=1+t*x^d`, use

```text
deg V=16,                         deg A=8.                (10)
```

For `1<=d<8`, the generic covariant has degree `d+15`.  At `d=8` its
nominal leading coefficient is proportional to `2d-16` and cancels, so
`deg J<=22`.  Substitution in the twelve terms of `(8)` gives the exact
maximum-degree ledger

```text
d=0,1,2,3,4,5,6,7,8:     deg K=96,
d=9:                     deg K=99.                       (11)
```

For `d<=8`, the degree-`96` face is the unique `64V^6` term.  At `d=9`,
the terms `-4AB^3J^2V` and `8B^3JV^3` reach degree `99`.  Their combined
coefficient is generically nonzero, as Section 4 computes.  Since `(3)` has
degree `44`, degree nine is therefore the first `d` in this stated monomial
family that changes the saturated infinity degree from `52` to `55`.

This is a minimality statement only inside `B=1+t*x^d`; it does not classify
arbitrary nonlinear `B,C_0`, or `E_0` deformations.

## 4. The primitive infinity wall

Write the leading response jets as

```text
V=v x^16+v_15 x^15+v_14 x^14+...,
A=a x^8+a_7 x^7+a_6 x^6+....                             (12)
```

From `(1)`,

```text
v=4/Gamma^2,                   a=2/Gamma,
a v_15=2v a_7,
a_6=a(4vv_14-v_15^2)/(8v^2).                             (13)
```

Substitution in `(8)` gives the first two infinity coefficients

```text
[x^99]K_t=-16t^4v^3(at-v),
[x^98]K_t=-72t^4v^2v_15(at-v).                           (14)
```

The factor `t^4` is support degeneration: at `t=0`, the degree-nine term of
`B_t` disappears.  The genuine primitive wall is

```text
t=t_*:=v/a=2/Gamma.                                      (15)
```

At this wall both coefficients in `(14)` vanish.  Define

```text
kappa_infinity=4vv_14-3v_15^2.                           (16)
```

The next coefficient is

```text
[x^97]K_(t_*)=-2v^6 kappa_infinity/a^4
             =-512 kappa_infinity/Gamma^8.              (17)
```

The exact accessory-field norms are

```text
Norm_(4111)(kappa_infinity)
 =494424620106921/3276800000000,

Norm_(3211)(kappa_infinity)
 =-102215864475620375014841556549072265625
   /12116574790945106558976.                             (18)
```

They are nonzero, so `(17)` is nonzero under every characteristic-zero
embedding.  Because the boundary `(3)` is monic, `(14)--(17)` imply the
exact degree drop `55 -> 53` in `(4)`.

The root-free primitive representatives of `Norm(Gamma*t-2)` are

```text
4111:
1250000000t^3+5131065625t^2+7626867738t+4689453125,

3211:
23887872t^3+293981184t^2-339674272t-14068359375.         (19)
```

Thus the wall can be recovered without selecting an accessory embedding.

## 5. Reciprocal square-root escape

Set

```text
w=1/x,                         delta=t-t_*,
epsilon=Gamma*t-2=Gamma*delta,
Hhat(t,w)=w^55H_t(1/w).                                  (20)
```

If `h_j(t)=[x^j]H_t`, exact quotient arithmetic gives

```text
h_55(t)=2048t^4(2-Gamma*t)/Gamma^8,
h_55(t_*)=h_54(t_*)=0,
h_53(t_*)=-512kappa_infinity/Gamma^8,
h_55'(t_*)=-32768/Gamma^11.                              (21)
```

Consequently, in the completed local ring at `(delta,w)=(0,0)`,

```text
Hhat=
 -32768 delta/Gamma^11
 -512 kappa_infinity w^2/Gamma^8
 +O(delta^2,delta*w,w^3).                                (22)
```

The nonzero coefficients in `(18),(22)` give two Puiseux branches over an
algebraic closure,

```text
w^2~-64 delta/(Gamma^3 kappa_infinity)
   =-64 epsilon/(Gamma^4 kappa_infinity).                (23)
```

Thus exactly two units of resultant intersection multiplicity move to
`x=infinity` at `(15)`.  Equation `(23)` is a ramified square-root escape,
not the one-root affine boundary loss of THM-3225.

## 6. FINITE-EXACT controls and characteristic-zero consequences

The exact companion uses the same good accessory reductions as THM-3225:

```text
passport   (p,u)       t control       deg H   gcd(B,g) gcd(H,g) gcd(H,H')
4111       (113,85)    2               55      0        0        0
4111       (113,85)    2/Gamma=85      53      0        0        0
3211       (101,64)    2               55      0        0        0
3211       (101,64)    2/Gamma=89      53      0        0        0.          (24)
```

The reductions of `kappa_infinity` are respectively `56 mod 113` and
`16 mod 101`.  The accessory roots are simple, `Gamma` is a unit, and all
displayed degrees are preserved.  Therefore a nonconstant characteristic-zero
gcd would reduce to a nonconstant gcd in `(24)`.  The rows certify that the
exact `t=2` and `t=t_*` residuals are squarefree and disjoint from `g`, while
their clutches are units modulo `g`.

The leading `y`-coefficient of `R_1` is the constant `4`.  Hence these simple,
boundary-disjoint resultant roots support respectively `55` and `53`
reduced transverse critical points, equivalently Morse critical points of
`P_t`, on the current four covers.  The baseline `t=0` recovers THM-3225's
degree-`52` constant clutch.

The finite-field enumerations in `(24)` are **FINITE-EXACT**.  They are not
an empirical sample over accessory roots: good reduction is used only for
the stated gcd nonvanishing implications.  No assertion is made about every
exceptional parameter outside the two named controls.

## 7. Failure anatomy and scope

The affine no-go left the degree-`96` term frozen.  Degree nine genuinely
changes that mechanism: it creates a degree-`99` face and a new primitive
wall where two roots escape quadratically.  The hoped-for repair still fails
at the first decisive consequence:

```text
new infinity channel -> disappearance of the critical resultant          FALSE.
```

At the wall, `53` reduced critical points remain.  Therefore the wall control
has no polynomial constant-Jacobian mate.  The generic control `t=2` likewise
has none.  The calculation does not construct a second coordinate, a marked
inverse pair, or an input to THM-3172's owner filtration.  It proves neither
`JC(2)` nor `DC(2)`, and it does not classify other nonlinear source
deformations.

The strongest survivor is the exact second-order infinity carry `(16),(23)`:
it is a new deformation coordinate and a reproducible obstruction object,
not a Jacobian-conjecture solution.

## 8. Exact reproduction

Run

```text
python3 04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py
python3 -O 04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py
```

and compare LF-normalized bytes with the declared output.  The companion is
self-contained, uses exact rational/accessory polynomial arithmetic, and has
no `assert` or `__debug__` gate.  Ordinary and optimized modes match the
stored transcript byte-for-byte.  The independent audit also reconstructed
the two wall polynomials and carry norms by a separate symbolic path.

QED.
