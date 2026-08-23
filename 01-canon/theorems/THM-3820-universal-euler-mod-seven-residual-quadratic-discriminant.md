---
id: THM-3820
title: "Universal Euler/mod-seven residual quadratic and source-P discriminant identity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  coefficient-free critical resultant for the
  pure-r nodal carrier has only three Euler layers.  After
  G=e^3Y, D=e^2(Y+Z), and t=e^7, its residual is an explicit quadratic F(t).
  The discriminant of F is 3^12 times the discriminant of the normalized
  source quadratic P times one exact square W^2.  Generically W is precisely
  the collision slope of the rational map from the two source-u sheets to
  their residual t-values.  This is a structural lemma only: it does not
  close degree-at-least-six carriers or any Jacobian-conjecture case.
source: jc_sparse_direct_search / universal residual discriminant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit rederived the Hamiltonian signs from the Poisson law, recovered
  P from the z-component and Q independently from the z^2/z^3 compatibility,
  and checked that G,D are genuinely coefficient-free variables before the
  profile specialization.  It then audited both Sylvester paths, the
  Euler/mod-seven substitution without division at e=0, the two source and
  residual discriminants, and the quadratic-algebra proof that W is exactly
  the sheet-collision slope.  The Y=0, Z=2, Y+2Z=0, R=0, and W=0 seams and
  the splitting-field quantifiers were checked separately.  Raw hashes
  match, and normal and optimized executions byte-match the frozen output.
  The deterministic companion recomputes the
  universal Sylvester resultant, all twenty terms of H, the three
  mod-seven Euler layers, the normalized quadratic, both discriminants, the
  square factorization, the direct normalized resultant, the generic
  quadratic-algebra reduction of the t-map, and sharp Y=0, Z=2, W=0, and R=0
  controls.  Normal and optimized runs byte-match the frozen transcript.
  Independent hostile audit is pending.
depends_on: []
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3817-quintic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_universal_euler_mod7_residual_quadratic_thm3820.py
output: 05-knowledge/results/jc2_universal_euler_mod7_residual_quadratic_thm3820.out
script_sha256: 0db08ef96e17c8e861b01b460d73054c121a7229dd23f75c8c81a64351eb0861
output_sha256: 5a84b2bbe0c9e5c973c3b0105ba4acac2ee916bafdc97e0bf675c4e726aa2942
semantic_sha256: 69696307c3a9767468ba2e738b92b423e9ac23af514d28433b3da79f1cc38d00
hash_basis: raw LF bytes
---

# THM-3820 -- the universal residual remembers the source double cover

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over a
field `k` of characteristic zero.
Let `e,G,D,u` be independent variables and put

```text
P = G u^2-(1+2u)(2e^3+u eD),
Q = e^2(1+2u)^3-729G^3u^2(1+u)^2.                    (1)
```

There is a polynomial `H in Z[e,G,D]` such that

```text
Res_u(P,Q)=G^3e^4H.                                   (2)
```

Under the polynomial substitution

```text
G=e^3Y,             D=e^2(Y+Z),             t=e^7,   (3)
```

the residual has exactly the three `e`-layers `3,10,17` and

```text
H(e,e^3Y,e^2(Y+Z))=e^3F(t,Y,Z),                      (4)

F=-(Y+2Z)-729t(Z-4)K+2125764t^2Y^3(Z-2)^2,
K=Z(Y+Z)^2+4(2-Z)(Y+2Z).                             (5)
```

Define

```text
R=(Y+Z)^2-8Z+16,
W=(Y+Z)(Z-4)^2+4Y(Z-2).                              (6)
```

Then the exact discriminant identity is

```text
Disc_t(F)=3^12 R W^2.                                 (7)
```

Moreover, after `(3)`, the normalized source equation is

```text
p=P/e^3=-(Y+2Z)u^2-(4+Y+Z)u-2,
Disc_u(p)=R.                                          (8)
```

Equivalently, before normalization,

```text
Disc_u(P)=e^2((D-4e^2)^2+8eG),
Disc_u(P)|_(3)=e^6R,
e^6 Disc_t(F)=3^12W^2 Disc_u(P)|_(3).                 (9)
```

The last display is a polynomial identity; it makes no localization claim at
`e=0`.

On the generic open set

```text
Y(Z-2)(Y+2Z)RW != 0,                                  (10)
```

the two roots of `F` generate the same quadratic extension as the two roots
of `p`.  More intrinsically, `W` records the information lost by sending a
source root `u` to its compatibility value

```text
phi(u)=(1+2u)^3/(729Y^3u^2(1+u)^2).                  (11)
```

Indeed, in the generic quadratic algebra `k(Y,Z)[u]/(p)`, the affine
coefficient of `phi` is

```text
[u] phi = (Y+2Z)W/(2916Y^3(Z-2)^2).                  (12)
```

Thus `R=0` merges the two source sheets, whereas, away from the degree-drop
divisors, `W=0` can keep the source sheets distinct while identifying their
two `t`-values.

## 1. Origin of the two source equations

The identities are coefficient-free, but their geometric source explains
why these particular polynomials matter.  On

```text
S=Spec k[r,z,e]/(r^2e-z^3+r)
```

use

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.     (13)
```

For the pure-r carrier

```text
A=e^2-z/3+r g(e),                                    (14)
```

its Hamiltonian components are

```text
{A,r}=r^2-9z^2(2e+r g'),
{A,z}=3gr^2-3(1+2re)(2e+r g'),
{A,e}=9gz^2-(1+2re).                                 (15)
```

Set `u=re`, then replace `g,g'` temporarily by independent symbols `G,D`.
The first equation in `(1)` satisfies

```text
P/e^2={A,z}/3.                                        (16)
```

The other two equations in `(15)`, together with the surface law, give

```text
z^2=(1+2u)/(9G),            z^3=u(1+u)/e.             (17)
```

Eliminating `z` from `(17)` gives `Q=0`.  This explains `(1)` without
assuming that `G` and `D` already come from a polynomial profile.

## 2. The coefficient-free Sylvester calculation

Direct expansion of the Sylvester determinant in `(2)` gives

```text
H=
 -729D^4e^2+1458D^3Ge+8748D^3e^4
 +2125764D^2G^3e^4-729D^2G^2-17496D^2Ge^3
 -34992D^2e^6-4251528DG^4e^3-8503056DG^3e^6
 +11664DG^2e^2+52488DGe^5+46656De^8-2De
 +2125764G^5e^2+8503056G^4e^5+8503056G^3e^8
 -2916G^3e-17496G^2e^4-23328Ge^7+G.                 (18)
```

Substituting `(3)` into these twenty terms shows that every `e`-exponent
is congruent to `3 mod 7`; the only exponents are `3,10,17`.  Collecting those
three layers is exactly `(4)`--`(5)`.  In particular, no coefficient profile,
degree bound, root condition, or division by `e,G,D` enters the quadratic
identity.

There is also a second check which preserves the elimination meaning.  Under
`(3)`, divide `P` by `e^3` and `Q` by `e^2`, replacing `e^7` by the independent
variable `t`.  This gives `p` from `(8)` and

```text
q=(1+2u)^3-729tY^3u^2(1+u)^2.                        (19)
```

Their direct resultant is

```text
Res_u(p,q)=Y^3F.                                      (20)
```

Thus, where `Y(Z-2)(Y+2Z)!=0`, the roots of `F` are precisely the
compatibility values of the roots of the source quadratic, counted with
multiplicity and with the degree-drop seams audited below.

## 3. Discriminants and the collision square

Write `F=at^2+bt+c`.  Since `2125764=4*3^12`, one has

```text
a=4*3^12Y^3(Z-2)^2,
b=-729(Z-4)K,
c=-(Y+2Z).                                            (21)
```

The whole factorization `(7)` reduces to the polynomial identity

```text
(Z-4)^2K^2+16Y^3(Z-2)^2(Y+2Z)=RW^2.                 (22)
```

Expanding either side proves `(22)`, hence `(7)`.  Separately, expanding `P`
from `(1)` as a quadratic in `u` gives

```text
P=(G-2eD)u^2-(4e^3+eD)u-2e^3,                        (23)
```

and its discriminant is exactly the first formula in `(9)`.  Substitution
`(3)` then proves `(8)`--`(9)`.

There is a reason the remaining factor in `(7)` is a square.  In the generic
quadratic algebra defined by `p`, the denominator of `(11)` is invertible:
`p(0)=-2`, while `p(-1)=2-Z`.  Euclidean inversion and reduction modulo `p`
give an affine function of `u`; its slope is `(12)`.  If `u_1,u_2` are the
two source roots, then

```text
phi(u_1)-phi(u_2)
 = (Y+2Z)W/(2916Y^3(Z-2)^2) * (u_1-u_2).             (24)
```

Squaring `(24)` is the conceptual content of `(7)`: the residual
discriminant is the source-sheet discriminant times the square of the
sheet-separation map.

## 4. Degree drops and hostile seams

The generic interpretation has four sharp boundary mechanisms.

1. `Y=0`.  The normalized direct resultant `(20)` has already lost its
   factor `Y^3`.  The divided residual nevertheless has the regular arm law

   ```text
   F(0,Z,t)=-Z(2+729t(Z-4)^3).                        (25)
   ```

   Before normalization this is

   ```text
   H(e,0,D)=De(-729e(D-4e^2)^3-2).                   (26)
   ```

   Thus `Y=0` is a divided special fibre, not a two-sheet `t`-image.

2. `Z=2`.  One source root reaches a pole of `(11)`, the coefficient of
   `t^2` vanishes, and

   ```text
   F(t,Y,2)=-(Y+4)+2916t(Y+2)^2.                      (27)
   ```

3. `Y+2Z=0`.  The source equation `p` drops from quadratic to linear.  The
   splitting-field statement therefore deliberately excludes this divisor.

4. `R=0` versus `W=0`.  For example,

   ```text
   (Y,Z)=(-1/2,5/2): R=0, W!=0,
   (Y,Z)=(-9/5,1):   W=0, R=216/25.                  (28)
   ```

   The first is a genuine source double root.  The second has two distinct
   source roots but a repeated residual value.  It is the canonical hostile
   example against treating `(7)` as an equivalence between the two branch
   divisors.

## 5. Exact scope and the higher-degree frontier

For an actual polynomial profile, one specializes

```text
G=g(e),              D=g'(e),
Y=g(e)/e^3,          Z=(e g'(e)-g(e))/e^3.            (29)
```

The theorem exposes a rigid double-cover skeleton behind every such profile:
the source discriminant `R` and the residual discriminant have the same
nonsquare part, while `W` measures loss under the compatibility map.  This is
the likely reason the quartic and quintic logarithmic-remainder systems admit
strong mod-seven normalizations.

It does **not** prove that a specialization of `H` has a root away from
`e g(e)=0`.  In `(29)`, `Y,Z,t` are linked rational functions of the same
variable `e`, not independent coordinates; the collision divisor `W=0`, the
source discriminant `R=0`, and boundary multiplicities can interact.  Hence
degree at least six, mixed carrier corrections, and the existence of a
polynomial Darboux mate all remain **OPEN**.

## 6. Exact replay

Run

```bash
python3 04-computation/jc2_universal_euler_mod7_residual_quadratic_thm3820.py
python3 -O 04-computation/jc2_universal_euler_mod7_residual_quadratic_thm3820.py
```

Both executions must byte-match

```text
05-knowledge/results/jc2_universal_euler_mod7_residual_quadratic_thm3820.out
```

The companion uses no inactive Python `assert` and freezes the full identity
packet and semantic packet.  It reports `CHECKS=35` and `RESULT=PASS`.
