---
id: THM-4291
title: "Lambda-zero genus-five tail and degree-forty-two equivariant hostile"
status: >
  PROVED FORMAL-LOCAL RELATIVE TO THM-4103/4230 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On one allowed W=0, Lambda=0 specialization, the
  genuine weighted exceptional face normalizes to a smooth genus-five tail
  carrying an abstract deck-equivariant degree-42 map to E0. The Keller
  invariant differential nevertheless vanishes there to vertical order
  eight, so an actual Keller extension is constant on that tail. A general
  positive-order statement is proved only for generically separable
  divisorial faces; repeated-root refinements, full degree persistence,
  exact-M=12 exclusion, and JC(2) remain open.
source: root/jc2-signal-20260830
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4288-a23-partial-normalization-relative-differential-and-etale-base-change-obstruction
related:
  - THM-4284-a23-conductor-defect-and-degree-shell-first-character-nondescent
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
script: 04-computation/jc23_lambda_zero_genus5_degree42_tail_thm4291.py
output: 05-knowledge/results/jc23_lambda_zero_genus5_degree42_tail_thm4291.out
script_sha256: 505ee506f08fbffb77452981d4b7f69ed6edc7c32421e095fb9e790b04ae526c
output_sha256: 46bc1f5f167dc2cf24b68b7496e3a6c04885315f5423db65b632d0db730ed561
independent_audit_script: 04-computation/jc23_lambda_zero_genus5_degree42_tail_independent_audit_thm4291.py
independent_audit_output: 05-knowledge/results/jc23_lambda_zero_genus5_degree42_tail_independent_audit_thm4291.out
independent_audit_script_sha256: 21d678579e6591e60e74b59ce1bbc2a8be8404c53f896543dc93894f9eee4ded
independent_audit_output_sha256: 99883743241680266e105be9f7b0d2bbd4e23b46907042c2cb63aad3d0678a5a
modpoly_audit_script: 04-computation/jc23_lambda_zero_genus5_degree42_tail_modpoly_audit_thm4291.gp
modpoly_audit_output: 05-knowledge/results/jc23_lambda_zero_genus5_degree42_tail_modpoly_audit_thm4291.out
modpoly_audit_script_sha256: deb5bdd0382810b6ae2d2a9fd94424abf1e38e73bca30e91ba0432c17acf6100
modpoly_audit_output_sha256: 7b7c88594aada79fe7fab13086f671d1ad25f0f5b7df15ebbce26b55bcf700dc
audit: >
  PASS. A primary SymPy path reconstructs the literal THM-4230 infinity
  chart, extracts its complete weight-24 face, normalizes the discriminant,
  checks smoothness, deck characters, the reciprocal elliptic quotient, the
  degree-42 eigenmap arithmetic, and the Keller form's vertical order. A
  standard-library sparse-polynomial and Q(sqrt(21)) audit independently
  rebuilds the face, quotient, special coefficient, and degree ledger. PARI
  independently computes the exact Phi_7(0,J) factorization. Normal,
  optimized, and fixed-hash-seed Python runs agree with the frozen outputs.
external: >
  Standard characteristic-zero facts used below are Weierstrass preparation,
  the genus and regular differentials of a squarefree hyperelliptic curve,
  the modular-polynomial criterion for a cyclic isogeny, the Hodge/Rosati
  degree pairing, and the fact that a morphism of smooth characteristic-zero
  curves with zero differential is constant. None supplies the exact local
  equation or arithmetic certificate.
---

# THM-4291 -- Lambda-zero genus-five tail and degree-forty-two equivariant hostile

**PROVED FORMAL-LOCAL RELATIVE TO THM-4103/4230 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE ABSTRACT TAIL MAP IS NOT A KELLER
SPECIALIZATION. REPEATED-ROOT REFINEMENTS, DEGREE PERSISTENCE, `M=12`, AND
`JC(2)` REMAIN OPEN.**

## 1. Statement and scope

Work over `C` in THM-4230's reduced `(2,3)`, exact-weight-twelve family. On
the `W=0`, `Lambda=U+Z=0` wall put `Z=-U` and impose the allowed lower-row
specialization

```text
Delta=Phi=Theta=eta=zeta_3=upsilon_5=xi_10=alpha_11=beta_11=0.       (1)
```

Retain the forced weights `2,4,6`, and set

```text
a=7168/135.                                                          (2)
```

> **Theorem.** At the infinity contact of the exact integral source:
>
> 1. The balanced weighted blowup has weights
>    `wt(sigma,b,q)=(1,1,12)`. Its exceptional face normalizes to
>
>    ```text
>    T_5: Z^2=(X^6-a)^2+4U.                                         (3)
>    ```
>
>    If `U!=0` and `18225U+12845056!=0`, this is a smooth genus-five
>    curve. The deck action is
>
>    ```text
>    delta:(X,Z)->(zeta_12^(-2)X,Z),                                (4)
>    ```
>
>    and its differential characters are `10,8,6,4,2 mod 12`; the target
>    character `10` occurs.
> 2. For two explicit algebraic values of `U`, the abstract curve with action
>    `(T_5,delta)` admits a nonconstant equivariant map to
>    `E_0:Y^2=X^3+1` of degree exactly `42`. Thus genus, deck character,
>    isogeny class, and degree alone do not exclude tail transfer.
> 3. The actual Keller pullback of the good-target invariant differential is
>
>    ```text
>    omega_0=sigma^8(-X^4 dX/Z+O(sigma)).                            (5)
>    ```
>
>    Hence it vanishes to exact vertical order eight on this tail. On every
>    resolution where the Keller map extends, its restriction to the strict
>    transform of `T_5` is constant. In particular, the degree-`42` map in
>    part 2 is an exact hostile to curve/action-only reasoning but is **not**
>    the Keller specialization.
> 4. More generally, every divisorial face above this contact whose prepared
>    quadratic initial form is generically separable has strictly positive
>    Keller-differential order. Repeated initial roots and their subsequent
>    Newton faces are not covered.

The closest proved mechanism is THM-4103's Keller residue identity. The
canonical hostile is the degree-`42` map in part 2. The corrected near miss is
the replacement of a hypothetical genus-eleven smoothing by the literal
equation below. The least-used decisive sidecar is not another character: it
is the **vertical order of the Keller differential**.

This theorem proves no Keller point exists on the displayed coefficient
specialization only if degree persistence to the remaining component is
separately supplied. No such result is claimed here. It neither promotes the
reserved THM-4290 nor closes exact `M=12` or `JC(2)`.

## 2. The literal infinity equation

Use THM-4230's exact integral variables `(sigma,S,P)` and put

```text
b=1/S,             r=P/S^2,             t=sigma*b.                   (6)
```

Multiplying THM-4230 equation (9) by `b^14` gives exactly

```text
F=(1-r)(b^12-Hhat(r,t))-sigma^12*b^12/2.                             (7)
```

Under `(1)`, with `K=2848/45`, the complete specialized polynomial is

```text
Hhat=U(r^6-r^4)
     -(1376/135)t^6r^3+(2848/45)t^6r^2
     +(8/3)t^8r^2-3t^10r.                                           (8)
```

The factor `b^12` in the smoothing term is load-bearing. Omitting it produces
the tempting but false local model `q(q-cb^12)=sigma^12*unit` and a spurious
genus-eleven tail at the wrong scale.

Set `q=r-1` and make the balanced substitution

```text
b=sigma X,                  q=sigma^12 Q.                            (9)
```

Every term of order below `24` vanishes. Dividing by `sigma^24` and setting
`sigma=0` gives

```text
F_0=2UQ^2-Q(X^12-aX^6)-X^12/2.                                     (10)
```

The fixed weight-six sum in `(8)` is exactly

```text
2848/45-1376/135=7168/135=a,                                       (11)
```

so it cannot be discarded as a unit correction.

As a quadratic in `Q`, `(10)` has discriminant

```text
Disc_Q(F_0)=X^12((X^6-a)^2+4U).                                    (12)
```

Writing

```text
4UQ-(X^12-aX^6)=X^6Z                                                (13)
```

normalizes the nonreduced `X^12` factor and gives `(3)`. Its squarefree
resultant is a nonzero rational unit times

```text
U^6(18225U+12845056)^5.                                             (14)
```

Thus the conditions in part 1 are necessary and sufficient for the displayed
degree-twelve branch polynomial to be squarefree; the standard hyperelliptic
genus formula gives `g(T_5)=5`.

## 3. Deck action and the character warning

The Kummer deck generator acts by

```text
sigma -> zeta_12 sigma,
(S,P) -> (zeta_12 S,zeta_12^2 P).                                  (15)
```

Therefore `b->zeta_12^(-1)b`, while `r` and `Q=q/sigma^12` are fixed.
Equation `(9)` then gives `(4)`. The standard basis

```text
X^j dX/Z,                    0<=j<=4,                               (16)
```

has character exponents

```text
-2(j+1)=10,8,6,4,2 mod 12.                                         (17)
```

The good target action is

```text
alpha:(x_E,y_E)->(zeta_12^4 x_E,-y_E),
alpha^*(dx_E/y_E)=zeta_12^(-2) dx_E/y_E.                            (18)
```

Thus the target character occurs already on `dX/Z`. Character absence cannot
kill this genuine face.

## 4. An exact equivariant degree-42 curve map

Put `d=a^2+4U`, choose `rho^12=d`, and scale `(3)` to

```text
C_c: eta^2=xi^12-2c xi^6+1,              c=a/rho^6.                 (19)
```

The reciprocal quotient

```text
f:C_c -> E_c,
w=xi^2+xi^(-2),                 v=eta/xi^3,                         (20)
E_c: v^2=w^3-3w-2c
```

is a degree-four morphism. Indeed the equation in `xi` over a generic `w` is
`xi^4-wxi^2+1=0`, and direct substitution verifies the elliptic equation.
Moreover

```text
f^*(dw/v)=2(xi^4-1)d xi/eta.                                       (21)
```

The elliptic `j`-invariant is

```text
j(E_c)=1728/(1-c^2)=1728+432a^2/U.                                 (22)
```

An independent PARI computation gives

```text
Phi_7(0,J)=J^2
 (J^2+34848505552896000J+11356800389480448000000)^3.                (23)
```

The nonzero roots are

```text
J=-17424252776448000
  +/-3802283679744000 sqrt(21).                                    (24)
```

For either root set

```text
U=432a^2/(J-1728).                                                  (25)
```

Then `j(E_c)=J`, so over a finite extension of `C(J)` there is a cyclic
degree-seven isogeny `psi:E_c->E_0` (use the dual if the modular polynomial
chooses the opposite direction). Also

```text
d=a^2J/(J-1728)!=0,                                                 (26)
```

so the tail is smooth.

Let `lambda=zeta_12^(-2)`, let `delta(xi,eta)=(lambda xi,eta)`, and put
`g=psi o f`. Base every map at either point at infinity. If
`e_k=g o delta^k`, equation `(21)` gives

```text
e_2=e_1-e_0.                                                        (27)
```

The pullback differentials agree, and the common base value removes the
possible translation constant. Since `alpha^2-alpha+1=0`, the integral map

```text
m=e_0-alpha e_1                                                   (28)
```

satisfies

```text
m o delta=alpha o m.                                               (29)
```

If `eta_E` is the invariant differential on `E_0`, then for a nonzero scalar
`s`

```text
g^*eta_E=s(xi^4-1)d xi/eta,
m^*eta_E=s(lambda^2-1)d xi/eta!=0.                                 (30)
```

The two character summands in `g^*eta_E` are orthogonal and have equal Hodge
energy: the reciprocal involution exchanges them. Since `deg(g)=4D` for
`D=deg(psi)`, each energy is `2D`, while

```text
|lambda^2-1|^2=3.                                                   (31)
```

Consequently

```text
deg(m)=6D=42.                                                       (32)
```

This map is integrally nontrivial. It is not divisible by the norm-three
prime `2-alpha`. Such a quotient would be an equivariant degree-`14` map with
pullback differential proportional to `d xi/eta`. The four deck-fixed points
`(0,+/-1)` and the two infinities must all map to the unique alpha-fixed point
`O`. Their fibre multiplicities are `1,1,5,5`, because `d xi/eta` has zero
order four at each infinity. Every other point in that fibre occurs with equal
multiplicity in a free six-element deck orbit. Hence every such map has degree
`12+6k`. This excludes both `14` and `34`; equation `(32)` realizes `42`.

Equations `(19)--(32)` are geometric over `C` after adjoining the displayed
roots and isogeny data. They do not assert descent to the unextended rational
seam field.

## 5. The Keller differential remembers the vertical scale

The abstract map above is a hostile only until the Keller residue is restored.
THM-4103 gives

```text
dA/(2C)=sigma^12 ds/(F_Q)_p.                                      (33)
```

For THM-4230's good target and integral source scalings

```text
A=sigma^(-4)x_E,       C=sigma^(-6)y_E,
s=sigma^(-1)S,         p=sigma^(-2)P,                              (34)
```

equation `(33)` becomes

```text
omega_0=dx_E/(2y_E)=sigma^9 dS/G_P.                                (35)
```

Here `G` is THM-4230 equation (9). Since `F=b^14G`, `P=r/b^2`, and
`S=1/b`, this is exactly, up to the global residue orientation,

```text
omega_0=-sigma^9 b^10 db/F_r.                                     (36)
```

On the weighted chart `(9)`, relative to the `sigma`-base,

```text
db=sigma dX,
F_r=sigma^12(F_(0,Q)+O(sigma))
   =sigma^12(X^6Z+O(sigma)).                                       (37)
```

Substituting `(37)` into `(36)` proves `(5)`:

```text
omega_0=sigma^8(-X^4 dX/Z+O(sigma)).                               (38)
```

The leading form `X^4dX/Z` is a regular nonzero differential on `T_5`, so
the vertical order is exactly eight, not merely positive. But restriction of
the full relative differential to `sigma=0` is zero. If the resolved Keller
map restricted nonconstantly to this smooth characteristic-zero tail, the
pullback of the nowhere-zero invariant differential on `E_0` would be
nonzero. Therefore that restriction is constant.

Notice the repaired character ledger: `X^4dX/Z` has exponent `2`, while the
factor `sigma^8` contributes exponent `8`; together they give the target
exponent `10`. Discarding vertical order creates the false impression that
the special differential should lie in the exponent-`10` tail line.

## 6. A generically separable divisorial-order lemma

Weierstrass preparation at the same `A_23` contact writes the exact equation,
up to a unit, as

```text
q^2+A(sigma,b)q+B(sigma,b)=0,
v(B)=12(v(sigma)+v(b)).                                             (39)
```

Let `v=v_E` be a divisorial valuation centered over `sigma=b=q=0`, and put

```text
s=v(sigma)>0,       beta=v(b)>0,       gamma=v(q)>0,
a_0=v(A),
N=min(2gamma,a_0+gamma,12(s+beta)).                                 (40)
```

At least two entries in the minimum agree because the prepared equation
vanishes. Assume its initial quadratic is generically separable along `E`,
equivalently

```text
v(2q+A)=N-gamma.                                                    (41)
```

The relative differential satisfies `v(db)>=beta`. Equations `(36)` and
`(41)` give

```text
v_E(omega_0)>=9s+11beta+gamma-N.                                  (42)
```

There are three possible cancelling pairs:

```text
2gamma=12(s+beta)       => gamma=6(s+beta),
a_0+gamma=12(s+beta)    => gamma>=6(s+beta),
2gamma=a_0+gamma        => gamma<=6(s+beta).                        (43)
```

In all three cases `(42)` yields

```text
v_E(omega_0)>=3s+5beta>0.                                          (44)
```

For the explicit genus-five face, `(s,beta,gamma)=(1,1,12)` and `(42)` is
the exact order eight computed in `(38)`.

If the initial quadratic has a repeated root, `(41)` can fail because
`v(2q+A)` increases. A translation and a subsequent Newton face are then
required, and the prepared constant term no longer carries enough information
by itself. This theorem makes no assertion about those refinements. That is
the precise remaining obstruction to a uniform exceptional-tail extinction
theorem.

## 7. Consequences and failure boundary

The source/target map exposed here is

```text
source:     exact M=12 infinity chart or its normalized exceptional tail;
target:     the good elliptic model E_0;
map:        Keller graph extension, versus the abstract map m;
preserved:  deck equivariance and total curve-map degree;
destroyed:  vertical sigma-order when one keeps only the special curve;
sidecar:    v_E(omega_0);
test:       equation (38), then repeated-root Newton refinement.             (45)
```

Thus:

1. A tail-free proof based only on genus, characters, isogeny class, or degree
   is false: `(32)` is a minimal exact hostile at the required degree.
   On this tail equivariance excludes degree `34` but permits degree `42`.
2. For the displayed genuine tail, the Keller differential repairs the false
   inference and proves constancy.
3. For every generically separable exceptional face, `(44)` supplies the same
   repair.
4. Repeated-root faces remain open. Until they are controlled, one cannot put
   all degree `34/42` on the genus-seven face and cannot invoke THM-4249's
   degree-divisible-by-four deck eigenspace as a complete exclusion.
5. The raw `A_23` residue at `sigma=0` is dualizing and may have poles; this
   theorem does not turn it into the ambient-regular Kähler form required by
   THM-4284/4288. There is no contradiction with their descent firewall.

## 8. Reproduction

Run

```bash
python3 -B 04-computation/jc23_lambda_zero_genus5_degree42_tail_thm4291.py
python3 -B 04-computation/jc23_lambda_zero_genus5_degree42_tail_independent_audit_thm4291.py
gp -q -s 200000000 04-computation/jc23_lambda_zero_genus5_degree42_tail_modpoly_audit_thm4291.gp
```

The first path uses exact SymPy algebra. The second is a clean-room
standard-library sparse-polynomial and `Q(sqrt(21))` reconstruction. The third
asks PARI for the classical level-seven modular polynomial rather than trusting
the frozen quadratic factor. **QED.**
