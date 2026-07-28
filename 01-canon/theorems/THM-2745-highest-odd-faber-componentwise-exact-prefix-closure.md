---
id: THM-2745
title: "Highest-odd Faber componentwise exact-prefix closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every member of
  the full chosen-sheet split degree-22 exact-square-prefix response family
  with at least one nonzero odd Faber coefficient is physically empty,
  including reducible and nonreduced members.  Its complete infinity divisor
  is five transverse points plus the highest-odd G2 branches.  A physical
  image component would contain exactly one response-pole point; pole order
  greater than one forces the unique source pole to a finite zero of U, where
  polynomial regularity of the exact prefix supplies two local lift equations.
  Their smooth-point resultant is nonzero and their G2 leading remainder is
  32q_*^6, excluding every component.  Combined with THM-2725, only the
  all-even zero-first-flux boundary remains in this degree-22 family.  The
  broader split branch, JC(2), and DC(2) remain open.
source: root/highest-odd-componentwise-exact-prefix-closure-2026-07-28
audit: componentwise-hostile-2026-07-28 (independent multinomial Faber reconstruction, infinity/component/DVR audit, corrected tangent-leading proof, and immutable replay)
depends_on:
  - THM-2719-full-split-odd-faber-generic-normalization-genus-four-hundred-nineteen
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2741-highest-odd-faber-response-pole-capacity-closure
related:
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
  - THM-2747-highest-odd-reduced-boundary-divisor-and-one-ended-factorization-atlas
script: 04-computation/jc2_degree22_highest_odd_componentwise_closure_20260728.py
output: 05-knowledge/results/jc2_degree22_highest_odd_componentwise_closure_20260728.out
script_sha256: 863923b4886f6ec7dee44e28e152109acc51fcfd778d550d391ee02fc577b48a
output_sha256: 619d2e49acdf60f5e2aa3ba9076957e2e38e9078363c833f7484f38311df3bb9
hash_basis: LF-normalized bytes
---

# THM-2745 -- the exact-square prefix closes every highest-odd component

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2741 makes the third response transverse on every highest-odd `G2`
branch and closes geometrically integral response members by a pole-count
contradiction.  A reducible member can distribute those branches among
components, so pole count alone stops.  The missing information is already
present on the source: once its sole response pole is forced to be finite,
the polynomial exact-square prefix has two regular coefficients whose local
lift equations exclude every possible infinity point one at a time.

## 1. Statement and response family

Work over `C` with the full chosen-sheet split degree-22 family

```text
F_23=Phi_22+sum_(j in J) a_j h^(22-j)Phi_j-lambda h^23,
G_24=Psi_22+sum_(j in J) a_j h^(22-j)Psi_j-W h^24       (1)
```

in `P(1,2,3,4)_[h,d,q,s]`, where

```text
J={1,2,3,5,6,7,9,10,11,13,14,15,17,19,21}.            (2)
```

Suppose at least one odd `a_j` is nonzero.  If `j` is the largest such
index, put

```text
r=22-j in {1,3,...,21},                g=gcd(r,6).     (3)
```

Then no reduced irreducible component of `(1)` can carry the response map of
a polynomial split exact-square-prefix Keller trajectory.  Since a reduced
source kills target nilpotents, the same conclusion holds for arbitrary
reducible or nonreduced members.

## 2. The complete infinity divisor

After removing nonzero scalar factors, the top degree-22 observables are

```text
Phi_22=unit*q*P,                 Psi_22=unit*Q,
R_22=unit*q*T,                                           (4)

P=-84d^2q^4s+3dq^6+280dq^2s^3-21q^4s^2-84s^5,

Q=-224d^3q^6+3360d^2q^4s^2-336dq^6s-3360dq^2s^4
  +3q^8+560q^4s^3+224s^6,

T=-84d^3q^4s+9d^2q^6+280d^2q^2s^3-105dq^4s^2
  -84ds^5+3q^6s+70q^2s^4.                              (5)
```

If `q=0`, equation `Q=224s^6` leaves only

```text
P_infty=[d:q:s]=[1:0:0].                               (6)
```

On the `q=1` index chart, pass to the coarse `mu_3` invariants

```text
z=s^3,                         rho=d/s^2.              (7)
```

The top equations become

```text
Pbar=3rho-21+z(-84rho^2+280rho-84),

Qbar=3+z(-336rho+560)
       +z^2(-224rho^3+3360rho^2-3360rho+224).          (8)
```

Their lexicographic basis is one linear equation in `rho` and

```text
p5(z)=20141047808z^5-14386462720z^4+1089822720z^3
      -21288960z^2-35910z+81.                          (9)
```

The polynomial `(9)` is squarefree; the ideals obtained by adjoining `z`,
`rho`, or the Jacobian of `(8)` are units.  Hence there are exactly five
distinct transverse coarse points, with

```text
q!=0,                         s!=0,                    d!=0,
ord(h)=1.                                                (10)
```

At those points the stripped response is

```text
Tbar=3+z(9rho^2-105rho+70)
       +z^2(-84rho^3+280rho^2-84rho),                  (11)
```

and `(Pbar,Qbar,Tbar)` is the unit ideal.  Thus `R_22` is nonzero and

```text
pole_order(R_25/h^25)=25                               (12)
```

at each smooth point.

At `(6)`, THM-2741 gives `3g` coarse normalization branches with

```text
ord(h)=6/g,
pole_order(R_25/h^25)=(150-6r)/g>=8.                  (13)
```

The invoice

```text
5*1+(3g)*(6/g)=23                                     (14)
```

is the full weighted hyperplane degree, so no infinity point or branch is
missing.  The exact gcd `gcd(qP,Q)=1` excludes a component inside `h=0` and,
by homogeneity, any common surface factor of `(1)`.

## 3. Every physical component would have one infinity point

Let `Y` be the normalization of the reduced irreducible component containing
the generic physical image.  The third-flux identity

```text
U R_source'=kappa,                         kappa!=0     (15)
```

makes the source map nonconstant.  It extends to a finite surjective morphism

```text
gamma:P1_x -> Y.                                       (16)
```

Every projective curve component meets the ample divisor `h=0`, and `(14)`
lists all its normalization points.  Every one is an `R_aff=R_25/h^25` pole.

THM-2723 classifies the source primitive as exactly one of

```text
U in C*,             R_source affine-linear,

U=u_0(x-a)^m,        R_source=s_0+s_1(x-a)^(1-m),
                     m>=2.                             (17)
```

In either case `R_source` has exactly one pole.  Distinct target poles have
disjoint nonempty inverse fibres under `(16)`, so `Y` has at most one infinity
normalization point and therefore exactly one.

If `e` is the source ramification index above it and `P_R` its response-pole
order, pullback gives source pole order `eP_R>=8e>1`.  This excludes the
constant-`U` line of `(17)`.  The unique source pole is consequently the
finite point `x=a`, and

```text
m-1=eP_R.                                               (18)
```

The divisibility boundary reserved originally under this theorem ID is thus
a proved intermediate statement.  Its real force is that all coefficients of
the original polynomial exact prefix are regular at `a`.

## 4. The two exact-prefix lift equations

Write the polynomial exact prefix as

```text
P_source=H^2+L,
H=U^2 z_source^2+Bz_source+C,
L=Az_source+E.                                         (19)
```

On the split sheet,

```text
w=Uz_source+B/(2U),
H=w^2+d_aff,                       L=q_aff w-s_aff.    (20)
```

Pull the weighted projective coordinates to the source DVR above `a`, passing
to its finite local index cover if necessary.  In this section `d,q,s` denote
homogeneous coordinates, so the affine coefficients in `(20)` are
`d/h^2,q/h^3,s/h^4`.  Put

```text
h w=(hU)z_source+beta,                beta=hB/(2U).    (21)
```

Coefficient comparison in `(19)--(21)` gives the exact identities

```text
beta^2+d=h^2 C,
q beta-s=h^4 E.                                          (22)
```

Here `hU` is regular and vanishes.  Since `d` and `C` are DVR-regular, the
first equation proves that `beta` is regular.  Let `omega_0` be its residue.
These equations are necessary physical lift conditions; the abstract response
curve forgot them.

## 5. The five smooth points are impossible

At a smooth infinity point, `(10)` and `(22)` give

```text
omega_0^2+d_0=0,                  q_0 omega_0-s_0=0,   (23)
```

with `q_0d_0s_0!=0`.  Substitute `d=-omega_0^2` and
`s=q omega_0` in the two top forms.  They factor as

```text
P=-8omega_0^2q^5(56omega_0^3+3q),

Q=q^6(7168omega_0^6+896omega_0^3q+3q^2).              (24)
```

After removing nonzero factors, the resultant in `q` is

```text
-76608 omega_0^6.                                     (25)
```

It cannot vanish because `d_0=-omega_0^2!=0`.  No physical component meets
one of the five smooth infinity points.

## 6. Every `G2` branch is impossible

At `P_infty`, the projective residues of `q,s` vanish; their **leading tangent
coefficients** must be used instead.  On the pulled-back `d=1` index cover,

```text
ord(h)=6e/g,                  ord(q)=ord(s)=re/g,
q=q_* t^(re/g)+...,           s=s_* t^(re/g)+...,
q_*s_*!=0.                                               (26)
```

The first equation of `(22)` gives `omega_0^2=-1`.  Since `E` is regular and
`r<=21<24`, the right side of the second has order at least `24e/g`, strictly
larger than `re/g`.  Hence its leading terms give

```text
s_*=q_* omega_0.                                       (27)
```

Every branch has `G2` tangent equation

```text
(q_*^2-s_*^2)(q_*^4-14q_*^2s_*^2+s_*^4)=0.            (28)
```

Substitution from `(27)` and `omega_0^2=-1` turns its left side into

```text
32 q_*^6!=0,                                           (29)
```

a contradiction.  Sections 5 and 6 exclude every possible unique infinity
point of `Y`, proving the theorem.

## 7. Reducibility, nilpotents, and exact residual

No integrality of the whole response member was used.  The generic image of
`P1` lies on one reduced irreducible component, and a morphism from a reduced
scheme kills target nilpotents and factors through the target reduction.
Embedded or generically nonreduced structure cannot change `(15)`, create a
second source pole, or erase `(22)`.

Combined with THM-2725, the exact residual inside `(1)` is

```text
split polynomial exact-square prefix;
reduced degree 22;
all odd Faber coefficients zero;
chosen-sheet first flux lambda=0.                      (30)
```

This does not derive `(1)` for every planar Keller pair, treat another reduced
degree or the broader split branch, or prove `JC(2)` or `DC(2)`.

## 8. Exact reproduction and hostile audit

Run

```bash
python3 04-computation/jc2_degree22_highest_odd_componentwise_closure_20260728.py
python3 -O 04-computation/jc2_degree22_highest_odd_componentwise_closure_20260728.py
```

against the declared output.  Normal and optimized modes byte-match the
stored transcript, both hashes agree, the companion compiles, and it contains
no optimization-sensitive `assert`.  It reconstructs all Faber faces, the
five-point Groebner basis, squarefreeness, nonvanishing and transversality
ideals, the complete hyperplane invoice, all eleven highest-odd pole rows,
the resultant `(25)`, and the tangent remainder `(29)`.

The hostile audit used generalized multinomial extraction rather than the
companion's recurrence, independently recovered the five-point quotient and
all response faces, and found the first draft's invalid substitution of the
zero projective `q,s` residues at `P_infty`.  After repair to the valuation
argument `(26)--(29)`, it rechecked weighted index covers, ramification,
component purity, finite-pole regularity, reducibility, and nilpotents.  It
found no remaining hypothesis gap.

QED.
