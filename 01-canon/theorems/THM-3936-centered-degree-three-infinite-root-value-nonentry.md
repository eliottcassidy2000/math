---
id: THM-3936
title: "Centered degree-three infinite-root-value maps have a non-unibranch ramification component"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the centered
  trace-zero linear-color binary-cubic
  grammar, assume a one-place discriminant component has polynomial
  normalization A1_u, deg A=3, and a degree-three repeated-root map with
  t(infinity)=infinity. Local trace and the finite Riemann--Hurwitz budget
  leave four collision-free pole rows. Exact color divisibility leaves one
  family in each row. Their implicit quintic or septic occurs to exponent
  one in the order discriminant, hence is genuine tame transposition
  ramification in the maximal cubic completion. Every surviving family has
  at least two normalization addresses over one target point; rank-three
  finite flatness forces those addresses onto one source ramification point.
  That source component is non-unibranch, so THM-3920 forbids an affine-plane
  Keller open. Together with THM-3933 this closes all values of t(infinity)
  in the stated centered degree-three grammar, not arbitrary root gauges,
  higher root degree, or JC(2).
source: jc_zero_debt_lift / complementary t(infinity)=infinity stratum after THM-3933, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root and zero_debt_lift child,
  2026-08-23). The audit independently reconstructed the exhaustive pole
  partition from completed-local trace and the finite Riemann--Hurwitz budget,
  including all four diagonal normal-form gauges. It checked every primitive
  incidence row, color resultant, cancellation seam, and coefficient-ideal
  boundary; the degree-one-in-color irreducibility argument; the exponent-one
  order-to-maximal discriminant bridge; and every multiple-address fibre. It
  also verified that rank-three finite flatness coalesces those addresses on
  the source ramification prime, so the conclusion is source non-unibranchness
  rather than only a branch-image collision. An imprecise early
  "equivalently" was repaired: trace zero and the prime degree, not target
  normalization alone, force the repeated root to generate. The assertion-free 65-gate
  companion byte-matches in normal and optimized mode; raw and semantic hashes
  and documentation checks pass.
depends_on:
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3930-two-pole-linear-color-aligned-one-place-branch-packet
  - THM-3931-degree-two-pole-cubic-principal-ramification-no-atlas
  - THM-3938-centered-degree-four-root-map-nonentry
  - THM-3941-all-degree-centered-cubic-pole-carrier-routing
script: 04-computation/jc2_centered_degree_three_infinite_root_value_thm3936.py
output: 05-knowledge/results/jc2_centered_degree_three_infinite_root_value_thm3936.out
script_sha256: c77b4a28a9b3a87cc57b364801c8d48dff88cf8bc25d2dd800306a1f1ed5c37d
output_sha256: df140a5f3e275ce18dafcc11ef55bc7217bdb42b80f5c6640b7613ef41db3197
semantic_sha256: 12141212df57e3d3e67e9a206682b64e0c2cbf64198096c4c57c4196f5c37adb
hash_basis: raw LF bytes
---

# THM-3936 -- infinity leaves four rows, and every row folds its ramification arm

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
characteristic zero. Consider the linear-color binary cubic

```text
Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3.                       (1)
```

Assume that an irreducible one-place discriminant component has normalization
`A1_u`, with

```text
A=A(u), C=C(u) in k[u],              deg A=3,               (2)
```

and that its repeated root `t=U/V in k(u)` satisfies

```text
deg(t:P1_u -> P1_t)=3,               t(infinity)=infinity.  (3)
```

We use the same centered root gauge as THM-3933: the primitive incidence
equation is

```text
a(A)t^3-c(A)t-2d(A)=0,                                  (4)
```

so `Tr_(k(u)/k(A))(t)=0`. We impose the coefficient-ideal gate

```text
(a,C,c,d)=k[A,C].                                         (5)
```

The normalization hypothesis includes `k(A,C)=k(u)`. Separately, trace zero
and nonconstancy give `t notin k(A)`: an element of the base field would have
trace `3t`. Since `[k(u):k(A)]=3` is prime, it follows that
`k(A,t)=k(u)`, and `(4)` is its minimal polynomial after removal of
coefficient content.

The conclusion is that none of these data can be the finite completion of a
Keller `A2` open. Together with THM-3933, it closes both possible values of
`t(infinity)` inside this fixed centered, degree-three, linear-color grammar.
It does not say that an arbitrary binary cubic can be centered by a root
translation without destroying its linear-color shape.

## 1. The complete pole and ramification ledger

Put `pi_A:P1_u -> P1_A`. Since `A(u)` is a cubic polynomial, infinity is
totally ramified and contributes two of the four Riemann--Hurwitz units.
Exactly two finite ramification units remain.

At a point of ramification index `e`, choose tame completed coordinates
`z=s^e`. The trace of a Laurent series in `s` keeps precisely the exponents
divisible by `e`, multiplying their coefficients by `e`. Two consequences
will be used repeatedly:

```text
an isolated simple finite pole requires e=2 or 3;
an isolated double finite pole requires e=3.               (6)
```

At infinity `e=3`. A pole of exact order three has a nonzero `s^(-3)`
term, which survives trace, so trace zero forbids it. Hence

```text
m_infinity=ord_infinity(t)_infinity is 1 or 2.              (7)
```

The remaining finite pole degree is respectively two or one. Suppose first
that two distinct finite pole supports have the same `A`-value. At either
pole the repeated homogeneous root is `[U:V]=[1:0]`. In the `U=1` chart,

```text
Phi=a+Cw+cw^2+dw^3,                  w=V/U,                 (8)
```

so multiplicity at `w=0` forces `a=C=0`. Thus the two supports are two
normalization addresses over the same target point. Section 2 proves that
this already makes the genuine source ramification component non-unibranch.

We may therefore classify the remaining rows assuming distinct pole
supports have distinct `A`-values. Combining `(6)`, `(7)`, and the two-unit
finite ramification budget gives exactly the following table. Affine changes
of `u,A`, a nonzero diagonal rescaling of `[U:V]`, and the compensating scalar
rescaling of `C` give the displayed normal forms.

```text
row  m_inf  finite poles       finite A-ramification  normalized A,t

A      2    one simple        e=2                    A=u^3+u^2,
                                                     t=u^2+pu+(p-1)/3+k/u

B      2    one simple        e=3                    A=u^3,
                                                     t=u^2+pu+k/u

C      1    one double        e=3                    A=u^3,
                                                     t=u+p/u+k/u^2

D      1    two simple        e=2 and e=2            A=u^3-3u,
                                                     t=u+q/(u-1)+k/(u+1).
                                                                  (9)
```

Here `k!=0` in rows A--C and `qk!=0` in row D; otherwise the stated finite
pole disappears. The constants in `(9)` are not guesses. For
`A=u^3+u^2`,

```text
Tr(u)=-1, Tr(u^2)=1, Tr(1/u)=0,
```

which gives the constant `(p-1)/3` in row A. For `A=u^3`, all displayed
nonconstant monomials have trace zero. For `A=u^3-3u`, `Tr(u)=0`, while
the residues lie at the two critical points and

```text
Tr(1/(u-1))=Tr(1/(u+1))=0.
```

This also proves that the table contains every trace-zero normal form, not
only four sample families.

## 2. Shared pole values already fold the source arm

The generic repeated root on the chosen component is exactly double. If it
were triple, its incidence derivative `3at^2-c` would vanish, contradicting
separability of the minimal polynomial `(4)`. At a generic finite nonzero
root, varying `C` changes the cubic value by `t^2`, so the universal
discriminant is transverse to the linear-color plane. The chosen component
is reduced and its odd discriminant valuation cannot be erased by the square
of an order-to-maximal-order index. It is genuine tame `(2,1)` ramification.

Let `B` be the maximal normalization of `k[A,C]` in the cubic field and let
`E` be the residue-degree-one ramification prime over the component. The
ring `B` is a normal finite module over the regular surface `k[A,C]`; normal
surfaces are Cohen--Macaulay, so miracle flatness makes `B` finite flat of
rank three. Every closed point of `E` is non-etale and consumes geometric
fibre length at least two. A length-three fibre cannot contain two distinct
such points.

Consequently the two addresses obtained from `(8)` land at the same point
of `E`. Since `E^nu` is finite birational to the component normalization
`A1_u`, they are two branches of `E` at that point. THM-3920 says every
irreducible boundary curve of an affine-plane open in a normal affine
surface is unibranch. A Keller open must delete `E`, so every shared-pole
row is impossible.

## 3. Exact color division in the four collision-free rows

For each row, eliminate `u` from `A-A(u)` and the denominator-cleared
equation `T-t(u)`. Its primitive trace-zero resultant supplies `a,c,d` in
`(4)`. The double-root derivative then determines the color:

```text
C(u)=-(3a(A(u))t(u)^2+c(A(u)))/(2t(u)).                   (10)
```

Polynomiality of `(10)` is rigid. After cancelling known monomial factors,
the exact numerator/denominator resultants are

```text
row A: -243 k^3(p-1)^6(27k+6p-2)^3,
row B: -27 k^6p^6,
row C: -27 k^6p^6,
row D: 1728 k^3q^3(k+q+2)^6.                              (11)
```

All boundary cases in `(11)` matter:

* In row A, `k!=0`. The branch `p=1` divides exactly. On the other resultant
  seam `k=(2-6p)/27`, polynomial division leaves the exact remainder

  ```text
  -(p-1)^3(3p-1)^3(3u-1)(3u+2)/27.                       (12)
  ```

  Thus only `p=1`, or `p=1/3,k=0`, can divide. The latter deletes the finite
  pole and belongs to a lower-degree stratum.
* In rows B and C, `k!=0`, so `(11)` forces `p=0`; direct division is exact.
* In row D, `qk!=0`, so `(11)` forces `q+k+2=0`; again direct division is
  exact. The excluded factors `q=0` and `k=0` are exactly the two finite-pole
  cancellations.

Write `delta=2q+2` in row D. The complete survivor table is

```text
row  A(u)       t(u)                         a,c,d
                C(u)

A    u^3+u^2    u^2+u+k/u=(A+k)/u           -A, -(A+k)^2, -(A+k)^3/2
                u(3u+4)(A+k)/2

B    u^3        u^2+k/u                     -A, 0, -(A+k)^3/2
                3u^2(A+k)/2

C    u^3        u+k/u^2                     -A^2, 0, -(A+k)^3/2
                3u^4(A+k)/2

D    u^3-3u     u+q/(u-1)-(q+2)/(u+1)       -(A^2-4),
                3(u^2-1)(u^2-5)(A+delta)/2   3(A+delta)^2,
                                              -(A+delta)^3/2.       (13)
```

Rows A--C satisfy `(5)` exactly when `k!=0`: `a` is a power of `A`, while
`d(0)=-k^3/2`. In row D,

```text
Res_A(A^2-4,A+delta)=delta^2-4=4q(q+2),                  (14)
```

so the unit ideal holds exactly when both finite poles are genuine. The
four rows therefore survive every coefficient-ideal test. The obstruction
below is geometric, not a hidden failure of `(5)`.

## 4. Four genuine one-place ramification curves

Let `H_i(A,C)` be the implicit equation of the parametrized component in
row `i`. Exact elimination and discriminant calculation give

```text
Res_u(A-A(u),C-C(u))=-H_i/8,                              (15)

Disc_T(Phi_A)=-(A+k)^3 H_A/4,
Disc_T(Phi_B)=-(A+k)^3 H_B/4,
Disc_T(Phi_C)=-(A+k)^3 H_C/4,
Disc_T(Phi_D)=-(A+delta)^3 H_D/4.                         (16)
```

The first three equations are

```text
H_A=27A^5+81A^4k+16A^4+36A^3C+81A^3k^2+48A^3k
    +72A^2Ck+27A^2k^3+48A^2k^2-4AC^2+36ACk^2+16Ak^3
    -8C^3-4C^2k,

H_B=27A^5+81A^4k+81A^3k^2+27A^2k^3-8C^3,

H_C=27A^7+81A^6k+81A^5k^2+27A^4k^3-8C^3.                (17)
```

For the last row,

```text
H_D=27A^7+81A^6delta+81A^5delta^2-648A^5-108A^4C
    +27A^4delta^3-1944A^4delta-216A^3Cdelta
    -1944A^3delta^2+2160A^3-108A^2Cdelta^2+432A^2C
    -648A^2delta^3+6480A^2delta-36AC^2+864ACdelta
    +6480Adelta^2-8C^3-36C^2delta+432Cdelta^2
    +2160delta^3.                                         (18)
```

Thus A and B are quintics, while C and D are septics; their normalization
has one polynomial infinity. Each `H_i` occurs to exponent **one** in the
order discriminant. Since passage to the maximal order changes a
discriminant valuation by an even nonnegative index valuation, exponent one
forces index zero at the generic point of `H_i` and maximal discriminant
valuation one. Therefore `H_i` is not an order artefact: it supports a
genuine tame residue-degree-one ramification prime `E_i`.

For completeness, every `H_i` is irreducible and `(13)` is its birational
normalization. Over `k(A)`, the cubic in `(1)` is linear in `C`. In any
factorization over `k(A)[C,T]`, one factor is `C`-independent and divides
both `T^2` and `aT^3+cT+d`; since `d!=0`, that factor is a unit. Hence the
generic cubic is irreducible. Moreover `t notin k(A)`: at the indicated
finite pole its pole order is respectively

```text
row A: 1 mod e=2;  row B: 1 mod e=3;
row C: 2 mod e=3;  row D: 1 mod e=2.                      (19)
```

Because the extension degree is prime, `k(A,t)=k(u)`. At a generic exact
double root, Euclid's algorithm recovers `t` rationally from the cubic and
its derivative, so `t in k(A,C)`. Thus `k(A,C)=k(u)`, and `(15)` is the
irreducible implicit equation rather than a power.

## 5. The surviving rows all identify too many addresses

The factor visible in `(13)` gives a single target point with multiple
normalization addresses:

```text
row A: A=-k, C=0; roots of u^3+u^2+k=0;
row B: A=-k, C=0; roots of u^3+k=0;
row C: A=-k, C=0; roots of u^3+k=0;
row D: A=-delta, C=0; roots of u^3-3u+delta=0.             (20)
```

Rows B and C have three distinct addresses because `k!=0`. Row D has three
distinct addresses because `delta^2!=4`, precisely `(14)`. Row
A has three distinct addresses except at `k=-4/27`, where it has two; it
never has only one because `k!=0` and a cubic `u^3+u^2+k` cannot be a cube
of a linear polynomial. Thus every row has at least two distinct points of
`H_i^nu=A1_u` over the point in `(20)`.

Apply the rank-three fibre argument of Section 2 to the genuine prime `E_i`.
If two normalization addresses landed at two different points of `E_i`, the
fibre would contain two non-etale points, each of length at least two,
contradicting total length three. Hence all the addresses in `(20)` land at
one point of `E_i`. They are distinct branches there. The **source
ramification curve itself** is non-unibranch; this is stronger than merely
saying that its branch image has a self-collision.

Every Keller `A2` open in the maximal completion would have to delete
`E_i`, but THM-3920 forbids a non-unibranch irreducible boundary curve.
Therefore none of the four survivor rows admits an affine-plane Keller
atlas.

## 6. Exact scope and the paired closure

THM-3933 proves the same nonentry when the centered degree-three root map is
finite at normalization infinity; the trace there forces the value to zero.
This theorem proves the complementary infinite-value case. Hence the pair
closes

```text
centered trace-zero incidence + fixed linear color + deg A=deg t=3. (21)
```

Still open are noncentered root gauges that cannot be translated without
mixing the color coefficient, root degree at least four, other binary-cubic
coefficient planes, alternative nonmonogenic cubic orders, and JC(2).

## Reproduction

```bash
python3 04-computation/jc2_centered_degree_three_infinite_root_value_thm3936.py
python3 -O 04-computation/jc2_centered_degree_three_infinite_root_value_thm3936.py
```

The frozen transcript is
`05-knowledge/results/jc2_centered_degree_three_infinite_root_value_thm3936.out`.
The companion is assertion-free and checks the local trace ledger, every
primitive resultant, all exceptional division seams, the four survivor
identities, coefficient-ideal boundaries, implicit equations, discriminant
factorizations, and address counts.
