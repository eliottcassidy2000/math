---
source: codex-2026-07-25-jc2-central-wall-r0
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; SUPERSEDED AS
  A CLOSURE STATEMENT BY THM-2345.
  In the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, the central-wall sublocus
  D=25B^2/126, R=20BC+21W=0, S!=0 is empty. Generic weighted ratios give a
  smooth genus-one blow-up cubic. Exactly two quadratic ratios give nodal
  rational cubics; the retained first-flux square Z=T^2 raises each to a
  genus-five double cover and excludes it. The later audited
  `jc2-degree18-central-twall-closure-opus-20260725` packet combines this
  face with S=0 and T=0 to close the full central ratio. Neither result
  proves JC(2).
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2313-degree-eighteen-bd-linear-ratio-closure
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-central-twall-closure-opus-20260725
---

# The degree-eighteen central `R=0` wall is empty

> **CANONICAL SUPERSESSION.** THM-2345 closes the entire common-root
> wall `126D=25B^2`.  This file is retained for its independent
> smooth/nodal normalization and Kummer genus-five atlas.

> **DIRECT SUBRESULTANT CLOSURE.** Promoted THM-2347 (`754346d86`)
> independently closes `20BC+21W=0`: the factorization `F=H S^2`,
> `deg H<=4`, forces `deg gcd(F,F')>=4`; after localization the selected
> principal-subresultant coefficients generate an ideal containing
> `(126D-25)^2`, so the branch lands on THM-2345's closed wall.  The
> remaining degree-eighteen debt is therefore genuinely off this wall.
> Its sharp next algebraic probe is the localized global
> principal-subresultant stratification
> `deg gcd(H_2 S_5^2,(H_2 S_5^2)')>=5` and
> `deg gcd(H_4 S_4^2,(H_4 S_4^2)')>=4`, away from
> `C(126D-25)(20C+21W)=0`.

> **NEXT PROVED REDUCTION, NOT A CLOSURE.** THM-2357
> (`THM-2357-degree-eighteen-h2-moving-root-reduction.md`, pushed candidate
> `b16cad056`) now attacks the first off-wall stratum.  The unique finite
> `e=3` branch of `H_2 S_5^2` gives a common root `P(r)=Q(r)=0`;
> `r=0` returns to THM-2345, while `r!=0` scales to `r=1`.  Smoothness on
> the normalization does **not** imply `Q'(r)!=0` in the polynomial-order
> plane model.  In the normalized chart
>
> ```text
> K=245+2835B+26244C,       Q'(1)=14K.
> ```
>
> Here `P,Q` are THM-2332's cleared uppercase covariants.  The older
> lowercase depressed-cubic coefficients differ by fixed rational scalars;
> their derivative and discriminant constants must not be compared without
> those scalars.  THM-2332 also supplies `H_2` squarefree, so
> `Disc(H_2)!=0` is inherited rather than a missing boundary.
>
> Dividing by `(y-1)^2` gives
>
> ```text
> R_10=4(y-1)p_3^3+49q_5^2=H_2 S_4^2.
> ```
>
> Direct top-down factorization reduces the identity to six equations
> `E_5,...,E_0` in `B,C,s_3,s_2`.  The first is linear in `C`, with pivot
>
> ```text
> L_piv=-10260B+1771s_2-1771.
> ```
>
> Hence there are two exact lanes: eliminate `C` when `L_piv!=0`, or impose
> `L_piv=0` together with the displayed irreducible degree-five residual
> `Phi`.  Neither lane is empty yet.  The `K=0` singular-order
> total-branch chart is explicitly open, as are `C=0`, the pivot wall,
> the localized five-equation main lane, and the `H_4 S_4^2` stratum.
> The theorem passes an independent algebra/scope audit as a reduction;
> it proves no new degree-eighteen emptiness statement and does not prove
> JC(2).

## 1. Inheritance and target

On THM-2262/2297's central ratio

```text
D=25B^2/126,
```

the repeated-branch equation splits into named weighted factors. This note
closes the whole factor

```text
R=20BC+21W=0,                 S!=0.                 (1)
```

The proof corrects an attractive but false shortcut: `S!=0` does **not**
imply that the original central tangent cone has three distinct lines.
Blowing up the central point is the faithful operation.

Throughout this note `B!=0`; the coordinate boundary belongs to the
already-routed lower-support strata.

## 2. The blow-up cubic

On (1), the depressed trigonal equation has

```text
p=y^2(245y^2+1890B),

q=y^2(539y^3+11340By+183708C).
```

Putting `U=yV` and removing the common exceptional factor gives the strict
transform

```text
H(V,y)
 =V^3+735y^2V+5670BV
   -3773y^3-79380By-1285956C.                       (2)
```

The map `(V,y) -> (yV,y)` is birational from (2) to the original spectral
curve away from the central exceptional divisor, so their normalizations
agree.

The homogeneous cubic is

```text
Hbar
 =V^3+735y^2V+5670BVw^2
   -3773y^3-79380Byw^2-1285956Cw^3.                (3)
```

At infinity, `r=V/y` obeys

```text
r^3+735r-3773=0,

Disc=-1972620783!=0.                                (4)
```

Thus all three points at infinity are distinct and smooth, uniformly in
the parameters. At `y=0`,

```text
H_y(V,0)=-79380B!=0,                                (5)
```

so the strict transform is smooth even when the original tangent cubic
has a repeated root.

## 3. Exact affine singular locus

An affine singular point has `y!=0`. Put

```text
r=V/y,                 z=B/y^2,                 c=C/y^3.
```

The equations `H_V=H_y=H=0` reduce exactly to

```text
2r^2+70r-49=0,

z=(10r-77)/540,

c=(7-3r)/972.                                      (6)
```

Hence

```text
r_plus/minus=(-35 plus/minus 21 sqrt(3))/2,

z_plus/minus=-7/15 plus/minus 7 sqrt(3)/36,

c_plus/minus=119/1944 minus/plus 7 sqrt(3)/216.     (7)
```

The corresponding weighted ratios are

```text
rho_plus/minus=C^2/B^3
 =-250/13041 plus/minus (500/117369)sqrt(3).        (8)
```

Conversely each ratio in (8) produces the indicated singularity, after
choosing the compatible square root `y` of `B/z`. Therefore every other
weighted ratio gives a smooth projective plane cubic, hence genus one.
A rational Keller trajectory into that normalization is constant; the
inherited genuine-deck contradiction then excludes it.

At the two exceptional points, the affine Hessian determinant is

```text
1166886 y^2(5 minus/plus 4sqrt(3)).                 (9)
```

Its norm is nonzero because

```text
(5-4sqrt(3))(5+4sqrt(3))=-23.
```

Each singularity is therefore an ordinary node. For a fixed ratio there
is exactly one compatible affine singularity, and infinity is smooth. A
reducible plane cubic cannot have exactly one transverse node and no other
intersection singularity, so both exceptional cubics are irreducible
nodal cubics with rational normalization.

The ratios genuinely lie on `S!=0`. On `R=0`,

```text
S=B^5(2888+244944rho),

S/B^5=(-41576 plus/minus 24000sqrt(3))/23,          (10)
```

whose norm is `24512/23`.

## 4. The retained square sidecar raises genus

Scale a nodal curve so its singular point has `y=1`, and write the line
through the node as

```text
V=r_plus/minus+lambda(y-1).
```

Define

```text
Q2_plus/minus(lambda)
 =3r_plus/minus lambda^2+1470lambda
   +(1470r_plus/minus-22638)/2,

Q3(lambda)=lambda^3+735lambda-3773.                 (11)
```

Substitution into (2) gives

```text
H=(y-1)^2[Q2_plus/minus+Q3(y-1)],
```

and hence the normalization

```text
y=1-Q2_plus/minus/Q3,

V=r_plus/minus-lambda Q2_plus/minus/Q3.             (12)
```

The first flux retains

```text
u=(35y^2+1080B-4yV)/1701,

Z=-2N_2/(5103y)=T^2.                               (13)
```

After (12),

```text
Z_plus/minus(lambda)
 =(8/729) N_plus/minus(lambda)/Q3(lambda)^3,        (14)
```

where `N_minus` is the `sqrt(3)`-conjugate of `N_plus` and `N_plus` has
degree nine. Exact reduction modulo five in

```text
F_25=F_5[s]/(s^2-3)
```

gives the low-to-high coefficient list

```text
(4+s), 3s, (1+2s), (2+4s), s,
(1+3s), (2+2s), 2s, 4s, (3+2s).                    (15)
```

Also

```text
Q3bar=lambda^3+2.
```

The exact Euclidean certificates are

```text
gcd(Nplusbar,Nplusbar')=1,

gcd(Nplusbar,Q3bar)=1,

gcd(Q3bar,Q3bar')=1.                                (16)
```

Consequently `N_plus/minus` has nine simple roots, `Q3` has three simple
roots, and the two sets are disjoint. Equation (14) has nine simple zeros
and three odd-order poles. Its value at infinity is finite and nonzero
because numerator and denominator both have degree nine and

```text
-22 plus/minus 12sqrt(3)!=0.                        (17)
```

The connected double cover

```text
tau^2=Z_plus/minus(lambda)                          (18)
```

therefore has twelve branch points. Riemann--Hurwitz gives

```text
2g-2=2(-2)+12=8,

g=5.                                                (19)
```

A Keller trajectory on the nodal spectral cubic supplies both
`lambda in C(x)` and the retained square root `tau=T in C(x)`, hence a
rational map from `P^1` to (18). It is constant by (19). Thus
`lambda,y,V,u,Z,T` are constant, and the genuine-deck contradiction closes
both nodal ratios.

## 5. The corrected tangent boundary

The original central tangent cubic is

```text
V^3+5670BV-1285956C.
```

It has repeated roots at

```text
C^2/B^3=-250/15309,
```

where nevertheless

```text
S/B^5=-1112!=0.                                     (20)
```

So `S!=0 -> ordinary triple point` is false. This ratio is not one of
(8), and (5) shows that its strict transform is smooth. It belongs to the
generic genus-one closure, not to a missing singular case.

## 6. Consequence and remaining wall

Within the exact scope of THM-2262/2297:

```text
D=25B^2/126, R=0, S!=0

has no survivor.                                    (21)
```

The independently hostile-audited companion
`jc2-degree18-central-wall-s0-closure-opus-20260725.md` closes the `S=0`
face. Hence the whole `R=0` wall is empty and every remaining central-wall
survivor must satisfy

```text
R!=0,                  S!=0,                  T=0,  (22)
```

where the final `T` is the named primitive weight-thirty repeated-branch
factor, not the first-flux square root in (13).

This result closes neither (22) nor the other degree-eighteen discriminant
components, split deck, even-leading descent, other Newton edges, JC(2),
or DC(2).

## 7. Information ledger

```text
source:
  the central repeated-branch factor R=0 and the first-flux square Z=T^2;

map:
  blow up the central spectral point, classify singular plane-cubic
  ratios, then retain the square sidecar over the nodal normalizations;

preserved:
  weighted ratio, normalization component, first flux, genuine deck, and
  the rational Keller trajectory;

destroyed by the blow-up alone:
  the square root T and its branch divisor;

restoring sidecar:
  the genus-five double cover tau^2=Z;

hostile correction:
  S!=0 does not imply distinct central tangent lines;

next decisive test:
  normalize the remaining primitive T=0 factor while retaining both the
  first-flux square and Keller one-form.                              (23)
```

## 8. Exact reproduction

The theorem-neutral standard-library companion is

```text
04-computation/jc2_degree18_central_r0_closure_probe.py
```

with stored transcript

```text
05-knowledge/results/jc2_degree18_central_r0_closure_probe.out.
```

Its LF byte hashes are

```text
script:
  23c7731ab0eb3497e25d35578f755a4489c27d7087733d4be09cf875b524cdaa;

output:
  8aebf3512d652979395b12f82c3fc42f0c8f1b28527fc34daf33f27825cd9ffe.
```

Normal and optimized runs are byte-identical to the stored output. The
companion checks (2)--(17), including the blow-up identity, infinity
discriminant, both singular ratios, Hessians, nonzero `S`, corrected
tangent warning, line-through-node parameterization, exact degree-nine
numerator, conjugation, all three finite-field gcd certificates, twelve
branch points, and genus five. The rational-curve and genuine-deck
conclusions are the mathematical proof above rather than computer
assumptions.
