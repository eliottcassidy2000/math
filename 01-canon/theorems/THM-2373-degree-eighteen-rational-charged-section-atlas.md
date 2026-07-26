---
id: THM-2373
title: "Degree-eighteen rational charged-section atlas"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. After the
  proved closure of every support of size at most two, each live
  degree-eighteen spectral survivor lies in at least one of the three
  adjacent charts BC!=0, CD!=0, or DW!=0. The Laurent ratios C/B, D/C,
  and W/D have weighted charge one. On their charts they therefore give
  unique rational, root-free scaling sections B'=C', C'=D', and D'=W',
  with explicit formulas and weight-zero transition functions. The
  sections preserve the weighted spectral equation and repeated-branch
  locus when the scale is retained. They are alternative gauges to
  THM-2357's moving-root normalization, not simultaneous extra
  normalizations, and do not by themselves preserve a forgotten Keller
  one-form or close H2, H4, JC(2), or DC(2).
source: codex-2026-07-25-jc-charged-section-atlas
depends_on:
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2328-degree-eighteen-bw-ratio-bank-closure
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
related:
  - THM-2357-degree-eighteen-h2-moving-root-reduction
  - THM-2371-degree-eighteen-h2-common-root-elimination
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
script: 04-computation/degree18_rational_charged_section_atlas_thm2373.py
output: 05-knowledge/results/degree18_rational_charged_section_atlas_thm2373.out
script_sha256: 25e08ca7783bf91e511e95532d294db402193959e05cdf463a018a6e28099310
output_sha256: 72a65082658fc4a1e6b40a8ba4e1196197f7ae8078a1b582dd4419a2deb36f66
hash_basis: working-tree bytes (LF)
---

# THM-2373 -- a root-free scaling atlas on the live degree-eighteen cone

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2297 exposes the weighted spectral cone

```text
G_0(u,y;B,C,D,W)=0,

wt(y,u,B,C,D,W)=(1,2,2,3,4,5).                     (1)
```

The usual weighted-projective quotient hides a practical field issue:
normalizing a coordinate of weight two, three, four, or five may require
extracting a root. The live higher-support locus has a better atlas. Three
adjacent Laurent ratios have charge exactly one, so their normalization
is rational over the original field.

## 1. The weighted action and the live support

Write the scaling action as

```text
rho.(y,u;B,C,D,W)
 =(rho y,rho^2 u;
   rho^2 B,rho^3 C,rho^4 D,rho^5 W).               (2)
```

THM-2297 proves

```text
G_0(rho^2u,rho y;
    rho^2B,rho^3C,rho^4D,rho^5W)
 =rho^6G_0(u,y;B,C,D,W).                           (3)
```

THM-2311 through THM-2328 close every support of size at most two among
`(B,C,D,W)`. THM-2332 therefore restricts every live repeated-branch
survivor to support size three or four.

Every three-element subset of the path

```text
B -- C -- D -- W                                   (4)
```

contains an adjacent pair. Hence the live locus is covered by

```text
U_BC={BC!=0},

U_CD={CD!=0},

U_DW={DW!=0}.                                      (5)
```

This is a cover of the current spectral survivor locus, not of the whole
weighted cone.

## 2. Three weight-one covariants

On the charts (5), define

```text
tau_BC=C/B,

tau_CD=D/C,

tau_DW=W/D.                                       (6)
```

Under (2),

```text
tau(rho.x)=rho tau(x).                              (7)
```

Thus each `tau` is a nonvanishing weight-one Laurent covariant. The
equation

```text
tau(rho.x)=1
```

has the unique solution

```text
rho=tau(x)^(-1).                                   (8)
```

No algebraic root or field extension is used.

## 3. Explicit rational sections

On `U_BC`, choose `rho=B/C`. The section is

```text
y'=By/C,

u'=B^2u/C^2,

B'=C'=B^3/C^2,

D'=B^4D/C^4,

W'=B^5W/C^5.                                      (9)
```

On `U_CD`, choose `rho=C/D`:

```text
y'=Cy/D,

u'=C^2u/D^2,

B'=C^2B/D^2,

C'=D'=C^4/D^3,

W'=C^5W/D^5.                                      (10)
```

On `U_DW`, choose `rho=D/W`:

```text
y'=Dy/W,

u'=D^2u/W^2,

B'=D^2B/W^2,

C'=D^3C/W^3,

D'=W'=D^5/W^4.                                    (11)
```

Equations (7)--(8) prove uniqueness **inside each chosen chart**. On an
overlap, two charts generally select different representatives of the
same weighted orbit.

By (3), each section preserves the zero locus of `G_0`; the repeated-
branch discriminant locus is preserved by the corresponding rescaling
of `y`. Recording the chart and `rho` reconstructs the original point
exactly.

## 4. Weight-zero transition functions

The ratios of the three charged covariants are invariant:

```text
I_BCD=tau_CD/tau_BC=BD/C^2,

I_CDW=tau_DW/tau_CD=CW/D^2,

I_BDW=tau_DW/tau_BC=BW/(CD).                      (12)
```

They satisfy

```text
I_BCD I_CDW=I_BDW.                                 (13)
```

If `s_i` and `s_j` are two section representatives, the scale carrying
`s_i(x)` to `s_j(x)` is

```text
rho_j/rho_i=tau_i/tau_j,                           (14)
```

the inverse of the corresponding invariant in (12). Thus the three
charts glue by rational weight-zero transitions.

## 5. What the atlas changes in the open reductions

The atlas replaces the root-bearing weighted-projective instruction

```text
"make one weighted coordinate equal to one"
```

by one of the root-free affine equalities

```text
B=C,                 C=D,                 D=W.     (15)
```

This is immediately suitable for exact elimination over `Q` or a finite
field whenever the relevant chart denominators are inverted. In
particular, the still-open `H_4 S_4^2` coefficient system may be split
into the three explicit lanes (15), with overlap consistency controlled
by (12).

For the `H_2` lane, THM-2357 already spends the same one-dimensional
scaling action to move its distinguished branch root to `y=1`. Equations
(9)--(11) are alternative gauges:

```text
moving root y=r to 1

or

making one adjacent parameter pair equal.          (16)
```

They cannot in general be imposed simultaneously. One may retain the
root as a variable and use this atlas, or use THM-2357's moving-root
system; the theorem supplies no second scaling degree of freedom.

The atlas preserves the spectral curve only together with its scale
sidecar. THM-2297 explicitly warns that discarding that sidecar is not a
target-preserving quotient of the Keller one-form, the third flux, or
the whole-polynomial Faber data. No emptiness or Jacobian consequence
follows from (15) alone.

## 6. Sharp boundary: the two-sparse even chart

On the `B,D`-only support, every Laurent monomial has weight

```text
2a+4b,
```

which can never equal one for integers `a,b`. The element `rho=-1`
fixes both nonzero coordinates, giving a genuine `mu_2` stabilizer. None
of the three adjacent charts (5) is live there.

This is the sharp obstruction to a global weight-one section. It does
not affect the current survivor locus because the entire two-sparse bank,
including the `B,D` lane, is already closed by the proved degree-eighteen
ratio theorems.

## 7. Consequence and scope

Every live degree-eighteen spectral survivor now has:

```text
a chart in {BC,CD,DW};

a unique rational section representative in that chart;

a rational inverse scale;

weight-zero overlap coordinates.                   (17)
```

The atlas preserves the field of definition, weighted spectral
equations, support labels, and branch multiplicities. It loses nothing
when the chart and scale are retained. If the scale is forgotten, it
loses the absolute Keller normalization and may not be used as a
Jacobian-conjecture quotient.

The theorem does not close the localized `H_2` ideals, the `H_4` stratum,
THM-2371's reserved common-root lane, `JC(2)`, or `DC(2)`.

## 8. Exact companion

The dependency-free `Fraction` companion:

- verifies that all five support-three/four patterns meet a chart;
- checks the full displayed `G_0` covariance on `73` rational samples;
- checks `219` chart normalizations and `219` transitions;
- checks `219` weight-zero invariants and their cocycle identity;
- verifies the closed formulas (9)--(11); and
- checks the `B,D` `mu_2`/no-weight-one hostile.

Run

```bash
python3 04-computation/degree18_rational_charged_section_atlas_thm2373.py
python3 -O 04-computation/degree18_rational_charged_section_atlas_thm2373.py
```

Both transcripts must match

```text
05-knowledge/results/degree18_rational_charged_section_atlas_thm2373.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit of the theorem text is pending. QED.
