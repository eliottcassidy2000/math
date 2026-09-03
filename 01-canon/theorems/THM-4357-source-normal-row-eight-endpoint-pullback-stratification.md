---
id: THM-4357
title: "Source-normal row-eight endpoint pullback stratification"
status: >
  PROVED FINITE-ROW COROLLARY RELATIVE TO THM-4308 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On THM-4308's five-dimensional source-parameter
  image, zeta_3=0 makes U,W,Z affine functions of xi_10 with three distinct
  rational roots. The late gates THM-4327-U/Z, 4334, 4337, 4339, and 4340
  have nonempty pullbacks of dimensions 4,4,3,3,2,4, while THM-4342, 4344,
  4350, 4351, and 4353--4356 have empty pullback. The smallest surviving
  slice among the named gates is the squarefree affine plane S_4339. These
  are finite row-eight
  statements only; no all-row lift, seam entry, Keller pair, gate extinction,
  JC(2), or DC(2) is asserted.
source: root + jc-row8-endpoint-pullback scout + independent gate auditor / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
related:
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
  - THM-4340-u-zero-repeated-cusp-planar-jacobian-extinction
  - THM-4356-double-zero-endpoint-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-539
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_row8_endpoint_pullback_stratification_thm4357.py
primary_output: 05-knowledge/results/jc2_source_normal_row8_endpoint_pullback_stratification_thm4357.out
primary_script_sha256: 0ebd463abcfd1cf25e1728da9d78faaf410001f5b6c47e39341102afabbfde87
primary_output_sha256: a139988460b2dc339210ed9d146c03b7e0d4eb5ffbbd7bec49e3de3cadd185c2
referee_script: 04-computation/jc2_source_normal_row8_endpoint_pullback_stratification_independent_referee_thm4357.py
referee_output: 05-knowledge/results/jc2_source_normal_row8_endpoint_pullback_stratification_independent_referee_thm4357.out
referee_script_sha256: af9b9c62b2ea7d44ae78648ca9aaeaab4a1c890057da2382957966d4da3c295c
referee_output_sha256: 3a9d6b1757d3cd3ef861ad0616da767fb557b972689c81e076863f16d269a92b
hash_basis: raw LF bytes
audit: >
  PASS. The 61-check rational primary and independently written 61-check
  symbolic referee reconstruct THM-4308's response, factor the zeta-zero
  trident, prove two unit ideals by explicit Bezout identities and Groebner
  bases, verify every equality and strict owner inequality in the late gate
  ledger, check exact witnesses, Jacobian ranks, source and full-fibre
  dimensions, and the S_4339 cubic discriminant. Normal and optimized runs
  byte-match both frozen outputs.
---

# THM-4357 -- Source-normal row-eight endpoint pullback stratification

**PROVED FINITE-ROW COROLLARY RELATIVE TO THM-4308 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS IS NOT AN ALL-ROW LIFT, SEAM-ENTRY, KELLER-PAIR,
GATE-EXTINCTION, `JC(2)`, OR `DC(2)` THEOREM.**

## 1. Source response and scope

Use THM-4308's source coordinates

```text
(Phi,eta,xi,alpha,beta) in A^5,       xi=xi_10,
alpha=alpha_11,                       beta=beta_11.       (1)
```

Its exact projected row-eight depth conditions fix

```text
Delta=896/15,       Theta=512/75,       K=-32/5,
upsilon_5=-731648/2025,  zeta_3=-3Phi/2,                 (2)
```

In quoted THM-4327 conditions below, `u:=upsilon_5`, `a:=alpha`,
`d:=Delta`, and `e=-1376/135`; this `u` is not THM-4308's source-chart
coordinate `u=x^2t`.

and give the response

```text
U=(475515904-109350xi)/200475,

W=-(4343625Phi^2-17172000xi+143826305024)/4009500,

Z=(12506118074368-173745000Phi^2-195463125Phi eta
   -926883000xi)/108256500.                              (3)
```

Every dimension below is a dimension inside the five-dimensional parameter
space `(1)`. Over each such point, THM-4308 has a separate affine
seven-dimensional terminal tangent fibre. A nonempty pullback below means an
exact finite row-eight jet, not an infinite formal solution or Keller pair.

The inheritance pass is:

- closest proved mechanism: THM-4308's three-parameter top response and exact
  seven-dimensional terminal fibre;
- canonical hostile: the THM-4339 gate is globally extinct but has a whole
  affine plane of row-eight finite jets;
- corrected near miss: extinction of a target gate does not make its finite
  source pullback empty, and an empty finite pullback does not prove general
  seam entry;
- least-used sidecar: the source-to-gate map before quotienting to a Boolean
  owner mask.

The live board was

```text
finite source image | coefficient wall | strict owner inequality
unit ideal | fibre dimension | squarefree carrier | all-row firewall.        (4)
```

## 2. The zeta-zero endpoint trident

Since the characteristic is zero, `zeta_3=0` is equivalent to `Phi=0`.
On that hyperplane, `(3)` factors as

```text
U=-(6/11)(xi-xi_U),
W= (424/99)(xi-xi_W),
Z=-(22886/2673)(xi-xi_Z),                               (5)

xi_U=237757952/54675,
xi_W=4494572032/536625,
xi_Z=1563264759296/115860375.                            (6)
```

The three pairwise differences are

```text
xi_U-xi_W=-58347587584/14488875,
xi_U-xi_Z=-28604827277312/3128230125,
xi_W-xi_Z=-3491293831168/682288875,                      (7)
```

so the roots are distinct. At the three roots the surviving pairs are

```text
xi=xi_U: W=-42434609152/2460375,
         Z= 5200877686784/66430125;

xi=xi_W: U=-10608652288/4829625,
         Z= 634780696576/14488875;

xi=xi_Z: U=-5200877686784/1042743375,
         W= 2539122786304/115860375.                     (8)
```

In particular, no two of `U,W,Z` vanish on `zeta_3=0`. Two useful stronger
nonfaces are the unit ideals

```text
<zeta_3,xi,W>=(1),              <Z,zeta_3,W>=(1).         (9)
```

The independent referee supplies the explicit Bezout certificates

```text
1=[W-(13Phi/18)zeta_3-(424/99)xi]
  /[-35956576256/1002375],

1=[Z+(11443/5724)W-13(6641Phi+3180eta)zeta_3/34344]
  /[634780696576/14488875].                              (10)
```

The first three-wall nonface is minimal: `zeta_3=xi=0` and
`zeta_3=W=0` are individually nonempty, and `W=xi=0` is nonempty over the
algebraically closed field at

```text
Phi^2=-143826305024/4343625 !=0.                         (11)
```

Thus independent Boolean switching of these owner walls strictly overcounts
the realizable row-eight supports.

## 3. Complete late-gate pullback ledger

The following table includes every equality and strict inequality used from
the cited gate. `dim` again excludes the affine seven-dimensional terminal
tangent fibre.

| gate | complete row-eight pullback | status/dim |
|---|---|---:|
| THM-4327-U: `U=0,WZ!=0`; `u!=0 => alpha^2-4Wu!=0`; `u=alpha=Delta=0 => eta^2-4eW!=0` | `xi=xi_U`; delete `WZ=0` and `alpha^2-4Wu=0`; the second antecedent is impossible because `u,Delta!=0` | nonempty, 4 |
| THM-4327-Z: `Z=0,U*beta*zeta_3!=0` | `Phi!=0,xi!=xi_U,beta!=0`, `eta=(12506118074368-173745000Phi^2-926883000xi)/(195463125Phi)` | nonempty, 4 |
| THM-4334: `Z=beta=0,U*W*zeta_3!=0` | the preceding graph, with `beta=0` and `W!=0` | nonempty, 3 |
| THM-4337: `Z=zeta_3=0,U*beta!=0` | `Phi=0,xi=xi_Z,beta!=0`; `eta,alpha` free | nonempty, 3 |
| THM-4339: `Z=beta=zeta_3=0,U*K*W*(U+W)!=0` | `Phi=0,xi=xi_Z,beta=0`; `eta,alpha` free; all inequalities automatic | nonempty, 2 |
| THM-4340: `U=0,WZ!=0` | `xi=xi_U`; delete `WZ=0` | nonempty, 4 |
| THM-4342: `Z=beta=zeta_3=K=0,U*W*(U+W)!=0` | `K=-32/5` | empty |
| THM-4344: `Z=beta=zeta_3=W=0,U*K*xi!=0` | `<Z,zeta_3,W>=(1)` | empty |
| THM-4350: `Z=beta=zeta_3=W=xi=0,U*K!=0,(alpha,Theta)!=(0,0)` | `<zeta_3,xi,W>=(1)` | empty |
| THM-4351: `Z=beta=zeta_3=W=xi=alpha=Theta=0,U*K!=0` | already `Theta=512/75!=0` | empty |
| THM-4353: `Z=beta=zeta_3=W=xi=K=0,U!=0` | `K=-32/5` | empty |
| THM-4354: `Z=beta=zeta_3=W=xi=U=0,K*alpha!=0` | `<zeta_3,xi,W>=(1)` | empty |
| THM-4355: `Z=beta=zeta_3=W=xi=U=alpha=0,K!=0` | `<zeta_3,xi,W>=(1)` | empty |
| THM-4356: `Z=beta=zeta_3=W=xi=U=K=0,Delta=5696/105` | `K=-32/5` and `Delta=896/15` | empty |

THM-4352 is a formal-local cusp law, not a coefficient gate, so it has no
row-eight pullback status.

The source-parameter dimensions in the table become respectively
`11,11,10,10,9,11` after adjoining the full terminal tangent fibre. None of
those dimensions is a claim about the moduli of actual Keller maps.

## 4. Exact nonempty controls and dimension proof

For THM-4327-U and THM-4340 take

```text
(Phi,eta,xi,alpha,beta)=(0,0,xi_U,0,0).                  (12)
```

Then `W,Z,upsilon_5` are nonzero and
`alpha^2-4W upsilon_5=-4W upsilon_5!=0`. For THM-4327-Z take

```text
Phi=1, xi=0, beta=1, alpha=0,
eta=12505944329368/195463125.                            (13)
```

This has `Z=0` and every required owner nonzero; changing `beta` to zero
gives the THM-4334 control. The THM-4337 and THM-4339 controls are

```text
(Phi,eta,xi,alpha,beta)=(0,0,xi_Z,0,1),
(Phi,eta,xi,alpha,beta)=(0,0,xi_Z,0,0),                  (14)
```

respectively.

For dimension, `U=0` is one independent linear equation. On `Phi!=0`,
`Z=0` is an irreducible graph because

```text
dZ/deta=-195463125Phi/108256500 !=0.                     (15)
```

Adding `beta=0` lowers its dimension by one. On `zeta_3=0`, first `Phi=0`,
and then `Z=0` is the independent equation `xi=xi_Z`; adding `beta=0` again
lowers dimension by one. The referee checks Jacobian ranks `1,1,2,2,3` and
all strict inequalities at `(12)--(14)`.

## 5. The smallest survivor `S_4339`

The complete two-dimensional pullback of THM-4339 is

```text
S_4339:
Phi=0,      xi=1563264759296/115860375,      beta=0,
eta,alpha arbitrary.                                      (16)
```

Uniformly on this affine plane,

```text
U=-5200877686784/1042743375,
W= 2539122786304/115860375,
K=-32/5,
U+W=17651227389952/1042743375.                           (17)
```

All four quantities are nonzero. Moreover, the cubic

```text
A(P)=K+Theta P+xi P^2+W P^3                              (18)
```

has the fixed discriminant

```text
483125535259306642688385993911761371136
----------------------------------------------------- !=0. (19)
7776331997934642451171875
```

Thus the whole plane lies in THM-4339's squarefree cubic stratum. It is the
unique smallest nonempty pullback among the named gates in the ledger. The exact next
higher-row test is therefore the row-nine bracket/depth condition restricted
to `(16)`, with only `(eta,alpha)` left as source coefficients. Quoting
THM-4339 cannot perform that test: its extinction conclusion is conditional
on an actual exact-seam Keller lift, whereas `(16)` consists only of finite
jets.

## 6. Audit and boundary of the result

The rational primary and independent symbolic referee each perform `61`
checks. They separately reconstruct `(2)--(3)`, the three roots and
factorizations, both unit ideals, every gate equality and inequality, the
positive controls, dimension ranks, full-fibre dimensions, and the cubic
discriminant. The referee imports no primary code and verifies `(9)` by both
the explicit identities `(10)` and reduced Groebner bases `[1]`. Normal and
optimized executions byte-match both frozen outputs in the frontmatter.

This theorem classifies intersections inside one finite projected response
image. It neither proves nor assumes that those jets extend to an all-row
`B_2` solution, terminate polynomially, arise from a Keller pair, or receive
an arbitrary hypothetical counterexample. Exact-seam entry, `JC(2)`, and
`DC(2)` remain open.
