---
id: THM-3249
title: "Cross-support upset atlas local sections and no constant gauge"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  The fixed
  twenty-two lawful upset templates of THM-3238 have strictly positive
  unique-reset sections on both complete support-(1,2), bank-I2 and
  support-(1,3), bank-I2 physical selector faces.  On the 239-state small
  face, six rows already suffice and an explicit positive perturbation uses
  all twenty-two.  Nevertheless an exact positive nineteen-state Farkas
  mixture has strictly positive expectation in every chart coordinate, so
  no one nonzero nonnegative chart weighting can expose both resets at once.
  This is a two-face transition-gauge obstruction, not arbitrary-support
  product-Gamma positivity or the Gaussian Moment Conjecture.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-03
audit: >
  The cache-free exact companion reconstructs all 239 small-face states and
  the sixteen additional full-face witness states from the pinned
  product-Gamma formula.  It verifies the six-row primitive certificate, the
  all-twenty-two positive perturbation, the failure of the transported
  THM-3238 weights, and every coordinate of the nineteen-state Farkas
  mixture using rational arithmetic.  Normal and optimized runs byte-match
  the stored transcript; independent proof and replay audit are pending.
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
related:
  - THM-3222-universal-product-gamma-reset-upper-filter-collar
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
script: 04-computation/gmc_cross_support_upset_atlas_holonomy_thm3249.py
output: 05-knowledge/results/gmc_cross_support_upset_atlas_holonomy_thm3249.out
script_sha256: 831f9086cfa34b85e63eb7095243ae5de398eed626f24b2feadd91fff25b83b3
output_sha256: bab45455dff1080caf47303615abbfd32f538b32b165cd67ad6373e0aeaaf816
hash_basis: LF-normalized bytes
---

# THM-3249 -- cross-support upset atlas local sections and no constant gauge

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-3238 closes one complete physical product-Gamma face with twenty-two
lawful coarsening-upset charts.  Transporting an atlas and transporting one
fixed positive weighting are different questions.  The present theorem
separates them exactly: the same charts work on the neighboring support
`(1,2)`, but their positive cones admit no constant common section across the
two faces.

## 1. Two faces and one chart atlas

For `f in {12,13}`, let `r_i^f(sigma)`, `1<=i<=22`, be the response of the
`i`-th degree/upset template in the table of THM-3238, evaluated in bank `I2`
at support `(1,2)` or `(1,3)` respectively.  Every `r_i^f` is lawful by
THM-3127.  On the small face the reduced pole alphabet and reset are

```text
P_12=(5,4,3,3,2,2,2,1,1,1,1),
Q_12=(1,2,2,3,4,5).                                  (1)
```

The complete nonempty submultiset bank has

```text
|S(P_12)|=(4+1)(3+1)(2+1)(1+1)(1+1)-1=239,          (2)
```

with depth census

```text
(5,13,24,35,42,42,35,24,13,5,1).                    (3)
```

For either face and every chart,

```text
r_i^f(Q_f)=0,                                        (4)
```

because the quotient alphabet at the reset is empty.

## 2. An exact positive local section on the small face

Only rows `5,8,10,15,16,21` are needed initially.  In the table below,
`M_i` is the exact maximum absolute response of row `i` on the 239 states,
`a_i` is a positive rational factor, and `c_i` is the primitive integer
coefficient obtained by clearing the six denominators of `a_i/M_i`.

| row `i` | `M_i` | `a_i` | primitive `c_i` |
|---:|---:|---:|---:|
|5|1538222247466205184|`213915797/685946076`|31540805104781016214807762327098361741908455758236214857697597793166077169575373070006205|
|8|10668668469065531392|`81328689/164203297`|7222551733565666973704472289597930205521464117990467362419953714461056410853911375394360|
|10|28568269398144|`6210895/131921019`|256386550655705646640018746246846057992771457973407640357382724455640369346859861323130803200|
|15|1484288|`255973/944937251`|28393046540266691701145009302241688233214029700961546492881520901040914942725957247847268249416960|
|16|283951084704|`36321374/498402635`|39927984611766798370362293023501771795759571792906966553729869583513560353170793810845831929856|
|21|386978585480040|`55150150/759378773`|29197138787760557148002850343923059397376740901033387386489565737312725443107887357993779200|

Put

```text
H_0^12(sigma)=sum_(i in {5,8,10,15,16,21}) c_i r_i^12(sigma). (5)
```

Exact evaluation gives

```text
H_0^12(Q_12)=0,
H_0^12(sigma)<0                     for every sigma!=Q_12. (6)
```

The negative/zero/positive census is `238/1/0`.  The closest negative state is
the singleton `(3)`, where

```text
H_0^12(3)
=-7485041558321314973110552874859175529199365055985729203025401847758119404861486274681270934476877696.
                                                               (7)
```

This six-row section lies on the boundary of the nonnegative chart orthant.
To obtain an honest positive twenty-two-chart section, keep the six displayed
coefficients and assign coefficient `1` to each of the other sixteen rows:

```text
H_+^12=H_0^12+sum_(i notin {5,8,10,15,16,21}) r_i^12. (8)
```

The exact census remains `238/1/0`; `(3)` is again closest, now with

```text
H_+^12(3)
=-7485041558321314973110552874859175529199365055985729203025401847758119404861486285486714204420184096.
                                                               (9)
```

Thus the open positive section cone on the small face is nonempty.  The same
expectation argument as THM-3238 proves

```text
C_D^(physical)(support (1,2),I2)={delta_(Q_12)}
                                      for every D>=14.          (10)
```

This is a second complete physical face, not merely the principal collar of
THM-3222.

The coefficients themselves do not transport.  Evaluating the primitive
THM-3238 weights on `(1,2),I2` gives

```text
negative/zero/positive=115/1/123,                    (11)
```

and its largest positive state is `(3,3)`.  The failure is therefore a gauge
change, not a failure of the twenty-two-chart atlas.

## 3. Exact nineteen-state common-gauge obstruction

Write `12:sigma` and `13:sigma` for a state in the corresponding face.  The
following positive rational mixture is the full Farkas witness.

| face:state | weight |
|---|---:|
|`12:(2,3)`|`110021555/210299496`|
|`12:(2,2)`|`244638624/855590729`|
|`12:(2,2,2)`|`13822588/86007735`|
|`13:(1)`|`981/345013099`|
|`13:(2,2)`|`740/175674781`|
|`13:(1,1,1,1)`|`1988/947209091`|
|`13:(1,1,1,2)`|`94/29209079`|
|`13:(2,2,2,3,3)`|`2741/436074315`|
|`13:(1,1,1,1,2,2,2)`|`1886/428306147`|
|`13:(4,4,5,5,6,7,8)`|`4040735/402653491`|
|`13:(3,3,4,4,5,5,7,8)`|`1681/222909716`|
|`13:(1,1,1,1,5,5,6,7,8)`|`495494/846338163`|
|`13:(1,1,1,1,2,2,2,3,3,4)`|`3301/978297006`|
|`13:(1,1,1,1,4,4,5,6,7,8)`|`498441/407384384`|
|`13:(2,2,3,3,4,4,5,5,6,7)`|`1/223456666`|
|`13:(2,3,3,4,4,5,5,6,7,8)`|`1338659/314691983`|
|`13:(2,2,2,3,3,4,4,5,5,6,7)`|`17/863033027`|
|`13:(1,1,1,1,2,2,2,5,5,6,7,8)`|`621161/731987841`|
|`13:(1,1,1,1,2,2,2,4,5,5,6,7,8)`|`11562097/875223752`|

Call these weights `w_(f,sigma)`, and for each chart put

```text
E_i=sum_(f,sigma) w_(f,sigma) r_i^f(sigma).          (12)
```

Direct exact rational evaluation gives

```text
E_i>0                              for every 1<=i<=22. (13)
```

The witness-table digest and the ordered expectation-vector digest are

```text
89c1b11710a1c1c0793dfa0b0e8b385b9fb1fb8a265d324b4bde463d436a38ae,
1d26244164b903e6d506d61962a4c49b580a73f0ce99a3d241fd884309616286. (14)
```

Suppose that one nonzero vector `theta in R_>=0^22` exposed both resets
strictly.  Every state in the table is nonreset, so

```text
sum_(f,sigma) w_(f,sigma)
  sum_i theta_i r_i^f(sigma)<0.                       (15)
```

Reversing the finite sums and using `(12)--(13)` makes the same quantity

```text
sum_i theta_i E_i>0,                                  (16)
```

a contradiction.  Hence no nonzero nonnegative weighting of these fixed
twenty-two templates is a strict common section for both faces.

The mixture in the table is a dual certificate for a chart-cone intersection.
It is not asserted to be a selector-feasible probability law on either face.

## 4. Holotopy interpretation and boundary

For a face `f`, define its positive exposure cone

```text
K_f={theta in R_>0^22:
     sum_i theta_i r_i^f(sigma)<0 for every sigma!=Q_f}. (17)
```

THM-3238 and `(8)--(9)` say `K_13` and `K_12` are both nonempty.  Equations
`(12)--(16)` say

```text
K_12 intersection K_13=empty.                         (18)
```

Thus the chart atlas transports but a constant positive gauge does not.  In
holotopy language, there are local positive sections and a nontrivial
transition requirement.  This is not yet a loop or a nonzero Cech class: a
third face and pairwise-compatible transition sections would be required for
that stronger conclusion.

THM-2267 and THM-2292 are useful abstract precedents for static charts whose
transition gluing carries extra information.  Their owner-cut and signed
gain-graph invariants are different types; neither proves `(18)`.  Here the
obstruction is exactly convex and is witnessed by the positive Farkas mixture.

## 5. Scope

The theorem concerns only bank `I2`, supports `(1,2)` and `(1,3)`, their
complete physical submultiset banks, and the twenty-two THM-3238 templates.
It proves neither that the atlas works on support `(1,4)` nor that some larger
atlas lacks a common gauge.  It does not establish arbitrary-support
product-Gamma positivity, arbitrary radial coefficients, NC2, or the Gaussian
Moment Conjecture.  Its positive conclusion is the complete `(1,2),I2` face;
its negative conclusion is precisely the constant-gauge transport wall.

QED.
