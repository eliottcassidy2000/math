---
id: THM-2608
title: "Alternative-rail future-clock marginal, lawful-clock collapse, and the missing transition index"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Replacing THM-2592's priority selector by the root-six theta-zero rail
  whenever available makes 1,012 of 1,092 target sections positive on 79 of
  84 rails.  After the future owner clock is marginalized, every displacement
  s != 6 gives seven dense 13 by 13 root cospans whose formal product has
  support 156 and twelve returns at every root displacement.  This prism is
  not composable: for every s and every nonzero C7-equivariant clock step,
  the seven clock-matched product is exactly zero (72/72 cases).  THM-2601's
  scalar coordinate gives exactly thirteen equivariant deep-root-to-target
  reference maps, but selects none.  Retaining all thirteen references and
  all future clocks produces a gauge-invariant positive 1,183-state formal
  return (612 paths per nonzero twist for s != 6), yet still not one physical
  seven-edge ancestry carrier.  The missing datum is a same-carrier next
  target index q' (or physical r-to-q' map) together with adjacent-clock
  support.  No scalar row is excluded; LRC(14) remains open.
source: deep-energy-audit-2026-07-28-clock-matched-cospan
depends_on:
  - THM-2592-fallback-rail-digit-diagonal-pullback-and-primitive-bockstein
  - THM-2601-linear-bockstein-sheet-coordinate-and-nonlinear-target-monodromy
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
related:
  - THM-2586-depth-five-arrival-to-future-root-diagonal
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2594-realized-theta-slaved-contraction-at-the-r5-window
script: 04-computation/lrc14_alternative_rail_clock_collapse_thm2608.py
output: 05-knowledge/results/lrc14_alternative_rail_clock_collapse_thm2608.out
script_sha256: ef50a685b1a1b3fa9d823d5336afbde0a58f53630b237f4a93be7a55062c37af
output_sha256: fb3d4793d755dd8152af2cca37dcceb3c2181b7c8f19f75d4647d7a1f2258052
hash_basis: LF-normalized bytes
---

# THM-2608 -- the alternative rail fills the marginal and misses the clock

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2592 found a positive common-`x` attachment only on three fallback
theta-zero rails.  Its final audit also located a sharp boundary: the other
theta-zero rail `(v,t)=(6,12)` is available on almost the entire bank, and
changing to it changes physical support rather than merely reselecting a
future digit.  This theorem performs that alternative selection and retains
the future owner clock before testing transition composition.

The answer has two opposite faces:

```text
future-clock marginal:       dense formal seven-edge root prism;
lawful C7 clock matching:    every seven-edge product is zero.            (1)
```

This is not a cancellation phenomenon.  The tensors are nonnegative and the
zero products contain literal zero edges.  The apparent prism is created by
forgetting which future clock supplied each positive root edge.

## 1. The alternative common-x tensor

Use THM-2592's canonical typed row, physical coordinate, clocks, present
packet, deep probe, and delayed word without modification.  Thus its exact
common-carrier tensor is

```text
J_(s,ell,q;ell',r)
 =13 int U(x)V(x-s/13) 1_(ell,v,t)(x)
       F_(ell',-q)(x) Delta_r(x) Q^+_(ell',v)(13^6 x) dx.  (2)
```

The indices have distinct types:

```text
s       : nonzero route-two displacement;
ell     : THM-2586/current owner clock;
q       : THM-2585 target-shift section s_5=-q;
ell'    : THM-2585/future owner clock;
r       : deepest-root translate.                          (3)
```

Choose the theta-zero rail in the order

```text
(v,t)=(6,12), then (0,0),                                 (4)
```

using the second only when the first is unavailable.  This is the reverse of
THM-2592's proved priority choice.  It is a new physical subcarrier, not an
operation on THM-2592's old 39 positive slices.

The exact integer rebuild gives

```text
positive (s,ell,q) sections:       1012/1092;
positive (s,ell) rails:              79/84;
global nonzero numerator content: 4,244,240.                (5)
```

All statements below use only positivity of these exact integer entries, so
one global content reduction or a common positive rescaling leaves them
unchanged.

## 2. The dense future-clock marginal

First forget `ell'` and define the nonnegative root cospan

```text
M_(s,ell)(q,r)=sum_(ell' in F_7) J_(s,ell,q;ell',r).       (6)
```

Write `supp M` for its Boolean support and multiply supports as relations.
For each

```text
s in F_13^* \ {6},                                       (7)
```

the seven edge-support sizes in the order `ell=0,...,6` are

```text
(156,156,144,156,156,156,156).                            (8)
```

Their formal product has support `156/169`.  More strongly, for every
`delta in F_13`, exactly twelve starting roots have an endpoint displaced by
`delta`:

```text
#{q : (q,q+delta) lies in supp(M_(s,0)...M_(s,6))}=12.    (9)
```

At `s=6`, the edge supports are

```text
(132,132,0,0,0,0,0),                                    (10)
```

and the product vanishes.

Equations (8)--(9) are an exact positive control for THM-2602's twisted
return criterion.  They are not yet a transition theorem: the summation in
(6) lets a different future clock `ell'` witness each edge and has erased the
clock needed to glue consecutive edges.

## 3. THM-2601 supplies a thirteen-choice reference torsor

THM-2601 gives the separating scalar values

```text
t_q=(9,3,10,1,8,0,7,2,12,11,4,6,5),                    (11)
```

and the transported target successor `S(t_q)=t_(q+1)`, whose value table is

```text
S=(7,8,12,10,6,9,5,2,0,3,1,4,11).                      (12)
```

The deep-root successor is `r -> r+1`.  Therefore an equivariant map from
the deep-root alphabet into the scalar target alphabet must satisfy

```text
phi(r+1)=S(phi(r)).                                      (13)
```

For each `c in F_13` there is exactly one solution,

```text
phi_c(r)=S^r(c).                                         (14)
```

Conversely every solution is determined by `c=phi(0)`.  Hence there are
exactly thirteen equivariant identifications.  If `P` is THM-2601's inverse
map `P(t_q)=q` and `q_0=P(c)`, then

```text
phi_c(r)=t_(q_0+r).                                      (15)
```

Thus feeding the output root into another target section amounts, in the
original `q` labels, to the reference choice

```text
q_next=q_0+r.                                            (16)
```

THM-2601 deliberately calls its separating functional noncanonical and does
not identify the target atlas with a physical deep-root deck.  Equations
(14)--(16) classify the missing choice; they do not make one choice physical.
The thirteen `(c,q_0)` pairs are

```text
(0,5),(1,3),(2,7),(3,1),(4,10),(5,12),(6,11),
(7,6),(8,4),(9,0),(10,2),(11,9),(12,8).                 (17)
```

For every one of these references, (8)--(9) remain true for all eleven good
displacements.  There are therefore `11*13=143` dense formal marginal
prisms and thirteen zero ones at `s=6`.  A choice of root gauge is not the
reason the lawful construction fails.

## 4. Every C7-equivariant clock matching collapses

A nonzero equivariant step on the regular owner-clock torsor has the form

```text
ell' = ell+d,                  d in F_7^*.                (18)
```

Traverse its clock cycle in the order

```text
0,d,2d,...,6d.                                            (19)
```

Before marginalizing `ell'`, define the clock-matched root cospan

```text
M^(d)_(s,ell)(q,r)=J_(s,ell,q;ell+d,r).                  (20)
```

For every `s in F_13^*` and every `d in F_7^*`, the sevenfold Boolean
product in the order (19) is zero.  This is all

```text
12*6=72                                                    (21)
```

lawful nonzero clock steps.

For each `s!=6`, the six edge-support signatures, indexed by `d=1,...,6`,
are

```text
d=1: (156,132,132,  0,  0,  0,144)
d=2: (  0,  0,  0,  0,144,  0,132)
d=3: (132,  0,144,  0,144,  0,132)
d=4: (  0,132,  0,  0,  0,132,132)
d=5: (  0,156,  0,  0,  0,144,  0)
d=6: (  0,  0,  0,144,132,  0,  0).                    (22)
```

At `s=6`, four signatures are wholly zero and the other two have only one
nonzero edge, of size `132`.  Every product therefore contains a zero
factor.  Reindexing `r` through any `phi_c` only permutes columns and cannot
repair (22).

This proves the promised clock-collapse statement.  It does not say that
the separate positive entries of (2) are false.  It says they cannot be
placed around one equivariant seven-clock cycle using the clock labels that
the tensor actually retained.

## 5. Retaining all clocks and all references gives a formal return

The strongest coefficient-level object available without imposing (18) is
the relation on `F_7 x F_13`

```text
C_(s,c)((ell,q),(ell',q'))>0

iff some r has J_(s,ell,q;ell',r)>0 and q'=P(c)+r.        (23)
```

This keeps all future clocks but uses one of the formal reference maps (16).
For every `s!=6` and every `c`, its exact support data are

```text
#supp C_(s,c)       =2760,
#supp C_(s,c)^7     =4320.                               (24)
```

For each nonzero base twist `a`, the number of same-clock paths

```text
(ell,q) -> (ell,q-7a)                                    (25)
```

in the seventh support is `47` or `48`, hence is always positive.  Across
the `11*13*12=1716` good `(s,c,a)` cases the exact histogram is

```text
47:1584,       48:132.                                   (26)
```

At `s=6`, all thirteen relations have support `264` and seventh power zero.

One need not select a reference to get a formal invariant.  Retain `c` as a
thirteenth state coordinate and take the block-diagonal union of all thirteen
relations (23).  This is an `13*7*13=1183`-state coefficient relation.  For
every `s!=6` it has

```text
#edge support          =13*2760=35880,
#seventh support       =13*4320=56160,
twisted return trace   =612                              (27)
```

for every `a!=0`.  Equation (27) is reference-gauge invariant and answers the
natural holotopy control positively: the complete coefficient torsor has
formal twisted returns.

It is still not a physical transition current.  A matrix product in (23)
may choose unrelated common-`x` witnesses at successive edges, and it may
change `ell` by a different amount at every step.  No theorem supplies a
single sevenfold fibre product realizing one of those 47/48 paths.  The
contrast between (22) and (27) is precisely the lost datum.

## 6. The minimum missing index

The exact type of THM-2592's tensor is a span

```text
(ell,q)  <-  common-x edge carrier  ->  (ell',r).         (28)
```

Its two legs have different meanings.  A composable ordered edge requires
one further target-section index on the same physical carrier, for example

```text
L_ell(q,r,q') >= 0,                                      (29)
```

with the following load-bearing properties:

```text
clock:      the output lies at the adjacent retained clock ell+1;
root map:   q' is the literal next target section, not a chosen relabelling;
ancestry:   consecutive L_ell share the intermediate physical state;
positivity: a sevenfold path is one nonnegative fibre product;
twist:      some path has q_7=q_0-7a.                    (30)
```

Equivalently, one may supply a physical map `iota_ell:r -> q'` satisfying
(30).  Summing `r` in (29) then gives THM-2602's ordered kernel
`K_ell(q,q')`.

Neither a scalar `t_q` nor the numerical equality of two copies of `F_13`
supplies `iota`.  Without the adjacent-clock clause, (23) gives a false
positive; without the next-target index, the thirteen maps (14) remain an
unresolved reference torsor.  These are independent obstructions.

## 7. Verification and exact scope

Run

```text
python 04-computation/lrc14_alternative_rail_clock_collapse_thm2608.py
python -O 04-computation/lrc14_alternative_rail_clock_collapse_thm2608.py
```

The companion rebuilds the alternative tensor from THM-2592's exact
integer/Fraction routines and checks:

1. all 84 rail masses and the census/content (5);
2. all 1,092 marginal sections and the support/return laws (8)--(10);
3. all thirteen equivariant maps (14), including successor covariance;
4. all 72 clock-matched products and the full signature atlas (22);
5. all 156 formal 91-state relations and their seventh powers; and
6. all 1,872 twisted traces, including the 1,183-state sums (27).

Normal and optimized runs byte-match the stored transcript.  No floating
point, random sample, or post hoc coupling is used.

The theorem is confined to the canonical typed row and to positivity support
of the alternative THM-2592 carrier.  It does not construct (29), identify a
semantic old head, preserve one named Perron ancestry sheet through seven
edges, recover a relation current, apply to the other 164 rows, exclude a
scalar profile, or prove LRC(14).  The ledger remains `165`.

QED (candidate; independent hostile audit pending).
