---
id: THM-4403
title: "Source-normal two-channel weight-eighteen row-twelve affine continuation"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4395 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Adjoining only the two late channels lambda18*p^3*y^4 and
  nu18*y^6 to the complete fixed residual-weight-at-most-fourteen source
  turns rows eleven and twelve into a global constant-pivot graph. The
  source projection through row twelve is A^4 over
  (Phi,eta,alpha11,c51), and every terminal fibre is A^9. The complete
  exact-t-valuation-twelve source diagonal has response rank one; its odd
  channels vanish, and no combination of its seven channels pays the second
  row-twelve condition. This is a restricted two-channel tail, not the
  complete weight-at-most-eighteen family, and it proves no chart/seam entry,
  all-row lift, Keller pair, JC(2), or DC(2) claim.
source: root + independent referee / JC2 and arXiv:2608.23777 continuation session, 2026-09-03
audit: >
  PASS. The primary certificate reuses the frozen row-eleven construction,
  computes both row-twelve source conditions, audits every exact-valuation-
  twelve monomial, and verifies the full bracket/depth graph at dense and
  Phi=eta=0 controls. A separate clean-room certificate imports only the
  audited THM-4308/4315 operators: it reconstructs the complete weight-at-
  most-fourteen state through row ten and the row-eleven late response before
  deriving row twelve. Its 123 exact checks reproduce every matrix, pivot,
  condition, response ratio, graph, and control. Normal, optimized, and
  hash-seeded clean-room streams byte-match the frozen LF transcript. No
  finite-field inference is used.
depends_on:
  - THM-4399-source-normal-row-eleven-late-weight-eighteen-response-absorption
  - THM-4395-source-normal-weight-fourteen-row-ten-global-affine-absorption
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension
primary_script: 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403.py
primary_output: 05-knowledge/results/jc2_source_normal_weight18_row12_affine_continuation_thm4403.out
primary_script_sha256: 76cd8f00c703f9d8093ae6f7cbbf9ae6ba0a3e2c144d566f98b016b109574e43
primary_output_sha256: a0fa7ba7b0e29dbff31340fa14785530ebcf91a75427fd55ae12cce653250055
independent_audit_script: 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.out
independent_audit_script_sha256: fdc4811c4a1ace2dee4512a62b0d8211b697c92d81bcf08e1e9230b0c6a61742
independent_audit_output_sha256: 7ed3ade1614d67458041206da70ae05f33b9dfc62d0c58dcca7376e45ee9a920
hash_basis: raw LF bytes
---

# THM-4403 -- source-normal two-channel weight-eighteen row-twelve affine continuation

**PROVED FINITE-ROW RELATIVE TO THM-4395 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THIS THEOREM ADJOINS ONLY TWO SELECTED WEIGHT-EIGHTEEN CHANNELS. IT
DOES NOT TREAT THE COMPLETE WEIGHT-FIFTEEN-THROUGH-EIGHTEEN SOURCE, CHART OR
SEAM ENTRY, AN ALL-ROW LIFT, POLYNOMIAL TERMINATION, A KELLER PAIR, `JC(2)`,
OR `DC(2)`.**

## 1. Universe and entry rows

Work over a characteristic-zero field in THM-4395's fixed source-normal
chart. Its complete residual-weight-at-most-fourteen source is

```text
H_14=H_12+c51*p^5*y+c23*p^2*y^3
           +c70*p^7+c42*p^4*y^2+c14*p*y^4.              (1)
```

Adjoin exactly

```text
H_tail=lambda18*p^3*y^4+nu18*y^6.                       (2)
```

Here `p=t(1+x^2*t)` and `y=x*t*p`. Consequently

```text
[t^r](p^3*y^4)=0  for r<11,       [t^11](p^3*y^4)=x^4,
[t^r](y^6)=0      for r<12,       [t^12](y^6)=x^6.      (3)
```

Thus `(2)` changes none of THM-4395's rows through ten. Before `(2)`, its
source graph is `A^6` over

```text
(Phi,eta,alpha11,c51,c23,c70),                           (4)
```

with solved coordinates `(beta11,c42,xi10,c14)` and terminal row-ten fibre
`A^8`.

The late-channel idea was prompted by the Hamiltonian correction in the
rank-two Poisson construction associated with arXiv:2608.23777: preserve the
already solved prefix and let a first-visible primitive pay the next
compatibility. That is an idea transfer only; no implication from the
Poisson counterexample is used here.

## 2. Row eleven: the first late response

The degree-capped row-eleven bracket system has shape `12 x 8`, rank eight,
and pivots `(0,...,7)`. Before `lambda18` it leaves the principal quartic
obstruction `F=0`, where

```text
F=
 712092531095885035687500*Phi^4
+821011731481169865000000*Phi^3*alpha11
+139278775876269887812500*Phi^3*c51
-316295855547474501562500*Phi^3*eta
+236647708617429900000000*Phi^2*alpha11^2
+80291186852342287500000*Phi^2*alpha11*c51
-43058909665599049687500*Phi^2*alpha11*eta
+6810413170511176171875*Phi^2*c51^2
-30932285940138480468750*Phi^2*c51*eta
+3733889920481087640000000*Phi^2*c70
+445628783011794805546875*Phi^2*eta^2
-76479730269188998356067968000*Phi^2
+80291186852342287500000*Phi*alpha11^2*eta
+13620826341022352343750*Phi*alpha11*c51*eta
+2152506377265652800000000*Phi*alpha11*c70
+205715422677291419531250*Phi*alpha11*eta^2
-44128098153451944185391360000*Phi*alpha11
+21736146783278091456000000*Phi*c23
+365157331857566100000000*Phi*c51*c70
+40145593426171143750000*Phi*c51*eta^2
-7301988813385081335237120000*Phi*c51
-829255929072250500000000*Phi*c70*eta
-91168842770934468750000*Phi*eta^3
+16987951185108221412403200000*Phi*eta
+6810413170511176171875*alpha11^2*eta^2
+12679418956912220016000000*alpha11^2
+365157331857566100000000*alpha11*c70*eta
+40145593426171143750000*alpha11*eta^3
-7301988813385081335237120000*alpha11*eta
+25358837913824440032000000*c51*eta
+4894705859649350400000000*c70^2
+1076253188632826400000000*c70*eta^2
-200365448846185933658357760000*c70
+59161927154357475000000*eta^4
-22073830342778447233850880000*eta^2
+2077158990491529775184308229636096.                  (5)
```

It is primitive and linear in `c23`, with

```text
dF/dc23=21736146783278091456000000*Phi.                 (6)
```

The late channel changes the exact equation to

```text
F+6864046352614134144000000*lambda18=0,                 (7)
```

so it gives the global graph

```text
lambda18=-F/6864046352614134144000000.                  (8)
```

No source localization is used. After `(8)`, the row-eleven projected-depth
universes and terminal system are

```text
pi_11(P_2): 102 x 228, rank 77,
pi_11(P_3): 114 x 361, rank 94,
joint:       45 x 12, rank 4, pivots (8,9,10,11).       (9)
```

They impose no further source equation and leave terminal fibre `A^8`.

## 3. Row twelve: two bracket conditions

After `(8)` and the row-eleven depth selection, the row-twelve bracket matrix
has shape `13 x 8`, rank eight, and pivots `(0,...,7)`. Its entire source
ideal is generated by

```text
E_nu=N+29652680243293059502080000000*nu18,
L=-23839*Phi-675*alpha11-675*c23+4266*eta,               (10)
```

where the exact 38-term constant part of the first generator is

```text
N=
 16207504788526268679854437500*Phi^4
+19972502596938878685405000000*Phi^3*alpha11
+3388192404837845491274062500*Phi^3*c51
-7694433044124931768682812500*Phi^3*eta
+6127519119231112400700000000*Phi^2*alpha11^2
+2078979701167698850237500000*Phi^2*alpha11*c51
-1333077286897766907413437500*Phi^2*alpha11*eta
+176342028224045884618359375*Phi^2*c51^2
-800929679848005674777343750*Phi^2*c51*eta
+128176893600549718154832000000*Phi^2*c70
+10895688995372876585505234375*Phi^2*eta^2
-1290637505065488046138160750592000*Phi^2
+2078979701167698850237500000*Phi*alpha11^2*eta
+352684056448091769236718750*Phi*alpha11*c51*eta
+77262672929094354925440000000*Phi*alpha11*c70
+5326589439383106725922656250*Phi*alpha11*eta^2
-812365382482794238187049984000000*Phi*alpha11
+116853525106903019667456000000*Phi*c23
+13107060586185649496280000000*Phi*c51*c70
+1039489850583849425118750000*Phi*c51*eta^2
-139393159906770648160569446400000*Phi*c51
-29765546945236416497400000000*Phi*c70*eta
-2360634845867806199343750000*Phi*eta^3
+317839241970920124304409640960000*Phi*eta
+176342028224045884618359375*alpha11^2*eta^2
+230461118960836511010816000000*alpha11^2
+83988471170586545385984000000*alpha11*c51
+13107060586185649496280000000*alpha11*c70*eta
+1039489850583849425118750000*alpha11*eta^3
-139373440874408858276000563200000*alpha11*eta
+65730107872632948562944000000*c23*eta
+460922237921673022021632000000*c51*eta
+224645295214590804368640000000*c70^2
+38631336464547177462720000000*c70*eta^2
-5739033021174331192331917590528000*c70
+1531879779807778100175000000*eta^4
-406434656654908878729682944000000*eta^2
+22658044532920649950285301653912420352.              (11)
```

Both pivots are constant. In triangular order, the two bracket graphs are

```text
c23_bar=-23839*Phi/675-alpha11+158*eta/25,
nu18_bar=-N(Phi,eta,alpha11,c51,c23_bar,c70)
          /29652680243293059502080000000.                (12)
```

Substitution into every one of the thirteen literal row-twelve bracket
coefficients gives zero.

## 4. The complete exact-valuation-twelve response line

All source monomials of exact `t`-valuation twelve are

```text
p^12, p^10*y, p^8*y^2, p^6*y^3, p^4*y^4, p^2*y^5, y^6,
```

with residual weights `(24,23,22,21,20,19,18)` and leading row `(1,x,...,x^6)`.
In the exact source-cokernel basis `(E_nu,L)`, their response matrix is

```text
[83988471170586545385984000000, 0,
 21910035957544316187648000000, 0,
 18780030820752271017984000000, 0,
 29652680243293059502080000000]
[0,                             0,
 0,                             0,
 0,                             0,
 0].                                                            (13)
```

It has rank one. Relative to the `y^6` response, the ratios are

```text
(3059/1080, 0, 133/180, 0, 19/30, 0, 1).                (14)
```

These numbers are not accidental. For THM-4315's row-twelve Student law

```text
mu_12(dx)=constant*(x^2+6)^(-13) dx,
```

the ratios in `(14)` are exactly

```text
E[1]/E[X^6], E[X]/E[X^6], ..., E[X^6]/E[X^6].           (15)
```

Thus parity kills every odd channel, while the four even channels all pay
the same Student--Stein quotient direction. The second generator `L` is a
different sidecar left by the prior depth-selected truncation. Every column
of `(13)` has zero `L` component, so **no combination of any exact-valuation-
twelve source monomials can pay `L`**. For example,

```text
(Phi,eta,alpha11,c23)=(0,0,1,0)
```

gives `L=-675`; arbitrary coefficients on all seven same-row monomials leave
that failure unchanged. The continuation works only because the inherited
lower-row coordinate `c23` remains free and pays `L` by `(12)`.

## 5. Projected depth and the global affine graph

After the bracket graphs `(12)`, append the full row-twelve tangent. The
exact projected-depth universes are

```text
pi_12(P_2): 117 x 267, rank 87,
pi_12(P_3): 130 x 424, rank 105.                         (16)
```

Their joint terminal system has shape `55 x 13`, rank four, and pivots
`(9,10,11,12)`. It leaves exactly one source condition:

```text
G=
-95023209310875*Phi^2
-157540468455000*Phi*alpha11
-21622415197500*Phi*c51
-4375612912500*Phi*eta
-21622415197500*alpha11*eta
-368130186720000*c70
-78770234227500*eta^2
+178582220348850176=0.                                  (17)
```

The constant `c70` pivot gives

```text
c70_bar=(178582220348850176
         -95023209310875*Phi^2
         -157540468455000*Phi*alpha11
         -21622415197500*Phi*c51
         -4375612912500*Phi*eta
         -21622415197500*alpha11*eta
         -78770234227500*eta^2)
        /368130186720000.                                (18)
```

Finally substitute `c70_bar` into `c23_bar`, `lambda18_bar`, `nu18_bar`, and
THM-4395's inherited four graphs. Every denominator is a nonzero rational
integer. Therefore the source projection is the global affine graph

```text
A^4_(Phi,eta,alpha11,c51),                               (19)
```

and the thirteen row-twelve tangent coordinates modulo the four depth pivots
give terminal fibre `A^9` over every source point.

The certificates evaluate every original bracket and depth coefficient at a
dense point `(Phi,eta,alpha11,c51)=(1,2,3,5)` and at the boundary point
`(0,0,1,2)`. All vanish, all four newly solved coordinates are nonzero, and
no division by `Phi`, `eta`, or another source parameter occurs.

## 6. Scope and mechanism

The positive mechanism has two distinct stages:

1. the first-visible late channel `p^3*y^4` pays the row-eleven principal
   obstruction without disturbing the solved prefix;
2. at row twelve, the Student-visible response line pays `E_nu`, but the
   Student-invisible condition `L` must be paid by the inherited coordinate
   `c23`; projected depth then spends `c70`.

This is sharper than a generic "add more coefficients" heuristic: the full
same-row source diagonal has rank one and provably cannot repair `L`.
Successful continuation depends on retained lower-row memory.

The universe remains strictly

```text
complete residual weight <=14
plus lambda18*p^3*y^4 plus nu18*y^6.                     (20)
```

It excludes every other residual-weight-fifteen-through-eighteen monomial,
as well as all higher weights. It supplies no converse from projected depth
to full `B_2` membership, no chart or seam entry, no all-row compatible
system, no polynomial termination, and no conclusion about `JC(2)` or
`DC(2)`.

## 7. Reproduction

Artifacts:

```text
04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403.py
05-knowledge/results/jc2_source_normal_weight18_row12_affine_continuation_thm4403.out
04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.py
05-knowledge/results/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.out
```

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403.py
python3 -B -O 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403.py
python3 -B 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_weight18_row12_affine_continuation_thm4403_independent_audit.py
```

The independent audit imports neither THM-4395 nor either late-response
implementation. It reconstructs the complete source from THM-4308/4315's
audited row operators and checks normal, optimized, and nondefault hash-seed
streams byte for byte.
