---
id: THM-4415
title: "Source-normal row-fourteen same-row response-rank obstruction"
status: >
  PROVED FINITE-ROW OBSTRUCTION RELATIVE TO THM-4410 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On THM-4410's global A^4 restricted source, the
  row-fourteen bracket has two active quartic conditions. Every polynomial
  Hamiltonian perturbation that leaves rows zero through thirteen unchanged
  has row-fourteen response in one Student-visible line, so one exact quartic
  J14 is untouched. The bracket-solvable base locus is exactly V(J14); the
  other condition has a constant pivot from p^2*y^6, the least-weight visible
  same-row channel. This is not the complete weight-at-most-twenty-two family
  and includes no projected-depth, chart/seam entry, all-row lift, Keller
  pair, JC(2), or DC(2) claim.
source: root + row13_late_scout / JC2 continuation session, 2026-09-05
audit: >
  PASS. The primary certificate uses the inherited elimination routine and
  audits the full valuation filtration. A separate certificate imports only
  the proved THM-4410 terminal state, builds the 15x9 coefficient matrix
  directly, and computes its six-dimensional left nullspace without reading
  the row-fourteen scout. Both obtain two active quartics, response rank one,
  identical parity and moment ratios, the same constant p^2*y^6 pivot, and
  exact expression hashes 878ff1... for the paid condition and d15ce1... for
  J14. The audit evaluates J14 nontrivially at the origin and a dense point.
  Normal, optimized, and fixed-hash-seed streams byte-match both frozen LF
  outputs. No finite field is used.
depends_on:
  - THM-4410-source-normal-least-weight-twenty-row-thirteen-affine-continuation
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4403-source-normal-two-channel-weight-eighteen-row-twelve-affine-continuation
  - THM-4395-source-normal-weight-fourteen-row-ten-global-affine-absorption
primary_script: 04-computation/jc2_source_normal_weight22_row14_conditional_scout_s616.py
primary_output: 05-knowledge/results/jc2_source_normal_weight22_row14_conditional_scout_s616.out
primary_script_sha256: 69aea05cf8169e3693ebfaecadccf93a577725beb269c235d45311a76e3c0166
primary_output_sha256: 6169be875c35859986a3490e76d74f45592a6497436178de6dafced8f4f01676
independent_audit_script: 04-computation/jc2_source_normal_row14_same_row_obstruction_independent_audit_s616.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_row14_same_row_obstruction_independent_audit_s616.out
independent_audit_script_sha256: 868f212a97f6c42469aa2cdee841cbff71e6713fbbd9d814cc6b9a2930a3d50a
independent_audit_output_sha256: e3acb38f275ac36166df3c2678f52b39a1578d43f9e306cd2913963c87da5a3e
hash_basis: raw LF bytes
---

# THM-4415 -- row fourteen has two active conditions but only one same-row response

**PROVED FINITE-ROW OBSTRUCTION RELATIVE TO THM-4410 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS IS A BRACKET-LEVEL THEOREM FOR ONE RESTRICTED
TAIL. IT DOES NOT TREAT THE COMPLETE RESIDUAL-WEIGHT-AT-MOST-TWENTY-TWO
SOURCE, PROJECTED DEPTH AT ROW FOURTEEN, CHART OR SEAM ENTRY, AN ALL-ROW LIFT,
POLYNOMIAL TERMINATION, A KELLER PAIR, `JC(2)`, OR `DC(2)`.**

## 1. Inherited restricted source

Work over a characteristic-zero field in THM-4410's fixed source-normal
chart. The Hamiltonian universe through row thirteen is

```text
H_14+lambda18*p^3*y^4+nu18*y^6+kappa20*p*y^6.          (1)
```

The five coordinates

```text
(c23,c70,lambda18,nu18,kappa20)                        (2)
```

are constant-denominator functions of

```text
(Phi,eta,alpha11,c51),                                 (3)
```

the source projection is the global `A^4` on `(3)`, and the row-thirteen
terminal fibre is `A^9`.

Substitute all graphs `(2)` and the row-thirteen depth selector. The
row-fourteen bracket coefficient matrix has

```text
shape 15 x 9,
rank 9,
column pivots (0,1,2,3,4,5,6,7,8),
raw cokernel dimension 6.                              (4)
```

Four raw cokernel coordinates vanish on the inherited source. Exactly two
primitive quartic conditions remain.

## 2. Exhaustive same-row Hamiltonian filtration

For a monomial in the polynomial Hamiltonian ring `K[p,y]`,

```text
ord_t(p^a*y^b)=a+2b,
[t^(a+2b)](p^a*y^b)=x^b.                              (5)
```

If a nonzero finite linear combination has least valuation `v`, its leading
coefficient is

```text
sum_(a+2b=v) c_(a,b) x^b.                              (6)
```

The exponents `b` in `(6)` are distinct, so the polynomial cannot vanish.
Consequently a Hamiltonian perturbation leaves every row below fourteen
unchanged iff all of its monomials have valuation at least fourteen. Modulo
terms of valuation greater than fourteen, every such perturbation is a unique
linear combination of

```text
p^14, p^12*y, p^10*y^2, p^8*y^3,
p^6*y^4, p^4*y^5, p^2*y^6, y^7.                       (7)
```

Thus `(7)` is not a sample or degree cutoff. It is the complete associated-
graded Hamiltonian diagonal relevant to row fourteen.

## 3. Rank-one response and the unpaid quartic

Let `rho_b` multiply the term of `(7)` with leading coefficient `x^b`. In a
primitive basis for the two active conditions, the first equation is

```text
P_14
+7363809000324505922074560000000000*rho_0
+1636402000072112427127680000000000*rho_2
+1178209440051920947531929600000000*rho_4
+1536794921806853409824256000000000*rho_6=0,           (8)
```

while the second equation is

```text
J_14=0,                                                (9)
```

with no `rho_b` term at all. Here `P_14` is the exact 29-term quartic frozen
in both companion outputs. The exact 28-term unpaid quartic is

```text
J_14=
-66170293788961831190426578125*Phi^4
-179894512500147625637441250000*Phi^3*alpha11
-24657789972417821823870000000*Phi^3*c51
-5413628051753727156750000000*Phi^3*eta
+492009926651486077837275000000*Phi^2*alpha11^2
+111708261200909110733100000000*Phi^2*alpha11*c51
+300376069086381618078630000000*Phi^2*alpha11*eta
+6286399947710619715420312500*Phi^2*c51^2
+38283138595085344778671875000*Phi^2*c51*eta
-45092607213843940505087812500*Phi^2*eta^2
+6658817276318586869516875333632000*Phi^2
+111708261200909110733100000000*Phi*alpha11^2*eta
+12572799895421239430840625000*Phi*alpha11*c51*eta
+530293065246571422615946875000*Phi*alpha11*eta^2
+54998920884002181487057780899840000*Phi*alpha11
+55854130600454555366550000000*Phi*c51*eta^2
+6624954319110052779115035770880000*Phi*c51
+162516929529399719951250000000*Phi*eta^3
+12735360517907219957620980940800000*Phi*eta
+6286399947710619715420312500*alpha11^2*eta^2
+7982620169970669746503680000000*alpha11^2
+2909159111239310981560320000000*alpha11*c51
+55854130600454555366550000000*alpha11*eta^3
+6623360605857808634838180986880000*alpha11*eta
+15965240339941339493007360000000*c51*eta
+123002481662871519459318750000*eta^4
+27507145681461912505930612408320000*eta^2
+689373657131467757187366756679487062016.             (10)
```

The full response matrix of `(7)` against the raw six-dimensional cokernel
has rank one. In the Student-visible coordinate its response vector is

```text
(1671525,0,371450,0,267444,0,348840,0),                (11)
```

and relative to `p^2*y^6` the ratios are

```text
(115/24,0,115/108,0,23/30,0,1,0).                     (12)
```

Parity kills every odd channel. The weight-21 term `y^7` is therefore
invisible, and

```text
p^2*y^6,                    residual weight 22,        (13)
```

is the least-weight visible row-fourteen channel. Its coefficient in `(8)`
is the nonzero constant

```text
1536794921806853409824256000000000.                   (14)
```

## 4. Exact bracket-solvable locus

Equation `(14)` lets `rho_6` solve `(8)` over every base point, with no
localization in `(3)`. Equation `(9)` is untouched by **every** Hamiltonian
of `t`-valuation at least fourteen, by the exhaustive argument of Section 2.
Therefore the row-fourteen bracket-solvable base locus, after allowing all
prefix-preserving Hamiltonian perturbations, is exactly

```text
boxed: V(J_14) subset A^4_(Phi,eta,alpha11,c51).       (15)
```

This is a proper hypersurface. At the origin,

```text
J_14(0,0,0,0)=689373657131467757187366756679487062016
              !=0,                                    (16)
```

and at the dense control `(Phi,eta,alpha11,c51)=(1,2,3,5)` its value is

```text
689753974312897399008554429060173258391 !=0.           (17)
```

Thus the global `A^4` continuation mechanism of rows eleven--thirteen
genuinely stops at row fourteen if one insists on changing only terms that
preserve the solved prefix. This is stronger than failure of one guessed
monomial: it is an associated-graded no-go for the whole Hamiltonian ring.

## 5. What the obstruction asks for next

The theorem does not say that the complete weight-at-most-twenty-two family
is empty. It identifies the missing resource. Any generic continuation must
use at least one coefficient of valuation below fourteen whose earlier-row
effects cancel against other omitted weight-15-through-20 terms, thereby
retaining a lower-row **memory sidecar** that reaches the `J_14` direction.
The alternative is to work on the proper hypersurface `(15)` and then audit
projected depth there.

This sharpens the late-channel strategy:

```text
rows 11--13: one newly visible Student direction suffices after inherited
             lower-row memory is spent at row 12;
row 14:      two active directions appear, but the new diagonal still has
             only the one Student direction.                          (18)
```

No statement is made about whether the complete omitted coefficient family
provides the required memory, whether `(15)` survives projected depth, or
whether either route enters a planar chart.

## 6. Independent reconstruction and reproduction

The primary certificate obtains `(4),(8)--(14)` with the inherited
elimination routine. The independent audit reads no row-fourteen scout: it
forms the `15 x 9` coefficient matrix directly, computes the left nullspace,
and contracts the eight basis vectors `(7)` into that cokernel. After
primitive normalization, both paths give the exact expression hashes

```text
paid condition with rho variables:
  878ff13b4b70f01dde023c73c4b44b0e4f310320e56fce87750fdf41c8bacf70
unpaid J_14:
  d15ce156dea6b7883378afebf23262d7950e351b845b0894986cbff39938884f. (19)
```

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight22_row14_conditional_scout_s616.py
python3 -B -O 04-computation/jc2_source_normal_weight22_row14_conditional_scout_s616.py
PYTHONHASHSEED=271828 python3 -B 04-computation/jc2_source_normal_weight22_row14_conditional_scout_s616.py
python3 -B 04-computation/jc2_source_normal_row14_same_row_obstruction_independent_audit_s616.py
python3 -B -O 04-computation/jc2_source_normal_row14_same_row_obstruction_independent_audit_s616.py
PYTHONHASHSEED=161803 python3 -B 04-computation/jc2_source_normal_row14_same_row_obstruction_independent_audit_s616.py
```

The primary and independent certificates perform `42` and `44` dynamic exact
checks respectively. All three modes byte-match each frozen LF output.
