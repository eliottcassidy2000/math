---
id: THM-3288
title: "Maximizing-witness lifted-walk rational series"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The static maximizing-witness relations of THM-3287 define three finite
  symmetric decorated graphs over the phase backbone, selected tree and full
  twenty-two-edge core.  Their length-n walk sequences have respectively
  18/32, 23/40 and 26/68 active-vertex/directed-arrow censuses, exact minimal
  prefix recurrence orders 10, 14 and 15, and the rational generating
  functions displayed below.  In the full core the even degree-fourteen
  polynomial annihilates only after one adjacency step: its initial scalar
  residual is -1392.  Thus the rational tail denominator has degree fourteen
  but the minimal prefix-valid scalar order and Hankel rank are fifteen.  The
  walks are static witness matchings, not chronological response dynamics.
source: root/creative-synthesis-recover/2026-08-03
audit: >
  The primary exact companion pins and replays the frozen THM-3287 relation
  companion, constructs each active decorated graph, enumerates a long exact
  walk prefix, derives the rational recurrence over Q, checks withheld terms,
  proves the leading and extended Hankel ranks, and factors every recurrence
  and generating denominator.  The independent hostile audit imports or
  executes neither companion: it reconstructs the eleven-by-twenty-two
  response bank and all maximizing relations from THM-3238/3254/3277/3278,
  builds the three symmetric adjacency matrices directly, and derives the
  shortest recurrences by exact Hankel linear algebra rather than
  Berlekamp--Massey.  It also verifies q10(A)1=0 and q14(A)1=0 for the first
  two families, while in the full core 1^T q14(A)1=-1392 and
  A q14(A)1=0.  Both normal and optimized modes byte-match their stored
  outputs and contain no assertion node or floating literal.
depends_on:
  - THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut
related:
  - THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
script: 04-computation/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.py
output: 05-knowledge/results/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.out
script_sha256: f442ff0dfd999fc314f420aacfc7102ba2ba17e033d3fd13e91407727af8a8da
output_sha256: 44a85f317bf9cb35c437e468633a88014de268cc3fb452836960655c4b764951
independent_audit_script: 04-computation/gmc_backbone_maximizing_witness_section_independent_audit_20260803.py
independent_audit_output: 05-knowledge/results/gmc_backbone_maximizing_witness_section_independent_audit_20260803.out
independent_audit_script_sha256: bab55f108a2dd1695465b4195b6565b15c8b8646028ebe756972ea0230621a98
independent_audit_output_sha256: 123ddb9ed988907fd758e6470c1474eaa36c5d52e394ea6efb42d9d114cdf332
hash_basis: LF-normalized bytes
---

# THM-3288 -- maximizing-witness lifted-walk rational series

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The finite relations proved in
[THM-3287](THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md)
admit exact transfer matrices.  Their rational series efficiently compute
arbitrarily long static witness-walk sequences.  The word “walk” here refers
only to the decorated relation graph; Section 7 records the chronological
no-go.

## 1. Decorated graph and walk sequence

Let `E` be one of the following undirected row-edge families:

```text
B = THM-3277's seven-edge phase-geodesic backbone,
T = THM-3277's eleven-edge selected tree T_*,
C = the full twenty-two-edge response core.             (1)
```

Choose either orientation `u -> v` of each row edge.  THM-3287 gives
`R_(v->u)=R_(u->v)^op`, so the symmetric decorated graph below is independent
of that choice.  For every ordered maximizing-witness pair

```text
(s,t) in R_{u->v},                                      (2)
```

put a decorated vertex at `(u,s)` and `(v,t)` and insert both directed
arrows

```text
(u,s) -> (v,t),             (v,t) -> (u,s).             (3)
```

Only decorated vertices incident to an arrow are retained.  Let `A_E` be
the resulting zero-one adjacency matrix and let `1_E` be its all-ones
column.  Define

```text
a_n(E)=1_E^T A_E^n 1_E,                    n>=0.         (4)
```

Thus `a_n(E)` counts length-`n` directed walks in the active decorated
graph, with all possible initial and terminal vertices.  In particular,
`a_0` is the active-vertex count and `a_1` is the directed-arrow count.  The
three exact graph censuses are

```text
family       active vertices       directed arrows
B                   18                    32
T                   23                    40
C                   26                    68.            (5)
```

The matrices are symmetric.  Since every row edge crosses THM-3278's
selector-origin bipartition, the decorated graphs are bipartite; this is the
source of the even-lag recurrence polynomials below.

## 2. Backbone series

For the backbone, the first fifteen coefficients are

```text
(a_0,...,a_14)=
(18,32,84,160,412,798,2060,4012,10426,20362,
 53286,104272,274628,538258,1425644).                   (6)
```

The minimal prefix-valid scalar recurrence has order ten:

```text
a_n = 16a_(n-2) - 94a_(n-4) + 244a_(n-6)
      -263a_(n-8) + 93a_(n-10),             n>=10.      (7)
```

Its recurrence characteristic polynomial is

```text
q_B(z)=z^10-16z^8+94z^6-244z^4+263z^2-93
      =(z^4-5z^2+3)(z^6-11z^4+36z^2-31).               (8)
```

At the vector level,

```text
q_B(A_B)1_B=0.                                          (9)
```

Consequently, for `G_B(x)=sum_(n>=0) a_n(B)x^n`, one has

```text
G_B(x)=N_B(x)/D_B(x),                                   (10)

D_B(x)=(1-5x^2+3x^4)(1-11x^2+36x^4-31x^6),

N_B(x)=2(279x^9+216x^8-762x^7-514x^6+623x^5
           +380x^4-176x^3-102x^2+16x+9).
```

The exact Hankel rank is ten, so no smaller scalar recurrence valid from its
own order index generates the complete prefix `(6)` and its continuation.

## 3. Selected-tree series

For the selected tree,

```text
(a_0,...,a_14)=
(23,40,104,212,536,1174,2950,6700,16850,39058,
 98502,231232,584846,1384502,3510318).                  (11)
```

The minimal prefix-valid scalar recurrence has order fourteen:

```text
a_n = 20a_(n-2) - 158a_(n-4) + 630a_(n-6)
      -1343a_(n-8) + 1493a_(n-10) - 774a_(n-12)
      +144a_(n-14),                         n>=14.      (12)
```

Its recurrence characteristic polynomial factors as

```text
q_T(z)=z^14-20z^12+158z^10-630z^8+1343z^6
             -1493z^4+774z^2-144
      =(z^2-2)(z^4-5z^2+3)
        (z^8-13z^6+54z^4-77z^2+24),                    (13)
```

and the vector identity is

```text
q_T(A_T)1_T=0.                                          (14)
```

The rational generating function is

```text
G_T(x)=N_T(x)/D_T(x),                                   (15)

D_T(x)=(1-2x^2)(1-5x^2+3x^4)
        (1-13x^2+54x^4-77x^6+24x^8),

N_T(x)=1152x^13+984x^12-5952x^11-4745x^10
       +10710x^9+7907x^8-8484x^7-5828x^6
       +3254x^5+2090x^4-588x^3-356x^2+40x+23.
```

The exact Hankel rank is fourteen, proving the stated prefix minimality.

## 4. Full-core series and the one-step initial atom

For the full core,

```text
(a_0,...,a_14)=
(26,68,274,1018,4134,16400,67270,270718,1114054,
 4498562,18529428,74885666,308530344,1247186562,
 5138801640).                                            (16)
```

Define the even degree-fourteen polynomial

```text
q_C(z)=z^14-34z^12+404z^10-2285z^8+6701z^6
             -10096z^4+7156z^2-1808
      =(z-2)(z+2)
        (z^12-30z^10+284z^8-1149z^6
         +2105z^4-1676z^2+452).                         (17)
```

Unlike `(9)` and `(14)`, this polynomial does not annihilate the initial
vector.  The exact load-bearing identities are

```text
1_C^T q_C(A_C)1_C = -1392,
q_C(A_C)1_C       != 0,
A_C q_C(A_C)1_C    = 0.                                 (18)
```

Therefore the degree-fifteen polynomial

```text
z q_C(z)
 =z(z-2)(z+2)
   (z^12-30z^10+284z^8-1149z^6
    +2105z^4-1676z^2+452)                              (19)
```

annihilates `1_C`.  Equivalently, the minimal recurrence valid from its own
order index is the order-fifteen law

```text
a_n = 34a_(n-2) - 404a_(n-4) + 2285a_(n-6)
      -6701a_(n-8) + 10096a_(n-10) - 7156a_(n-12)
      +1808a_(n-14) + 0a_(n-15),            n>=15.      (20)
```

The same even relation fails at `n=14` by exactly `-1392`.  This is not a
cosmetic trailing zero: the exact Hankel rank is fifteen, so an order-
fourteen prefix-valid recurrence is impossible.

Nevertheless the eventual recurrence uses only fourteen lags, and the
rational generating function has degree-fourteen tail denominator

```text
G_C(x)=N_C(x)/D_C(x),                                   (21)

D_C(x)=(1-4x^2)
        (1-30x^2+284x^4-1149x^6
         +2105x^8-1676x^10+452x^12),

N_C(x)=-2(696x^14-3808x^13-8672x^12+17140x^11
          +21970x^10-24644x^9-22573x^8+15495x^7
          +11000x^6-4630x^5-2661x^4+647x^3
          +305x^2-34x-13).
```

In particular, the coefficient of `x^14` in `N_C` is `-1392`, exactly the
initial scalar residual in `(18)`.  Multiplication by `D_C` vanishes from
degree fifteen onward, not from degree fourteen.  Equations `(18)--(21)`
reconcile the degree-fourteen rational tail with the minimal prefix order
and Hankel rank fifteen.

## 5. Exact mechanism and minimality

For any of the three families, the finite transfer identity `(4)` makes the
series rational.  The primary companion performs the following exact steps:

1. replay the frozen THM-3287 response and relation companion and require a
   byte match with its stored transcript;
2. build the active decorated nodes and symmetric arrow set without
   multiplicity;
3. compute more than four times the active-vertex count of exact terms;
4. derive a rational Berlekamp--Massey recurrence, verify every withheld
   term, and compute the leading and extended exact Hankel ranks; and
5. recover and factor the recurrence characteristic polynomial, generating
   denominator and generating numerator coefficientwise.

The independent audit follows a separate path.  It reconstructs the
reset-link response bank and every relation without importing the primary
relation companion, forms each symmetric adjacency matrix directly, and
solves the shortest exact Hankel linear system instead of using
Berlekamp--Massey.  It then checks the vector annihilator identities `(9)`,
`(14)` and `(18)`.  The vector identities prove the recurrences for all
future indices; the full-rank leading Hankel matrices prove minimality.  The
initial coefficients then determine the numerators in `(10)`, `(15)` and
`(21)` exactly.

For any fixed family, companion-matrix powering evaluates `a_n` from the
displayed recurrence in `O(r^3 log n)` exact ring-arithmetic operations, with
`r` equal to `10`, `14`, or `15`; this is an arithmetic-operation statement
and makes no claim about coefficient bit complexity.

## 6. Reproduction

Run the primary companion in both modes:

```text
python3 04-computation/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.py
python3 -O 04-computation/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.out.
```

Its source/output hashes are

```text
f442ff0dfd999fc314f420aacfc7102ba2ba17e033d3fd13e91407727af8a8da
44a85f317bf9cb35c437e468633a88014de268cc3fb452836960655c4b764951.
```

Run the independent reconstruction in both modes:

```text
python3 04-computation/gmc_backbone_maximizing_witness_section_independent_audit_20260803.py
python3 -O 04-computation/gmc_backbone_maximizing_witness_section_independent_audit_20260803.py
```

and compare with its declared transcript and frontmatter hashes.  The audit
pins both frozen primary companions but deliberately does not pin either
theorem file, so promotion creates no theorem--audit hash cycle.

## 7. Static-walk boundary

The decorated arrows in `(3)` are not physical state transitions.  THM-3287
proves that every dominance arrow has count-vector `L1` distance two, while
a physical one-pole arrow has distance one and every Q-directed reset-link
arrow ends at `Q`.  A decorated walk may also revisit rows and witnesses;
it is not one of THM-3287's vertex-simple row-path lifts.

Accordingly, `(4)` counts sequences of compatible static maximizing-witness
matchings only.  The rational forms do not compose response updates, produce
a chronological trajectory, continue after reset, or prove a positive
current.  They imply no Gaussian-moment, factorial, Jacobian, tournament, or
`LRC(14)` statement.
