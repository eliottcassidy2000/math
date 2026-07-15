---
id: THM-815
title: Scale-one Hamming-four safe-component closure
status: PROVED by two independent reductions (sharp interval-comb component ladder; collar-cycle/doubling box), with a finite-recursion theorem through radius six + FINITE-EXACT (666,705 nested Hamming-four rows, two higher-radius initial censuses, two 35,640-row component/endpoint replays, and an independent 768,735-row C++ collar certificate)
source: codex-2026-07-15-S10 Hamming-four continuation
depends_on: [LRC(<=13), THM-795, THM-800, THM-804, THM-806, THM-810, THM-816]
related: [THM-770, THM-800, THM-804, THM-810, THM-816, THM-820, THM-845, HYP-6820]
verification:
  - 04-computation/lrc13_h4_scale_one_component_ladder_codex_S10.py
  - 05-knowledge/results/lrc13_h4_scale_one_component_ladder_codex_S10.out
  - 04-computation/lrc13_scale_one_hamming_four_collar_closure_codex_S10.cpp
  - 05-knowledge/results/lrc13_scale_one_hamming_four_collar_closure_codex_S10.out
proof_companion:
  - 07-reflections/lrc13-hamming-four-independent-collar-doubling-proof-codex-S10.md
---

# THM-815 — Scale-one Hamming-four safe-component closure

Put

```text
delta=1/13,                         [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

### A. Every proper scale-one four-lift is loose

Let `R` be any four-element subset of `[12]`, and for every `r in R` let

```text
u_r=r+13h_r,                       h_r>=1.               (1)
```

Then

```text
B=([12] minus R) union {u_r:r in R}
```

satisfies

```text
M(B)>1/13.                                               (2)
```

Thus the scale-one Hamming-four chart left by THM-810 is empty.

### B. Uniform full-residue rigidity through Hamming radius four

Together, THM-795, THM-800, THM-804/806, THM-810, the present theorem, and
THM-816 show the following.  Let `13` not divide `c`, start with `c[12]`, and
properly replace at most four coordinates by positive speeds in the same
nonzero residue classes modulo `13`.  Every nontrivial resulting packet is
loose at `1/13`.

Indeed THM-810 says that a hypothetical tight four-coordinate packet either
has common scale or lies on its order-three coset interface.  The latter is
empty by THM-816.  In the common-scale case divide by `c`; four genuine lifts
give Part A, while any zero-height coordinate reduces to the already closed
Hamming-at-most-three cases.  Multiplication by `c` is onto on the circle, so
division preserves `M`.

### C. The scale-one charts through radius six are recursively finite

More generally, fix `m<=6`, choose `R subset [12]` with `|R|=m`, and replace
every `r in R` by a proper lift `r+13h_r`.  Whether the resulting scale-one
Hamming-`m` packet is tight is decidable by a finite exact recursion.  Thus
Parts A--B close radius four, while radii five and six are finite residuals;
the first shallow radius at which this discrepancy recursion supplies no
speed bound is seven.

Indeed, order the lifted speeds `x_1<...<x_m`.  After adjoining the first
`j<m` lifts to `[12] minus R`, the prefix has at most eleven speeds.  Settled
`LRC(<=13)` gives it maximin at least `1/12>1/13`, so its strict-safe set has
a positive-length component.  Suppose inductively that the first `j` lifts
range over a finite set.  Among those finitely many prefixes let `L_j>0` be
the minimum length of a longest strict-safe component.  If the completed
packet were tight, its `m-j` remaining danger combs would cover that component.
Formula (8) therefore gives the uniform next-speed bound

```text
x_(j+1) <= floor(22(m-j)/[13(13-2(m-j))L_j]).
```

The denominator is positive exactly because `m-j<=6`.  This proves the
induction and leaves finitely many full packets, all testable by (29)--(30).
The exact first step is already modest enough to record:

```text
m=5: min L_0=1/52,  uniquely at R=(3,4,5,11,12),
     component=(9/26,19/52),                  x_1<=146;
m=6: min L_0=31/1430, uniquely at R=(2,3,5,6,8,12),
     component=(53/143,51/130),               x_1<=468.
```

These bounds leave respectively `40,590` and `194,040` labelled first-prefix
states. THM-820 supplies an independent collar dichotomy at `m=5`; intersecting
its branch inequalities with `x_1<=146` gives the explicit joint boxes
`(146,292,584,1168,2336)` and
`x<=146,v<=1986,max(v,w,y,z)<=7944`. For `m=7`, the initial mean danger
density is `14/13`; the coefficient
`13-2m` changes sign, so this union-bound recursion becomes noncoercive before
the first lift.  No claim of radius-five or radius-six emptiness is made here,
and the density barrier is a limitation of this method, not a proof that the
radius-seven chart is intrinsically infinite.

## 1. Strict-safe components

For a finite positive speed set `Q`, write

```text
E_Q={t in R/Z:min_(q in Q)||qt||>delta}.                (3)
```

For one speed `q`, its strict-safe bands in `(0,1)` are

```text
((13k+1)/(13q),(13(k+1)-1)/(13q)),       0<=k<q.        (4)
```

Intersecting the finitely many ordered unions (4) gives every component of
`E_Q` with exact rational endpoints.  Let

```text
L(Q)=maximum length of a component of E_Q.              (5)
```

All cores below have nonempty strict-safe sets; this is also certified
directly by the endpoint computations.

## 2. Sharp one-interval danger discrepancy

For a positive speed `u`, put

```text
D_u={t:||ut||<=delta}.
```

For every interval `I` of length `L`,

```text
|I intersect D_u| <= 2L/13 + 22/(169u).                (6)
```

### Proof

After scaling by `u`, remove all complete unit periods.  On the remaining
circle interval of length `s`, the danger arc has length `a=2/13`, so its
overlap is at most `min(a,s)`.  Its excess over mean density is at most

```text
max_(0<=s<=1)(min(a,s)-as)=a(1-a)=22/169.               (7)
```

Divide by `u`.  Endpoints have measure zero, so the closed danger convention
does not affect (6).  This is the one-component instance of the discrepancy
lemma independently used in THM-816. ∎

If `m` danger combs, all of speeds at least `y`, cover an interval of length
`L`, summing (6) gives

```text
y <= 22m/[13(13-2m)L].                                  (8)
```

This inequality makes the four-lift chart recursively finite without using
the collar cycle or a lower-runner LRC theorem.

## 3. The exact safe-component ladder

Numerically order the four replacements:

```text
x<v<w<z.                                                 (9)
```

They are distinct because their residues modulo `13` are distinct.  Put
`P=[12] minus R`.

### 3.1 The first speed satisfies `x<=105`

The exact census of all `C(12,4)=495` eight-speed cores gives

```text
min_R L(P)=1/78.                                        (10)
```

The minimum is unique at

```text
R=(3,4,5,12),       P=(1,2,6,7,8,9,10,11),
component=(4/13,25/78).                                 (11)
```

If `B` were tight, the four replacement combs would cover this component.
Apply (8) with `m=4` and `y=x`:

```text
x <= floor(88/[13*5*(1/78)])=105.                      (12)
```

### 3.2 The second speed satisfies `v<=118`

Range over every label of `x` and every proper lift `x<=105`.  There are
exactly `14,025` nine-speed cores `P union {x}`.  Their exact minimum is

```text
min L(P union {x})=7/1144,                              (13)
```

uniquely at

```text
R=(3,4,5,12),       x=51 of label 12,
component=(27/104,38/143).                              (14)
```

The remaining three combs cover that component under tightness.  Formula
(8) with `m=3` gives

```text
v <= floor(66/[13*7*(7/1144)])=118.                    (15)
```

### 3.3 The third speed satisfies `w<=83`

Over every ordered labelled pair `x<v` satisfying (12)--(15), the exact
`191,070`-row ten-core census gives

```text
min L(P union {x,v})=23/5434,                           (16)
```

uniquely at

```text
R=(4,6,10,12),      x=95 of label 4,  v=110 of label 6,
component=(196/1235,233/1430).                          (17)
```

The last two combs and (8) first give

```text
w<=floor(44/[13*9*(23/5434)])=88.                      (18)
```

Consequently `x<v<88`.  Repeating the exact ten-core minimum only on that
necessary `98,145`-row domain sharpens it to

```text
min L(P union {x,v})=1/221,                             (19)
```

uniquely at

```text
R=(3,4,5,12),       x=51 of label 12,  v=57 of label 5,
component=(4/13,69/221).                                (20)
```

Thus

```text
w<=floor(44/[13*9*(1/221)])=83.                        (21)
```

The minimizer in (20) remains in the resulting smaller domain, so this is a
fixed-point refinement rather than an iterative heuristic.

### 3.4 The final speed satisfies `z<=50`

The complete `313,965`-row census of ordered labelled triples

```text
x<v<w<=83                                                (22)
```

gives

```text
min L(P union {x,v,w})=1/325.                           (23)
```

There are exactly three minimizers:

```text
R=(4,6,10,12) in every row,

(x_label,x; v_label,v; w_label,w; component)
=(4,30; 12,38; 10,75; (157/975,32/195)),
=(6,19;  4,30; 10,75; (11/65,56/325)),
=(10,23; 12,25; 4,30; (2/13,51/325)).                  (24)
```

The last closed danger comb must contain the connected component in (23),
so it must contain it in a single tooth of length `2/(13z)`.  Hence

```text
z<=floor((2/13)/(1/325))=50.                            (25)
```

Now (9) implies `x<v<w<51`.  On the resulting `49,005` eleven-core rows the
minimum remains `1/325`, uniquely at

```text
R=(4,6,10,12),      x=23 of label 10,
v=25 of label 12,   w=30 of label 4,
component=(2/13,51/325).                                (26)
```

Thus (25) is again stable under its own consequence.

The six nested censuses contain `666,705` exact rows in total.  Each domain
is defined only by a previously proved necessary bound and the numerical
ordering (9); no sampled height cutoff enters the reduction.

## 4. Exact closure of the `z<=50` box

For labels `1,...,11`, the allowed positive lifts at most `50` are
`r+13,r+26,r+39`; label `12` has `25,38`.  Therefore the full residual has

```text
C(11,4)*3^4 + C(11,3)*2*3^3 = 35,640                  (27)
```

rows.

For a row, form every component `(l,h)` of the eleven-speed core
`P union {x,v,w}`.  Put

```text
c=(l+h)/2,                       eta=(h-l)/2.            (28)
```

Because `D_z` is closed,

```text
(l,h) subset D_z
 iff ||zc||+z eta<=1/13.                                (29)
```

If `c=C/N` and `eta=H/N`, (29) is the integer comparison

```text
13*(dist(zC,N)+zH)<=N.                                 (30)
```

The verifier computes every component and tests (30).  The result is

```text
rows                         35,640
rows satisfying (29) on every component   0.            (31)
```

The digest of every row and its first failed component is

```text
118e9413d8e9b4daf3a240b96a6f70d4760ae0771485cec44a1cdb3af8f704cf.
```

The closest first-failure record has

```text
R=(1,9,11,12), labelled lifts=(40,22,37,25),
(x,v,w,z)=(22,25,37,40),
positive cleared surplus=3718/137566.                   (32)
```

Thus every row has a nonempty strict-safe interval, proving Part A.
Because the packet is a complete nonzero residue transversal modulo `13`, it
already has `M(B)>=1/13` at thirteenth-grid points; the strict interval makes
the inequality strict. ∎

### Independent endpoint-cell replay

As a separate final-box audit, the verifier does not choose a last speed or
use (29).  For each of all `35,640` full twelve-speed packets, it sorts every
threshold endpoint from (4), tests exact midpoints of the resulting open
cells, and finds a strict witness.  This replay again has zero failures and
has digest

```text
82823bf934a438e4dcbcff2724f322da095ae23b54a333a55e91e8eb54face8c.
```

The smallest margin of its deterministically first selected witnesses is
`1/1274`, on `R=(1,3,4,10)` with labelled lifts `(14,16,17,49)` and witness
`99/2548`.  This statistic is not asserted to be the packet's exact maximum;
the audit only needs one strict endpoint-cell witness per row.

### Independent collar/doubling proof

A separately derived proof is retained in the companion listed above.  It
starts from the same owner collar but uses it as the main reduction rather
than as telemetry.  If every lift exceeds `24`, the handoff digraph forces the
unique cyclic band word `(2,2,2,5)` and missing-label packet
`a{1,2,4,8}`.  Ratio bands bound the spread by four, and an eight-runner safe
interval bounds the least lift by `245`.  Exact containment rejects all
`626,962` rows in this all-large box.

Outside that branch the collar forces a least lift `14<=x<=24` and the
recursive doubling box

```text
x<v<=2x,             v<w<=2v,             w<z<=2w.
```

The independent C++ verifier rejects all `141,773` rows in this anchored box,
for `768,735` exact rows and zero tight packets in total.  Its branch digests
are

```text
all-large  27c45d31f19370b8b3c30e79f378b5b3ed9b1f9538062ac2f80e7dd056a6a64e
anchored   07594ab0e69196583fdf667b4d54c8a048a1b4d2b2a87924d26a7da4d8bc7542
```

and the stored-output SHA-256 is
`f098acc358f534f4edf75e1affcaa03ff0bf9cda83f058d5fc86cfc984d2dca0`.
This proof shares the exact component-containment terminal predicate but has a
different unbounded-to-finite reduction, providing a genuine independent
cross-check of Part A.

## 5. The collar four-cycle is a structural sidecar

The THM-806 owner collar still explains why Hamming four is the first new
shallow combinatorial arity.  At the inverse thirteenth for a missing owner,
the retained eight-speed core has a uniform left-safe collar of length
`1/156`.  If all four replacements exceed `24`, tightness forces a cross
handoff at every own-tooth exit.  Directed two- and three-cycles are
impossible by the exact THM-806 ratio argument, so a directed four-cycle must
occur.

For an arrow with speed ratio `lambda` and residue ratio `z`, the half-open
handoff condition is

```text
-1<z-2lambda-13m<=1.                                   (33)
```

Let `a=z-13m` be the positive integer band centre.  Distinct labels give
`a>=2`.  Around a four-cycle, `product(lambda)=1`; hence

```text
product(a_i-1)<=16<product(a_i+1),
product(a_i)=1 (mod 13).                                (34)
```

Some ratio is strictly below one, so some `a_i=2`; then the left inequality
in (34) gives every `a_i<=17`.  The exact finite integer lemma on
`{2,...,17}^4` has only the four rotations of

```text
(2,2,2,5).                                               (35)
```

Therefore the four missing labels must be

```text
R=a{1,2,4,8} in F_13^*.                                 (36)
```

Off these twelve scaled label packets, tightness would force a height-one
replacement `u<=24`.  This classification is sharp only as a collar
predicate: the row

```text
R={1,2,4,8},       (u_1,u_2,u_4,u_8)=(79,54,30,34)     (37)
```

realizes the live cycle

```text
1 -> 8 -> 4 -> 2 -> 1
```

but is loose with exact `M=3/19`, witnessed at `1/19`.  The safe-component
ladder, not (36), is what closes the exceptional cosets.

## 6. Tournament Analysis and assumption challenge

At the collar layer, vertices are the four **owner-exit obligations**.  The
pair observable is the exact half-open predicate (33); silent pairs follow
the numerical tie Hamiltonian path.  Reversing the four live arrows is the
switch.  For (37), both gauges have

```text
score sequence            (2,2,1,1)
directed triangles        2
SCC sizes                 (4)
Hamiltonian paths         5
edge flips                4.                            (38)
```

The strongly connected shadow detects the quartic cycle but cannot decide
looseness: its method-limit row has margin far above `1/13`.

The assumption that vertices should remain runners or residues was therefore
challenged at the recursive stage.  Runner vertices lose component length;
residues lose lift magnitude; gap or fixed-section vertices lose moving comb
scale; tooth vertices lose which core-safe component they must cover; Fourier
modes retain average density but not the connected interval used in (25).
The predicate-preserving carrier is instead the bipartite incidence between
strict-safe **components** and remaining danger **combs**, decorated by exact
width and endpoint ownership.  It is not naturally antisymmetric, so forcing
it into a tournament destroys the cover predicate.  Tournament Analysis is
faithful telemetry at the collar sidecar and deliberately not substituted
for the component proof.

## Reproduction

```bash
python3 04-computation/lrc13_h4_scale_one_component_ladder_codex_S10.py
```

The replay uses only Python integers and `fractions.Fraction`.  It proves the
sharp discrepancy constants, checks all six nested Hamming-four census minima
and their digests, verifies the first radius-five/six component minima and
caps, closes the final containment box, runs the independent endpoint-cell
audit, verifies the four-cycle band lemma and method-limit maximum, and records
both tournament gauges.
