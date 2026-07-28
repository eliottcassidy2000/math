---
id: THM-2648
title: "Two-rainbow thinning full-carry cover and noncanonical selector boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.  For every
  pair S=F_13\A, T=F_13\B with |A|=|B|=2, the complete two-edge relation
  S x T contains two eleven-edge rainbow matchings whose carry-colour sets
  have disjoint two-point holes and hence cover all thirteen carries.  If the
  two missing-pair steps differ up to sign, both matchings are the two affine
  bijections A->B.  In the step-matched/high-energy THM-2645 class, exactly
  one affine orientation is rainbow; an explicit nonlinear normal-form
  matching supplies the second disjoint-hole chart.  Exhaustively there are
  11,154 affine and 1,014 nonlinear charts, and every two-chart cover has
  colour multiplicities 1 on four carries and 2 on nine.  Two charts are
  minimal.  This is an abstract rainbow thinning, not a physical LRC
  selector: current eleven-sheet rows do not provide a same-base positive
  product relation or a lawful measurable restriction to the selected
  edges.
source: deep-energy-audit-2026-07-28-two-rainbow-thinning
depends_on:
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
  - THM-2645-eleven-sheet-multiplicity-full-character-spectrum-and-energy-split
related:
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - HYP-2233-missed-problem-frontier-carrier-atlas
script: 04-computation/lrc14_two_rainbow_full_carry_cover_thm2648.py
output: 05-knowledge/results/lrc14_two_rainbow_full_carry_cover_thm2648.out
script_sha256: 963bc4996deaf3eeec5b78268cfbf592736b35f1d8d7b05dd8410df70282bd1c
output_sha256: 3438805ecf31bc835f9f7bfcd487450f14f748e030011c13d6426e732069c499
hash_basis: LF-normalized bytes
---

# THM-2648 -- two rainbow charts thin every eleven-sheet pair

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.**

THM-2642 shows that an eleven-by-eleven relation pair is maximally thick at
the level of carry support.  THM-2645 shows that its exact multiplicity still
has every charged character, but is far too impure for a fixed-branch gate.
The missing operation is thinning.  A rainbow-transversal mechanism, an
underused lane flagged abstractly in HYP-2233, gives an exact answer here:
two thin charts suffice to cover every carry.

## 1. Rainbow matchings

Let

```text
G=F_13,
A={a_0,a_1},       B={b_0,b_1},
S=G\A,             T=G\B.                               (1)
```

The abstract two-edge relation is the complete product `S x T`; an edge
`(s,t)` has carry colour

```text
chi(s,t)=s+t.                                             (2)
```

A **rainbow matching** is the graph of a bijection

```text
f:S->T                                                     (3)
```

for which `s->s+f(s)` is injective.  It has eleven edges, eleven distinct
carry colours, and therefore a two-point hole set

```text
H_f=G\{s+f(s):s in S}.                                    (4)
```

Two rainbow matchings form a full-carry cover exactly when their hole sets
are disjoint.

## 2. The two affine charts off the matched-step wall

Orient `A,B` as in (1), set

```text
d_A=a_1-a_0,             d_B=b_1-b_0,
t=d_B/d_A.                                                 (5)
```

There are exactly two affine bijections of `G` carrying `A` to `B`:

```text
f_+(x)=b_0+t(x-a_0),
f_-(x)=b_1-t(x-a_0).                                      (6)
```

They restrict to bijections `S->T`.  Their colour maps are

```text
x+f_+(x)=(1+t)x+b_0-t a_0,
x+f_-(x)=(1-t)x+b_1+t a_0.                               (7)
```

Thus `f_+` is rainbow iff `t!=-1`, while `f_-` is rainbow iff `t!=1`.
Whenever

```text
t notin {1,-1},                                           (8)
```

both are rainbow.  Because each colour map in (7) is then a bijection on all
of `G`, restricting away from `A` gives the exact hole formulas

```text
H_+={a_0+b_0,a_1+b_1},
H_-={a_0+b_1,a_1+b_0}.                                   (9)
```

The four entries in (9) are distinct under (8): every cross-equality would
force either `a_0=a_1`, `b_0=b_1`, or `d_B=+-d_A`.  Hence

```text
H_+ intersect H_-=empty                                  (10)
```

and the two affine matchings cover every carry.

This argument works over every odd field.  The remaining matched-step wall
`t=+-1` is exactly THM-2645's higher-energy class.

## 3. The matched wall needs one nonlinear chart

On `t=1`, `f_+` is rainbow and `f_-` has constant colour.  On `t=-1`, the
roles reverse.  The affine failure is real, but its two missing carries are
not structural.  Over `F_13` there are explicit nonlinear second charts.

### 3.1 Parallel normal form

Take

```text
A=B={0,1},                    S=T={2,3,...,12}.            (11)
```

The affine rainbow is `f(x)=x`; its colour holes are `{0,2}`.  On the ordered
domain `x=2,...,12`, define the second target vector

```text
(f(x))=(11,12,5,2,6,3,8,9,4,10,7).                      (12)
```

It is a permutation of `T`, its eleven sums are distinct, and

```text
H_f={4,11}.                                               (13)
```

The holes (13) are disjoint from `{0,2}`.

### 3.2 Antiparallel normal form

Take

```text
A={0,1},             B={0,12},
S={2,...,12},        T={1,...,11}.                        (14)
```

The affine rainbow is `f(x)=x-1`; its holes are `{1,12}`.  On the same
ordered source domain, define

```text
(f(x))=(10,11,4,1,5,2,7,8,3,9,6).                       (15)
```

This is a permutation of `T`, has eleven distinct sums, and has

```text
H_f={3,10},                                               (16)
```

again disjoint from the affine holes.  Both (12) and (15) are genuinely
nonlinear: the affine map determined by their first two edges fails on later
edges.

### 3.3 Transport to arbitrary origins

For any matched pair, orient it so `d_B=d_A` or `d_B=-d_A`.  The change of
coordinates

```text
x=a_0+d_A u,                   y=b_0+d_A v               (17)
```

sends it to (11) or (14).  It sends carry colour `u+v` to

```text
a_0+b_0+d_A(u+v),                                      (18)
```

an affine bijection of `G`.  Therefore it preserves rainbow injectivity and
disjointness of hole sets.  Transporting (12) or (15) supplies the required
second chart for every matched pair.

## 4. Complete atlas and sharp chart count

Among the `78^2=6,084` ordered relation pairs, THM-2645's step census is

```text
step-distinct: 5,070,          step-matched: 1,014.       (19)
```

The construction above gives

```text
affine rainbow charts:       2*5,070+1,014=11,154,
nonlinear rainbow charts:                       1,014.   (20)
```

For every pair the two hole sets are disjoint.  Hence the chart-colour union
is all of `G`; exactly the four hole colours occur in one chart and the other
nine colours occur in both:

```text
cover multiplicity:       1^4 2^9.                       (21)
```

Across the complete bank this is the exact census

```text
multiplicity one:     6,084*4=24,336,
multiplicity two:     6,084*9=54,756.                    (22)
```

The two charts need not be edge-disjoint.  They share exactly one selected
edge in the `5,070` step-distinct pairs and exactly three in the `1,014`
matched pairs.  Their union therefore retains respectively

```text
21 of 121 edges       or       19 of 121 edges.           (23)
```

Two charts are minimal: one eleven-edge rainbow matching has only eleven
colours and cannot cover thirteen.  Thus (21) is a sharp two-chart thinning
of the full `121`-edge relation.

## 5. Holotopy meaning and noncanonical boundary

Each chart is thin over its eleven-colour domain.  Together the charts form a
two-open-set analogue of a relation-valued local system: they cover the carry
torsor and have a nine-colour overlap on which two branches coexist.  This
explains exactly how a maximally saturated relation can contain private local
branches without having one global holonomy permutation.

The construction does not contradict THM-2644's purity/return no-go for the
full multiplicity.  By (23) it first discards `100` or `102` of the `121`
increment pairs.  Nor is the selector canonical.  The affine charts require the labelled missing
pairs; on the matched wall the nonlinear templates additionally choose a
normalization and break the residual affine symmetry.  Different choices can
give different thin branches with the same original relation.

Most importantly, current LRC rows do not instantiate the product `S x T` as
a same-base positive physical transition table.  A static coefficient row or
a signed tomographic character does not authorize selecting the eleven
specific physical pair events in (3), (12), or (15).  The physical
predecessor carry is target-neutral, and the formal affine lift changes other
factor phases.  A lawful application needs a positive joint packet plus a
measurable selector preserving its word, clock, owner, and endpoint current.

No such selector is proved here.  No scalar row is excluded and the LRC
ledger remains `165`.

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_two_rainbow_full_carry_cover_thm2648.py
python3 -O 04-computation/lrc14_two_rainbow_full_carry_cover_thm2648.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_two_rainbow_full_carry_cover_thm2648.out.
```

The dependency-free referee uses explicit optimization-safe guards.  It

1. verifies both nonlinear normal-form templates, their target bijections,
   colour injectivity, holes, and nonlinearity;
2. constructs both charts on all `6,084` relation pairs;
3. checks all `12,168` rainbow matchings and `133,848` selected edges;
4. verifies every affine and transported nonlinear hole formula, all
   `1,014` degenerate affine controls, and the chart census (20); and
5. proves disjoint holes, full thirteen-colour union, and (21) on every pair,
   recovering the global census (22), and checks the exact one/three-edge
   chart overlaps and retained-edge counts (23).

The script and stored output have the LF-normalized SHA-256 hashes declared
in the frontmatter.

QED (candidate; independent audit pending).
