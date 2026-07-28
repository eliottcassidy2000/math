---
id: THM-2648
title: "Two-rainbow thinning full-carry cover and noncanonical selector boundary"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  For every
  pair S=F_13\A, T=F_13\B with |A|=|B|=2, the complete two-edge relation
  S x T contains two eleven-edge rainbow matchings whose carry-colour sets
  have disjoint two-point holes and hence cover all thirteen carries.  If the
  two missing-pair steps differ up to sign, both matchings are the two affine
  bijections A->B.  In the step-matched/high-energy THM-2645 class, exactly
  one affine orientation is rainbow; an explicit three-cycle flip supplies
  the second disjoint-hole chart.  Across all 6,084 relation pairs there are
  exactly 11,154 affine rainbow charts, while the displayed atlas chooses
  1,014 nonlinear charts, one for each matched-step pair.  Every two-chart
  cover has colour multiplicities 1 on four carries and 2 on nine.  Both its
  even and chart-swap sectors retain all twelve charged carry characters,
  with energies 9/169 and 1/13 on C_2 x C_13.  Two charts are minimal, and
  the matched-wall union is sharply only fourteen edges.  This is an abstract
  rainbow thinning, not a physical LRC
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
script_sha256: 3a81c3ee39738ac695dd505f15752c4bd8bf561185af57881cc93330aaf2b4f5
output_sha256: e9ee94f42936a9f3d5a2075d97f1f361998bfdba464b1b5c185cc65705e36756
hash_basis: LF-normalized bytes
---

# THM-2648 -- two rainbow charts thin every eleven-sheet pair

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

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

### 3.1 Parallel normal form: a ternary flip

Take

```text
A=B={0,1},                    S=T={2,3,...,12}.            (11)
```

The affine rainbow is `f(x)=x`; its colour holes are `{0,2}`.  On the ordered
domain `x=2,...,12`, rotate the three targets over `x=6,7,8` and leave the
other eight edges fixed:

```text
(f(x))=(2,3,4,5,7,8,6,9,10,11,12).                       (12)
```

It is a permutation of `T`, its eleven sums are distinct, and

```text
H_f={3,12}.                                               (13)
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
(f(x))=(1,2,3,4,6,7,5,8,9,10,11).                        (15)
```

This is a permutation of `T`, has eleven distinct sums, and has

```text
H_f={2,11},                                               (16)
```

again disjoint from the affine holes.  Both (12) and (15) are genuinely
nonlinear: the affine map determined by their first two edges fails on later
edges.  Each is an alternating six-cycle flip between two perfect matchings:
three old edges are replaced by three new edges.  Of the three old carry
colours, one is retained and the other two are exchanged for the affine
chart's two missing colours.  This is a literal `3`-move repairing a
`2`-point defect.

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

### 3.4 Sharpness and the local `C_2*C_3` frame

No second rainbow chart can differ from the matched-wall affine chart in only
two edges.  A two-edge change in a bijection is a transposition.  In the
parallel form, transposing the targets at sources `x,y` makes both new carry
colours equal `x+y`; in the antiparallel form both equal `x+y-1`.  Either way
the result is not rainbow.  One moved edge is impossible for a permutation.
The charts (12) and (15) move exactly three edges, so they are sharp.  They
share eight of eleven edges with the affine chart and their union has

```text
11+11-8=14 edges.                                         (18a)
```

On the active triples the affine chart and the two inverse three-cycle
repairs are exactly the three one-factors of `K_(3,3)`: after cyclically
labelling the triples by `C_3`, they are `v=u+r` for `r=0,+1,-1`.  The field
reflection fixes `r=0` and swaps `r=+-1`, producing the local symmetry

```text
C_3 semidirect C_2 = S_3.                                 (18b)
```

Keeping one nonlinear repair gives the sharp two-chart/fourteen-edge cover
but chooses an orientation.  Keeping both gives a reflection-stable
three-chart atlas: the eight outside edges are common and the active block is
all nine edges of `K_(3,3)`, for seventeen union edges.  Its carry
multiplicity is `1^2 2^2 3^9`, with centered `C_13` energy `94/169`; each
nontrivial `C_3` chart character retains every nonzero carry mode and has
energy `4/117`.  This is an exact co-occurrence of the binary reflection and
ternary matching grammars, but it still requires an occurrence-level chart
label.

## 4. A complete two-chart atlas and its chart-type census

Among the `78^2=6,084` ordered relation pairs, THM-2645's step census is

```text
step-distinct: 5,070,          step-matched: 1,014.       (19)
```

The construction above gives

```text
affine rainbow charts (all):       2*5,070+1,014=11,154,
nonlinear charts chosen in this atlas:              1,014.   (20)
```

The second line is the size of the **selected atlas**, not a census of all
nonlinear rainbow matchings.  Already in the parallel normal form the two
target vectors

```text
(2,4,11,5,6,7,8,9,3,10,12),
(2,10,3,5,6,7,8,9,11,4,12)                              (20a)
```

are distinct nonlinear rainbow bijections with the same hole pair `{6,9}`.
Thus even the complete carry profile does not recover edge occurrences.

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

Let `w(c)` count the number of the two charts carrying colour `c`.  Then
`w` has the universal rational shape `1^4 2^9`.  A rational vector on
`C_13` with one vanishing nonzero Fourier mode is constant, by irreducibility
of `Phi_13`; hence this nonconstant profile has

```text
what(k)!=0                         for every k!=0.         (23)
```

With the normalized DFT, its exact charged energy is

```text
sum_(k!=0)|what(k)|^2
 =13^(-1)sum_c|w(c)-22/13|^2
 =36/169,                                                 (24)
```

so some chart-incidence mode has square at least `3/169`.  Thus the rainbow
thinning preserves THM-2645's full charged spectrum and regularizes both
dense energy classes to the lower value.

There is a stronger statement if the occurrence bit is retained.  Let
`h_0,h_1` be the two chart-colour indicators and put

```text
w=h_0+h_1,                    d=h_0-h_1.                 (24a)
```

The signed profile `d` has two entries `+1`, two entries `-1`, and nine zero
entries.  It too is rational and nonconstant, so `Phi_13` forces every
nonzero Fourier mode of `d` to survive.  For the normalized transform on
`C_2 x C_13`,

```text
Hhat(epsilon,k)
 =1/26 sum_(j,c) h_j(c)(-1)^(epsilon j) zeta^(-kc),
sum_(k!=0)|Hhat(0,k)|^2=9/169,
sum_(k!=0)|Hhat(1,k)|^2=1/13.                            (24b)
```

Thus both the even and chart-swap charged sectors retain all twelve carry
characters.  Swapping the two charts changes only the sign of the odd sector,
but that sector exists only after adding an occurrence-level `C_2` chart bit.
It is not present in the current physical LRC packet.

The two charts need not be edge-disjoint.  They share exactly one selected
edge in the `5,070` step-distinct pairs and exactly eight in the `1,014`
matched pairs.  Their union therefore retains respectively

```text
21 of 121 edges       or       14 of 121 edges.           (25)
```

Two charts are minimal: one eleven-edge rainbow matching has only eleven
colours and cannot cover thirteen.  Thus (21) is a sharp two-chart thinning
of the full `121`-edge relation, while (18a) is sharp on the matched wall.

## 5. Holotopy meaning and noncanonical boundary

Each chart is thin over its eleven-colour domain.  Together the charts form a
two-open-set analogue of a relation-valued local system: they cover the carry
torsor and have a nine-colour overlap on which two branches coexist.  This
explains exactly how a maximally saturated relation can contain private local
branches without having one global holonomy permutation.

The construction does not contradict THM-2644's purity/return no-go for the
full multiplicity.  By (25) it first discards `100` or `107` of the `121`
increment pairs.  Nor is the selector canonical.  The formulas label the two
affine maps after orienting the missing pairs, although their unordered pair
is intrinsic.  On the matched wall the displayed templates choose one of the
two inverse three-cycle repairs.  Other nonlinear charts, including (20a),
can have the same holes.  Thus the physical noncanonicity is occurrence
selection, not a forced breaking of residual affine symmetry.

Perfect-matching inherited colouring can isolate a chosen chart only after
edge occurrences are duplicated and labelled by the chart bit: colour both
endpoints of an occurrence `(s,t,j)` by `j`, and a monochromatic matching
selects chart `j`.  The current LRC model contains no such chart-coloured
same-base function-valued kernel.  Since the charts overlap, one unduplicated
edge label cannot encode both occurrences.

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
   `1,014` degenerate affine controls, and the constructed-atlas chart-type
   census (20); and
5. proves disjoint holes, full thirteen-colour union, and (21) on every pair,
   recovering the global census (22), verifies both `w` and `d` in every
   charged sector and the energies (24)/(24b), proves the distance-two
   hostile and distance-three repair census in both normal forms, and checks
   the exact one/eight-edge chart overlaps and retained-edge counts (25).

The script and stored output have the LF-normalized SHA-256 hashes declared
in the frontmatter.

One independent immutable audit rederived both affine charts, the two
matched-wall nonlinear templates, their transported hole formulas, the
complete `11,154` affine census plus the chosen `1,014`-chart nonlinear
atlas, the one/eight-edge overlap split, the `21/14` retained-edge counts,
and the abstract-selector versus physical-packet boundary.  It also replayed
the exact companion in normal and optimized modes against the stored
transcript.

A second independent immutable audit rederived the strengthened incidence
claim from `w=1^4 2^9`: rational cyclotomic irreducibility forces every one
of the twelve charged modes to survive, while normalized Parseval gives
exactly `36/169` centered energy.  It retained the occurrence bit and proved
the chart-swap profile has every charged mode and energy `1/13` in the
normalized `C_2 x C_13` transform.  The same audit caught and repaired the
false reading of `1,014` as the total nonlinear census: it is only the chosen
atlas size, as witnessed by (20a).  Both audits independently reproduced the
character guards and the declared LF-normalized hashes.

QED.
