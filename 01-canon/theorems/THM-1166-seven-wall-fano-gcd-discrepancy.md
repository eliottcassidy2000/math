---
id: THM-1166
title: Seven-wall quadratic/Fano gcd discrepancy -- sharp global pair/tree credits, common-dilate closure, and exact covered-needle gcd-error laws
status: PROVED complementary S75/S76 packages plus sharp three-comb addendum.  At radius 1/14 every pair has overlap at least 1/91, every triple has pair sum at least 51/1183 and an edge at least 1/63, and every union of three distinct danger combs has measure at most 36/91, with equality only at scaled/permuted (1,12,13).  The original argument gives every seven-packet a 110/1183 global tree and every seven-comb packet global uncovered mass at least 1/12.  THM-1221 strictly supersedes those two seven-packet global constants by the strict-spectrum values 15/154.  The Fano, forest, and density-weighted gcd-error necessities here remain live; neither theorem alone supplies uniform localization, crown collapse, or LRC(14)
source: codex-2026-07-18-S75/S76; codex-2026-07-24-sharp-three-comb-union
depends_on: [LEM-042, LEM-043, THM-965, THM-1153, THM-1155]
related: [THM-576, THM-856, THM-1025, THM-1156, THM-1221, THM-1226, THM-1234, THM-2137, THM-2168, HYP-7678, HYP-7870]
script:
  - 04-computation/lrc14_seven_wall_fano_gcd_codex_20260718.py
  - 04-computation/lrc14_seven_wall_pair_tree_referee_codex_S76.py
  - 04-computation/lrc14_three_comb_union_referee_thm1166.py
output:
  - 05-knowledge/results/lrc14_seven_wall_fano_gcd_codex_20260718.out
  - 05-knowledge/results/lrc14_seven_wall_pair_tree_referee_codex_S76.out
  - 05-knowledge/results/lrc14_three_comb_union_referee_thm1166.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCSevenWallFanoGCD.lean
referee:
  - 04-computation/lrc14_fano_gcd_discrepancy_referee_codex_20260718.py
  - 05-knowledge/results/lrc14_fano_gcd_discrepancy_referee_codex_20260718.out
---

# THM-1166 -- seven-wall Fano/gcd discrepancy

Work at the actual LRC(14) radius.  For a positive integer `s`, put

```text
D_s={t in R/Z: ||st||<1/14}.
```

Let `s1,...,s7` be distinct positive speeds, let `I` be an interval of length
`L`, and suppose their union covers `I`.  In the LRC application, the
protected-interval construction and the already-known lower-case LRC(7)
input supply the separate hypothesis

```text
L>=1/(7m),                 m=max(retained core).         (1)
```

THM-1153's first-order coefficient is zero here.  The point of this theorem is
that the exact arithmetic pair mass is nevertheless large enough to impose
scale constraints once periods and Fano incidence are retained.

## 1. Exact pair and three-speed floors

Write

```text
rho(a,b)=measure(D_a intersect D_b).
```

After dividing by `gcd(a,b)`, take coprime `a<=b`, and put
`F(r)=r(14-r)` for the representative `0<=r<14`.  The folded form of
LEM-042's trapezoid identity is

```text
rho(a,b)=
 [4ab+F((a+b) mod 14)-F((b-a) mod 14)]/(196ab).         (2)
```

Since the folded difference is at least `-49`,

```text
rho(a,b)>=1/49-1/(4ab).                                 (3)
```

For `ab>=27`, the right side is greater than `1/91`.  The complete 36-row
coprime bank `ab<=26` has unique minimum

```text
rho(1,13)=1/91.                                         (4)
```

Thus (4) is the universal pair floor, with equality only at reduced ratio
`1:13`.

There is a stronger collective fact.

> **Three-speed lemma.**  For any three distinct speeds,
>
> ```text
> rho12+rho13+rho23>=1/24.                               (5)
> ```

**Proof.**  Suppose the sum were below `1/24`.  By (4), every individual edge
would then be below

```text
tau=1/24-2/91=43/2184.                                  (6)
```

LEM-043's defect bound

```text
rho(a,b)>=1/49-1/(7b')
```

for reduced `a'<b'` forces `b'<=198`, because

```text
1/49-tau=11/15288,       b'<2184/11<199.                (7)
```

Exact evaluation leaves 106 oriented reduced ratio types below `tau`.  For
ordered speeds `x<y<z`, all three rational ratios

```text
r=y/x,       s=z/x,       s/r=z/y
```

would have to lie in that list.  Its exact quotient closure has 173 triples.
Their minimum pair sum is

```text
11/252=1/24+1/504,                                      (8)
```

attained only by the two ratio words `(9/4,27,12)` and `(12,27,9/4)`.
This contradicts the assumption and proves (5).  The 106-type and 173-triple
enumerations are frozen, rational, and independently recompute (2). ∎

Summing (5) over the 35 triples of seven vertices counts every edge five
times.  Therefore

```text
R:=sum_(i<j)rho(si,sj)>=7/24.                            (9)
```

## 2. The optimal quadratic certificate and the j=4 chamber

At a point `t`, let `C(t)` be the number of active danger combs and put

```text
Q(C)=C-(2/7)binom(C,2)=C(8-C)/7.                        (10)
```

For `C=0,...,7`, its values are

```text
0, 1, 12/7, 15/7, 16/7, 15/7, 12/7, 1.                (11)
```

Thus `1_{C>=1}<=Q(C)`.  Integration gives

```text
measure(union_i D_si)<=1-(2/7)R<=11/12,                 (12)
measure(uncovered)>=2R/7>=1/12.                         (13)
```

The coefficient `2/7` is optimal within the normalized family
`C-alpha*binom(C,2)`: imposing the majorant at `C=7` forces
`7-21alpha>=1`, hence `alpha<=2/7`.

This also has a multiplicity-debt reading.  If `mu_k=measure{C=k}`, the mean
`integral C=1` gives

```text
mu_0=sum_(k>=2)(k-1)mu_k,
R=sum_(k>=2)binom(k,2)mu_k<=(7/2)mu_0.                  (14)
```

Equality in the last comparison can put nonprivate debt only at depth seven.
On a covered chamber,

```text
Q(C)-1=(C-1)(7-C)/7.                                   (15)
```

The slack vanishes at private depth one and full depth seven and peaks exactly
at `C=4`.  Hence the open `j=4` flood tail is not an unrelated side problem:
it is the central chamber of the unique symmetric quadratic certificate whose
two boundary depths are tight.

## 3. Common-dilate branch: a genuine positive seventh crown

Suppose a common integer `G` divides all seven deleted speeds.  Their union is
`1/G`-periodic.  By (13), every period contains uncovered mass at least
`1/(12G)`.  Hence a connected interval contained in the union has length at
most the total covered mass of one period:

```text
L<=(1-2R/7)/G<=11/(12G).                                (16)
```

Combining (16) with the lower-LRC needle (1) yields the sharper exact form

```text
G/m<=7-2R<=77/12.                                       (17)
```

This crosses the zero-density wall on the common-dilate branch.  It is a
local continuous Mirsky-Newman statement: global pair mass creates a hole in
every common period, and a covered protected needle cannot cross that hole.

## 4. Fano decomposition for mixed triple periods

Fix any labelled Fano plane on the seven comb vertices.  Its seven lines are
triples, every vertex lies on three lines, and every pair lies on exactly one
line.  For a line `ell`, define

```text
c_ell(t)=sum_(i in ell)1_{D_si}(t),
f_ell(t)=c_ell/3-(2/7)binom(c_ell,2),
R_ell=sum_({i,j} subset ell)rho(si,sj),
G_ell=gcd(si:i in ell).                                 (18)
```

The four values of `f_ell` at `c_ell=0,1,2,3` are

```text
0, 1/3, 8/21, 1/7.                                     (19)
```

Fano incidence gives the exact pointwise identity

```text
Q(C)=sum_(ell) f_ell.                                   (20)
```

The function `f_ell` has period `1/G_ell`, maximum `8/21`, and mean

```text
nu_ell=(1-2R_ell)/7.                                    (21)
```

Use the sharp periodic-remainder lemma: if `0<=f<=M` is `P`-periodic with
mean `nu`, then every interval `J` satisfies

```text
integral_J f<=nu|J|+nu(1-nu/M)P.                        (22)
```

Indeed, remove all full periods.  On a remainder of length `r`, its mass is
at most `min(Mr,nu P)`; relative to `nu r`, the maximum excess occurs at
`r=nu P/M` and equals `nu(1-nu/M)P`.

Put

```text
e_ell=nu_ell(1-21nu_ell/8).                             (23)
```

Since the seven combs cover `I`, (10), (20), and (22) give the exact necessary
metric-Fano inequality

```text
(2LR)/7<=sum_(ell) e_ell/G_ell.                         (24)
```

Every line is a three-speed packet, so (5) gives `R_ell>=1/24`.  Consequently
`0<=nu_ell<=11/84<4/21`, and the quadratic in (23) is increasing on this
range.  Hence

```text
e_ell<=11/128.                                          (25)
```

Insert (1), (9), and (25) into (24): for **every** labelled Fano plane,

```text
sum_(ell) m/G_ell>=32/231.                              (26)
```

In particular every labelled plane has a line satisfying

```text
G_ell/m<=1617/32.                                       (27)
```

Equivalently, if all seven triple gcds of some plane are at least `G`, then
`G/m<=1617/32`.  This gives the Fano/`chi_7` probe a concrete metric consumer:
`chi_7` may address or compare the 30 labelled planes, but it cannot discard
their numerical line periods and pair masses.

## 5. Adaptive forest certificate

The same periodic-remainder argument yields, on any interval of length `L`,

```text
measure(I intersect D_s)<=L/7+6/(49s),
measure(I intersect D_si intersect D_sj)
 >=L rho_ij-rho_ij(1-rho_ij)/gcd(si,sj).                (28)
```

For every forest `F` on the seven comb vertices, the pointwise forest-Hunter
inequality and the cover of `I` therefore force

```text
sum_({i,j} in F)
 [L rho_ij-rho_ij(1-rho_ij)/gcd(si,sj)]
 <=(6/49)sum_i 1/si.                                    (29)
```

This is an adaptive exact edge-gcd certificate.  Negative edge lower bounds
may be omitted; a proof search should maximize the positive forest weight,
not preselect a speed-order path.

THM-1226 later combines (29) with its disconnected-strict-spectrum projective
charge bound.  The resulting exact crown is
`H/L>=59625/1485836`; this is sharper than using the elementary `H/7` tree
cap and closes the gcd-period localization estimate on every disconnected
`G_gt={rho>1/63}` branch.

## 6. Why global pair credits still do not localize for free

There is an exact protected-needle guardrail.  Let the retained core contain
six distinct odd speeds at most an odd `m>=11`, and put

```text
I=[1/2-1/(14m),1/2+1/(14m)].                            (30)
```

This interval is safe for every core speed.  Choose an arbitrarily large `a`
and the six further deleted speeds

```text
a+1, 2a+1, 2a+3, 3a+1, 4a+1, 5a+1.                    (31)
```

If `b=qa+c` and both `a,b` were dangerous at `t`, then

```text
||ct||<(q+1)/14.                                        (32)
```

On (30), odd `c` instead gives

```text
||ct||>=1/2-c/(14m).
```

For `(q,c)=(1,1),(2,1),(2,3),(3,1),(4,1),(5,1)`, this contradicts (32).
Thus the six pair intersections with the central speed `a` are all empty on
`I`; they form an empty spanning star while all seven deleted speeds remain in
a bounded ratio cone and `a/m` is arbitrarily large.

This construction does **not** say that the seven combs cover `I`, and so it
does not refute (24), (26), or a cover-conditioned adaptive Fano theorem.  It
does refute any cover-free rule that inherits a fixed global Hunter tree on an
arbitrary protected needle.

## 7. Carrier and verification audit

Runner order, Fano line-mass `R_ell` order, and line-period order all give
transitive tournaments after ties are lexicographically resolved: score
histogram `0,...,6`, no directed cycles, seven singleton SCCs, and one
Hamiltonian path.
They lose which edges share a Fano line, the values `R_ell,G_ell`, and the
endpoint discrepancy in (24).  Candidate vertices were challenged as runners,
pairs, Fano lines, slow gaps, wall events, and proof obligations.  The faithful
current carrier is

```text
labelled Fano incidence hypergraph
  + (R_ell,G_ell)
  + protected-needle endpoint/owner data.               (33)
```

The exact script checks (2) against an independent trapezoid evaluator on
6,360 pairs, the 36-row pair tail, all 106 strict-low ratios and 173 compatible
quotient triples, all 30 labelled Fano planes and all 128 activity patterns,
every displayed rational constant, and the empty-star inequalities.  Normal
and optimized runs are byte-identical.  Frozen hashes are

```text
source  34ad825245f7288bd2253db1472b875f6303c2e604d62473be05312586c0ee0b
output  b2781fdf746bf3a9584508b9f8b0ab4f73369c5da5138cfa77323bf4fa37c802
```

The separate referee script independently checks the sharp period-cell
discrepancy, the full pair-density range `1/91<=rho<=1/14`, all 128 Fano
activation masks, and the weaker density-only selected-line consequence

```text
m[H_ell+(13/28)sum_(pairs in ell)1/gcd]>=3/91.
```

Its source/output hashes are respectively
`51516e7156ada06424af7cf4b594d70a73fb3e13d87db5db1f8170063615a340`
and `da0cea50f3eb8998c36c0e6a1d280262d6a50ca61b23b0c8be84848d2504d3a4`.
The resulting common-dilate constant `3523/36` is a cross-check only and is
strictly dominated by (17).

What has closed is the common-dilate seventh crown and a mixed-period Fano/gcd
necessary law.  What remains open is a uniform way to turn (24) or (29) into a
contradiction for arbitrary seven deleted speeds, together with the broader
`r=5`, `r=6`, crown-collapse, n=12 equality, and LRC(14) obligations.

## S76 independent pair-tree and periodic-discrepancy strengthening

Put

```text
D_s={t in R/Z : ||st||<1/14}
rho(a,b)=measure(D_a intersect D_b).
```

Open or closed danger arcs give the same measures below.  The point of this
theorem is to separate three statements that the seventh-deletion discussion
had blurred:

1. a strong second-order credit exists on the **whole circle**;
2. a covered interval still obeys a useful credit after an explicit gcd
   positioning error is charged;
3. a chosen global edge or tree need not carry any of its mass on an
   arbitrary child interval.

The third fact is a theorem-level obstruction, not a caveat added after the
calculation.

## 1. Exact pair floor, sharp triple bar, and triple pair sum

THM-965 gives, after reducing `(a,b)` to a coprime pair with `a<b`,

```text
rho(a,b)
 =1/49+[fold_14(a+b)-fold_14(b-a)]/(196ab),
fold_14(r)=(r mod 14)(14-(r mod 14)).                    (1)
```

In particular

```text
rho(a,b)>=1/49-1/(4ab).                                  (2)
```

> **Theorem A (pair and triple floors).**
>
> 1. For every two distinct positive integer speeds,
>
>    ```text
>    rho(a,b)>=1/91,                                      (3)
>    ```
>
>    with equality exactly when the reduced pair is `(1,13)`.
> 2. Among every three distinct speeds, some pair satisfies
>
>    ```text
>    rho(a,b)>=1/63.                                      (4)
>    ```
>
>    The constant is sharp: the triple `(1,12,27)` has pair masses
>    `1/84,1/63,1/63`.
> 3. Every three distinct speeds satisfy the aggregate inequality
>
>    ```text
>    rho(a,b)+rho(a,c)+rho(b,c)>=51/1183.                 (4a)
>    ```
>
>    This is sharp, with equality exactly for a scaled and permuted copy of
>    `(1,13,169)`.

**Proof.**  For (3), the right side of (2) is strictly greater than `1/91`
once `ab>=27`.  Exact evaluation of the coprime cells `ab<=26` has unique
minimum `1/91` at `(1,13)`.

For (4), (2) is strictly greater than `1/63` once `ab>=56`.  The complete
strict-low list in the remaining coprime cells is

```text
(1,10):1/70, (1,11):1/77, (1,12):1/84, (1,13):1/91,
(2,11):1/77, (3,10):1/70, (3,11):1/77.                  (5)
```

If `x<y<z` had all three edges strict-low, the reduced ratios `y/x` and
`z/y` would occur in (5), and their product `z/x` would also have to occur
there.  The exact `7*7=49` product table has no such product.  Hence the
strict-low graph is triangle-free.  Direct substitution proves sharpness at
`(1,12,27)`.

We first record the compact intermediate certificate behind (4a).  If the
sum were below `1/24`, then by (3) each individual edge would be strictly
below

```text
1/24-2/91=43/2184.                                      (5a)
```

The tail (2) is strictly above `43/2184` once `ab>=348`.  Exact evaluation
of every coprime cell `ab<=347` leaves precisely `106` reduced channels
below (5a).  Sort the speeds as `x<y<z`; if `r=y/x` and `u=z/y`, then the
three reduced ratios are `r,u,ru`.  Among the `106^2` ordered channel pairs,
exactly `173` have `ru` in the same list.  Their minimum total pair mass is

```text
11/252=1/24+1/504,                                      (5b)
```

attained within this candidate table at consecutive ratios `(12,9/4)` and
`(9/4,12)`.  This proves the intermediate sum floor `1/24`.

The sharp upgrade has a different finite reduction.  For an edge `e`, let
`p_e` be the product of the numerator and denominator of its reduced speed
ratio, and order these products as `p_1<=p_2<=p_3`.  Suppose the triple sum
is at most `S=51/1183`.  Summing (2) gives

```text
S >= 3/49-3/(4p_1),
```

so

```text
p_1 <= 3/[4(3/49-S)]=8281/200<42.                       (5c)
```

The smallest-product edge also has mass at least `1/91`, while each of the
other two edges has product at least `p_2`.  Hence

```text
S >= 1/91+2/49-1/(2p_2),
p_2 <= 1/[2(1/91+2/49-S)]=8281/144<58.                  (5d)
```

Thus `p_1<=41` and `p_2<=57`.  The two corresponding edges share a vertex.
There are `124` oriented nonunit reduced ratios of product at most `41` and
`184` of product at most `57`.  Normalize their common speed to one.  Exact
evaluation of the `124*184-124=22,692` distinct configurations `(1,r,u)`
has minimum

```text
51/1183,                                                 (5e)
```

and its only oriented minimizers are `(r,u)=(1/13,13)` and the reversal.
They are precisely the triple `(1,13,169)` after rescaling.  This proves
(4a), its sharpness, and its equality statement.  All finite evaluations
use integers and rational arithmetic; the companion referee replays them
independently from both the folded and tent formulas. ∎

The intermediate table is complete because (2) proves its tail.  The sharp
table is complete because (5c)--(5d) force two incident edge products into
the displayed finite banks.  Neither is a bounded search standing in for an
unproved asymptotic claim.

## 2. Seven speeds: a global tree and seven Fano obligations

> **Theorem B (global tree floor).**  Every seven distinct positive integer
> speeds admit a spanning tree `T` such that
>
> ```text
> sum_({i,j} in T) rho(s_i,s_j)
>   >=110/1183.                                           (6)
> ```

**Proof.**  Call an edge high when `rho>=1/63`.  Theorem A says the low graph
is triangle-free.  If the high graph is connected, a high spanning tree has
weight at least `6/63=2/21>110/1183`.  Otherwise the low graph contains every
edge between distinct high components.  Three high components would give a
low triangle, so there are exactly two.

Let `m<1/63` be the largest crossing weight and take high spanning trees
inside the two components, five internal edges in total.  For any internal
tree edge of weight `w`, adjoining any vertex of the other component and
using (4a) gives

```text
w+2m>=51/1183.
```

Thus every internal tree edge weighs at least
`max(1/63,51/1183-2m)`, and adding a crossing edge of weight `m` gives

```text
5 max(1/63,51/1183-2m)+m >=110/1183.                     (6a)
```

Indeed the two affine branches meet at
`m=(51/1183-1/63)/2=145/10647`, where their common value is
`110/1183`; each branch increases away from that point. ∎

Now label the seven speeds by the points of any Fano plane.  Each of its
seven lines is a triple and every pair belongs to exactly one line.

> **Corollary B1 (Fano line debt).**  Every Fano line has pair-mass sum at
> least
>
> ```text
> 51/1183,                                                (7)
> ```
>
> and consequently
>
> ```text
> sum_(i<j) rho(s_i,s_j)>=51/169.                         (8)
> ```

Indeed every line is a three-speed packet, so (7) is Theorem A(3).  This
strictly improves the coarser consequence `1/63+2/91=31/819` obtained from
Theorem A(1)--(2) alone.  One may also select a distinct high edge from every
Fano line, giving seven selected edges of total weight at least `1/9`;
connectivity of those seven selected edges is not asserted.  The
spanning-tree conclusion is (6), proved by the component argument rather
than by silently identifying a Fano transversal with a tree.

> **Corollary B2 (global safe-set floor).**  The complement of seven danger
> combs has measure at least
>
> ```text
> 110/1183.                                               (8a)
> ```

For any spanning tree `T` and any nonempty active vertex set `A`, the forest
induced by `A` has at most `|A|-1` edges.  Hence, pointwise,

```text
1_(union_i D_i)
 <=sum_i 1_(D_i)-sum_({i,j} in T)1_(D_i intersect D_j).
```

Integrate this Hunter inequality for the tree from Theorem B.  Since each
danger comb has measure `1/7`, the sum of the seven single masses is one;
the uncovered mass is therefore at least the tree weight `110/1183`.

**Later sharpening.**  THM-1221 resolves the strict/equality spectrum at the
`1/63` wall.  Its first strict rung `5/308` and complete closed-channel
component analysis improve both (6) and (8a) to `15/154`.  The present
`110/1183` proof is retained because its triple-sum/component interpolation is
independent and because the positioned gcd-error inequalities below remain
the relevant localization interface.

THM-1234 supplies a complementary complete-pair sharpening
`R_7>=22/65`.  Feeding that value through this theorem's quadratic and Fano
consumers gives global safe mass `44/455` and, for every Fano labeling,
`sum_ell m/G_ell>=512/3185`.  This scalar floor is weaker than THM-1221 by
`1/1430`, but the retained line-period ledger is stronger structural data.

## 3. The quadratic wall and the positioned gcd error

Let seven combs and a real interval `I` of length `L` be given.  Set

```text
C(t)=#{i:t in D_(s_i)},
H=sum_i 1/s_i,
P=sum_(i<j) measure(I intersect D_(s_i) intersect D_(s_j)).
```

The critical quadratic is

```text
Q(C)=C-(2/7)binom(C,2)=C(8-C)/7.                         (9)
```

For integer `1<=C<=7`, `Q(C)>=1`.  If `c_l` is the active count on Fano
line `l`, then the exact pointwise decomposition is

```text
21Q(C)=sum_(Fano lines l) [7c_l-6 binom(c_l,2)].          (10)
```

The line summands at `c_l=0,1,2,3` are `0,7,8,3`.  Thus the Fano plane is
not decorative: it partitions the complete pair ledger into seven triple
obligations while counting each vertex three times and each pair once.

The needed positioning statement is elementary and sharp at the level of
arbitrary periodic indicators.

> **Lemma C (density-weighted one-interval periodic discrepancy).**  If a
> measurable set `B` is periodic with period `1/g` and has Haar density
> `rho=measure(B)`, then every interval `I` satisfies
>
> ```text
> |measure(I intersect B)-L rho|<=rho(1-rho)/g.           (11)
> ```

**Proof.**  Put `f=1_B-rho`.  Over one period of length `1/g`, the positive
variation of its primitive is

```text
(1-rho) measure(B intersect one period)=rho(1-rho)/g,
```

and the negative variation has the same magnitude.  Remove all complete
periods from `I`.  The remaining integral of `f`, even if it wraps across a
period boundary, lies between minus the total negative variation and the
total positive variation.  This is (11).  A single interval of density
`rho` in each period shows the constant is sharp for arbitrary periodic
indicators. ∎

For `B=D_a intersect D_b`, one may take `g=gcd(a,b)`.  This improves the
generic two-boundary estimate when the carrier is one interval.

> **Theorem C (positioned tree and covered-interval Fano/gcd inequality).**
> Put
>
> ```text
> g_ij=gcd(s_i,s_j),
> rho_ij=rho(s_i,s_j),
> e_ij=rho_ij(1-rho_ij)/g_ij,
> E=sum_(i<j) e_ij.
> ```
>
> For every `I`, some spanning tree `T` on the seven combs satisfies
>
> ```text
> sum_({i,j} in T) measure(I intersect D_(s_i) intersect D_(s_j))
>   >=(102/1183)L-(2/7)E.                                (12)
> ```
>
> If the seven danger combs cover `I`, then
>
> ```text
> (51/169)L <= H/2+E.                                    (13)
> ```
>
> Moreover, for every Fano labeling there is a Fano line `l` such that
>
> ```text
> (51/1183)L <= H/14+sum_({i,j} subset l) e_ij.           (14)
> ```

**Proof.**  Lemma C gives, edge by edge,

```text
measure(I intersect D_(s_i) intersect D_(s_j))
  >=L rho_ij-e_ij.
```

Summing and using (8) gives

```text
P>=(51/169)L-E.                                          (15)
```

A uniformly random labelled spanning tree of `K_7` contains each edge with
probability `2/7`.  Its expected restricted weight is therefore `(2/7)P`,
so the maximum tree has at least that weight.  Combining this with (15)
proves (12).

Now assume coverage.  Fragmentation for one comb gives

```text
measure(I intersect D_(s_i))<=L/7+1/(7s_i).
```

Integrating `Q(C)>=1` over the covered interval gives

```text
L<=sum_i measure(I intersect D_(s_i))-(2/7)P,
```

and hence

```text
P<=H/2.                                                   (16)
```

Combining (15) and (16) proves (13).  The seven Fano line overlap sums add
to `P`, so some line contributes at most `P/7<=H/14`.  That line's global
projective mass is at least `51/1183` by (7); applying Lemma C to its three
edges proves (14). ∎

This is a genuine positive second-order statement at the zero-coefficient
wall.  In the notation of HYP-7678 it proves the requested positioned-tree
shape with explicit coefficient `eta=102/1183`, but with error `(2/7)E`
instead of `C H`.  Consequently HYP-7678 holds on every branch where
`E<=kappa H`, with `C=2kappa/7`.  The general limitation is equally exact:
for pairwise-coprime high speeds `E` remains order one while `H` tends to zero,
so no absolute `kappa` follows and (12) does not close the wall.

> **Corollary C1 (common-dilate scale bound).**  If `s_i=g a_i`, where the
> seven `a_i` are distinct positive integers, then every covered interval
> obeys
>
> ```text
> gL<=1073/1183=0.9070160... .                            (17)
> ```

Indeed the union of the danger combs and its safe complement are periodic
with period `1/g`.  By Corollary B2 the safe set has global measure at least
`110/1183`, hence safe mass at least `110/(1183g)` in every fundamental period.
A connected covered interval lies in one component of that period's danger
union, whose length is at most its total mass `1073/(1183g)`.  This proves
(17).

For comparison, the density-weighted error itself gives a weaker but direct
check: since `rho_ij<=1/7`, one has `e_ij<=6/(49g)`, so
`E<=18/(7g)`; also
`H<=(1/g)sum_(i=1)^7 1/i=363/(140g)`.  Substitution in (13) yields
`gL<=61009/4760`.  Periodicity plus the global tree is what sharpens this to
(17).

If this is THM-1153's protected needle for a six-speed core of maximum `m`,
so `L>=1/(7m)`, then

```text
g/m<=1073/169=6.3491124... .                             (18)
```

This bounds the common scale relative to the core maximum, not the smallest
deleted speed.  It sharply removes large-common-factor protected needles but
does not touch the pairwise-coprime high-speed branch.

THM-1221 later replaces the global safe mass `110/1183` used here by
`15/154`, sharpening (17)--(18) to `gL<=139/154` and `g/m<=139/22`.
THM-1226 further identifies the exact projective-height charge hidden in
`E`; in particular it proves that the pairwise-coprime limitation is an
unbounded obstruction for the gcd-period estimate, even though the underlying
localized overlap can be maximal.

## 4. Why the global tree cannot simply be restricted

For every `N>=1`, triangle inequality on the circle gives

```text
D_N intersect D_(N+1) subset {t:||t||<1/7}.              (19)
```

Therefore

```text
measure([1/7,6/7] intersect D_N intersect D_(N+1))=0      (20)
```

on an interval of length `5/7`, while the global pair mass is positive and
`1/N+1/(N+1)` tends to zero.  Consequently no constants `eta>0,C<infinity`
can make the per-edge statement

```text
measure(I intersect D_a intersect D_b)
  >=eta|I|-C(1/a+1/b)                                    (21)
```

hold for all pairs and intervals.

There is also a tree-level version of the warning.  On
`{N,N+1,...,N+6}`, the six-edge order path consists entirely of consecutive
pairs, so its restricted weight on `[1/7,6/7]` is zero.  Thus a tree chosen
from global projective weights cannot be frozen and inherited by an arbitrary
child interval.

This does **not** refute HYP-7678.  That target may choose its tree after
seeing the interval, and its interval is covered by all seven combs inside a
specific core-safe carrier.  What (20) refutes is the tempting proof that
localizes each edge of (6) independently.  Coverage, adaptive ownership, and
where the locally best tree changes are indispensable hypotheses.

## 5. Tournament and carrier audit

The theorem supplies a natural switch on the seven speed vertices:

```text
pair observable: rho(s_i,s_j),
switch: high iff rho>=1/63,
gauge: orient high edges up speed order and strict-low edges down it.
```

The strict-low undirected graph is triangle-free.  The switched tournament
need not be transitive: on `(1,2,3,4,5,6,10)` its score histogram is
`(1,2,2,3,3,5,5)`, it has seven directed triangles, and its tie Hamiltonian
path is `(10,1,2,3,4,5,6)` in speed labels (39 Hamiltonian paths total).

This tournament preserves the global low/high channel, the Fano line debt,
and the component proof of the MST floor.  It destroys exactly the data that
(20) shows the interval problem consumes: tooth phase, endpoint address,
owner chronology, and the position of the interval inside the gcd-period
cell.  Two better vertex sets are now explicit:

- the seven Fano lines as triple proof obligations, carrying (14); and
- wall-crossing/owner chambers on `I`, carrying the locally optimal graphic
  basis and its gcd-period address.

The challenged assumption is therefore resolved in a precise direction:
runner vertices suffice for the global theorem, while the remaining
localized theorem must use obligation or chamber vertices.

## 6. Verification and honest frontier

The dependency-free Fraction referee checks:

- `7,140` folded-versus-tent pair evaluations;
- the complete `ab<=26` pair-floor and `ab<=55` strict-low cells;
- all `49` oriented low-ratio products;
- the complete `ab<=347` intermediate channel bank (`106` channels and
  `173` compatible quotient triples);
- all `22,692` configurations in the sharp product-`41`/product-`57`
  two-edge reduction;
- all `31,824` seven-subsets of `{1,...,18}` against (6)--(8);
- all `128` active-set instances of the Fano decomposition (10);
- `1,740` exact interval-discrepancy samples; and
- the consecutive-path localization obstruction through scale `1000`.

Normal and optimized Python runs are byte-identical.  Frozen hashes are
stable because every certificate check uses an explicit `require` that raises
`RuntimeError`; the referee contains no optimization-sensitive `assert`.

```text
source  960df28c1c0be54bfd9f25856a9691cf02b5ab6299f16ec95dbc85c09b2f9b71
output  632c31bd652d62042d6b8984f2c3600f14deecc07c2ced7918a576922eef3307
```

What is closed is the global projective/Fano second order and the exact
gcd-positioned necessary inequality for a covered needle.  What remains open
is an error controlled by harmonic endpoint mass rather than the
density-weighted gcd debt `E`, or an owner-incidence argument that uses
coverage to defeat the localization counterexample.  HYP-7678, the seventh
ratio jump, and LRC(14) remain open.

## S77 sharp three-comb union addendum

The pair floor alone does not give the sharp union bound because the triple
intersection enters inclusion--exclusion with the opposite sign.  The exact
pair formula nevertheless leaves only a small finite cleanup.

> **Theorem D (sharp global three-comb union).** For any three distinct
> positive integers `a,b,c`,
>
> ```text
> measure(D_a union D_b union D_c)<=36/91.             (S77.1)
> ```
>
> Equality holds exactly for scaled and permuted copies of `(1,12,13)`.
> Equivalently, every three distinct danger combs have common safe mass at
> least `55/91`.

Put

```text
rho_ab=measure(D_a intersection D_b),
tau=measure(D_a intersection D_b intersection D_c),
Q=rho_ab+rho_ac+rho_bc-tau.                           (S77.2)
```

Since `tau` is at most the smallest of the three pair intersections, `Q` is
at least the weight `T` of a maximum spanning tree on the three speeds:

```text
Q>=T=the sum of the two largest pair overlaps.         (S77.3)
```

For an edge, let `p` be the product of numerator and denominator of its
reduced speed ratio. Order the three edge products as
`p_1<=p_2<=p_3`. The two edges with products `p_2,p_3` form a spanning tree.
The folded pair lower bound in the S76 package therefore gives

```text
T>=2(1/49-1/(4p_2)).                                  (S77.4)
```

If `p_2>=64`, the right side exceeds `3/91` by exactly

```text
2(1/49-1/(4*64))-3/91=3/81536.                        (S77.5)
```

It remains only to consider `p_2<=63`. The two smallest-product edges share
a vertex. Normalize that common speed to one. Their two oriented reduced
ratios belong to the complete bank

```text
R_63={r!=1: numerator(r) denominator(r)<=63},
|R_63|=208.                                           (S77.6)
```

The exact `208*207=43,056` ordered pairs of distinct ratios contain only
fourteen primitive triples with `T<=3/91`. Exact interval intersection and
an independent complete-breakpoint arrangement agree on every triple
intersection. Their values of `Q` are

```text
(1,12,13)       3/91
(1,12,27)       29/756
(1,12,40)       17/420
(1,12,66)       13/308
(1,13,169)      50/1183
(4,9,108)       8/189
(1,10,55)       3/70
(3,10,120)      3/70
(1,12,156)      47/1092
(1,13,156)      47/1092
(12,13,156)     47/1092
(2,11,132)      10/231
(1,12,144)      11/252
(2,11,110)      17/385.                              (S77.7)
```

Thus `Q>=3/91`, with equality only at `(1,12,13)` up to scaling and
permutation. Three-term inclusion--exclusion now gives

```text
measure(D_a union D_b union D_c)
 =3/7-Q
 <=3/7-3/91
 =36/91.                                              (S77.8)
```

This proves Theorem D. It also promotes the `j=3` value and minimizer in
THM-576 from bounded exact search to a global theorem. The aggregate sharp
pair-sum floor `51/1183` in Theorem A is compatible but is not by itself
sufficient for (S77.1): one must retain the triple intersection, and the
maximum-tree reduction is genuinely strict-low only after the fourteen-row
cleanup. For example, `(1,12,27)` has maximum-tree weight `2/63<3/91`, while
its net redundancy is `29/756>3/91`.

The independent referee is optimization-stable and its normal and `-O`
outputs are byte-identical after LF normalization:

```text
source  ad8121cd66929205f4929b03642d9a2a16e01a1469c0caff43b6045575dd1bb6
output  50775cbfcfcbb7cbee7ea6c051ce44e06c984e7f91ef8d6892394f4441ea5eda
```
