---
id: THM-4224
title: "Order-ten minimal plus-min two-bad-owner obstruction"
status: >
  PROVED / REFUTED / FINITE-EXACT / VERIFIED-EXACT / OPEN FIREWALL.
  Refutes B_ij<=End_i End_j+min(End_i,End_j) first at strong order ten and
  resolves every tail fiber of an infinite two-reversal family; proves no
  fixed multiple of either the min term or the two-owner successor-rook count
  can repair the product bound; and proves the family has a strictly positive
  exact HYP-9081 five-copy gap in every order. In the fixed transitive-order
  gauge, all 97 strong labelled order-eleven presentations with exactly two
  reversed arcs satisfy HYP-9081, 25B_ij<=27End_iEnd_j, and the coefficient-one
  successor-rook repair. This sphere is not an isomorphism-class census. At
  order eleven all three assertions remain OPEN outside it; globally the
  coefficient-one rook repair is refuted by the order-twelve X_7.
source: codex-five-copy-switching-continuation-20260826
depends_on:
  - THM-4219-no-sink-endpoint-energy-floor-and-near-ordinal-sharpness
  - THM-4223-cyclic-cut-cover-boolean-mobius-hierarchy-and-two-bad-owner-obstruction
  - THM-4225-bad-owner-upper-zeta-successor-rook-hierarchy
related:
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
script: 04-computation/tournament_plus_min_obstruction_thm4224.cpp
output: 05-knowledge/results/tournament_plus_min_obstruction_thm4224.out
script_sha256: 4c0b1ac8f80cf41d9091420c852e1bd7595fed94105ff577535d4939fc88afb6
output_sha256: fd7392cb630dbcd27bc4c8317fcc0a4c3c998cf321e4cb26e34eb88b43b2ab37
secondary_script: 04-computation/tournament_order11_tail_fiber_thm4224.py
secondary_output: 05-knowledge/results/tournament_order11_tail_fiber_thm4224.out
secondary_script_sha256: 3c9ad5de462dcc71d4312bba12b42086d4e1e43520b5209b4a807f762d059c18
secondary_output_sha256: 6f88b500029fa9b9ddaff00503e5b26676880073a6e34ef121e8cc1e05cef451
hash_basis: raw LF bytes
audit: >
  PASS. A literal scan of all 10! linear orders and all 9! cyclic orders of
  the order-ten witness independently counts H, every End_i and D_i, c,
  N_89, and B_89. A separate subset-DP engine checks every displayed family
  formula for m=1..12. Homebrew nauty 2.9.3 supplies every strong class
  through order nine for the exact ratio census. A standard-library companion
  checks the algebraic consequences of the proved interior-tail identities,
  independently checks X_m through m=7, and exhausts the fixed-gauge
  exactly-two-reversal order-eleven sphere:
  97 strong labelled presentations and 5,335 pair rows. Optimized/unoptimized,
  hash-seed, and sanitizer controls are frozen in the outputs.
---

# THM-4224 -- order-ten minimal plus-min two-bad-owner obstruction

**PROVED / REFUTED / FINITE-EXACT / VERIFIED-EXACT / OPEN FIREWALL.**

THM-4223 proved the owner-refined cyclic cut-cover hierarchy and found that
the pair-product bound

```text
B_ij<=End_i End_j                                            (1)
```

first fails in order nine. Its complete census found no failure through
order nine of the explicitly **OPEN** enlargement

```text
B_ij<=End_i End_j+min(End_i,End_j).                          (2)
```

The equality row behind that extrapolation extends to an infinite family,
but in the wrong direction: order nine is exactly its last nonfailure. This
theorem supplies the correction lineage. It does not record a mistake,
because THM-4223 never promoted `(2)` beyond FINITE-EXACT / OPEN.

## 1. The two-reversal family

For `m>=1`, let `X_m` have ordered vertex set

```text
0<1<2<3<x_1<...<x_m<z.                                   (3)
```

Orient every edge from the lower to the higher vertex except

```text
3->0,                         z->3.                       (4)
```

Put `r=x_m` and `n=m+5`. The cyclic word

```text
0,1,2,x_1,...,x_m,z,3                                    (5)
```

is a directed Hamilton cycle, so every `X_m` is strong.

Retain THM-4223's counts: `c` is the number of directed Hamilton cycles
modulo rotation, `End_i` counts Hamilton paths ending at `i`, `N_ij` counts
spanning two-path covers with terminal set `{i,j}`, and

```text
B_ij=N_ij-End_i-End_j+c                                   (6)
```

counts cyclic orders whose exact bad-owner set is `{i,j}`.

> **Theorem 1 (exact endpoint and two-bad-owner profile).** For every
> `m>=1`,
>
> ```text
> c=1,                         End_z=5,
> End_r=5*2^(m-1),            N_rz=16*2^m-7,
> B_rz=27*2^(m-1)-11.                                    (7)
> ```

### The four-vertex core

Let `L={0,1,2,3}` and `T=X_m[L]`. Its five Hamilton paths are

```text
0123,       1230,       1302,       2301,       3012.     (8)
```

This small list controls all the counts in `(7)`. For later use, direct
enumeration of the sixteen subsets `A subseteq L`, with `B=L\A`, gives

```text
A:          -  0  1  01 2  02 12 012 3  03 13 013 23 023 123 0123
H(T[A])H(T[B]):
            5  1  3   1 3   1  1   1 1   1  1   3  1   3   1    5
H(T[A])End_3(T[B]):
            1  1  1   1 1   1  0   1 0   0  0   0  0   0   0    0.
                                                                    (9)
```

Here `H(empty)=1` and `End_3(empty)=0`. The two row sums are `32` and `7`.

### Proof of Theorem 1

Vertex `z` has the unique outneighbor `3`. In a directed Hamilton cycle it
must therefore be followed by `3`. Once a directed word enters the middle
set `{x_1,...,x_m}`, it can only move upward inside that set or move to `z`.
Rotating a Hamilton cycle after `z` consequently forces

```text
3,0,1,2,x_1,...,x_m,z.                                  (10)
```

The core prefix `3012` is the unique path in `(8)` beginning at `3`, so
`c=1`.

If a Hamilton path ends at `z`, delete its last vertex. Since `m>=1`, the
remaining path consists of one of the five core paths `(8)`, followed by all
middle vertices in increasing order. Conversely each such word appends to
`z`. Hence `End_z=5`.

Now consider a Hamilton path ending at `r`. It contains the forced adjacency
`z,3`. Choose a subset `A` of `{x_1,...,x_(m-1)}`. In one of the five core
paths `(8)`, insert the increasing word on `A`, followed by `z`, immediately
before its occurrence of `3`; then append the complementary middle word and
`r`. Every adjacency is directed. Conversely, the middle vertices before
`z` must form exactly such an increasing block, and after the forced `z,3`
adjacency the remaining middle vertices must form the final increasing
suffix ending at `r`. This is a bijection, proving
`End_r=5*2^(m-1)`.

It remains to count `N_rz`. Label the two components by their distinct
endpoints. Partition `L` into `A` for the path ending at `r` and `B=L\A`
for the path ending at `z`; partition the other `m-1` middle vertices
arbitrarily between the two paths. Each component is a Hamilton path of its
assigned core, followed by its assigned middle vertices in increasing order,
then its terminal `r` or `z`. This gives

```text
2^(m-1) sum_(A subseteq L) H(T[A])H(T[L\A])              (11)
```

objects, except when the `z`-component receives no middle vertex and its
core path ends at `3`, because `3->z` is the reversed edge. The number of
excluded objects is

```text
sum_(A subseteq L) H(T[A])End_3(T[L\A]).                  (12)
```

Using the two sums `32,7` in `(9)` gives

```text
N_rz=32*2^(m-1)-7=16*2^m-7.                              (13)
```

Substitution into `(6)` gives

```text
B_rz=(16*2^m-7)+1-5*2^(m-1)-5
    =27*2^(m-1)-11.                                      (14)
```

This proves `(7)`. QED.

## 2. Minimal refutation and failure anatomy

Since `min(End_r,End_z)=5`, Theorem 1 gives the exact repaired excess

```text
B_rz-End_r End_z-min(End_r,End_z)=2^m-16.                (15)
```

Thus `X_4`, of order nine, is an equality case for `(2)`, while every
`X_m`, `m>=5`, refutes `(2)`. For `m=5`, use upper-triangle row-major
encoding with bit `1` meaning `i->j` for `i<j`. The order-ten word is

```text
110111111111111111111111111110111111111111111            (16)
```

and its exact data are

```text
H=165,                  c=1,
End=(1,1,2,1,5,10,20,40,80,5),
N_89=505,               B_89=421>80*5+5=405.             (17)
```

Swapping `r,z` in `X_4` gives exactly THM-4223's labeled order-nine equality
witness. THM-4223's exhaustive strong-class census proves that `(2)` has no
failure at any order at most nine. Therefore `(16)` is a **minimal-order
counterexample**, with the lower-order part FINITE-EXACT and the order-ten
witness literal.

The `+min` term at order nine has no stable switching-fiber interpretation.
There it equals `End_z=5`, the five core paths `(8)`, merely because

```text
B_rz-End_r End_z=2^m-11=5                 when m=4.        (18)
```

At the next scale the binary middle assignments double while the five core
paths do not. The first failed implication was treating an equality at the
last enumerated scale as evidence for one exceptional object per smaller
endpoint path.

> **Corollary 2 (no fixed min repair).** No universal constant `K` satisfies
>
> ```text
> B_ij<=End_i End_j+K min(End_i,End_j)                    (19)
> ```
>
> for all strong tournaments.

Indeed, on `X_m`,

```text
B_rz-End_r End_z=2^m-11,             min(End_r,End_z)=5.  (20)
```

The obstruction is therefore unbounded relative to the proposed sidecar.

There is one deliberately open multiplicative diagnostic:

```text
25B_ij<=27End_i End_j.                                    (21)
```

On `X_m`,

```text
B_rz/(End_r End_z)=27/25-11/(25*2^(m-1)),                 (22)
```

so any universal multiplicative constant must be at least `27/25`. An exact
strong-class census through order nine finds no failure of `(21)`; its
maximum ratios by order `3,...,9` are

```text
0, 1, 1, 1, 1, 1, 26/25.                                 (23)
```

This is only FINITE-EXACT evidence plus one proved family. Inequality `(21)`
remains **OPEN**, and section 4 shows why even its proof would not establish
HYP-9081.

## 3. Global five-copy compensation on the same family

The failure of every fixed local min repair does not threaten HYP-9081 on
this family. Let `D_i` count one-defect permutations owned by `i`, as in
THM-4219/4223. For `1<=t<m`, a direct cut after the unique bad adjacency,
followed by the same core/middle partition used above, gives

```text
H=5*2^m+5,

D_0=0,
D_1=2^m+1,                 D_2=5*2^m+4,
D_3=3*2^m+2,               D_z=27*2^m-12,

D_(x_t)=17*3^(t-1)2^(m-t+1)-5*2^m+2^(t+4)-5,             (24)

D_r=34*3^(m-1)+11*2^m-12.                                (25)
```

For completeness, the exact enumeration behind `(24)--(25)` is as follows.
Cutting a one-defect word after its bad owner produces two directed paths,
the first ending at that owner and the second starting at an inneighbor of
the owner. Conversely every such ordered cover concatenates to a unique
one-defect word. Thus

```text
D_i=sum_(S proper subset V, i in S) End_i(X_m[S])
       *sum_(a in V\S, a->i) Start_a(X_m[V\S]).            (26)
```

Every directed component word is numerically increasing except that it may
contain either or both of the good descent blocks `z3` and `30`. Assigning
the middle vertices to the increasing runs in `(26)` gives the following
finite-core table:

```text
owner       exact contribution
0           0
1           2^m+1
2           5*2^m+4
3           3*2^m+2
z           27*2^m-12
x_t,t<m     17*3^(t-1)2^(m-t+1)-5*2^m+2^(t+4)-5
r           34*3^(m-1)+11*2^m-12.                         (27)
```

The factors `2` and `3` record assignment to two or three increasing runs;
the constants remove empty runs whose exposed boundary is `03` or `3z`.
This exhausts the possibilities because there are no other directed
descents. Equation `(26)` makes the table a disjoint count rather than a
recurrence fit. The independent audit also reconstructs every row from
subset Hamilton-path DP for `m=1,...,12`.

Summing `(24)--(25)` yields

```text
q=sum_iD_i
 =34*3^m+(34-5m)2^m-5m-44.                              (28)
```

Substitution into THM-4219's exact one-defect form gives the five-copy gap

```text
G_m=q^2+2(n-4)Hq+n(n-3)H^2-5sum_iD_i^2

   =1632*6^m-(8174/3)4^m-(476/3)3^m
      +2168*2^m+1298/3.                                  (29)
```

> **Theorem 3 (all-order global compensation).** Every `X_m`, `m>=1`,
> satisfies HYP-9081 strictly: `G_m>0`.

For `m=1`, `(29)` gives `G_1=3186`. For `m>=2`, multiply `(29)` by three and
divide its possibly negative part by `6^m`. The desired comparison is

```text
4896>8174(2/3)^m+476(1/2)^m.                             (30)
```

It already holds at `m=2`, where the right side is less than `3752`, and
the right side decreases with `m`. The remaining terms
`6504*2^m+1298` are positive. This proves `G_m>0`. QED.

The mechanism is global compensation: the two-bad mass at `{r,z}` outruns
every fixed min correction, while defects at all middle owners disperse the
total collision energy enough to leave the positive gap `(29)`. A proof of
HYP-9081 must exploit such cross-owner structure rather than bound each
`B_ij` independently.

## 4. Moment firewall for local sidecars

Even the open multiplicative diagnostic `(21)` would not imply HYP-9081.
Here is an exact abstract Boolean tensor showing the missing coordinate.
Take `n=6` and set

```text
C_empty=C_V=1,
C_0=C_1=5,                         C_01=25,
C_(V\{0})=C_(V\{1})=5,            C_(V\{0,1})=25,
C_012=12,                          C_345=36,               (31)
```

with every other `C_S=0`. Its layer totals are

```text
(sum_(|S|=k) C_S)_(k=0)^6=(1,10,25,48,25,10,1).          (32)
```

Hence it has the required total `120=5!`, reversal symmetry `C_k=C_(6-k)`,
and first moment `sum_S |S|C_S=360=6!/2`. Its owner marginals are

```text
(48,48,48,72,72,72)=24*(2,2,2,3,3,3),                   (33)
```

matching the universal identity

```text
sum_(S containing i)C_S=(n-2)! d_i^-                     (34)
```

for the feasible strong indegree sequence `(2,2,2,3,3,3)`; its Landau
inequalities are strict before the total sum.

Nevertheless, the low layers give

```text
End_0=End_1=6,                 B_01=25<=36,
H=16,                          D=(50,50,0,0,0,0),          (35)
```

and therefore

```text
q^2+2(n-4)Hq+n(n-3)H^2-5sum_iD_i^2
 =21008-25000=-3992.                                      (36)
```

This is not asserted to be a tournament-realizable tensor. It proves that
Boolean nonnegativity, total cyclic mass, reversal moments, owner marginals,
a feasible strong score sequence, and even `(21)` do not logically force
HYP-9081. Actual successor-compatible switching data remain indispensable.

## 5. Interior tail fibers and the order-eleven two-reversal sphere

### 5.1 Every interior tail fiber

Fix `1<=t<m` and put `u=2^(t-1)`. The endpoint `x_t`, rather than only the
terminal `r=x_m`, has a closed profile:

```text
End_(x_t)=5u,                 N_(x_t,z)=32u,
B_(x_t,z)=27u-4.                                      (37)
```

To end a Hamilton path at `x_t`, every later tail vertex
`x_(t+1),...,x_m` is forced into the increasing block before `z`; after the
unique descent `z->3`, an arbitrary subset of the earlier tail lies in that
block and its complement forms the final increasing suffix to `x_t`. The
five core paths `(8)` and the `2^(t-1)` earlier-tail choices give the first
formula in `(37)`.

For a two-path cover ending at `{x_t,z}`, no later tail vertex can lie in the
`x_t`-component: without `z` it has no descent to a smaller tail vertex.
Consequently every later tail vertex lies in the `z`-component. This component
is nonempty beyond its core, so the seven terminal exclusions `(12)` disappear.
The full core partition sum `32` and the same earlier-tail choices give
`N_(x_t,z)=32u`. Substitution into `(6)`, using `c=1` and `End_z=5`, gives the
last formula in `(37)`.

With `P_(t,z)=End_(x_t)End_z=25u`, the three exact diagnostics are

```text
B_(x_t,z)-P_(t,z)=2^t-4,
B_(x_t,z)-P_(t,z)-min(End_(x_t),End_z)=2^t-9,
27P_(t,z)-25B_(x_t,z)=100.                              (38)
```

Thus the product bound fails on the interior fiber exactly for `t>=3`, the
`+min` repair fails exactly for `t>=4`, and the open `27/25` candidate has
constant positive cross-multiplied slack `100`. At `t=3`, the order-nine row
has `B/P=104/100`; at `t=4`, the order-ten row has
`B=212>200+5`. The mechanism is the forced nonempty later tail: it removes
the seven exclusions present at the terminal fiber, rather than changing the
core partition count.

### 5.2 Successor rooks do not give an all-order local repair

Use THM-4225's two-owner successor-rook count

```text
rho_ij=d_i^- d_j^- - kappa_ij.                           (39)
```

On an interior tail pair and on the terminal pair respectively,

```text
rho_(x_t,z)=(t+3)(m+3)-(t+2),
rho_(r,z)=m^2+5m+7.                                     (40)
```

Indeed `d_(x_t)^-=t+3`, `d_z^-=m+3`, and their common-inneighbour count is
`t+2`; put `t=m` for the terminal formula. Since

```text
B_(r,z)-End_r End_z=2^m-11,                             (41)
```

while both `rho_(r,z)` and `rho_(r,z)+min(End_r,End_z)` grow only
quadratically, no fixed `K>=0` can make either

```text
B_ij<=End_i End_j+K rho_ij,
B_ij<=End_i End_j+K(rho_ij+min(End_i,End_j))             (42)
```

valid for all strong tournaments.

The coefficient-one boundary is exact:

```text
End_r End_z+rho_(r,z)-B_(r,z)=m^2+5m+18-2^m,
End_r End_z+min(End_r,End_z)+rho_(r,z)-B_(r,z)
                                      =m^2+5m+23-2^m.   (43)
```

Both margins are nonnegative exactly for `m<=6`. Direct evaluation handles
that finite range; at `m=7` both are negative, and the forward difference
`2m+6-2^m` is then strictly decreasing. Thus the order-eleven `X_6` is the
last terminal survivor,

```text
(B,P,min,rho)=(853,800,5,73),       margins=(20,25),     (44)
```

whereas order twelve `X_7` is the first terminal failure,

```text
(B,P,min,rho)=(1717,1600,5,91),     margins=(-26,-21).   (45)
```

Thus no fixed additive multiple of the normalized two-rook count `rho`, with
or without `min`, captures the weighted successor completions lost by that
local pair statistic. This does not refute the raw upper-zeta bound
`B_ij<=Z_ij=(n-3)!rho_ij`, which is tautological.

### 5.3 Exact fixed-gauge order-eleven sphere

Fix the labelled transitive order `0<...<10` and reverse exactly two of its
`55` arcs. Among the `binom(55,2)=1485` resulting labelled presentations,
exactly `97` are strong. The exact reachability classification is the disjoint
union

```text
{ {(0,10),e}: e notin {(0,10),(0,1),(9,10)} }
union
{ {(0,b),(c,10)}: 1<=c<=b<=9 },                         (46)
```

of sizes `52` and `45`. Their total `97` equals `n^2-2n-2` at `n=11`; the
companion compares `(46)` with direct strong-connectivity tests on all `1485`
candidates.

> **Finite-exact order-eleven sphere.** On all `97` strong presentations and
> all `97*55=5,335` unordered owner pairs, the five-copy gap is strictly
> positive and
>
> ```text
> 25B_ij<=27End_iEnd_j,              B_ij<=End_iEnd_j+rho_ij.  (47)
> ```

There are four product failures and three `+min` failures inside the sphere.
The largest local ratio is `107/100`, uniquely at the interior pair
`(x_5,z)` of `X_6`, where `(B,P)=(428,400)`. The largest positive normalized
rook excess is `53/73`, at the terminal pair of `X_6`.

The smallest five-copy gap is `27,894,642`, attained by three labelled
presentations of one tournament: a transitive ten-vertex prefix plus a final
vertex beating exactly prefix vertices `0` and `8`. There the exact target
and energy are `1,134,462,087` and `221,313,489`, so their reduced ratio is
`378154029/73771163>5`, while the largest local `B/P` is only `1`. Local
pair severity and global five-copy closeness are therefore genuinely
different coordinates even in this small shell.

This sphere is presentation-dependent. In the fixed gauge `X_6` has reversal
set `{(0,3),(3,10)}`, but swapping labels `0,1` gives the isomorphic
three-reversal presentation `{(0,1),(1,3),(3,10)}`. Hence `(47)` is neither an
isomorphism-class census nor a radius-two ball. For arbitrary strong order
eleven, both inequalities in `(47)` and HYP-9081 remain **OPEN** outside this
sphere. All-order `27/25`, HYP-9081, and `(OS+)` also remain **OPEN**; the
coefficient-one rook repair is globally **REFUTED** by `(45)`.

## 6. Audit and exact scope

The audit has four independent lanes:

1. a literal order-ten engine scans all `10!` linear permutations and the
   `9!` cyclic representatives beginning at vertex zero, with no subset DP;
2. a subset Hamilton-path DP checks `(7)`, `(24)--(29)` for `m=1,...,12`;
3. `gentourng -c` supplies every strong isomorphism class through order nine
   for the exact ratio census `(23)`;
4. the standard-library companion checks the algebraic consequences of
   `(37)--(45)`, independently verifies every interior tail pair for
   `m=2,...,7`, and exhausts the exact sphere `(46)`, including its 5,335 pair
   rows and gauge controls. The all-`m` quantifier comes from the prose
   forced-tail bijection, not this finite replay.

Replay the new lane with

```bash
python3 -B 04-computation/tournament_order11_tail_fiber_thm4224.py
python3 -B -O 04-computation/tournament_order11_tail_fiber_thm4224.py
PYTHONHASHSEED=4224 python3 -B \
  04-computation/tournament_order11_tail_fiber_thm4224.py
```

All streams byte-match the frozen output. The stored outputs freeze commands,
values, compiler controls, and hashes.
The all-order family formulas and positivity are proved above; the literal
and finite DP ranges are audits, not the source of their quantifiers.

HYP-9081 remains **OPEN** for arbitrary strong tournaments. Inequality `(21)`
also remains **OPEN**. No local product, fixed-min repair, scalar reversal
moment, or score-compatible Boolean tensor is promoted into a sufficient
criterion.
