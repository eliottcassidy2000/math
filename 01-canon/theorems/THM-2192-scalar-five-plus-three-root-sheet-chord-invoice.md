---
id: THM-2192
title: "Scalar five-plus-three root-sheet chord invoice"
status: >
  PROVED + VERIFIED-EXACT. In any scalar five-unit/three-deep candidate,
  every generic phase safe for the three divided deep teeth forces an exact
  thirteen-root ownership law: over the guard-danger base arc the five unit
  doubletons partition the ten guard-safe sheets, while over the guard-safe
  base arc there is either a singleton-plus-four-chord partition or a
  five-chord cover with exactly one incidence defect. This gives a signed
  three-deep-union overlap invoice. On the guard-safe base arc, the induced
  anchored chord diagram on F_13 excludes eight of the 252 possible
  multisets of five unit residue lengths, uniformly at every 13-adic depth.
  On the guard-danger base arc, a ten-vertex perfect-matching carrier realizes
  only 216 of the 252 profiles. Each of the other 36 either belongs to the
  already-empty eight or forces the three divided blockers to cover the whole
  fat guard comb almost everywhere; the odd guard and distinct positive
  blocker coefficients make that fat-guard cover impossible in the actual
  LRC scalar lane. Thus all 36 missing guard-danger profiles are uniformly
  empty there. Three additive endpoint-potential inequalities detect 17 of
  those missing profiles, but the other 19 lie inside the convex hull of the
  feasible length-count vectors. They are genuine holes in the coloured
  perfect-matching/Hafnian support, invisible to every affine-linear
  length-count separator. The unit-annulus law also eliminates
  every branch in which the maximum valuation among the three actual blockers
  is repeated; a survivor has a unique deepest blocker. The scalar
  five-plus-three tail is narrowed but not emptied, so this is not a proof of
  LRC(14).
source: codex-2026-07-24-scalar-five-plus-three-root-sheets
depends_on:
  - THM-1155
  - THM-2138
  - THM-2148
  - THM-2168
related:
  - THM-2186
  - THM-2190
  - THM-2193
script: 04-computation/lrc14_scalar_five_plus_three_root_sheet_chords_thm2192.py
output: 05-knowledge/results/lrc14_scalar_five_plus_three_root_sheet_chords_thm2192.out
script_sha256: 6c8796b6aeeb7539157974c9ec2b6eb91d32263e1d72cb569e7dbd970657b50e
output_sha256: 7a2416e8fedca07a04c0c00cc82b39e73a4eb7e997cf017c5189fc8cead83104
hash_basis: working-tree bytes (LF)
---

# THM-2192 -- the scalar `5+3` root-sheet chord invoice

For a nonzero integer `a`, put

```text
D_a={t in R/Z:||at||<1/14},
E_a={t in R/Z:||at||<1/7},
C_a={t in R/Z:||at||>1/7}.                           (1)
```

Let

```text
H,q_1,...,q_5,v_1,v_2,v_3
```

be nonzero integers such that `H,q_1,...,q_5` are units modulo thirteen.
Suppose that, outside a null set,

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(13v_j).          (2)
```

This is the exact one-dimensional cover in the fully scalar `(5,3)` branch
of THM-2168.  No distinctness assumption is needed for the local theorem,
although the LRC application has positive distinct terminal values.
For the fat-guard exclusion in Section 4.2, we use the additional facts from
THM-2168's actual scalar lane that `H` is positive and odd and the three
`v_j` are distinct and positive.

Define

```text
G={y:||v_j y||>1/14 for j=1,2,3},
S(y)=#{i:y in D_(q_i)}.                              (3)
```

Away from the finitely many band boundaries and the image under
multiplication by thirteen of the null exceptional set in (2), one has

```text
y in G intersection E_H  implies S(y)=0,
y in G intersection C_H  implies S(y)<=1.            (4)
```

There is a sharper sheetwise classification behind (4).  For

```text
R_y={x in R/Z:13x=y},
B(y)={x in R_y:||Hx||>1/7},
A_i(y)={x in R_y:||q_i x||<1/14},                    (5)
```

the following statements hold.

1. If `y in G intersection E_H`, then the five sets `A_i(y)` are disjoint
   doubletons and

   ```text
   B(y)=disjoint_union_(i=1)^5 A_i(y),       |B(y)|=10. (6)
   ```

2. If `y in G intersection C_H` and `S(y)=1`, the unique active lower
   label gives a singleton, the other four sets are doubletons, and all five
   sets partition the nine guard-safe sheets:

   ```text
   B(y)=disjoint_union_(i=1)^5 A_i(y),       |B(y)|=9.  (7)
   ```

3. If `y in G intersection C_H` and `S(y)=0`, all five sets are
   doubletons.  Put

   ```text
   w(y)=sum_i |A_i(y)\B(y)|,
   e(y)=sum_(x in B(y)) (#{i:x in A_i(y)}-1).         (8)
   ```

   Then every term in the second sum is nonnegative and

   ```text
   w(y)+e(y)=1.                                       (9)
   ```

   Thus exactly one of two things happens: one chord has one endpoint on a
   guard-unsafe sheet and the other four chords match the remaining eight
   safe sheets, or all five chords lie on the safe sheets and exactly one
   safe sheet is double-owned.

All assertions are insensitive to changing open sets to closed sets at the
finite boundary set.

## 1. The thirteen-root count

Fix a generic `y` and enumerate its thirteen roots.  Because `H` is a
thirteen-unit, their `H`-values form the complete translated
thirteen-grid.  Direct interval counting gives

```text
#{x in R_y:||Hx||<1/7}
 =3 if y in E_H,
 =4 if y in C_H.                                     (10)
```

Hence `|B(y)|` is respectively ten or nine.  Similarly, for every unit
`q_i`,

```text
|A_i(y)|
 =1 if y in D_(q_i),
 =2 if y notin closure(D_(q_i)).                     (11)
```

Every deep tooth is constant on the root fibre:

```text
13v_j x=v_j y                 for every x in R_y.    (12)
```

When `y in G`, all three deep teeth are therefore safe on all thirteen
sheets.  The five unit masks must cover `B(y)`.  Their total root capacity
is

```text
sum_i |A_i(y)|=10-S(y).                              (13)
```

For `y in E_H`, equations (10) and (13) give

```text
10<=10-S(y),
```

so `S(y)=0`; equality of total capacity with the size of `B(y)` proves the
disjoint partition (6).  For `y in C_H`, they give

```text
9<=10-S(y),
```

so `S(y)<=1`.  If `S(y)=1`, equality gives (7).  If `S(y)=0`, coverage
gives nonnegative excess multiplicities on `B(y)`, and counting all ten
incidences gives

```text
10=w(y)+9+e(y),
```

which is (9).  This proves (4)--(9).

The almost-everywhere hypothesis causes no loss.  Multiplication by
thirteen maps null sets to null sets.  A violation of one of the strict
conclusions persists on an open phase interval, in which one can choose
`y` off that image and off every band boundary.

## 2. The signed three-deep-union invoice

Let

```text
V=D_(v_1) union D_(v_2) union D_(v_3).               (14)
```

Equation (4) is the pointwise almost-everywhere inequality

```text
sum_(i=1)^5 1_(D_(q_i))(y)
       <=1_(C_H)(y)                    on (R/Z)\V.    (15)
```

Every nonzero danger comb has measure `1/7`, while `C_H` has measure
`5/7`.  Integrating (15) and subtracting from these equal total masses
gives the necessary signed overlap invoice

```text
sum_(i=1)^5 measure(D_(q_i) intersection V)
              >=measure(C_H intersection V).         (16)
```

Thus a prospective closure can work entirely on the divided-deep union:
it is enough to prove the strict reverse of (16).  Unlike an unsigned
capacity estimate, (16) retains whether the deep masks absorb unit danger
or guard safety.

The guard-safe case is never vacuous.  The union bound gives

```text
measure(E_H union V)<=2/7+3/7=5/7,                   (17)
```

so `G intersection C_H` has positive measure.  A generic phase in this set
always supplies one of the two chord carriers in (7)--(9).

## 3. A repeated deepest valuation is impossible

Let the valuations of the three **actual** blocker coefficients be

```text
lambda_j=nu_13(13v_j)>=1,
beta=max_j lambda_j.                                 (17a)
```

Suppose that `beta` is attained at least twice.  Put

```text
m=beta+1,                    N=13^m.                 (17b)
```

Normalize the unit guard coefficient to one modulo `N` and restrict (2) to
THM-2138's unit annulus `U_N`.  If an actual blocker has depth `beta`,
write it as `13^beta a` with `a` a thirteen-unit.  At a primitive annulus
point `t=z/N`,

```text
13^beta a t=az/13.                                  (17c)
```

This is a nonzero thirteenth root, so its circle norm is at least
`1/13>1/14`.  Every depth-`beta` blocker is therefore safe throughout
`U_N`.

If exactly two blockers have depth `beta`, the proposed cover of `U_N`
uses at most five unit masks and one positive-valuation mask.  For `m>=3`,
let `M_m` be THM-2138's maximum unit-annulus mask size.  At odd depth,

```text
5M_m+M_m=|U_N|-10.                                  (17d)
```

At even depth, THM-2138's positive-valuation deficit gives

```text
5M_m+(M_m-20)=|U_N|-10.                             (17e)
```

Both bounds are strictly below the universe size.  If all three blockers
have depth `beta`, there are only five unit masks.  At odd `m` this is
smaller than the left side of (17d).  At even `m>=4`, THM-2138 gives
`|U_N|=6M_m-10` and `M_m>10`, so again `5M_m<|U_N|`.  At the sole omitted
base `m=2`, necessarily `beta=1`, so all three positive-depth blockers are
deepest and safe; THM-2138's exact mod-169 base has `|U_169|=110` and
unit-mask maximum `20`, whence `5*20<110`.

No torsion endpoint occurs because a power of thirteen is coprime to seven
and fourteen.  Any uncovered torsion point has all inequalities strict and
therefore thickens to an uncovered open interval, contradicting the
almost-everywhere cover (2).  After relabelling, every surviving profile
must therefore satisfy

```text
lambda_1<=lambda_2<lambda_3.                         (17f)
```

On the depth-`lambda_3+1` annulus, the unique deepest blocker is safe and
the two shallower positive-valuation masks remain active.  This isolates
the genuinely new two-active-mask obstruction beyond THM-2138's scalar
`5+2` theorem.

## 4. The anchored `F_13` chord carrier

Fix such a generic phase in `G intersection C_H`.  Label its roots by
`F_13` so that increasing the label by one increases the guard value by
`1/13`.  The four guard-unsafe sheets are then four consecutive vertices.

For each unit label define

```text
d_i=H q_i^(-1) mod 13,
ell_i=min(d_i mod 13,-d_i mod 13) in {1,2,3,4,5,6}.  (18)
```

When `A_i(y)` is a doubleton, its two endpoints differ by `+d_i` or `-d_i`.
Indeed, in the terminal coordinate the two dangerous roots are adjacent;
changing to guard-root coordinates multiplies their difference by
`Hq_i^(-1)`.  Thus `ell_i` is the unoriented chord length.  A singleton
still retains its coefficient label `ell_i`, although no chord is visible
on that fibre.

Equations (7)--(9) leave exactly three finite carrier types on the nine
safe vertices:

```text
M: one singleton and a perfect matching on the other eight vertices;
X: one safe-to-unsafe cross chord and a perfect matching on the
   other eight safe vertices;
D: five internal chords with degree sequence (2,1,1,1,1,1,1,1,1). (19)
```

The carrier is anchored: the omitted four-consecutive block is part of the
data.  Enumerating (19) gives `244` possible length multisets among the
`252` multisets of size five from `{1,...,6}`.  The eight impossible
profiles are

```text
(3,3,3,3,3),
(3,6,6,6,6),
(4,5,6,6,6),
(4,6,6,6,6),
(5,5,5,6,6),
(5,5,6,6,6),
(5,6,6,6,6),
(6,6,6,6,6).                                        (20)
```

Consequently every scalar `5+3` branch whose normalized unit residues have
one of the profiles (20) is empty, uniformly in the three deep
coefficients and at every 13-adic depth.

For transparency, the finite enumeration has only

```text
M: 9*7!!*6=5670 anchored vertex-labelled states,
X: 9*4*7!!=3780 anchored vertex-labelled states,
D: 9*C(8,2)*5!!=3780 anchored vertex-labelled states. (21)
```

The terminal labels are not assigned to the matching edges in this census;
only their length multiset is retained.  The companion checks these counts,
the `244/252` profile census, the exact list (20), and the root-step
arithmetic for all nonzero residue classes.  It uses integer and rational
arithmetic only and gives byte-identical normal and optimized-Python output.

### 4.1. The guard-danger matching fork

There is a second anchored carrier whenever `G intersection E_H` has
positive measure.  Choose a generic phase in that set.  Equation (6) says
that the three consecutive guard-unsafe vertices are omitted and the five
unit chords form a perfect matching of the remaining ten vertices.  There
are exactly

```text
9!!=945                                                (21a)
```

such anchored vertex-labelled matchings.  Their length multisets realize
only `216` of the `252` ambient profiles.  The `36` missing profiles are

```text
(1,1,1,1,2), (1,1,1,2,3), (1,1,2,2,2), (1,1,2,3,3),
(1,2,2,2,3), (1,2,3,3,3), (1,3,4,4,4), (1,4,4,4,5),
(1,5,5,5,5), (2,2,2,2,2), (2,2,2,3,3), (2,2,2,4,4),
(2,3,3,3,3), (2,3,5,5,5), (2,4,4,4,4), (2,4,5,5,5),
(2,5,5,5,6), (2,6,6,6,6), (3,3,3,3,3), (3,3,3,3,4),
(3,3,3,4,4), (3,4,4,4,4), (3,5,5,5,5), (3,5,6,6,6),
(3,6,6,6,6), (4,4,4,4,5), (4,4,6,6,6), (4,5,5,5,5),
(4,5,5,6,6), (4,5,6,6,6), (4,6,6,6,6), (5,5,5,5,6),
(5,5,5,6,6), (5,5,6,6,6), (5,6,6,6,6), (6,6,6,6,6).
                                                               (21b)
```

All eight profiles in (20) occur in (21b), as they must: the guard-safe
carrier already excludes them without a case split.  For any of the other
`28` profiles in (21b), the exact conclusion is the following fat-guard
alternative:

```text
E_H subset D_(v_1) union D_(v_2) union D_(v_3)
                                              almost everywhere. (21c)
```

Indeed, if (21c) failed on positive measure, then `G intersection E_H`
would have positive measure after discarding the finite deep-band
boundaries.  A generic phase there would supply the forbidden ten-vertex
matching.  Thus (21c) follows.  In the abstract local theorem, where `H`
need not be odd and the `v_j` need not be distinct, this remains a genuine
structural alternative.  The actual LRC scalar lane excludes it in the next
subsection.  The companion separately verifies the `945`, `216/252`, exact
`36`, and eight-profile subset assertions.

### 4.2. The odd-guard obstruction empties the fat-guard fork

Assume now the actual scalar-lane hypotheses: `H` is positive and odd, and
the `v_j` are distinct positive integers.  Then

```text
E_H is not contained almost everywhere in
             D_(v_1) union D_(v_2) union D_(v_3).     (21d)
```

The proof has a finite-kernel part and a one-interval part.  First replace
each open `D_a` by

```text
bar(D_a)={t:||at||<=1/14}.
```

An almost-everywhere open cover of the open set `E_H` implies the pointwise
closed cover

```text
E_H subset union_(j=1)^3 bar(D_(v_j)).                (21e)
```

Indeed, a point at which all three inequalities were strictly reversed
would thicken to an uncovered open interval.

Restrict (21e) to

```text
F=ker(H) isomorphic to Z/HZ.
```

THM-2148's finite three-colour lemma says that either some restricted
character is trivial or the three restrictions are the nonzero characters
of a `C_2 x C_2` quotient.  A cyclic group has no such quotient.  After
relabelling,

```text
H divides v_1.                                       (21f)
```

Write `v_1=Hk_1`.  Since `bar(D_(k_1))` has measure `1/7`, there is an open
phase interval in `E_1` on which it is inactive.  Over any phase `y` in
that interval, (21e) says that the two remaining translated bands cover
the `H`-root torsor.  A nontrivial restricted character has odd image order
`m>=3`.  Any translate of a closed radius-`1/14` band contains at most

```text
floor(m/7)+1
```

values of the `m`-grid, and

```text
2(floor(m/7)+1)<m                 for odd m>=3.       (21g)
```

This is THM-2148's translated-grid bound with its closed-endpoint
convention; no generic-position assumption is hidden here.
Thus two nontrivial restrictions cannot cover the torsor, so a second
restriction is trivial and `H|v_2`.  The two closed quotient combs
`bar(D_(k_1)),bar(D_(k_2))` overlap on an open neighbourhood of zero.
Their union therefore has measure strictly below `2/7=measure(E_1)`.
Choose an open phase interval in `E_1` where both are inactive.  The third
band must cover the whole root torsor there, which is impossible for a
nontrivial restriction.  Hence

```text
H divides v_1,v_2,v_3.                               (21h)
```

Put `v_j=Hk_j` and order the distinct positive quotients
`k_1<k_2<k_3`.  The cover would descend to

```text
E_1 subset D_(k_1) union D_(k_2) union D_(k_3)
                                               almost everywhere. (21i)
```

THM-1155's exact interval-fragmentation bound says that a comb `D_k`
meets an interval of length `L` in measure at most

```text
L/7+1/(7k).
```

Applying it to `E_1`, of length `2/7`, gives

```text
1/k_1+1/k_2+1/k_3>=8/7.                              (21j)
```

If `k_1>=2`, distinctness bounds the left side by
`1/2+1/3+1/4=13/12<8/7`.  Thus `k_1=1`.  On

```text
J=(1/14,1/7)
```

the first comb is absent, so the same bound for the other two gives

```text
1/k_2+1/k_3>=5/14.                                   (21k)
```

The pointwise closed consequence (21e), now in the quotient, can be tested
at `1/11` and `1/13`.  A nonzero residue on either grid has norm at least
`1/11` or `1/13`, both strictly greater than `1/14`.  Therefore one of
`k_2,k_3` is divisible by eleven and one is divisible by thirteen.

If different quotients carry the two divisibilities, their reciprocal sum
is at most

```text
1/11+1/13=24/143<5/14,
```

contrary to (21k).  If the same quotient `r` carries both, then `r>=143`;
(21k) forces the other quotient to be `2`.  But `D_2` is empty on `J`,
while the fragmentation bound gives

```text
measure(J intersection D_r)
 <=1/98+1/(7r)
 <=1/98+1/1001
 =157/14014<1/14=measure(J).                          (21l)
```

This is the final contradiction and proves (21d).  Consequently every one
of the `36` profiles missing from the guard-danger matching census is
uniformly empty in the actual scalar `5+3` lane; exactly `216/252` profiles
are not excluded by this carrier test.

Oddness is load-bearing in the translated finite-kernel step.  For even `H`, two
nontrivial restrictions can have the same order-two direction and occupy
complementary translated cosets.  The theorem neither asserts that this
abstract exception occurs globally nor classifies it.  Distinctness is
load-bearing in (21j), and sign is forgotten by `D_a=D_(-a)`; the actual
positive distinct coefficients supply the needed absolute-value
distinctness.  The first structural conclusion `H|v_j` for at least one
`j` uses neither oddness nor distinctness; oddness enters only when the
translated two-band cover is pruned.

### 4.3. Pairwise linear shadows and nineteen matching-support holes

The `36` exclusions have two mathematically different mechanisms.  Let

```text
c_r=#{i:ell_i=r},                       1<=r<=6.
```

Every perfect matching of

```text
S={3,4,...,12} subset F_13
```

obeys the three additive edge inequalities

```text
sum_(r=1)^6 r c_r<=25,
c_2<=4,
c_2+3c_3+2c_4<=12.                                 (21m)
```

All three have short endpoint-potential certificates.  For the first, put
`p_v=dist_(F_13)(v,8)`; the triangle inequality gives
`ell(u,v)<=p_u+p_v`, and `sum_(v in S)p_v=25`.  For the second, the four
vertices `{5,6,9,10}` meet every length-two edge of `S`.  For the third,
use edge weights

```text
w_2=1,       w_3=3,       w_4=2,       w_1=w_5=w_6=0
```

and the endpoint potentials, in the order `3,4,...,12`,

```text
(1,0,0,2,3,3,2,0,0,1).
```

Every edge satisfies `w_(ell(u,v))<=p_u+p_v`, and the potentials sum to
`12`.  Summing any of these pointwise inequalities over the five disjoint
matching edges proves (21m).

Exactly `17` of the `36` missing profiles violate at least one inequality
in (21m): the first detects the twelve profiles of total length above
`25`, the second detects `(2,2,2,2,2)`, and the third detects

```text
(2,3,3,3,3), (3,3,3,3,3), (3,3,3,3,4), (3,3,3,4,4).
```

The other `19` are

```text
11112, 11123, 11222, 11233, 12223, 12333, 13444,
14445, 15555, 22233, 22244, 23555, 24444, 24555,
25556, 34444, 35555, 44445, 45555,                    (21n)
```

where a word denotes the sorted length multiset.  These are not merely
missed by the three chosen potentials.  Every profile in (21n) lies in the
convex hull of feasible matching profiles.

Here is a compact exact certificate.  Name the feasible anchors

```text
A=11111, B=12222, C=13333, D=11333, E=22333,
F=44444, G=55555, H=22223, I=22224, J=22225,
K=33335, L=22666, M=33444.
```

One matching witness for each anchor, with hyphens separating multi-digit
vertices, is

```text
A:34,56,78,9-10,11-12;   B:34,57,68,9-11,10-12;
C:36,45,7-10,8-11,9-12; D:34,56,7-10,8-11,9-12;
E:35,46,7-10,8-11,9-12; F:3-12,48,59,6-10,7-11;
G:38,49,5-10,6-11,7-12; H:35,46,79,8-11,10-12;
I:3-12,46,57,8-10,9-11; J:35,46,7-12,8-10,9-11;
K:36,47,5-10,8-11,9-12; L:39,4-11,57,6-12,8-10;
M:3-12,47,59,6-10,8-11.
```

Interpreting each word as its vector `(c_1,...,c_6)`, the nineteen convex
identities are

```text
11112=3A/4+B/4;             11123=A/2+B/4+C/4;
11222=A/4+3B/4;             11233=A/4+B/4+C/2;
12223=3B/4+C/4;             12333=(D+E)/2;
13444=3A/20+C/4+3F/5;       14445=A/5+3F/5+G/5;
15555=A/5+4G/5;             22233=(H+E)/2;
22244=3I/4+F/4;             23555=J/4+K/4+G/2;
24444=I/4+3F/4;             24555=J/4+F/5+11G/20;
25556=J/12+L/3+7G/12;       34444=(M+F)/2;
35555=K/4+3G/4;             44445=4F/5+G/5;
45555=F/5+4G/5.                                      (21o)
```

It follows immediately that no affine-linear one-sided inequality in the
six length counts which holds on every feasible matching profile can exclude
any profile in (21n).

Equivalently, form the coloured matching polynomial

```text
P_S(x_1,...,x_6)
 =sum_(perfect matchings M of S) product_({u,v} in M) x_(ell(u,v)).
                                                               (21p)
```

This is the Hafnian of the symmetric adjacency matrix whose `uv` entry is
`x_(ell(u,v))`.  The profiles in (21n) are zero coefficients lying inside
the Newton polytope of `P_S`: genuine matching-support holes.  Thus a
length-additive binary relation sees only `17/36` exclusions.  The other
`19` require the global vertex-disjointness sidecar, and any recursive lift
must retain either the matching itself, the Hafnian support, or the
phase-by-phase matching movie.  The companion verifies all `135`
edge-potential inequalities, the exact `17+19` split, feasibility of every
anchor, and all nineteen rational identities in (21o).

## 5. Why the octagon `5+3` word is not the map

THM-2186's octagon target is the minimum of all pairwise distances among
eight phase points.  The present scalar target observes eight distances
from an anchored origin, together with a separate guard threshold and the
unit/deep partition.  Forgetting the anchor does not preserve danger.

An exact hostile control is the cyclic gap vector

```text
(1,1,1,1,1,7,8,8)/28.                               (22)
```

Place a vertex at the origin: a scalar terminal is dangerous.  Rotate the
same labelled cyclic configuration until the origin is the midpoint of an
`8/28` gap: every vertex is then at distance at least `1/7>1/14` from the
origin.  The cyclic order and the complete gap vector up to rotation are
unchanged, while the scalar predicate reverses.

The exact connection ledger is therefore

```text
source:       anchored scalar terminal phases;
map:          forget the origin and retain cyclic gaps;
preserved:    cyclic order and all adjacent gap magnitudes;
destroyed:    origin/guard phase, unit/deep type, coefficient winding;
needed sidecar:
              the anchored four-sheet guard block, residue chord labels,
              and the integer winding coefficients;
decisive test:
              the two rotations of (22).             (23)
```

The shared numerical phrase `5+3` is not an incidence-preserving map from
the octagon star order to the scalar Fano residual.

## 6. Exact scalarization kernel and the rank-seven no-go

The relation-rank interface can be stated without ambiguity.  Let

```text
A:Q^9 -> Gamma_Q
```

send the nine coordinate vectors to

```text
(g,c_1,...,c_5,u_1,u_2,u_3),
```

and let `h:Gamma_Q->Q` be the LRC evaluation.  Put

```text
L_char=ker(A),                 L_eval=ker(h A).       (24)
```

There is an exact sequence

```text
0 -> L_char -> L_eval -> im(A) intersection ker(h) -> 0. (25)
```

In a rational-rank-two lane, the dimensions are `7,8,1`.  The final
one-dimensional quotient is exactly the transverse evaluation-kernel
sidecar.  In the fully scalar lane, `rank(A)=1` and `h` is nonzero on its
image.  Therefore

```text
L_char=L_eval,                  dim L_char=8.         (26)
```

There is no transverse descended relation to force in precisely the branch
under study.

For comparison with THM-2190 and THM-2193, let the original
thirteen-speed relation space `Lambda(V)_Q` have dimension twelve.  Even
under the additional hypothesis of an injective lift

```text
T:Q^9 -> Q^13,              V T=h A,                 (27)
```

the scalar relation image `T(L_eval)` has dimension eight.  If a bounded
relation space `S subset Lambda(V)_Q` has dimension `r`, Grassmann's
inequality gives only

```text
dim(S intersection T(L_eval))>=r+8-12=r-4.           (28)
```

Thus rank seven guarantees only three bounded scalar relations.  Eight
independent bounded scalar relations are sufficient, by maximal minors, to
bound the primitive nine-coefficient scalar row; the guaranteed
intersection is short of that by five dimensions.  The precise missing
sidecar is either:

```text
dim(S intersection T(L_eval))=8
```

with a height-controlled inverse for `T`, or five additional independent
bounded relations inside `T(L_eval)` beyond the worst-case rank-seven
guarantee.  Neither (27) nor this internal-rank condition is supplied by
THM-2190 or THM-2193.  Thus even THM-2193's explicit height
`78*7^21` rank-seven harvest does not by itself close scalar `5+3`.

## 7. Scope

The theorem supplies:

1. an all-depth root-sheet equality law;
2. the signed overlap target (16);
3. the unique-deepest valuation reduction (17f);
4. the faithful anchored monomer-dimer carrier;
5. thirty-six uniformly empty residue profiles in the actual scalar lane;
6. the odd-guard theorem excluding the only fat-guard escape; and
7. the exact split of those exclusions into `17` endpoint-potential
   shadows and `19` nonlinear holes in the matching/Hafnian support.

It does not prove the strict reverse of (16) for every coefficient row,
classify the remaining `216` residue profiles, bound the three deep
valuations in the unique-deepest lane, control its two active
positive-valuation masks on the surviving `216` guard-danger profiles, or
prove LRC(14).  The next exact target is to propagate the anchored chord
ownership through multiplication by thirteen and show that the two shallower
divided-deep combs cannot shelter every ownership transition. QED.
