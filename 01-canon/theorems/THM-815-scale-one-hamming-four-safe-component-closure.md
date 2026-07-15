---
id: THM-815
title: Scale-one Hamming-four safe-component closure
status: PROVED by two independent reductions (sharp interval-comb component ladder; collar-cycle/doubling box), with a finite-recursion theorem through radius six, a nonprimitive H6 contraction addendum, complete closure of the primitive-core three-antipodal-pair H6 stratum, and closure of the first eleven two-antipodal-pair rows + FINITE-EXACT (666,705 nested Hamming-four rows, two higher-radius initial censuses, two 35,640-row component/endpoint replays, an independent 768,735-row C++ collar certificate, a 136,288-prefix nonprimitive H6 chamber certificate, ten primitive-core H6 slice trees with 3,699 total states, six residual H6 trees with 2,653,600 total states, and eleven `f=2` trees with 1,517,085 total states)
source: codex-2026-07-15-S10 Hamming-four continuation; codex-2026-07-15-S11 nonprimitive H6 addendum; codex-2026-07-15-S14 primitive AP-pin reduction, residual six-row exact closure, and first `f=2` cap-stratum closure
depends_on: [LRC(<=13), THM-795, THM-800, THM-804, THM-806, THM-810, THM-816]
related: [THM-770, THM-800, THM-804, THM-810, THM-816, THM-820, THM-845, HYP-6820]
verification:
  - 04-computation/lrc13_h4_scale_one_component_ladder_codex_S10.py
  - 05-knowledge/results/lrc13_h4_scale_one_component_ladder_codex_S10.out
  - 04-computation/lrc13_scale_one_hamming_four_collar_closure_codex_S10.cpp
  - 05-knowledge/results/lrc13_scale_one_hamming_four_collar_closure_codex_S10.out
  - 04-computation/lrc13_hamming_six_nonprimitive_contraction_scout_codex_S11.cpp
  - 05-knowledge/results/lrc13_hamming_six_nonprimitive_contraction_scout_codex_S11.out
  - 04-computation/lrc13_hamming_six_primitive_ap_pin_scout_codex_S14.py
  - 05-knowledge/results/lrc13_hamming_six_primitive_ap_pin_scout_codex_S14.out
  - 04-computation/lrc13_hamming_six_open_f3_exact_closure_codex_S14.cpp
  - 05-knowledge/results/lrc13_hamming_six_open_f3_exact_closure_codex_S14.out
  - 04-computation/lrc13_hamming_six_f2_rank_scout_codex_S14.cpp
  - 05-knowledge/results/lrc13_hamming_six_f2_rank_scout_codex_S14.out
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
the first lift.  Part C by itself made no claim of radius-five or radius-six
emptiness.  Downstream THM-845 closes radius five, and the addendum below closes
only the nonprimitive slice of radius six.  Primitive radius six remains open
outside the partial reduction in C.2 below.
The density barrier is a limitation of this method, not a proof that the
radius-seven chart is intrinsically infinite.

### C.1 Addendum: only `2[12]` survives in nonprimitive scale-one H6

Let

```text
W=([12] minus R) union {r+13h_r:r in R},
|R|=6, h_r>=1.                                         (C1)
```

If `gcd(W)>1`, then all six retained labels in `[12] minus R` are divisible
by that gcd.  Since `[12]` contains fewer than six multiples of every integer
at least three, the gcd is two and the retained set is exactly
`{2,4,6,8,10,12}`.  Thus

```text
R={1,3,5,7,9,11}.                                      (C2)
```

The retained label `2` makes the gcd exactly two.  For odd `r`, the lift
`r+13h_r` is even exactly when `h_r` is odd, so nonprimitivity is equivalent
to (C2) together with all six heights odd.  Write

```text
r=2i-1,                 h_r=2k_i+1,       1<=i<=6.
```

Division by two gives the exact contraction

```text
W/2={1,...,6} union {6+i+13k_i:1<=i<=6}.               (C3)
```

If no `k_i` is positive, then `W=2[12]`.  If between one and five are
positive, (C3) is a proper scale-one Hamming-`j` packet with `1<=j<=5`, hence
is loose by THM-795/800/804/806/815/845.  It remains only to close the single
top-half chamber in which all six `k_i` are positive:

```text
{1,...,6} union {r+13h_r:7<=r<=12, h_r>=1}.            (C4)
```

Apply Part C's longest-component cap separately at every numerically ordered
prefix of (C4).  The exact recursion has the following depth census; `edges`
is the number of legal next labelled lifts, and `dead` is the number of
prefixes with none:

```text
depth       0       1        2         3       4     5  6
nodes       1      54     3612    130515    2104     2  0
edges      54    3612   130515      2104       2     0  -
dead        0       0        0    129772    2103     2  -
max cap   132     430      683       608     315    28  -
```

There are zero covering prefixes at every depth.  The two depth-five leaves
are

```text
(7:33,10:36,9:48,8:60,11:63),  unused label 12,
    L=17/3120, cap=28, least legal next speed=64;
(7:33,10:36,9:48,8:60,12:64),  unused label 11,
    L=47/8580, cap=28, least legal next speed=76.        (C5)
```

Thus neither leaf has an outgoing edge, and no hypothetical tight completion
of (C4) exists.  Independently reconstructing the strict-safe sets of both
leaves as complements of the full closed-danger union gives exact endpoint-
for-endpoint agreement with the intersection recursion.  The complete tree
has `136,288` prefixes, trace
`919c6848d4e1187a2cef093e58982ae6`, and `313` cached speeds.

There is also a useful equality-side Kakeya ledger.  The strict-safe set
`E_{1,...,6}` has twelve components, total measure `27/65`, and longest
component `1/13`.  The six periodic danger combs `D_7,...,D_12` cover all
twelve components.  Their component--comb incidence graph has 34 edges and
component-degree histogram `1:4,2:2,3:2,4:2,6:2`.  Exactly four components
have zero overlap debt and a unique full owner:

```text
[7/39,12/65] and [53/65,32/39]       owner 11,
[14/65,3/13] and [10/13,51/65]       owner 9.           (C6)
```

This is a component/needle rigidity fingerprint of the AP equality, not a
tournament proof.  Pairwise comb overlap is symmetric; an arbitrary switch
would forget the union-cover predicate.  The faithful carrier is the weighted
bipartite incidence of safe components with danger combs.

The `-O3` and `-O0` builds produce byte-identical output, and an
AddressSanitizer/UndefinedBehaviorSanitizer build is clean.  The source
asserts every frozen depth counter and cap, both trace halves, the cache size,
and both deepest paths.  The frozen hashes are

```text
source  ee57510a4796e23da1408b383af1146478f067a5e0e98d2ad52220e55e2e8bf1
output  aa87c107ff7e5fdc6cb6a3803ecc0751ac41c267dd59540e4bd6dd764f7769ed
```

Consequently, among **nonprimitive** scale-one Hamming-six packets, `2[12]`
is the only possible tight packet.  This addendum makes no claim about
primitive scale-one H6 packets, arbitrary-scale H6 packets, or global `n=12`
sporadic emptiness. ∎

### C.2 Addendum: all primitive-core three-antipodal-pair H6 rows are loose

Continue with (C1), but now restrict to missing-label rows whose retained core
is primitive, `gcd(P)=1`.  Put

```text
P=[12] minus R,
f(R)=#{r in {1,...,6}:{r,13-r} subset R}.              (C7)
```

Exactly one of the `C(12,6)=924` missing-label rows has a nonprimitive retained
core, namely (C2).  The other 923 primitive-core rows split as follows:

```text
f                  0      1      2      3
number of rows    63    480    360     20.              (C8)
```

This is deliberately a statement about the retained core, not merely the gcd
of the completed packet.  Every packet over these 923 cores is primitive.
Conversely, on the exceptional row (C2), mixed height parity makes `W`
primitive even though `P` has gcd two.  Those mixed-parity packets are outside
this addendum and remain open; C.1 treated only the all-odd-height,
nonprimitive branch.

This statistic has a direct AP-grid meaning.  For `a in {1,...,12}`, the
point `a/13` belongs to the strict-safe set `E_P` exactly when the antipodal
pair `{a^(-1),-a^(-1)}` is contained in `R`.  Hence precisely `2f(R)` of the
twelve nonzero thirteenth-grid points lie in `E_P`.

#### The oriented germ handoff lemma

Fix a full pair and orient one of its labels as `r`.  Set `a=r^(-1) mod 13`,
and write `[x]` for the representative of `x mod 13` in
`{-6,...,-1,1,...,6}`.  Define

```text
c(z)=z-1 if z>0,                 c(z)=12+z if z<0,
ell_r(P)=min_(p in P) c([ap])/p, B_r(P)=2/ell_r(P).     (C9)
```

The left-hand germ of `E_P` at `t_0=a/13` has exact reach
`ell_r(P)/13`.  The danger tooth of the owner lift `u_r` has left reach
`2/(13u_r)`.  For every other missing label `s`, its first leftward danger
tooth begins at distance `c([as])/(13u_s)` from `t_0`.  Consequently, if the
completed packet is tight, then for every oriented owner `r` either

```text
u_r <= B_r(P),                                           (C10)
```

or there is an `s in R minus {r}` such that

```text
u_s >= c([as])u_r/2.                                    (C11)
```

Indeed, if (C10) fails, the owner tooth stops strictly inside the core-safe
germ.  Unless some other tooth has begun by that point, the interval between
the two endpoints is strict-safe for the completed packet.  Comparing the two
exact reaches gives (C11).  This argument is height-independent.

For an `f=3` row all six missing labels are oriented owners.  If no owner
uses (C10), choose one provider in (C11) for every owner.  The resulting
functional digraph has a directed cycle.  Multiplication around that cycle
forces the product of its edge weights `c([as])/2` to be at most one.  Exact
enumeration of the 409 simple directed cycles in each of the twenty complete
six-vertex handoff graphs gives

```text
minimum cycle product       1/16       1       3/2
number of rows                  6       2        12.     (C12)
```

In each of the two product-one rows, every edge on either nonexpanding cycle
has weight one.  Equality in the multiplied inequalities would force two
distinct-residue lifts to be equal, which is impossible.  Thus fourteen rows
cannot avoid a boundary exit.  In six of them (C10) admits no proper lift at
all, so the rows are already loose:

```text
{1,2,3,10,11,12}   {1,3,4,9,10,12}   {1,3,5,8,10,12}
{1,4,6,7,9,12}     {2,3,4,9,10,11}   {2,3,5,8,10,11}.  (C13)
```

For the other eight rows, enumerating the proper lifts `u_r=r+13h_r` below
the exact boundary caps in (C10) leaves only the following alternatives:

```text
R={1,2,5,8,11,12}     => u_5=18
R={1,4,5,8,9,12}      => u_5=18
R={1,5,6,7,8,12}      => u_5=18
R={2,3,6,7,10,11}     => u_6=19
R={2,4,6,7,9,11}      => u_6=19
R={2,5,6,7,8,11}      => u_5=18 or u_6=19
R={3,4,6,7,9,10}      => u_6=19
R={4,5,6,7,8,9}       => u_5=18 or u_6=19.             (C14)
```

The six rows not decided by the germ-cycle argument are exactly

```text
{1,2,4,9,11,12}   {1,2,6,7,11,12}   {1,3,6,7,10,12}
{2,4,5,8,9,11}    {3,4,5,8,9,10}    {3,5,6,7,8,10}.   (C15)
```

Their minimum cycle product is `1/16`; the germ-cycle quotient alone makes no
looseness claim for them.  The exact residual-component recursion below closes
all six.

#### Exact exhaustion of the ten boundary slices

There are ten fixed-coordinate slices in (C14), because each of the two
disjunctive rows contributes two.  In a fixed slice, adjoin the displayed
speed to `P`, then numerically order the five remaining lifts.  At a depth-`j`
prefix, let `L_j` be the longest component of its strict-safe set.  If a tight
completion existed, the next lift would obey (8) with `m=5-j`.  Enumerating
every unused residue class, every proper congruent lift above the preceding
one, and every lift up to this exact cap therefore includes every tight
completion exactly once.  The fixed speed need not be the smallest speed; it
is already part of the prefix, while the other five are sorted among
themselves.

The ten exact trees are:

```text
R                  fixed    root L      depth nodes          depth edges
12581112           u5=18    19/585      1,28,301,119,11,0    28,301,119,11,0
1458912            u5=18    101/2574    1,22,196,119,0,0     22,196,119,0,0
1567812            u5=18    5/117       1,21,194,95,1,0      21,194,95,1,0
23671011           u6=19    21/494      1,20,170,130,5,0     20,170,130,5,0
2467911            u6=19    11/247      1,20,166,66,0,0      20,166,66,0,0
3467910            u6=19    87/2717     1,30,334,166,4,0     30,334,166,4,0
456789             u5=18    7/156       1,20,171,196,14,0    20,171,196,14,0
456789             u6=19    11/247      1,20,172,80,3,0      20,172,80,3,0
2567811            u5=18    19/468      1,21,192,115,3,0     21,192,115,3,0
2567811            u6=19    4/117       1,26,280,150,8,0     26,280,150,8,0
aggregate                              10,228,2176,1236,49,0 228,2176,1236,49,0
```

All ten trees have zero covering prefixes and no depth-five state.  At every
one of their 3,699 states, an independent complement-of-closed-danger-union
construction reproduces the recursive strict-safe interval union exactly.
The frozen tree trace is

```text
f8c8465455fdbf1f21aec0438a8894c11f56accf53579dc84658a703aa5c227e. (C16)
```

Thus each single forced slice is empty, and the union of the two slices
exhausts each disjunctive row in (C14).  All eight rows in (C14) are loose.
Together with (C13), this closes fourteen of the twenty primitive-core `f=3`
rows by local reduction plus ten small slice trees.

#### Exact exhaustion of the six product-`1/16` rows

For each row in (C15), start from its six-speed retained core and numerically
order all six proper replacement lifts.  At a depth-`j` prefix, let `L_j` be
the longest component of the exact residual strict-safe set.  If a tight
completion existed, its `6-j` remaining danger combs would cover that
component, so (8) gives the exact next-speed cap

```text
x_(j+1) <= floor(22(6-j)/[13(13-2(6-j))L_j]).          (C17)
```

The denominator is positive through all six levels.  Enumerating every unused
residue label and every proper congruent lift above the preceding speed and at
most (C17) therefore contains every hypothetical tight completion exactly
once.  The six frozen recursion trees are

```text
R                  depth nodes                         trace128
12491112           1,64,4769,195705,7340,50,0          063c1b3f17d4520feaf265b370a1c09d
12671112           1,64,4761,195502,5875,12,0          3b468adfcff92e6bafc5cf50ea6b5000
13671012           1,97,10343,620068,17195,39,0        4844a2fbb06e7c9aad37d047eefeecab
2458911            1,93,9560,550797,10885,21,0         a60c4579e7f6cff7404360e8db9bc6f4
3458910            1,79,7032,349248,10523,10,0         144d8b3b67d164cdb04afbf0bafffa3c
3567810            1,97,10348,620350,22597,70,0        242306835e1bb0c1289ce8d030c29a7a
aggregate          6,494,46813,2531670,74415,202,0
```

The aggregate edge counts are

```text
494,46813,2531670,74415,202,0,                         (C18)
```

and the trees contain `2,653,600` states.  Every row has zero covering prefix,
and none reaches depth six.  Independently taking the complement of the full
closed-danger union reconstructs all six root safe sets and all 202 deepest
dead leaves endpoint-for-endpoint, for 208 independent checks.  The source
asserts every row's depth counts, dead counts, maximum caps,
minimum-longest-component values, trace halves, cache size, and crosscheck
count.

Thus all six rows in (C15) are loose.  Together with (C13)--(C14), **all twenty
primitive-core `f=3` rows are closed**, reducing the primitive-core
missing-label frontier from 923 rows to the 903 rows with `f<=2`.  The
exceptional odd-label row still has a separate primitive mixed-parity branch,
so 904 missing-label patterns retain open primitive assignments.  The packet
`2[12]` remains the nonprimitive AP equality from C.1.  This is not a proof of
primitive scale-one H6 emptiness or global `n=12` sporadic emptiness. ∎

### C.3 Addendum: the first two-antipodal-pair cap stratum is loose

Now take the 360 primitive-core missing-label rows with `f=2`.  For a row `R`
and retained core `P=[12] minus R`, let `L(P)` be the exact longest component
of `E_P`.  Before any replacement is placed, (C17) is

```text
C_1(R)=floor(132/(13L(P))).                              (C19)
```

The exact root census has `113<=C_1(R)<=396`.  Exactly eleven rows satisfy
`C_1(R)<=132`; in lexicographic atlas indexing they are

```text
5, 47, 62, 71, 74, 76, 77, 128, 158, 331, 346.          (C20)
```

Run the same numerically ordered, labelled longest-component recursion used
above, now on each entire row rather than on a fixed-coordinate slice.  The
frozen aggregate depth census is

```text
11, 599, 39415, 1414388, 62443, 229, 0,                 (C21)
```

or `1,517,085` states.  Every covering count is zero and no state reaches
depth six.  Per-row node counts, 128-bit traces, and cache sizes are frozen;
the complement of the full closed-danger union independently reconstructs
all eleven roots and the 229 deepest dead leaves, giving 240 crosschecks.

Consequently all eleven primitive-core `f=2` rows with `C_1<=132` are loose.
The primitive-core frontier falls from 903 to 892 rows.  Adding the exceptional
odd-row mixed-parity pattern gives 893 open scale-one H6 label patterns;
`2[12]` remains the separate nonprimitive equality.  This closes a root-cap
sub-stratum, not all `f=2` rows, arbitrary-scale H6, or the global sporadic
branch. ∎

The statistic `f` should not be mistaken for a geometric quotient.  It is
invariant not only under reflection `r -> 13-r`, but under every multiplication
by a unit modulo 13 and, for six-element rows, under `R <-> P` (full and empty
antipodal pairs exchange in equal numbers).  It preserves the number of
antipodal AP cusps and hence the number of forced germ owners, but not their
identities.  It also destroys the signs of the pairs, the retained speeds,
exact safe endpoints and longest components, divisor obligations, height
congruences, and Kakeya overlap.
These label operations need not even preserve the primitive-core domain, and
they are not common integer dilations of a normalized scale-one packet.

There is a precise tournament-related lesson.  The pairwise observable is the
provider start `c([as])/(13u_s)` versus the owner reach `2/(13u_r)`, the switch
is (C11), and the gauge is the oriented left germ at `a=r^(-1)`.  Thresholding
does not produce a tournament: a vertex pair can carry both arrows or neither
because its two orientations use different gauges.  Accordingly score
histograms, SCCs, and Hamiltonian-path counts are not proof invariants here;
the weighted directed-cycle products in (C12) carry the fourteen-row local
reduction but stop at (C15).  The increasing `(speed,label)` path used in the
exact recursions is a transitive enumeration tournament, with score sequence
`0,...,5`, no directed cycles, six singleton SCCs, and one Hamiltonian path.
It is not a quotient of the cover predicate.  The faithful global carrier that
closes (C15) and the first `f=2` cap stratum is the residual endpoint union
together with the labelled bank of unplaced danger combs.

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
