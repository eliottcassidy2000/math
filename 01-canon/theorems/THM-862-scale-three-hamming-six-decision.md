---
id: THM-862
title: Scale-three Hamming-six common-sheet classification and exact terminal closure
status: PROVED STRUCTURAL + FINITE-EXACT — the primitive c=3 common-sheet bank has 212 presentations and 1,504 unit contexts with exact affine toothpick codes; the complete no-height-cutoff metric recursion visits 950,540,566 candidate edges, finds zero covering terminals and two certified nonempty terminals, and therefore closes every primitive proper AP-centred c=3 Hamming-six packet as strictly loose. This theorem closes only that face; THM-957/958/960/962/963/969/970 subsequently close c=4,...,10, while c>=11, H5 ramification, non-AP-centred/deep-sheet branches, and global n=12 sporadic emptiness remain open.
source: codex-2026-07-15-S15/S16 c=3 transport audit; terminal closure codex-2026-07-17-S57
depends_on: [THM-810, THM-815, THM-823, THM-857, THM-859, THM-860, THM-861]
related: [THM-844, THM-847, THM-858, THM-957, THM-958, THM-960, THM-962, THM-963, THM-969, THM-970, HYP-6820]
verification:
  - 04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out
  - 04-computation/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.out
  - 04-computation/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp
  - 04-computation/lrc13_scale_three_hamming_six_depth_two_combine_codex_S16.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.out
  - 04-computation/lrc13_scale_three_hamming_six_depth_two_geometry_cache_crosscheck_codex_S16.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_depth_two_geometry_cache_crosscheck_codex_S16.out
  - 04-computation/lrc13_scale_three_hamming_six_depth_three_geometry_batch_codex_S16.cpp
  - 05-knowledge/results/lrc13_scale_three_hamming_six_depth_three_geometry_batch_codex_S16.out
  - 05-knowledge/results/lrc13_scale_three_hamming_six_depth_three_benchmark_codex_S16.md
  - 04-computation/lrc13_scale_three_hamming_six_terminal_combine_codex_S57.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_terminal_census_codex_S57.out
  - 04-computation/lrc13_scale_three_terminal_packet_certificate_codex_S57.py
  - 05-knowledge/results/lrc13_scale_three_terminal_packet_certificate_codex_S57.out
---

# THM-862 — the scale-three sheet stalk is classified

Let `R subset F_13^*` have six elements and put

```text
A=3([12] minus R) union {w_r:r in R},
w_r=3r (mod 13),                 w_r>0,       w_r!=3r.   (1)
```

Assume that `A` is primitive and `M(A)<=1/13`.  At common scale three the
effective order of a replacement is

```text
D_r=3/gcd(3,w_r) in {1,3}.                              (2)
```

For `D_r=3`, write `e_r=w_r mod 3 in {+1,-1}`.  An order-one unit is
trivial.  THM-810/823's germ argument makes common-sheet coverage a necessary
condition for (1).  This theorem classifies that condition completely and now
runs every resulting arithmetic continuation to an exact terminal verdict.

## Theorem

### A. Exact common-sheet bank

After excluding the all-order-one word, which would make every speed in (1)
divisible by three, the complete primitive `c=3` common-sheet bank is

| effective orders | presentations | unit contexts |
|---|---:|---:|
| `1^2 3^4` | 84 | 336 |
| `1 3^5` | 84 | 672 |
| `3^6` | 44 | 496 |
| **total** | **212** | **1,504** |

There are no rows with one, two, or three order-three colours.  In the
all-order-three stratum, `26` presentations have eight unit words and `18`
have sixteen, so

```text
26*8+18*16=496.                                        (3)
```

### B. The signed-provider law and toothpick codes

Put

```text
H={1,5,8,12},
X_-={4,5},                       X_+={8,9}.              (4)
```

At an order-three owner `o`, the owner itself supplies normalized sheet zero.
An off-owner order-three label `r` supplies a nonzero sheet exactly when

```text
r/o in X_- union X_+.                                  (5)
```

It supplies the normalized sign

```text
-e_r  if r/o in X_-,                 +e_r if r/o in X_+.(6)
```

Thus an order-three owner is covered on all three sheets exactly when both
signs occur among (6).  This is one signed not-all-equal provider obligation
per owner.  An order-one colour covers every sheet at its own owner and none
at a distinct replacement owner.

On every surviving presentation the apparently higher-arity obligations
collapse further: the valid unit words form an affine binary code cut out by
disjoint signed pair equations.  Write `e_i e_j=+` for equal signs and
`e_i e_j=-` for opposite signs.  The mixed rows normalize to

```text
C=H:             e_1 e_12=+,   e_5 e_8=+;              (7)
C=H union {2}:   e_1 e_12=+,   e_5 e_8=+,  e_2 free.   (8)
```

For six order-three labels there are exactly five multiplicative orbits:

| normalized `C` | orbit size | affine toothpicks | free labels | words | contexts |
|---|---:|---|---|---:|---:|
| `{1,2,3,4,9,10}` | 12 | `e1e2=e4e9=e3e10=+` | none | 8 | 96 |
| `{1,2,3,5,8,12}` | 12 | `e5e8=e1e12=+` | `2,3` | 16 | 192 |
| `{1,2,5,6,8,11}` | 12 | `e1e6=-`, `e5e8=e2e11=+` | none | 8 | 96 |
| `{1,2,5,8,11,12}` | 6 | `e5e8=e1e12=+` | `2,11` | 16 | 96 |
| `{1,3,4,9,10,12}` | 2 | `e4e9=e3e10=e1e12=+` | none | 8 | 16 |

The last row is the quadratic-residue six-set and its nonresidue mate.  The
pair equations in (7)--(8) and the table are exact descriptions of the full
unit fibres, not only necessary parity tests.  This matching-code form is a
scale-three **toothpick code**: every live sheet stalk is assembled from
disjoint signed two-vertex constraints and free stalks.  The name refers only
to this matching decomposition; it does not revive THM-841's false local
toothpick recurrence for the Farey ladder.

### C. Exact multiplicative and reflection orbits

Let `T_a`, `a in F_13^*`, multiply every replacement label by `a` while
leaving its order and attached unit fixed.  Let

```text
J:e_r -> -e_r                                             (9)
```

on every order-three colour.  Both actions preserve common-sheet coverage and
commute, so the sheet action is `F_13^* x <J>`, of order twenty-four.

The exact orbit census is:

| objects | `1^2 3^4` orbit sizes | `1 3^5` orbit sizes | `3^6` orbit sizes | total orbits |
|---|---|---|---|---:|
| presentations under `F_13^*` | `6x12, 2x6` | `7x12` | `3x12, 1x6, 1x2` | 20 |
| contexts under `F_13^*` | `24x12, 8x6` | `56x12` | `36x12, 10x6, 2x2` | 136 |
| contexts under `F_13^* x <J>` | `12x24, 4x12` | `28x24` | `18x24, 5x12, 1x4` | 68 |

The reflection `J` is free, and no multiplier composed with `J` stabilizes a
context.  Hence it pairs the `136` multiplier-only orbits into `68` sheet
orbits.

These `68` orbits have the following compact canonical schemes.  Unit words
use digits `1,2` in increasing order of the displayed order-three set.

```text
1^2 3^4:
  C={1,5,8,12},
  B in {{2,3},{2,4},{2,6},{2,7},{2,9},{2,11},{4,6},{4,9}},
  e in {1111,1221};                                     16 reps

1 3^5:
  C={1,2,5,8,12},
  B in {3,4,6,7,9,10,11},
  e in {11111,11221,12111,12221};                       28 reps

3^6:
  C={1,2,3,4,9,10},     e in {111111,111221,112112,112222};
  C={1,2,3,5,8,12},     e in {111111,111221,112111,112221,
                                121111,121221,122111,122221};
  C={1,2,5,6,8,11},     e in {111211,112221,121212,122222};
  C={1,2,5,8,11,12},    e in {111111,111121,112211,112221,
                                121121,122221};
  C={1,3,4,9,10,12},    e in {111111,112211}.           24 reps
```

The two `1^2 3^4` label pairs `{2,11}` and `{4,9}` have stabilizer
`{+/-1}`.  In the fourth all-order-three row, the words
`111111,112211,121121,122221` have the same stabilizer; its other two words
have trivial stabilizer.  On the quadratic-residue row, `111111` has the
six-element quadratic-residue stabilizer and `112211` has `{+/-1}`.  All
other displayed context representatives have trivial stabilizer.

Multiplicative inversion is **not** an action on this bank.  It maps back only

```text
86/212 presentations = all 84 rows of type 1^2 3^4 plus 2 QR rows,
352/1504 contexts     = all 336 contexts of that type plus 16 QR contexts.
                                                               (10)
```

It maps no `1 3^5` row back to the bank.  A dihedral multiplicative quotient
would therefore be false.

### D. Sheet symmetry is not metric symmetry

The `68` sheet representatives cannot replace the `1,504` arithmetic
contexts in an exact component recursion.  An explicit all-order-three
context supplies the guardrail.  For

```text
C={1,2,3,4,9,10},                     e=111111,          (11)
```

the least proper packet, its unit reflection, and its label-doubled image are

```text
Q ={1,4,15,16,18,19,21,22,24,25,33,36},
JQ={14,15,17,18,21,24,29,32,33,35,36,38},
T2Q={3,9,19,25,27,28,30,31,33,34,36,37}.               (12)
```

Exact breakpoint evaluation gives

```text
M(Q)=1/4,             M(JQ)=7/26,             M(T2Q)=14/65. (13)
```

The respective witness pairs are

```text
{7/20,13/20},         {1/52,51/52},           {32/65,33/65}. (14)
```

Thus both `T_a` and `J` are only sheet-equivariant.

### E. Exact first metric layer and terminal recursion

For a fixed context, the complete proper replacement rays are

```text
D_r=1:  w_r=3r+39k,                     k>=1,
D_r=3:  w_r=u_r+39k,                    k>=0,           (15)
```

where `u_r` is the CRT representative in `[1,38]` satisfying

```text
u_r=3r (mod 13),                    u_r=e_r (mod 3).     (16)
```

Numerically order the six replacements.  At a prefix with `j` replacements,
longest strict-safe component length `L_j`, and `s=6-j` later combs, every
tight continuation obeys THM-815's cap

```text
x_(j+1)<=floor(22s/[13(13-2s)L_j]).                    (17)
```

There are `110` distinct retained-core roots among the `1,504` contexts.  By
stratum there are `84,66,44` distinct roots; these sets overlap across
strata.  Direct exact root-component evaluation gives:

| type | contexts | first `D=1` edges | first `D=3` edges | all first edges | per-context range |
|---|---:|---:|---:|---:|---:|
| `1^2 3^4` | 336 | 10,440 | 22,268 | 32,708 | 53--161 |
| `1 3^5` | 672 | 10,200 | 54,396 | 64,596 | 66--152 |
| `3^6` | 496 | 0 | 49,608 | 49,608 | 68--183 |
| **total** | **1,504** | **20,640** | **126,272** | **146,912** |  |

The exact root real-cap range on this survivor bank is

```text
4752/13 <= B_root <= 1188.                              (18)
```

The `146,912` count and (18) are theorem outputs.  They are not sampled height
cutoffs.

For execution planning only, scale THM-861's frozen `c=2` ratios by this
first layer.  Its depth-two ratio predicts

```text
146912*(641866/6266)=47148908896/3133
                         about 15,049,125 depth-two nodes, (19)
```

while its logical nodes per root predict

```text
1504*(41882982/64)=984,250,077 logical nodes.            (20)
```

Equations (19)--(20) were workload estimates, not `c=3` censuses.  They led to
the sharded terminal run in Part H.  That run retains every arithmetic lane and
uses THM-857's complete-tooth and streaming-cap certificates; its exact total
is slightly below the root-based estimate.

Two exact representation reductions and one shortcut guardrail refine that
plan.  Before the first order-three insertion every speed is divisible by
three.  For `Q=3Q'`,

```text
E(Q)={t:3t mod 1 lies in E(Q')},             L(Q)=L(Q')/3. (20a0)
```

The THM-815 cap and the order-one ray `3(r+13k)` scale by the same factor, so
root and order-one-only prefixes may be computed in the scale-one quotient;
literal scale-three components must be materialized at the first order-three
insertion.  Also every residual is invariant under `t -> 1-t`, so exact
half-circle storage is legal provided the `1/2` cut and all open endpoint flags
are retained.  Neither reduction identifies arithmetic lanes.

The generic complete-safe-tooth implication from tooth containment alone is
automatic only once at most three combs remain.  Indeed, a full retained safe
tooth for the just-inserted speed `x` has length `11/(13x)`, so the next cap is
at most

```text
2s x/(13-2s)<x                         exactly when s<=3. (20a1)
```

At the first two insertion depths (`s=5,4`), containment alone must not trigger
the shortcut; an exact larger child component could still make the numerical
cap test pass.  The streaming-cap certificate remains valid there.  This
distinction is essential for a proof-facing depth-two scout.

### F. Exact first-child geometry cache and the affine Walsh projector

There is an exact batching map, but it is a cache map rather than a quotient
of proof states.  For a cap-admissible first edge let `R` be the sorted
missing-label six-set and `x_1` the inserted speed.  Then

```text
geometric key              G=(R,x_1),
logical lane key           Lambda=(R,x_1; five remaining step-39 rays). (20a)
```

The literal child `E(3([12] minus R) union {x_1})`, all its endpoints, its
longest component and therefore its next THM-815 cap depend only on `G`.
Across the `146,912` logical first edges there are exactly

```text
22,262 geometric keys,
fibre multiplicities {2:2454, 4:8753, 6:6479, 12:2375, 18:2201}. (20b)
```

Thus an implementation may compute and cache the expensive literal child
only once per `G`, a factor

```text
146912/22262=73456/11131 about 6.60                         (20c)
```

at the first layer.  It must still propagate every lane in the fibre: after
adjoining the exact five future progressions to `G`, all `146,912` lane keys
in (20a) are distinct.  This is not merely a label artefact.  A ray is the
actual set `{b+39k:k>=0}`; its base modulo thirteen recovers its owner label.
Consequently distinct remaining-base tuples are distinct arithmetic future
languages.  Recursively, `(R,x_1,...,x_j)` is always a valid component-cache
key, while the remaining labelled rays stay in the proof state.

The multiplicities in (20b) are the partition functions of Part B's matching
codes under one pinned local unit.  If `M` is a disjoint signed matching, put
`e_i in {+1,-1}` and give an edge `ij` sign `s_ij`.  Its exact code indicator
is

```text
1_C(e)=2^(-|M|) product_(ij in M)(1+s_ij e_i e_j).       (20d)
```

Hence, for the normalized Walsh transform `fhat`,

```text
E_(e in C) f(e)
 =sum_(A subset M) (product_(ij in A)s_ij) fhat(V(A)),    (20e)
```

where `V(A)` is the union of the endpoints of the chosen matching edges.
The Walsh support is exactly the unions of toothpicks.  A three-edge code has
degree histogram `1,3,3,1` in degrees `0,2,4,6`; a two-edge code has
`1,2,1` in degrees `0,2,4`.  This gives the exact Fourier meaning of the
`26` eight-word and `18` sixteen-word all-order-three presentations.

Conditioning is recursively self-similar on the **code** side.  Pinning one
unit halves every fibre.  Pinning both endpoints of one matching edge imposes
no second independent condition, whereas pins in different matching/free
components halve it again.  For one presentation this gives conditioned
sizes `2`, `4`, or `4/8`; the larger multiplicities `6,12,18` in (20b) are
sums of compatible order presentations over the same `(R,x_1)`.  This is the
rigorous toothpick recursion that survives.  It organizes conditional unit
enumeration but does not identify the metric lanes.

One useful corollary is that every separable one-ray statistic
`f(e)=sum_i f_i(e_i)` has its code average equal to its unrestricted uniform
average: (20e) has no nonconstant degree-zero-or-one contribution.  The
first-edge **aggregate count** is such a statistic for a fixed presentation.
The maximin, literal component erosion, and terminal cover indicator are not;
the Walsh identity therefore supplies no emptiness theorem by itself.

Finally, ordering the six least proper ray bases gives `1,151` distinct
root numeration clocks among the `1,504` contexts.  Orienting one ray toward
another when its least base is smaller produces only the transitive
tournament fingerprint

```text
scores {0,1,2,3,4,5}, triangles 0, SCCs 1^6,
Hamiltonian paths 1.                                    (20f)
```

It is the priority-queue clock for the numerical recursion, not a cover
carrier.  It forgets the base values, the literal components, and the five
future languages, so equal clocks cannot be merged.

### G. Exact depth-two metric census

The complete sharded depth-two scout has now replaced the workload estimate
in (19) by an exact census.  Across all `1,504` arithmetic contexts,

```text
depth-zero contexts                         1,504
cap-admissible first edges                146,912
cap-admissible second edges            14,992,263.       (20g)
```

The second layer splits by sheet stratum and effective order as follows:

| stratum | depth-two nodes | `D1->D1` | `D1->D3` | `D3->D1` | `D3->D3` |
|---|---:|---:|---:|---:|---:|
| `1^2 3^4` | 3,408,353 | 225,672 | 901,620 | 911,740 | 1,369,321 |
| `1 3^5` | 6,469,464 | 0 | 1,069,716 | 1,078,764 | 4,320,984 |
| `3^6` | 5,114,446 | 0 | 0 | 0 | 5,114,446 |
| **total** | **14,992,263** | **225,672** | **1,971,336** | **1,990,504** | **10,804,751** |

Equivalently the second insertion has order-one/order-three totals

```text
2,216,176 + 12,776,087 = 14,992,263.                     (20h)
```

No root or first-edge lane is dead, and no prefix through depth two covers.
Exactly one depth-two prefix has no cap-admissible third insertion.  In
context `1448`,

```text
R=(4,6,7,8,9,12),
D=(3,3,3,3,3,1),              e=(2,1,1,1,2,0),
(label,speed) insertions=(9,14),(4,38).                   (20i)
```

Its longest child component is

```text
(183/494,56/143),       length=115/5434.
```

With four combs remaining, (17) gives cap `63`, while the four least legal
future speeds are `70,73,76,75`.  Thus this lane is rigorously dead.  Every
other one of the `14,992,262` depth-two lanes reaches the depth-three
frontier.  The depth-two longest-component minimum is `11/51389`, with
multiplicity two, and the largest next cap is `6,324`.

Literal geometry still gives a cache, not a state quotient.  Grouping by

```text
G_2=(R,x_1,x_2)                                          (20j)
```

reduces the second layer to exactly `4,307,561` geometric children, with
fibre multiplicities

```text
{1:212990, 2:2123879, 3:535281, 4:337242,
 6:892407, 9:164606, 18:41156}.                           (20k)
```

The cache factor is `14992263/4307561`, about `3.48`.  Distinct lanes in a
fibre still carry different unused labelled step-39 rays and cannot be merged
for continuation.  The numerical-order tournaments remain transitive at all
`146,912` first prefixes, but conditioning flips `552,554` pair orientations;
this is scheduling telemetry, not cover data.

The primary C++ engine reconstructs all contexts independently from literal
CRT masks, materializes every depth-two child exactly, and hash-combines eight
canonical context shards.  A separate one-core Python implementation builds
the `110` roots as complements of merged closed danger combs, reconstructs
all `22,262` first children, and obtains (20g)--(20k) without using the C++
interval engine.  It also independently rebuilds the unique dead geometry
and its cap certificate.  The exact count is only `56,862` below the estimate
in (19), a relative error under `0.38%`.

This near-total depth-two survival was a negative computational result with a
positive design consequence: the full recursion had to retain arithmetic
lanes while exploiting only proof-preserving geometric shortcuts.  Part H now
records the completed run; nothing in the depth-two census alone implied that
terminal verdict.

### H. Exact terminal census and scale-three closure

The generic engine applies (17) at every expanded prefix and enumerates all six
labelled step-39 replacement rays in increasing numerical order.  Its complete
depth-zero through depth-six node vector is

```text
(1,504, 146,912, 14,992,263, 931,412,556,
 3,984,352, 4,481, 2),                                  (20l)
```

with `950,540,566` candidate edges.  The sheet-stratum split is

| stratum | contexts | nodes at depths `0,...,6` | candidate edges | covers | nonempty terminals |
|---|---:|---|---:|---:|---:|
| `1^2 3^4` | 336 | `336,32708,3408353,214874008,691319,758,0` | 219,007,146 | 0 | 0 |
| `1 3^5` | 672 | `672,64596,6469464,391549538,1609804,1240,0` | 399,694,642 | 0 | 0 |
| `3^6` | 496 | `496,49608,5114446,324989010,1683229,2483,2` | 331,838,778 | 0 | 2 |

Every shortcut has a one-sided proof meaning.  A contained complete safe band
after the third or later insertion leaves a component too long for the
remaining numerically larger combs; this certifies `879,373,305` depth-three
nodes.  The streaming version of (17) certifies the vector

```text
(0,0,1,50,593,411,3,980,598,4,479,2).                  (20m)
```

The single depth-two entry is exactly (20i).  At depth six the same exact
intersection routine detects a nonempty interval rather than materializing an
unneeded full list.  The two such terminals are contexts `888` and `1502`,
both in the all-order-three stratum.  Fully ungated replay of context `888`
reproduces its complete node/dead/edge vector and the same loose verdict;
context `1502` has a separate ungated referee in the S57 audit.

For a concrete terminal certificate, context `888` follows the labelled path
`9:1,10:4,2:19,3:22,12:23,4:25` and ends at

```text
{1,3,4,15,18,19,21,22,23,24,25,33}.
```

An independent closed-danger complement reconstruction gives exactly twenty
safe components, and its exact maximin is `1/9`, attained at
`4/27,8/27,19/27,23/27`.  Thus the rare terminal survivor is not a near-cover
or an accounting ambiguity; it is substantially looser than `1/13`.

Most importantly,

```text
covering terminals = 0.                                (20n)
```

The largest discrepancy cap encountered is only `6,324`, and the largest
enumerated candidate speed is `6,317`.  These extrema are recomputed per row,
max-aggregated by each shard and the independent combiner, and lie far below
the checked 32-bit carrier boundary.

The recursion has no lift-height cutoff and uses exact rational endpoints.
If a terminal residual is nonempty, it contains a time at which all twelve
speeds have strict clearance above `1/13`, so that packet is loose.  If it were
empty, it would be a covering/tight terminal.  THM-815's cap proves that every
tight continuation occurs among the enumerated rays.  Therefore (20n) proves:

> Every primitive proper AP-centred common-scale-three Hamming-six packet
> satisfying the necessary common-sheet classification has `M(A)>1/13`.

The four root-aligned shards cover contexts
`0:100`, `100:508`, `508:1004`, and `1004:1504`.  The independent combiner
rejects a duplicate or missing index, checks every row and shard sum, freezes
the full 1,504-row hash, and rederives (20l)--(20n).  The exact ray-order
tournament is transitive at each materialized prefix—zero directed cycles,
singleton SCCs, and one Hamiltonian path—but no tournament quotient is used.
The proof state remains the literal component union, unused labelled rays,
last speed, and shortcut ancestry.

## Proof

### 1. Derivation of the signed-provider law

Fix an order-three provider `r`, owner `o`, and put `a=r/o in F_13^*`.
Let `u` be the CRT class with

```text
u=3r (mod 13),                         u=e_r (mod 3).    (21)
```

On sheet `ell`, coverage asks that the centred residue

```text
z=<u(o^(-1)+13ell)>_(39)                               (22)
```

satisfy `-3<z<=3`.  Modulo thirteen, `z=3a`.  The only possibilities are

| `a` | 1 | 4 | 5 | 8 | 9 |
|---:|---:|---:|---:|---:|---:|
| `z` in `(-3,3]` | 3 | -1 | 2 | -2 | 1 |
| `z mod 3` | 0 | -1 | -1 | +1 | +1 |

The ratio `a=12` would give `z=-3` and is excluded by the left-open
orientation.  No other ratio hits the window.  Reducing (22) modulo three
gives

```text
z=e_r(o^(-1)+ell) (mod 3).
```

Since `e_r^(-1)=e_r`, the normalized owner sheet

```text
y=ell+o^(-1) (mod 3)                                  (23)
```

is precisely zero for the self-provider and (6) for an off-provider.  All
three sheets occur exactly under the signed NAE condition.

At order one, the oriented interval contains only its own nonzero residue.
After lifting to common scale three, that colour fills all three sheets at
its own owner and none at another owner.  This proves the symbolic law.

### 2. Classification and affine collapse

Let `C` be the order-three labels and `B=R minus C` the order-one labels.
If `C` is empty, every core and replacement speed in (1) is divisible by
three, contradicting primitivity.

If `|C|=4`, the four order-three colours must cover their four owners without
help from `B`.  THM-810's equality classification gives

```text
C=aH,                         e_r=e_(-r) on C.           (24)
```

There are three cosets, `C(8,2)=28` choices for `B`, and four unit words,
giving `84` presentations and `336` contexts.

If `|C|=5`, the five order-three colours form THM-823's all-order-three flag:

```text
C=aH union {b},                     b in 2aH.            (25)
```

The two antipodal equations on `aH` remain and `e_b` is free.  Choose the
coset, the forward flag, and the one order-one label outside `C`; this gives

```text
3*4*7=84 presentations,                     84*8=672 contexts. (26)
```

The signed law, evaluated over its finite unit words, rejects `|C|<=3`.
For `|C|=6`, evaluating (5)--(6) over all
`C(12,6)2^6=59,136` set/word pairs gives the five rows in Part B.  Direct CRT
masks and the signed NAE test agree.  Substituting the displayed pair
equations reproduces the complete direct-mask word set on every representative.
Multiplication supplies the stated orbit sizes, proving (3) and Part A.

The verifier checks the direct-mask/NAE identity on all

```text
sum_(k=1)^6 C(12,k)2^k=94,448                            (27)
```

order-three set/word states before decorating by `B`.  It then constructs all
presentations and contexts and independently tests every affine code.

### 3. What the sheet actions preserve

Let `q_o` be the standard integer inverse of `o mod 13`, reduced modulo three.
Equation (23) shows that `T_a` leaves every provider ratio and unit unchanged
but translates the literal mask at owner `o` by

```text
q_o-q_(ao) (mod 3).                                    (28)
```

The unit flip `J` reflects it by

```text
ell -> -ell-2q_o (mod 3).                              (29)
```

These owner-dependent affine relabellings preserve whether the three sheets
are covered, proving the group action in Part C.  Orbit-stabilizer calculation
gives the three exact rows there.  Inversion replaces provider ratios by their
inverses, but the support `{4,5,8,9}` is not inverse-stable; the exact failures
in (10) follow.

Because (28)--(29) depend on the owner, neither is induced by one circle map
on the continuous strict-safe set.  The exact values (13) prove this failure
without relying only on that observation.  The verifier obtains them from
the complete piecewise-linear candidate set with denominators
`2u`, `u+v`, and `|u-v|`.

### 4. Exact recursion and terminal certificates

The CRT gives (15)--(16).  Distinct labels have distinct residues modulo
thirteen, so every terminal has one numerical replacement order.  Before the
sixth insertion there are at most eleven speeds, so the lower-runner theorem
makes the strict-safe prefix nonempty.  THM-815 then gives (17), exactly as in
THM-857 and THM-861.

At the root, `E(3P)` is the threefold inverse image of `E(P)`, so its component
lengths are those of `E(P)` divided by three.  Enumerating (15) under the exact
root cap gives the table in Part E and the range (18).

There is one legitimate dynamic quotient.  While a numerical prefix contains
only core and order-one speeds, all its speeds remain divisible by three, so
the `Z/3Z` deck and the scale-one quotient are exact.  An order-three insertion
is not `Z/3Z`-deck-invariant, so it invalidates that quotient as a certified
continuation and must materialize literal components.  More generally,
translation by the reciprocal of the current prefix gcd, time reversal
`t->-t`, common dilation, and permutation of already inserted speeds are
genuine metric equivariances.  Multiplicative residue relabelling, unit reflection,
provider-graph isomorphism, and inversion are not.

At every prefix the next legal member of each unused ray is determined by the
last numerical speed.  Formula (17) bounds the finitely many admissible members.
If an exact child intersection is nonempty at depth six, the packet is strictly
loose.  If a child contains a component whose (17) cap is already below the
least possible next ray, no covering continuation exists; the complete-safe-
band shortcut is the special case (20a1).  These implications justify every
early verdict counted in (20m).  Induction on the six unused labels therefore
makes the depth-six search exhaustive, and the zero cover count (20n) proves
the scale-three closure.

### 5. Why the cache and projector are exact

For fixed `R` and `x_1`, the child speed set is literally
`3([12] minus R) union {x_1}`.  This proves that its strict-safe components
and every functional of those components descend to `G`.  Conversely, after
one label is used, each remaining replacement language is one progression
with modulus thirty-nine and its frozen least member from (15)--(16).
Enumerating the exact first edges from Part E and grouping by `G` proves
(20b); adjoining the sorted five-base tuple makes the grouping injective.
The companion verifier asserts the full eight-row stratum/order fibre ledger,
not only the coarser multiplicity histogram.

Equation (20d) is the product of the indicator identities
`1_(e_i e_j=s)=(1+s e_i e_j)/2`.  Since the matching edges are disjoint,
expansion gives one different Walsh character for every edge subset, proving
(20e) and the degree histograms.  The direct verifier evaluates every Walsh
coefficient on all seven canonical codes.  The pinning rule follows by
considering each two-vertex or free connected component separately.

For the numeration tournament, distinct owner residues modulo thirteen imply
distinct least ray bases.  Orienting by a scalar ranking is transitive and
has the fingerprint (20f).  Direct enumeration gives the `1,151` distinct
root Hamiltonian paths.

## Tournament analysis and challenged vertices

Use the order-three labels as vertices and draw the signed edge

```text
r -> o   iff r/o in {4,5,8,9},
sign=- on {4,5},              sign=+ on {8,9}.          (30)
```

For each canonical order-three set, record the number of edges, the number of
unordered pairs carrying zero/one/two arrows, the SCC sizes, and code size:

| `C` | edges | pairs `0/1/2` | SCC sizes | code size |
|---|---:|---|---|---:|
| `{1,5,8,12}` | 8 | `2/0/4` | `4` | 4 |
| `{1,2,5,8,12}` | 10 | `4/2/4` | `4,1` | 8 |
| `{1,2,3,4,9,10}` | 12 | `5/8/2` | `6` | 8 |
| `{1,2,3,5,8,12}` | 14 | `6/4/5` | `4,2` | 16 |
| `{1,2,5,6,8,11}` | 12 | `5/8/2` | `6` | 8 |
| `{1,2,5,8,11,12}` | 12 | `7/4/4` | `4,1,1` | 16 |
| `{1,3,4,9,10,12}` | 12 | `3/12/0` | `6` | 8 |

Every row has equally many positive and negative edges.  None is a tournament:
some pairs are absent or bidirected.  The quadratic-residue row is closest,
with one arrow on twelve pairs and three absent pairs.  Forcing the absent or
double pairs into a tournament destroys the NAE predicate.  Even the coarse
signed graph fingerprint fails to distinguish the third and fifth
all-order-three rows; their labelled signed incidence and affine codes are the
sharper carrier.

Challenge the runner-vertex assumption in stages:

```text
runners
 -> owner-sheet obligations
 -> signed provider incidences
 -> affine toothpick equations
 -> literal strict-safe components and remaining rays.  (31)
```

The matching code preserves the complete common-sheet predicate and destroys
continuous component geometry, lift height, and numerical comb order.  The
faithful object for the next computation is therefore an **orbit bundle**:

```text
(one of 68 signed sheet templates,
 every arithmetic lane in its 1,504-context fibre,
 literal strict-safe components,
 remaining labelled step-39 rays,
 last numerical speed,
 exact shortcut ancestry).                             (32)
```

This is also the Kakeya-needle guardrail.  The provider graph describes how
periodic danger needles meet the three local sheet germs; it does not decide
whether their full continuous combs cover the circle.

## Bounded depth-three geometry-batch prototype

The exact depth-two cache can be used without quotienting proof states.  The
prototype processes one retained-root block at a time.  It attaches every
logical lane—with context, labels, used mask, last speed, least future ray,
and shortcut ancestry—to its metric `(R,x_1,x_2)` parent, transposes proposed
third speeds, computes each `(R,x_1,x_2,x_3)` intersection once, and then
replays lane-specific streaming caps.  Complete-tooth detection is metric;
the least-future-speed test is not.

On contexts `20,...,27`, all sharing retained root mask `2199`, raw and
batched evaluators agree row-by-row and in two independent lane fingerprints:

```text
logical depths 1/2/3       568 / 40,680 / 1,685,358
depth-three verdicts       1,542,729 full-tooth
                             138,945 streaming-dead
                               3,684 materialized/live
                                   0 covering
metric calls depth 1       568 -> 143
metric calls depth 2       40,680 -> 17,447
metric calls depth 3       142,629 -> 86,316
CPU seconds                11.8168 -> 7.38807  (1.59944x)
```

The combined raw/batch process used about seven MB maximum RSS.  Timings are
only a benchmark; logical counts, call counts, and fingerprints are exact.
On context 103, `--literal-crosscheck` additionally compares every streamed
comb intersection endpoint-for-endpoint with the original materialized-band
routine; all `97,386` depth-three lanes agree.  Warning-clean `-O0`/`-O3`
builds agree after timing fields are removed, and the sanitizer pilot is
clean.

Nine earlier raw pilots project `568.6M--1.426B` ungated depth-three logical
nodes and `2.33--12.65` uncached CPU-hours, while three light full-depth pilots
gate over 99% at depth three and die by depth five.  These are workload
projections, not bounds or a census.  The prototype validates root-block
sharding and metric transposition; by itself it proved no additional lane
loose.  The later unbatched exact terminal census in Part H supersedes its open
planning status.

## Scope guardrail and reproduction

This theorem proves the complete primitive `c=3` **common-sheet**
classification, the affine unit codes, the exact sheet orbit bank, both cached
prefix layers, and the full terminal recursion over all `1,504` arithmetic
languages.  It closes the primitive proper AP-centred Hamming-six face at
common scale three.  THM-957, THM-958, THM-960, THM-962, THM-963,
THM-969, and THM-970 separately close `c=4,...,10`.  This theorem does
**not** close `c>=11`, the finite smooth-ramified Hamming-five bank, radius seven
and higher, non-AP-centred/deep-sheet packets, or global `n=12`
sporadic-branch emptiness.

Reproduce the frozen output with

```bash
python3 \
  04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py \
  > /tmp/thm862.out
cmp /tmp/thm862.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out

python3 -O \
  04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py \
  | cmp - \
  05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out

python3 \
  04-computation/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.py \
  > /tmp/thm862-cache.out
cmp /tmp/thm862-cache.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.out

python3 -O \
  04-computation/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.py \
  | cmp - \
  05-knowledge/results/lrc13_scale_three_hamming_six_geometry_cache_walsh_codex_S16.out

clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp \
  -o /tmp/thm862-depth2
for shard in 0 1 2 3 4 5 6 7; do
  /tmp/thm862-depth2 --context-start $((188*shard)) --context-limit 188 \
    > /tmp/thm862-depth2-$shard.out
done
python3 \
  04-computation/lrc13_scale_three_hamming_six_depth_two_combine_codex_S16.py \
  /tmp/thm862-depth2-{0,1,2,3,4,5,6,7}.out \
  > /tmp/thm862-depth2.out
cmp /tmp/thm862-depth2.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.out

python3 \
  04-computation/lrc13_scale_three_hamming_six_depth_two_geometry_cache_crosscheck_codex_S16.py \
  > /tmp/thm862-depth2-crosscheck.out
cmp /tmp/thm862-depth2-crosscheck.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_depth_two_geometry_cache_crosscheck_codex_S16.out

# Complete terminal recursion; these four cuts are retained-root aligned.
/tmp/thm862-depth2 --context-start 0 --context-limit 100 --depth 6 \
  > /tmp/thm862-terminal-0.out
/tmp/thm862-depth2 --context-start 100 --context-limit 408 --depth 6 \
  > /tmp/thm862-terminal-1.out
/tmp/thm862-depth2 --context-start 508 --context-limit 496 --depth 6 \
  > /tmp/thm862-terminal-2.out
/tmp/thm862-depth2 --context-start 1004 --context-limit 500 --depth 6 \
  > /tmp/thm862-terminal-3.out
python3 \
  04-computation/lrc13_scale_three_hamming_six_terminal_combine_codex_S57.py \
  /tmp/thm862-terminal-{0,1,2,3}.out \
  > /tmp/thm862-terminal.out
cmp /tmp/thm862-terminal.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_terminal_census_codex_S57.out

python3 \
  04-computation/lrc13_scale_three_terminal_packet_certificate_codex_S57.py \
  > /tmp/thm862-terminal-packet.out
cmp /tmp/thm862-terminal-packet.out \
  05-knowledge/results/lrc13_scale_three_terminal_packet_certificate_codex_S57.out
```

Frozen SHA-256 values are

```text
source  f87b0436c509313d5d90f5f9abe2d06b591e9ab639b5988bf206e4115cd7f39b
output  c8413de89655592b5009aad83596330750d9b5ca9cb407af692fd06f5e353ba8

cache source  ee7813878c7e589a30c32657c93c4d5ce106470ceaa6c125b5f09f7290797ec7
cache output  69632620030da8ec7f2d2584f37b9f450e1a3da28b6186ec48971f9b2032d72c

shared C++ source   b1b01dfdef4fba032d8641b65dad0747e43eed536821d4ab62c2137601008a9e
depth-two combiner  a946796336fb798f2e02e841a9839a6bd0ab2a33077ab352dc58ce2b98473edd
depth-two output    b3f2dc605e5a7f24e1ea796dd0721388657fac6cb0d2f59038b367b5c3e13b38

terminal combiner   a62579afde20a35e1c1ab3fafc8ff2306ca866a26154d865d85d4dce37d4c695
terminal output     ffc41d5eb38ae3c7338ea75dd0498ae7b7dc3d3e443584d34a8a59f5d554022e
packet referee      a7e7d1d4a07eea6f5abea76f16649392b9d653d4a87602071cf8f24f272bcdbd
packet output       bd20c2bc414e0dd9262c89dc04bdb3a06cf8fb246cf8b214d1c23af5119dc5f1

independent source  8cee05a32b863c369fff2f0f09fcb0f247648778fac27c6ea8eb3500cf272fa6
independent output  935cf7b809e178378230329ccdf1ffe027b20f7243105865fe9bc174aad0a58b

depth-three batch source  5b53b1cad43ac30a7c2b12c7a5fdb97513c449c10f832f7d38a190118823bc4f
depth-three batch output  016b9260477d6c1189edadff3980cda991eb384b07d5f2fe0a05fc39ed4f604c
depth-three planning note e9a5616d5efee35a4a1cf4f201bde56ad876cea43c59776e3cc9841f703506a7
```
