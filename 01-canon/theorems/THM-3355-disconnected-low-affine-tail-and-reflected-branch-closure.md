---
id: THM-3355
title: "Disconnected-low weighted horn-tree and reflected-branch closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  Every disconnected intrinsic-low
  six-level assignment closes on every one of the 649 upper-median bodies.
  A centered primitive-grid lemma gives the physical edge floor 1/294 except
  on four oriented scale-one horns; an exact finite head, coefficient-positive
  tail, complete-multipartite forest census, and debt-sensitive K1,5 repair
  close those horns without asserting a false uniform floor.  Together with
  THM-3352 this closes the current 561-body reflected six-distinct-level k=1
  branch.  It does not close arbitrary k<=1, projected k=2,3, the rung,
  physical entry, or LRC(14).
source: codex-major-frontiers-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1250-six-private-needles-force-fully-located-spanning-tree
  - THM-3350-connected-low-full-tree-atlas-dense-closure-and-uniform-tail
  - THM-3352-connected-low-all-head-universal-physical-forest-closure
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
script: 04-computation/lrc14_disconnected_weighted_horn_tree_closure_thm3355.py
output: 05-knowledge/results/lrc14_disconnected_weighted_horn_tree_closure_thm3355.out
script_sha256: 78975824a16471cd6aad7194443f292773a9af08c748d945733768d64c2116dd
output_sha256: 1f18288cd8c2886f85b0a1839e1cb0737d9141a8de0e543e68a28ef6d7e9c565
hash_basis: LF-normalized bytes
---

# THM-3355 -- disconnected-low weighted horn trees close the reflected branch

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

## 1. Statement

Use the reflected physical setup of THM-2941 and the notation of THM-3350.
There are six distinct positive levels `q_0,...,q_5`.  For each pair write

```text
(q_i,q_j)=g(P,Q),             gcd(P,Q)=1, P<Q.             (1)
```

Join `i,j` in the intrinsic low graph when `P+Q<=7`, and suppose this graph
is disconnected.  Fix any of the `649` upper-median bodies, its ruler `L`,
its selected whole body-safe cell, and any assignment of its six endpoint
labels to the levels.  Then a complete-graph Hunter tree has physical pair
credit strictly exceeding the exact singleton debt.  Hence every such
disconnected-low assignment closes.

THM-3352 proves the same conclusion when the intrinsic low graph is
connected.  The two cases are exhaustive, so their composition closes every
six-distinct-level reflected assignment in the current upper-median bank and,
through THM-2941's reflected reduction, all `561` residual bodies in this
declared `k=1` branch.

The proof does **not** establish the formerly sought uniform edge floor on
four oriented horns.  It keeps their location and pays for them with the
remaining forest and a smaller level-sensitive debt.

## 2. Inherited component and debt geometry

THM-3350 enumerates the connected intrinsic-low shapes on `r=1,...,5`
vertices:

```text
1, 8, 94, 1295, 19389.                                  (2)
```

Consequently the unordered disconnected six-vertex component profiles have
generating coefficient

```text
[x^6] product_(r=1)^5 (1-x^r)^(-c_r)=36520.             (3)
```

Every profile is realizable: after placing one primitive connected component,
scale the next by more than six times the current maximum.  Internal low
ratios are unchanged and every new cross ratio is high.  Thus this is not a
formal-partition overcount.

For any profile the cross-component graph is complete multipartite and hence
connected.  Every forest in it extends to a spanning tree of `K_6`; all added
physical overlaps are nonnegative.  It is therefore enough to find a cross
forest whose certified edge weights exceed the debt.

For ruler `L`, endpoint labels `e_i`, and distinct levels `q_i`, the exact
singleton debt is

```text
D(e,q)=sum_i e_i/[7(Lq_i-e_i)].                           (4)
```

The rank monotonicity and exact `649*6!` census in THM-3350 give

```text
D(e,q)<=Dmax=186636088362/11773143757375.                 (5)
```

The target edge weight is

```text
eta=1/294,
5 eta-Dmax=570672686921/494472037809750>0.                (6)
```

Thus five `eta`-edges close without any compatibility assumption on their
individual minimizing contexts.

## 3. The centered primitive-grid lemma

Fix a body-safe cell, ordered endpoint labels `e,f`, ruler `L`, and the pair
in `(1)`. Here `e` labels the lower level `gP`, while `f` labels the higher
level `gQ`; this orientation is part of the channel data. Put

```text
p=gP, q=gQ,       z=gLP-e, w=gLQ-f,       C=Qe-Pf.       (7)
```

The one-first-tooth response `f_0` from the full-turn decomposition has

```text
mean(f_0)=L/(49z),       osc(f_0)<=L/(7w),
Lip(f_0)<=L/w.                                             (8)
```

Its superlevel sets are unions of at most two circular intervals: after
removing the constant full-turn baseline, `f_0` is the circular convolution
of two interval indicators, hence a periodized trapezoid.

The physical overlap is the sum of `gP` samples of `f_0` with circular step
`w/z`.  Split them into `g` consecutive blocks of length `P`.  In one block
replace the true step by `Q/P`.  Since `gcd(P,Q)=1`, the replacement points
are one complete translated `P`-grid.  Layer-cake discrepancy for a union of
at most two circular intervals gives

```text
sum_(one complete P-grid) f_0
 >=P mean(f_0)-2 osc(f_0).                                (9)
```

Indeed each superlevel count differs from `P` times its length by at most two;
integrating the count inequality proves `(9)`.

The two steps differ exactly by

```text
w/z-Q/P=C/(Pz).                                           (10)
```

Center the comparison at a median grid index.  The total index displacement
in one block is

```text
min_t sum_(r=0)^(P-1)|r-t|=floor(P^2/4).                 (11)
```

Using `(8)`--`(11)` in all `g` blocks proves the uniform located bound

```text
I_(L,j,e,f)(gP,gQ)
 >= gLP/(49z)
    -2gL/(7w)
    -gL |Qe-Pf| floor(P^2/4)/(P wz).                     (12)
```

No phase is averaged away: the starting translate of each complete grid is
arbitrary.  Primitivity is load-bearing.  At

```text
(L,j,e,p,f,q)=(1680,870,5,1792,10,2688)                 (13)
```

the exact physical overlap is zero, while illegally inserting the
nonprimitive pair `(896,1344)` into `(12)` gives the positive value
`5457530976/271903091713`.

## 4. Endpoint reduction and the regular edge floor

For fixed `L,e,f,P,g`, the right side of `(12)` is fractional-linear in `Q`
on each side of the single kink `Q=Pf/e`.  Hence its minimum on the relaxed
interval `P+1<=Q<=8P` occurs at an endpoint or at the kink.  At the kink the
determinant loss vanishes and the remaining loss is no worse than at `P+1`.
It follows that every primitive `P<Q<8P` is bounded below by

```text
min(B(P,P+1),B(P,8P)),                                   (14)
```

where `B` denotes `(12)`.

The value `B(P,8P)` is only the right endpoint of the relaxed algebraic
envelope. The pair `(P,8P)` is generally nonprimitive, and the primitive-grid
lemma is not being applied to it as a physical channel.

For odd `P=421+2t` and even `P=422+2t`, clear the positive denominator in
`294(B-eta)`.  On every small-ruler ordered lane except four listed below at
`g=1`, and on every lane at `g=2,3`, both endpoint numerators have strictly
positive coefficients as polynomials in `t`. The companion also verifies
coefficientwise that `N'D-ND'>=0` for every cleared margin `N/D`, so the
minimum occurs at the base point of its parity. The exact bank contains

```text
g=1 nonhorn records   5200,
g=2,3 all records    10432,          total 15632.         (15)
```

The weakest tail margin is

```text
1687815/5745649440454>0                                  (16)
```

at `(g,L,e,f,P,Q)=(1,168,1,12,421,422)`.  Thus `(14)` proves
`I>=eta` for every declared regular lane from `P=421` onward.

For `P<=420`, `g=1,2,3`, all `2530` small-ruler contexts, all primitive
`P<Q<8P` with `P+Q>=8`, and `(P,Q)!=(3,5)`, an endpoint screen closes
`3,065,046` context/range rows.  The exact THM-3352 mass engine evaluates the
remaining

```text
6,144,244 physical masses,       failures 0.             (17)
```

The minimum is

```text
92/7645=1/294+19403/2247630                              (18)
```

at `(g,P,Q;L,j,e,f)=(1,4,5;168,90,12,6)`.  Eight fixed
context shards have separately frozen counts, minima, and semantic hashes.

The inherited exact companions complete the other regimes:

- every `L>=4592` lane has mass at least `eta` at every pair dilation;
- `q>=8p` has mass at least `23/4655>eta`;
- the primitive `3:5` lane has mass at least `eta` at every dilation; and
- every other moderate-ratio lane with `g>=4` has mass at least `eta`.

Consequently every high cross edge is **regular**, with weight at least
`eta`, except possibly

```text
L=168, j=90, e=12, f in {1,2,3,4},
g=1, P>=421, P<Q<8P.                                    (19)
```

These four oriented horns are honest.  The asymptotic values of the bound
`(12)` minus `eta` are respectively

```text
-25/37632, -3/6272, -11/37632, -1/9408.                 (20)
```

So the proof must not relabel them as regular by continuity or by a sampled
finite floor.

## 5. Deleting the horns leaves almost all of the tree

Identify the possible horn vertices with labels `1,2,3,4,12`; label `6` is
the sixth vertex forced by `L=168`.  For each nontrivial set partition of six
vertices, form its complete multipartite cross graph and delete the four
undirected pairs `{12,f}`, `f=1,2,3,4`.  Across all `B_6-1=202` partitions,
the resulting component census is

```text
components       1    2   5
partitions     150   51   1,                              (21)
```

so the maximum regular-forest sizes are `5,4,1`.  The unique five-component
case is

```text
{1,2,3,4,6} | {12}.                                     (22)
```

This census also has a direct graph proof.  All deleted edges meet `12`.
Unless every other vertex lies in one low component as in `(22)`, the
remaining complete multipartite graph away from `12` is connected, and at
most one additional component can occur.

If the graph in `(21)` is connected, five regular edges and `(6)` close.  If
it has two components, take four regular forest edges.  Should one deleted
cross edge actually have the reverse orientation or lie outside `(19)`, it
is regular and joins the components, reducing to the five-edge case.
Otherwise a genuine horn is present, so the level on label `12` is at least
`421`.  Rank monotonicity then gives the sharper exact debt bound

```text
D<=D421=443767487288/52278303328335,
3/294-D421=1255584224873/731896246596690>0.              (23)
```

Thus even three of the available four regular edges suffice.  This closes
all nonexceptional partitions.

## 6. The exceptional `K_(1,5)` star

It remains to treat `(22)`.  Its cross graph is the star from `12` to
`1,2,3,4,6`.  The edge `12--6` is always regular.  Let

```text
b=#{f in {1,2,3,4}: q_f<q_12}.                           (24)
```

Those `b` edges have the reverse horn orientation and are regular, so the
star contains at least `1+b` regular edges.  If `b=4`, all five edges are
regular and `(6)` applies.

Assume `0<=b<=3` and `q_12>=421`; the finite head already handles smaller
`q_12`.  The low component on `{1,2,3,4,6}` is connected, and every low ratio
is at most six. Any vertex below `q_12` crosses at most `b+1` low edges before
first reaching a level above `q_12`. Therefore each such lower level is strictly greater
than

```text
q_12/6^(b+1).                                             (25)
```

If label `6` is actually above `q_12`, moving it into the relaxed lower bank
only increases the debt bound; `(25)` remains conservative.  Debt decreases
with every level, so the worst row occurs at `q_12=421`, with the lower levels
the first distinct integers above `421/6^(b+1)` and all upper levels beginning
at `422`.  Exact permutation over the endpoint labels gives:

```text
b  regular edges  first lower level  maximum debt D_b
0       1                71          10171035532358753244424/87499329204988335395285245
1       2                12          65302219886882882438/90087463329358292024115
2       3                 2          140545706290782894/31861415435911397875
3       4                 1          7350964239952/883497095223435.          (26)
```

In every row the available regular credit is strictly larger:

```text
(1+b)/294-D_b =
12072720679782123134489227/3674971826609510086601980290,
3832763666951738490749/630612243305508044168805,
2583990888487810609/446059816102759570250,
32685830817806/6184479666564045,                         (27)
```

respectively.  Extend these regular star forests to complete spanning trees.
Nonnegative added edges cannot reduce the Hunter credit, so `(27)` closes the
exceptional partition.

## 7. Assembly and the affine near miss

Sections 4--6 prove every disconnected-low assignment.  THM-3352 proves every
connected-low assignment.  This closes the reflected six-distinct-level
physical branch, including independently scaled disconnected components; it
does not assume a common scale across those components.

The earlier affine-ray route is not a dependency.  Its primitive quotient
`22890 -> 14168` and strict carrier cutoff `14913` survive as an independent
sidecar, but MISTAKE-377 repairs an omitted `9|c|<=p` hypothesis in its
many-turn routing.  The corrected residual occurrence count is `8,079,264`,
not `8,013,156`.  The weighted horn-tree proof avoids that scan entirely.

## 8. Verification and scope

The companion pins its inherited structural, large-ruler, `3:5`, `g>=4`,
mass-engine, and output dependencies by LF-normalized SHA-256.
It contains no Python `assert` nodes.  Ordinary and optimized eight-worker
replays are byte-identical to the stored transcript.

Besides `(15)--(18)` and `(21)--(27)`, it performs:

- `100,000` seeded exact checks of `(12)` against physical mass;
- nine fast/reference/literal three-route positive controls;
- a physical low-edge zero control;
- the nonprimitive hostile `(13)`;
- all `467,280` body/rank debt rows; and
- all multipartite spanning-tree counts, whose exact min/max/sum are
  `1/1296/67392`.

This theorem closes the current `561` reflected bodies only in the declared
six-distinct-level `k=1` lane.  Projected `k=2,3`, arbitrary `k<=1`, the
six-body/seven-tail rung, semantic physical entry from projected states, and
LRC(14) remain **OPEN**.
