---
id: THM-3350
title: "Connected-low full-tree atlas, dense closure, and uniform high-forest tail"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  There are exactly 305,909
  unlabelled primitive six-level shapes with connected intrinsic low graph.
  Complete-graph Hunter trees have homogeneous credit exceeding singleton
  debt on every labelling and every one of the 649 upper-median bodies.  The
  two one-high-edge shapes close at every common dilation; every other shape
  closes uniformly from common scale 11 onward.  The remaining finite bank
  has 261,254 unlabelled shape-scale heads before body/label compression.
  This is not the whole reflected branch or LRC(14).
source: root/lrc-math-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3349-reflected-low-two-star-selected-half-limit-all-dilations
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
atlas_script: 04-computation/lrc14_connected_low_full_tree_atlas_thm3350.py
atlas_output: 05-knowledge/results/lrc14_connected_low_full_tree_atlas_thm3350.out
atlas_script_sha256: 8020b953084b19940845349980a83547a5d73bbd2e96f3cc45556402b015bb67
atlas_output_sha256: 75aac002408daa90c63f6706b909fab6aab0fb50f4d728626282d6ffef34607b
dense_script: 04-computation/lrc14_connected_low_dense_full_tree_all_dilations_thm3350.py
dense_output: 05-knowledge/results/lrc14_connected_low_dense_full_tree_all_dilations_thm3350.out
dense_script_sha256: a44755f9cbc4f9ae91cd465ea6773744eb964988f6adf9733d28ec19e6d643b3
dense_output_sha256: 22803a615cc75a836c359d6a6db3f19beef79642e344042b9bcdbe4bed893cd9
tail_script: 04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py
tail_output: 05-knowledge/results/lrc14_connected_low_uniform_high_forest_tail_thm3350.out
tail_script_sha256: 78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9
tail_output_sha256: 2d409089c5ceb6e9ab7d3ae8aa99bc21a711fc98944e72f1da5daa31874270fa
hash_basis: LF-normalized bytes
---

# THM-3350 -- connected-low full trees, dense closure, and a uniform tail

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

## 1. Statement

Let `q=(q_0,...,q_5)` have six distinct positive integer coordinates.  For
each pair write

```text
(q_i,q_j)=d(P,Q),                 gcd(P,Q)=1.             (1)
```

Call `ij` intrinsically low when `P+Q<=7`, and suppose that this low graph is
connected.  Quotient by common scaling, but retain all `6!` labellings.  On
each of the `649` upper-median bodies behind THM-2941, use its deterministic
upper-median body-safe cell and the complete-graph Hunter functional.

Then:

1. there are exactly `305,909` unlabelled primitive shapes, hence
   `220,254,480` labelled primitive rays;
2. on every body and labelling, the maximum homogeneous complete-graph tree
   credit is strictly greater than the singleton debt;
3. the two exceptional shapes

   ```text
   (1,2,3,4,6,12),        (2,3,4,6,8,12)                 (2)
   ```

   close at every common integer dilation `s>=1` on every body and labelling;
4. every other connected-low shape closes, on every body and labelling, for
   every common dilation `s>=11`.

The exact shape-dependent tail is sharper than the uniform value `11`.  Its
finite complement contains

```text
261,254 unlabelled (shape,scale) heads,                  (3)
```

before exploiting body, labelling, phase, or channel coincidences.  Those
heads remain open here.  Thus this theorem is a homogeneous closure plus a
dense all-scale theorem plus a uniform non-dense tail, not a proof of the
whole reflected branch or LRC(14).

## 2. The exact intrinsic atlas

The unoriented low-ratio alphabet is

```text
{4/3,3/2,2,5/2,3,4,5,6}.                               (4)
```

Normalize a connected projective set so that its least element is `1`.
Starting from `{1}`, adjoin a new vertex along any oriented ratio in `(4)`,
then clear denominators and divide by the global gcd at size six.  This is
complete: root a low spanning tree at the least vertex and insert vertices
in parent-before-child order.

The companion performs two independent exact enumerations:

- ordinary connected growth with global set deduplication;
- duplicate-free reverse search, whose canonical parent deletes the greatest
  removable non-root vertex.

Their layer counts agree exactly:

```text
1, 8, 94, 1,295, 19,389, 305,909.                       (5)
```

The primitive-row digest and the compact labelled-orbit digest are

```text
6d0b14b4608ba696bb107689afecf985a8ad70f5e5e3ada8698b70e5bede8681,
357926e866c6fd1b2705175a2caf57b30f44e04729e157f16ed66c54238b5a43. (6)
```

The exact low-edge histogram is

```text
edges  5       6      7      8     9    10  11  12  13  14
count  171119  88196  32244  10046 3131 894 210 58  9   2. (7)
```

The only fourteen-low-edge rows are `(2)`: their respective unique high
edges are `1--12` and `3--8`.  All other `305,907` rows have at least two
high edges.  The unique maximum primitive height is `7,776`, on
`(1,6,36,216,1296,7776)`.

## 3. Why the missing object is the complete tree

Any two distinct simple edges of `K_6` form a forest, since a simple cycle
has at least three edges.  Every forest in `K_6` extends to a spanning tree.
This observation needs neither connectivity of the high graph nor adjacency
of the selected high edges.

For a high edge, `P+Q>=8`.  THM-2941 `(25i4)` gives the pointwise projective
fibre floor

```text
Phi_(P,Q)(z)>=1/105.                                    (8)
```

Its homogeneous located overlap is an average of this fibre, so it is at
least `1/105`.  All other edge overlaps are nonnegative.  Thus every
non-dense shape has a complete-graph tree of homogeneous credit at least
`2/105`.

For body labels `e_i`, ruler `L`, and distinct levels `q_i`, the scale-one
singleton debt is

```text
D(e,q)=sum_i e_i/[7(q_i L-e_i)].                        (9)
```

If `r_i` is the rank of `q_i`, then `q_i>=r_i` and every summand decreases
with its level.  Hence `D(e,q)<=D(e,r)`.  The exact `649*6!=467,280` rank
census has maximum

```text
D_max=186636088362/11773143757375                       (10)
```

at body `(1,2,3,4,6,12)` and rank assignment `(6,5,4,3,2,1)`.  Exact
subtraction gives

```text
2/105-D_max=112842806764/35319431272125>0.              (11)
```

This proves the homogeneous gate for all non-dense rows.  The dense companion
checks the same strict gate for every labelling of `(2)`.  Therefore the
maximum homogeneous complete-tree credit beats debt on all `220,254,480`
labelled shapes and all `649` bodies.

Connectedness of the low graph is used to define and enumerate this finite
atlas.  It is not being smuggled into the forest-extension argument.

## 4. All-channel midpoint/shear transport

The low-channel restriction in THM-3349 is unnecessary for the transport
estimate.  Fix one edge, ruler `L`, safe cell `j`, endpoint labels `e,f`, and
raw levels `(1)`.  Put `h=sd`, `R=ej mod L`, `S=fj mod L`, and

```text
F(x,t)=chi(Px-(R+et)/L) chi(Qx-(S+ft)/L),
I_s=integral_0^1 F(ht,t)dt,
H(t)=integral_0^1 F(x,t)dx,                             (12)
```

where `chi` is the circular radius-`1/14` indicator.  Slabwise contraction
toward each midpoint gives exactly as in THM-3349

```text
|I_s-J_h| <= gamma_P eP/(hLP-e)+gamma_Q fQ/(hLQ-f),     (13)

gamma_k = 1/2                    if k is even,
          (k^2+1)/(2k^2)         if k is odd,            (14)
```

where `J_h` is the composite midpoint sum of `H`.

For `C=Qe-Pf`, the static fibre has the periodized-tent form

```text
H(t)=Phi(z_0+Ct/L),
Phi(z)=[T_A(z)-T_B(z)]/(PQ),
A=(P+Q)/14,                 B=|P-Q|/14.                 (15)
```

Distributionally on the circle,

```text
Phi''=(delta_(-A)+delta_A-delta_(-B)-delta_B)/(PQ),      (16)
```

with the two `B=0` atoms combined.  Consequently `TV(Phi')<=4/(PQ)`
per traversed period.  A path of length `|C|/L` crosses each phase atom at
most `floor(|C|/L)+1` times, whence

```text
TV(H') <= 4|C|(floor(|C|/L)+1)/(LPQ).                   (17)
```

Composite midpoint quadrature for a continuous piecewise-affine function is
bounded by `TV(H')/(8h^2)`: decompose into an affine function and slope-jump
hinges, each of which costs at most its jump divided by `8h^2`.  Combining
`(13)` and `(17)` yields the general, multiwrap-safe bound

```text
|I_s-I_infinity| <= B_e(s),                             (18)

B_e(s)=gamma_P eP/(sdLP-e)
      +gamma_Q fQ/(sdLQ-f)
      +|C|(floor(|C|/L)+1)/(2(sd)^2LPQ).                (19)
```

This is symmetric under exchanging endpoints.  When `C=0`, the last term
vanishes and the static fibre is constant; the exact companions evaluate the
limit directly rather than invoking a singular Bernoulli quotient.

## 5. The two dense shapes close at every dilation

For every body and every labelling of either shape in `(2)`, the dense
companion computes all fifteen exact homogeneous limits and selects a maximum
spanning tree by Kruskal, breaking exact ties by the lexicographically largest
edge.  It uses `(19)` on the five selected edges to find a tail threshold and
literal rational arc intersection for every earlier scale.

The complete census is

```text
assignments                 2*649*6! = 934,560
tail threshold histogram    1:934,335, 2:223, 3:2
literal finite heads        227
failures                    0.                            (20)
```

The weakest exact head margin is

```text
1359521303047/13562237063235>0.                         (21)
```

It occurs on body `(1,2,3,4,6,12)`, levels
`(4,3,2,1,12,6)`, scale one, and tree
`((1,4),(0,3),(3,4),(1,2),(0,5))`.  Thus both dense shapes close at every
common dilation.  Here the tail propagation uses a scaled inequality, not
merely ordinary monotonicity.  If `s>=S`, each shear summand in `(19)` and
the determinant summand satisfy

```text
B_e(s) <= (S/s) B_e(S),                                 (21a)
```

the determinant term in fact decaying quadratically.  Likewise the exact
singleton debt satisfies `D(s)<=(S/s)D(S)`.  Thus, for the fixed tree chosen
at its certified threshold `S`,

```text
I_infinity(T)-sum_(e in T) B_e(s)
 >= (S/s)(I_infinity(T)-sum_(e in T) B_e(S)),            (21b)
```

because the omitted `(1-S/s)I_infinity(T)` is nonnegative.  Strict closure
at `S` therefore persists at every `s>=S`; the literal checks cover exactly
the earlier scales.

## 6. A context-aware high-forest tail by scale eleven

For a non-dense shape and scale `s`, assign each high edge the certified gain

```text
g_e(s)=1/105-E_e(s),                                    (22)
```

where `E_e(s)` is the maximum of `(19)` over the `180` ordered endpoint-label
contexts that actually occur in the `649` bodies.  For each ordered label
pair the companion uses the least ruler of a body containing it.  This is
uniformly lawful: both shear terms decrease with `L`, and for fixed `c`,

```text
c(floor(c/L)+1)/L                                       (23)
```

decreases between quotient jumps and drops at every jump.  The exact context
digest is

```text
ee75c24771753b676cbc296a283c5f4efd87053f7266e716cf3ac36d8a4cedaf. (24)
```

Kruskal selects a maximum forest of positive gains `(22)` in the high graph.
Extend it to a spanning tree of the complete graph; added edges have
nonnegative physical overlap.  At scale `s`, singleton debt is at most
`D_max/s` (strictly less for `s>1`).  Therefore

```text
sum_(e in forest) g_e(s)>D_max/s                         (25)
```

is a uniform physical Hunter certificate for every body and labelling of
that shape.  The exact first-certified-scale histogram on the `305,907`
non-dense shapes is

```text
S       1       2       3      4     5    6   7   8  11
count 101280  159228   36631   7204  1046 175 329  7   7. (26)
```

To propagate `(25)`, freeze the positive high forest chosen at its first
certified scale `S`.  The contextwise maximum errors still obey
`E_e(s)<=(S/s)E_e(S)` for `s>=S` (a maximum of quantities with the same
bound has the same bound).  Consequently each frozen edge satisfies

```text
g_e(s) >= (S/s)g_e(S)+(1-S/s)/105 >= (S/s)g_e(S).        (26a)
```

Summing over the frozen forest and using `(25)` at `S` gives a strict lower
bound `(S/s)(D_max/S)=D_max/s`.  This scaled monotonicity, rather than the
mere fact that both error and debt decrease, proves that every row closes by
its displayed threshold and that all close for `s>=11`.  Summing `S-1` over
the shapes gives the residual count `(3)`.

## 7. Hostile boundary: low-only trees are false

The complete graph is essential.  On

```text
body=(1,2,3,4,6,12), L=168, upper-median cell j=90,
q=(1,3,9,27,108,432),                                   (27)
```

the low graph is the path `01,12,23,34,45`, with channels
`1:3,1:3,1:3,1:4,1:4`.  Every homogeneous low-edge overlap and every
scale-one physical low-edge overlap is exactly zero, while

```text
debt=5824887360027324/3056780423523449161>0.             (28)
```

At cell `12` the same low tree has limit `29/224`, so cell reselection can
repair this particular row.  At the fixed upper-median cell, however, only
the complete/high tree is a uniform object: its scale-one physical margin is

```text
362436495428246378/3056780423523449161>0.                (29)
```

Thus the discarded low-only compiler fails structurally, not because its
constant was slightly weak.

## 8. Exact verification and scope

The atlas companion uses only exact integer and rational arithmetic and two
complete enumerators.  The dense companion scans all `934,560` assignments,
checks the universal `467,280`-row rank-debt bank, and compares optimized and
literal rational intersection masses on the finite heads.  The uniform-tail
companion reconstructs the atlas and all `180` body contexts, then freezes
the threshold histogram and semantic digest.  Ordinary and `python3 -O`
transcripts byte-match for all three companions; their hashes are in the
front matter and results index.

This theorem assumes six distinct levels, connected intrinsic low graph,
common dilation, one of the `649` upper-median bodies, and its selected safe
cell.  It does not prove that every reflected assignment enters this family,
does not close the `261,254` finite heads, does not handle independently
scaled low components, and does not prove LRC(14).
