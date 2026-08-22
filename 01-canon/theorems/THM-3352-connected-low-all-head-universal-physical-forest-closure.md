---
id: THM-3352
title: "Connected-low all-head universal physical-forest closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  Every connected-low primitive
  six-level ray closes at every common integer dilation on all 649
  upper-median bodies and all 720 labellings.  The 261,254 finite heads left
  by THM-3350 admit universal exact physical high-forest certificates, with
  zero failures.  This is not the whole reflected branch or LRC(14).
source: root/lrc-math-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3350-connected-low-full-tree-atlas-dense-closure-and-uniform-tail
engine: 04-computation/lrc_general_reflected_pair_mass_thm3352.py
compiler: 04-computation/lrc14_connected_low_all_heads_universal_forest_thm3352.py
output: 05-knowledge/results/lrc14_connected_low_all_heads_universal_forest_thm3352.out
argmin_audit: 04-computation/lrc14_all_head_channel_argmin_reference_audit_thm3352.py
argmin_audit_output: 05-knowledge/results/lrc14_all_head_channel_argmin_reference_audit_thm3352.out
reference_engine: 04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py
reference_literal_audit: 04-computation/lrc_general_reflected_pair_mass_reference_literal_audit_thm3352.py
reference_literal_output: 05-knowledge/results/lrc_general_reflected_pair_mass_reference_literal_audit_thm3352.out
tail_script_sha256: 78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9
engine_sha256: afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b
compiler_sha256: aebfe98ab72f7eb0dc1718dfb17529e5b3f9288c6ed97d57f048bf3b29281f19
output_sha256: 048251e4aa79c3005eba677b1123514d64f474fefd00fe0d9c70d1da2692a961
argmin_audit_sha256: d0999f7b9d2a67454d9a3774bfac56ed535ea793f548feb603f1a8aa3416f0e0
argmin_audit_output_sha256: 0286c34ff0976c97d16510af6e9ffab35efdaed2aee8c7611d19c1327ecabd34
reference_engine_sha256: da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e
reference_literal_audit_sha256: a5091bb697e03bdf59fd6f31f61bea99d02a2d628e64d66851b0a8b1a8f02ba4
reference_literal_output_sha256: 1f9ec86a82695b21da2aabb27514460fe87a0cae7a2fea02c47a17b3a5f79f87
hash_basis: LF-normalized bytes
---

# THM-3352 -- connected-low all-head universal physical-forest closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

## 1. Statement

Use the notation and hypotheses of THM-3350: six distinct positive levels
whose reduced-ratio low graph `P+Q<=7` is connected, common integer dilation
`s>=1`, any of the `649` upper-median bodies, any labelling, and the body's
fixed upper-median safe cell.  Then the complete-graph Hunter certificate
closes every such ray at every scale.

THM-3350 already closed both dense one-high-edge shapes at every scale and
all other connected-low shapes from their exact thresholds, uniformly by
scale `11`.  Its finite complement was

```text
261,254 unlabelled (shape,scale) heads.                 (1)
```

Every head in `(1)` closes here.  Consequently all `305,909` unlabelled
primitive shapes, equivalently all `220,254,480` labelled primitive rays,
close under every common dilation on every one of the `649` bodies.

## 2. Exact arbitrary-ratio physical overlap

For a body ruler `L`, safe cell `j`, endpoint labels `e,f`, and levels `p,q`,
put

```text
z=Lp-e,   w=Lq-f,   r=ej mod L,   t=fj mod L.           (2)
```

The engine computes exactly the intersection mass in `[0,1]` of the two
reflected tooth unions.  It first swaps the clauses if necessary so that
`z<=w`; ordering merely by `p<=q` is false when `p=q` because the labels in
`(2)` can reverse the tooth counts.

For the nominal first-tooth indices `k=0,...,p-1`, the relative second-tooth
phase is the affine residue progression

```text
rho_k=(ak+b) mod z,
a=w mod z,
b=((rw-tz)/L) mod z.                                   (3)
```

The overlap of two unbounded full teeth is the periodized trapezoid

```text
([A-L|x|]_+-[B-L|x|]_+)/(zw),
A=(L/14)(z+w),   B=(L/14)(w-z).                         (4)
```

Write `R=floor((A-1)/L)=uz+v`, `0<=v<z`.  In the
periodized triangular sum, the `2u` complete lifts have cancelling residue
slopes and contribute

```text
2uA-Lzu^2.                                              (5)
```

Only the low prefix `rho<=v` and high suffix `rho>=z-v` remain.  Their counts
and residue sums are given by the Euclidean `floor_moments` recursion.  Apply
the same identity at peak `B` and subtract.  This proves the bulk formula at
arbitrary ratios, including any number of transverse lifts.

On the atlas domain

```text
L>=168,   14|L,   1<=e,f<=14,   p,q>=1,                (6)
```

only first-tooth indices `{-1,0,p-1,p}` can be extra, absent, or clipped.
Indeed an actual index satisfies

```text
-r/L-1/14 < k < p-(r+e)/L+1/14,                        (7)
```

and `r+e+L/14<2L`.  Every nominal tooth with `1<=k<=p-2` is therefore full;
only `k=0,p-1` can clip, while `-1,p` are the only possible extra teeth.  The
engine removes each affected nominal unbounded tooth and restores its exact
clipped intersection using the periodic primitive

```text
G(y)=floor(y)/7+min({y},1/14)+max(0,{y}-13/14),         (8)
```

whose derivative almost everywhere is the period-one indicator
`1_{||y||<1/14}`.  After the affine substitution `y=(wx-t)/L`, `(8)` gives
the exact intersection of a clipped first tooth with the entire second
union.  Thus the engine is an exact rational evaluator, not a limiting or
floating-point approximation.

The floor-moment routine is a Euclidean descent: after normalization its
modulus/slope pair is replaced by the usual reciprocal remainder pair, while
the sample count is reduced by the associated floor quotient.  Thus it has
the standard logarithmic recursion depth of a Euclidean floor sum.  This
complexity fact is useful computationally but is not used in the closure
argument.

## 3. Universal physical channel minima

Across the `649` bodies, their upper-median cells, and ordered endpoint
labels, duplicate removal leaves exactly

```text
4,044 feasible contexts (L,j,e,f).                      (9)
```

Their sorted exact digest is
`9f7c7f24f81c409c09becfa921aa53387e5093710657bd3e5a0e935ecf4ea6c2`.

For each channel-scale type `c=(s,d,P,Q)` occurring on a high edge of a head,
define the exact physical floor

```text
m_c=min_(the 4,044 contexts) I_(L,j,e,f;sdP,sdQ).       (10)
```

There are exactly `4,148` types in `(10)`.  This independent minimization is
lawful because it is used only as a one-sided lower bound.  In any concrete
body and labelling, every high edge of type `c` has physical weight at least
`m_c`, even though different minima may occur in incompatible contexts.

For each head, assign every high edge its nonnegative floor `(10)` and run
exact Kruskal, retaining a maximum spanning forest of the high graph.  Freeze
that forest.  Its physical credit in every actual context is at least the sum
of its floor weights.  Since any
forest in `K_6` extends to a spanning tree and all added intersections are
nonnegative, this is a valid Hunter lower certificate.  Comparing it with
the universal singleton-debt envelope `D_max/s` from THM-3350 proves closure.

This explains why incompatible minimizing contexts cause no problem: the
argument never claims that their minima occur simultaneously.  It compares
the actual weights of one fixed acyclic set edge by edge with universal lower
bounds.

## 4. Exact census

The head scales are

```text
s       1       2      3     4    5   6   7  8  9  10
count 204627   45399   8768  1564 518 343  14  7  7   7. (11)
```

The selected universal-forest sizes are

```text
edges       2    3     4       5
heads      48   576   6,110   254,520.                  (12)
```

All `261,254` exact comparisons are strict; failures are zero.  The weakest
margin is

```text
1017482276391320181594/73150827741949428345875 > 0,     (13)
```

on shape `(2,3,4,6,12,24)` at scale one.  Its universal forest uses the two
edges

```text
(1,5), channel 3*(1,8), floor 3276/183455,
(0,5), channel 2*(1,12), floor 2016/169343,              (14)
```

with total credit `924612948/31066820065`.

The exact semantic digests are

```text
contexts         9f7c7f24f81c409c09becfa921aa53387e5093710657bd3e5a0e935ecf4ea6c2,
channel minima   40472e196369baf0db300d96c028c2b7f78836bf8a06373b2170d4769ee9eff4,
head outcomes    c84b550ace3ce3f0433227ad3587c92d674037a1b5f37531f27075742e038a9d. (15)
```

The compiler freezes the engine hash, all atlas counts, `(11)--(13)`, and
all three digests in `(15)`.  Ordinary and `python3 -O` replays are required to
produce byte-identical output.

A separate companion recomputes all `4,148*4,044=16,774,512` fast physical
context masses and then evaluates every minimizing context with the slower
reference engine.  All `4,148` argmins agree and reproduce the channel digest
in `(15)`.  The reference engine is in turn checked against a literal rational
interval merger on `1,044,591` hostile, exhaustive equal-level, and seeded
arbitrary-ratio cases, with zero mismatches.  These independent layers audit
every numerical lower bound used by the closure, not merely the weakest row.

## 5. Scope

This closes the entire connected-low common-dilation atlas.  It does not say
that every reflected assignment has a connected intrinsic low graph; it does
not handle independently scaled disconnected low components, and it does not
prove the full reflected branch or LRC(14).  The next structural target is
therefore the disconnected-low component geometry, not another connected-low
finite head.
