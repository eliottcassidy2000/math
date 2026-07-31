---
id: THM-2892
title: "Eight-body five-slot heavy-triangle closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.
  Every eight-speed body in {1,...,14} has positive good-set measure after
  adjoining any five distinct speeds at least 15.  The heavy-triangle census
  closes all 10202 honest residual first-apex carriers and all 1444 residual
  active roots; THM-885 supplies the disjoint 495 low bodies.  With THM-741
  this closes the complete v8<=14 sector, not unrestricted LRC(14): the
  v8>=15 / at-most-seven-in-window sector remains open.
source: codex/lrc-j5-recursion-2026-07-29; independent mac-mini/root audits
depends_on:
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-741-near-AP-four-slot-closure-all-2002-bodies-in-1-14
  - THM-885-covering-case-decomposition-j56-sweep
  - THM-2885-eight-body-global-top-fifteen-and-top-ten-hitting-gate
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-2893-complement-cap-finite-core-flag-lemma
related:
  - THM-2885-eight-body-global-top-fifteen-and-top-ten-hitting-gate
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-2893-complement-cap-finite-core-flag-lemma
verification:
  - 04-computation/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.out
  - 04-computation/lrc14_j5_singleton_component_seal_audit_codex_20260729.py
  - 05-knowledge/results/lrc14_j5_singleton_component_seal_audit_codex_20260729.out
---

# THM-2892 -- eight-body five-slot heavy-triangle closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.**

The ordinary and optimized locked replays are byte-identical.  Independent
audits accepted the heavy-triangle logic, forced-edge dichotomy, equality
boundaries, root recomposition, and exact scope.

## 1. Statement

Let

```text
E in C({1,...,14},8),                    G_E=the good set of E,
Q subset {15,16,...},                    |Q|=5.             (1)
```

Then

```text
|G_(E union Q)|>0.                                      (2)
```

for every such `E,Q`.  Together with THM-885 on the `495` bodies contained
in `{1,...,12}` and THM-741 whenever a ninth speed enters `{1,...,14}`, this
closes the complete sector in which at least eight normalized speeds lie in
`{1,...,14}`.

It does **not** prove unrestricted LRC(14).  The sector

```text
v_8>=15,
```

equivalently at most seven normalized speeds in `{1,...,14}`, remains
separate.

## 2. Inherited exact reduction

THM-2885 assigns ten globally ranked first apices to every eight-body root.
THM-2888 then closes `1064/2508` active roots meeting `{13,14}` using only
genuine terminal apices.  It leaves exactly

```text
1444 active roots,
10202 honest nonterminal first-apex carriers.           (3)
```

No merely finite rank-three branch is terminal in `(3)`.

For a residual apex carrier write

```text
C=G_E minus D_a,       h=|C|,       r=#components(C),
U(x,y)=|C intersect (D_x union D_y)|.                   (4)
```

THM-2888 supplies a global pair cap

```text
U(x,y)<=B<5h/7                                           (5)
```

on every target carrier, subject to the forced-edge normalization below.

## 3. The twenty-six forced-edge normalizations

Exactly `26` of the target carriers have the raw THM-2888 cutoff above
`2500`.  Their largest raw cutoff is

```text
130827.                                                  (6)
```

For each, THM-2888 selects its exact maximizing edge `e={x,y}` and proves
that the literal residual `C minus (D_x union D_y)` cannot be covered by
two further dangers.  Thus every four-set containing both endpoints of
`e` closes.

Delete exactly `e` from the allowed pair maximum and recompute the global
cap `B'`.  Every four-set avoiding the already closed edge has all six of
its pairs governed by `B'`.  The largest deleted-edge cutoff is only

```text
1497.                                                    (7)
```

and every unnormalized target has cutoff at most

```text
2409.                                                    (8)
```

The target-normalization census pays `282` exact pair checks on the forced
residuals and `321` on the deleted maxima.  Direct full-family
reconstruction verifies every forced residual.  One deletion suffices on
all `26` targets; no deleted cap is already pair-direct.

## 4. Heavy triangles and literal residuals

For an unnormalized carrier use `B` in `(5)`; for a normalized carrier use
the deleted cap `B'` on the branch avoiding its forced edge.  In either case
write

```text
theta=h-B_* > 2h/7,
H={w : |C intersect D_w|>=theta/2}.                     (9)
```

THM-2893 with `(k,s,ell)=(4,2,3)` proves:

- any hypothetical four-cover makes every one of its six pairs
  `theta`-heavy;
- at least three of its vertices lie in the finite core `H`; and
- those three form a heavy triangle whose literal residual must be covered
  by the fourth speed.

The strict good-set discrepancy estimate supplied by THM-735 through the
canonical THM-732 tail bound seals `H` globally.  The exact census finds

```text
finite H vertices             55159       (maximum 65/carrier),
heavy edges                  140846       (maximum 1021/carrier),
heavy triangles              372209       (maximum 7814/carrier),
carriers with a triangle       4578,
H-only K4s                  1007764.                      (10)
```

The H-only `K4` count is diagnostic only: THM-2893 forces three, not
necessarily four, cover vertices into `H`.  The consequence object is the
literal residual behind every heavy triangle.

For such a triangle `L`, put

```text
R=C minus union_(w in L) D_w,          m=|R|,    s=#components(R).   (11)
```

Every residual has `m>0`.  Scan all allowed fourth speeds through

```text
W=max(14,floor((99/70)s/(6m))).                         (12)
```

For every unscanned `w>=W+1`, the strict discrepancy estimate gives

```text
|R intersect D_w|
 < m/7+(99/70)s/(7w)
 <= m/7+(99/70)s/(7(W+1))
 < m.                                                   (13)
```

The finite scan also has maximum coverage strictly below `m` on every row.
Consequently none of the `372209` triangles extends to a four-cover.
There are

```text
372209/372209 one-comb closures,
0 extendible triangles,
0 empty triangle residuals,
maximum W=476.                                          (14)
```

The smallest margin of the locked discrepancy certificate occurs at

```text
E=(1,5,6,7,8,9,11,13),       apex a=19 at rank 2,
L=(37,121,130),
m=1099999/24444420,           s=42,       W=219.         (15)
```

Its finite-head maximizer is speed `25`, with coverage

```text
7111/459800,
```

but the global cap is tail-dominated:

```text
m/7+(99/70)42/(7*220)=19249981/427777350,
m-cap=1/285184900>0.                                   (16)
```

The first unscanned speed is `220`, so the strict first inequality in `(13)`
is load-bearing for this particular discrepancy certificate.  This number
is not a near-covering margin.

THM-2893's longest-component singleton seal gives an independent geometric
audit.  If `[a,b]` is a residual component, coverage up to measure zero by
one danger comb requires containment in the closure of one tooth.  Equivalently
there must be an integer `k` with

```text
w*b-1/14 <= k <= w*a+1/14,
```

so noncontainment is exactly

```text
ceil(w*b-1/14)>floor(w*a+1/14).                         (16a)
```

In particular, if `ell=b-a`, a singleton cover requires

```text
w<=floor(1/(7*ell)).                                    (16b)
```

For the row `(15)`, the longest component is

```text
[939/1820,881/1694],       ell=911/220220,
1/(7*ell)=31460/911,       floor(1/(7*ell))=34.         (16c)
```

Only speed `27` among the allowed speeds `15,...,34` passes the necessary
closed-tooth containment test.  A full exact scan through `34` still has
maximum coverage `7111/459800` at speed `25`, hence the exact finite
coverage margin is

```text
m-7111/459800=137171599/4644439800>0.                   (16d)
```

For `w>=35`, `(16b)` excludes coverage without discrepancy.  Taking tooth
closures in `(16a)` deliberately retains endpoint equality, so this
independent seal cannot manufacture a positive-measure conclusion by
discarding boundary points.  Thus the singleton closure is not intrinsically
dependent on the strictness in `(13)`.  The exact-cover five-set hypergraph
is empty.

## 5. Carrier and root recomposition

Of the `10202` target carriers,

```text
5624 are triangle-free,
4578 have triangles and close by (13),
10202/10202 close,
0 remain.                                               (17)
```

Return these newly proved terminal apices to the THM-2888 weighted
top-fifteen gate.  Every initially closed root stays closed, and every one
of the `1444` residual roots now has positive complement margin.  Hence

```text
2508/2508 active roots close.                            (18)
```

The minimum final weighted margin is

```text
3973/89664120
```

at

```text
E=(2,5,7,8,10,12,13,14).                                (19)
```

For that root the initial terminal set is

```text
{19,22,38,44,55,63}
```

and its retained leading five speeds are

```text
(18,27,17,46,34).
```

Together with the `495` THM-885 low bodies, all `3003` eight-body roots
close.

The project decomposition strata sharpen `(18)`:

```text
exactly one of {13,14} in E:
  roots 1584/1584, target apices 6032, triangles 237601;

both 13 and 14 in E:
  roots 924/924, target apices 4170, triangles 134608.   (20)
```

## 6. Verification ledger

The locked verifier and its transcript are

```text
04-computation/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.py
SHA-256 0a386864dee44144130d25060301ed0f6d8c3cd02136b6aecf58efa3ae3a790d

05-knowledge/results/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.out
SHA-256 0b595ca0cdf6c6fca6f4e255ec81d0aadc2e804f0b7bb2237925ad5575898447
```

It hash-pins THM-2888 and its transcript, reproduces the complete
`30030`-carrier atlas, restricts to the honest `10202` post-gate targets,
replays every forced-edge normalization, and uses exact rational arithmetic
throughout.  Its controls include:

```text
20404 H vector/scalar comparisons,
4578 sequential/simultaneous triple-subtraction comparisons,
744418 triangle-residual vector/scalar comparisons,
direct full-family reconstruction of every one of 372209 residuals. (21)
```

Ordinary `python3` and optimized `python3 -O` replays are byte-for-byte
identical and both match the stored transcript.  Their wall times were
`2244.24` and `2247.15` seconds.  Every discovery count, digest, extremum,
and equality convention is hard-pinned in the source.

The independent longest-component audit is

```text
04-computation/lrc14_j5_singleton_component_seal_audit_codex_20260729.py
SHA-256 770e293ee62f312803ba508e87fcae03d245a6860c6f83c76c7c9808cc326f87

05-knowledge/results/lrc14_j5_singleton_component_seal_audit_codex_20260729.out
SHA-256 bacf919f0260703c2a0ce5f43959dad741a4c9aed03df517902cdb3d6dc72061
```

Its ordinary and optimized executions are also byte-identical.  All four
hashes are SHA-256 hashes of the raw stored bytes.

## 7. Exact scope

This theorem closes the complete eight-body/five-slot rung.
It does not change the distinguished runner, does not count finite heads as
positivity, and does not prove unrestricted LRC(14).  The live complementary
sector has at most seven normalized speeds in `{1,...,14}`, equivalently
`v_8>=15`.
