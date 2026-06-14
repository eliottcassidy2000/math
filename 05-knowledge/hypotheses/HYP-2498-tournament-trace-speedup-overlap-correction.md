# HYP-2498 - Tournament trace speedups hit their first overlap correction at c6

**Status:** OPEN speedup/proof program with one proved correction identity and exact finite audits.
**Source:** codex-2026-06-13.
**Extends:** THM-118, THM-498, HYP-2492, HYP-2493, HYP-2487, T815.
**Artifacts:** `04-computation/tournament_structure_speedup_patterns_codex.py`, `05-knowledge/results/tournament_structure_speedup_patterns_codex.out`, `07-reflections/tournament-trace-speedups-and-overlap-corrections.md`.
**Rebase integration:** kind-pasteur-2026-06-13-S5 / THM-499 proves the sharper
`n<=6` OCF formula `H=1+2(c3+c5)+4D`, where `D` is the number of
vertex-disjoint directed-triangle pairs.  This HYP is the trace-side mirror:
`p33_meet = C(c3,2)-D`, so the first `c6` trace correction records the
complementary intersecting-pair geometry.

## Claim

For tournament structure computations, the trace method has a clean first
boundary:

```text
c_k(T) = tr(A^k)/k  for k=3,4,5,
```

because tournaments have no loops and no directed `2`-cycles. Any non-simple
closed walk must split into at least two directed cycles, so it has length at
least `3+3=6`.

At length `6`, the trace is still efficiently correctable, but the correction
is not a scalar rooted-triangle term. It is the first cycle-overlap term:

```text
tr(A^6) = 6*c6 + 3*c3 + 6*p33_meet,
```

where `p33_meet` is the number of unordered pairs of distinct directed
triangles with nonempty vertex intersection. Equivalently,

```text
c6 = (tr(A^6) - 3*c3 - 6*p33_meet)/6.
```

The naive half-return correction

```text
sum_v tau_v^2,  tau_v=(A^3)_{vv},
```

is insufficient, because many non-simple length-6 closed walks are rotations
of a two-triangle figure-eight that do not return to their starting vertex at
time `3`.

## Proof Sketch

A closed directed walk of length at most `5` in a tournament cannot repeat a
vertex. If it did, the repeated vertex would split the walk into shorter
closed directed walks, each of length at least `3`, impossible below length
`6`.

For length `6`, every non-simple closed walk splits into two directed
triangles.

There are two cases:

1. The two triangles are the same triangle. Traversing a directed triangle
   twice contributes `3` rooted closed walks, one from each vertex. This gives
   `3*c3`.
2. The two triangles are distinct and intersect. Each unordered intersecting
   pair gives one cyclic figure-eight word and all `6` of its rotations as
   rooted closed walks. This gives `6*p33_meet`.

All remaining length-6 closed walks are simple directed `6`-cycles, each
contributing `6` rotations. This proves the identity.

## Computed Evidence

The script validates the formulas by brute force:

```text
exhaustive n=3..6 checked=33864, mismatches={}
random n=7..9 samples per n=80, mismatches={}
```

Benchmark on a fixed random `n=14` tournament:

```text
c5 trace = 1816, brute = 1816, speedup about 13.4x
c6 corrected = 7113, brute = 7113, speedup about 106.1x

trace6 = 59949
3*c3 = 309
6*p33_meet = 16962
6*c6 = 42678
```

The failed scalar diagnostic is visible in the same row:

```text
sum_v tau_v^2 = 6855,
```

which is not the full non-simple correction `3*c3 + 6*p33_meet = 17271`.

For Hamiltonian paths, subset DP beats permutation brute force on a fixed
random `n=8` tournament:

```text
H = 235
DP speedup about 40.2x
```

## Information Tournament

The scout exhausts all `32768` labelled tournaments on `n=6` vertices and
builds a Tournament Analysis over invariants:

```text
vertices = score, c3, c4, c5, c6, H
observable = U(Y|X) = 1 - entropy(Y|X)/entropy(Y)
gauge = X -> Y when X explains a larger fraction of Y than Y explains X
tie path = score -> c3 -> c4 -> c5 -> c6 -> H
```

Result:

```text
outdegrees: H=5, score=4, c4=3, c5=2, c3=1, c6=0
directed_3_cycles: 0
SCC sizes: six singleton SCCs
Hamiltonian paths: 1
champion order: H > score > c4 > c5 > c3 > c6
```

This tournament is an information order, not a computation order. `H` is most
information-rich, while `c6` has low entropy and is often explained by richer
invariants.

The structural surprise is the bucket audit:

```text
keys=(score) -> H has 9 mixed buckets
keys=(score,c5) -> H has 7 mixed buckets
keys=(score,c5,c6) -> H has 1 mixed bucket
keys=(c3,c4,c5,c6) -> H has 0 mixed buckets
```

Thus, for `n=6`, the low corrected cycle vector `(c3,c4,c5,c6)` determines
`H`.  THM-499 explains why: the primitive non-spectral ingredient is
`D=alpha_2`, the disjoint directed-triangle-pair count.  Since
`p33_meet=C(c3,2)-D`, the `c6` correction is seeing the complementary
intersecting-pair side of the same boundary.

## Program

The next speed frontier is not "use traces blindly." It is:

1. Build a pattern-correction engine for `tr(A^k)` for `k>=7`, where the
   correction terms are closed-walk support types.
2. Track which corrections are scalar counts and which require placement data
   such as intersecting or disjoint cycle pairs.
3. Test whether corrected trace vectors `(c3,...,ck)` continue to compress or
   determine `H` at `n=7`, and compare each trace correction with the OCF
   conflict ingredients `alpha_2, alpha_3, ...`.
4. Use the information tournament as a guide for which invariants carry unique
   structure and which are shadows of richer ones.

The deeper pattern is that efficient computation becomes proof-shaped exactly
when the speedup records the obstruction it had to subtract. For `c6`, that
obstruction is intersecting triangle pairs.
