---
id: THM-809
title: The lower-first B3/S2 metagraph address is injective on all n=8 complement lines
status: FINITE-EXACT (all 1,048,576 n=8 complement lines)
source: codex-2026-07-15-S13
depends_on: [THM-553, THM-796, THM-801]
related: [THM-785, THM-790, HYP-6880]
verification:
  - 04-computation/mobius_cech_n8_frontier_codex_S13.cpp
  - 05-knowledge/results/mobius_cech_n8_frontier_codex_S13.out
---

# THM-809 — lower-face recursion already decides `n=8`

For an apex-zero `n=8` complement line `e`, let

```text
(a,a'), (b,b'), (c,c')
```

be the ordered merged-node pairs of its high, gap, and low `n=7` B3 faces.
Let `UABC` be the four-role colour word and let `S2_tau` be THM-801's raw
four-state mirror-pair counts on crossing layer `tau`; on the fixed layer keep
the zero/one counts.  Define the lower-only key

```text
Lambda(e)=((a,a'),(b,b'),(c,c'),UABC,(S2_tau)_(3<=tau<=8)). (1)
```

Then `Lambda` is injective on all `1,048,576` `n=8` complement lines.
Consequently THM-801's `Omega+B2` address is also injective at `n=8`, because
it contains every field of (1) and additionally retains the ordered upper
node pair.

This is stronger than constructing an `n=8` atlas and using its canonical
codes to repair lower collisions: no upper tournament classification occurs.

## Exact refinement ladder

The line partition evolves as follows:

| carrier | cells | excess | collision cells | maximum multiplicity |
|---|---:|---:|---:|---:|
| B3 lower node pairs + `UABC` | 1,048,158 | 418 | 418 | 2 |
| plus `tau=3` | 1,048,324 | 252 | 252 | 2 |
| plus `tau=4` | 1,048,428 | 148 | 148 | 2 |
| plus `tau=5` | 1,048,502 | 74 | 74 | 2 |
| plus `tau=6` | 1,048,524 | 52 | 52 | 2 |
| plus `tau=7` | 1,048,576 | 0 | 0 | 1 |
| plus fixed `tau=8` | 1,048,576 | 0 | 0 | 1 |

Thus every intermediate ambiguity is a double collision.  The new
size-three crossing layer `tau=7` is exactly the first complete separator;
the fixed layer adds no information after it.  This localizes the `n=8`
frontier much more sharply than the yes/no statement: within-layer crossing
position, not an upper isomorphism code, is the decisive new datum.

## Why the computation is exact

The verifier reads the previously certified `n=7` tiling-to-node array, whose
ranks lie in `0,...,271`, and enumerates the dense `2^20` complement-line
indices.  For each line it selects the unique apex-zero endpoint, applies the
three literal face maps, and records the ordered face/complement node ranks.

The nonfixed crossing-layer sizes at `n=8` are

```text
1,1,2,2,3                    for tau=3,...,7,
```

and the fixed layer has size three.  A four-state count vector on a layer of
size `s` has `binom(s+3,3)` possibilities.  Hence the exact mixed radix of S2
is

```text
4*4*10*10*20*4=128000 < 2^17.                            (2)
```

Six nine-bit node ranks, four colour bits, and (2) fit in 75 bits.  The
program sorts these exact integers; it uses no probabilistic hash and retains
literal line indices for every collision.  No full-key collisions occur.

Independent regressions give:

- all `2^20` literal B3 face-line triples are distinct;
- the closed ten-atom `UABC` law has zero failures;
- the B2 skew-depth histogram is exactly `2^11 binom(9,k)`;
- the final radix and every per-layer composition sum are checked literally.

## Positional moment threshold

For a mirror layer of at most three positions, a subset is determined by its
cardinality and the sum of its position labels.  Therefore adjoining the first
position moment separately to each of the four S2 states reconstructs every
raw layer assignment at `n=8`; the analogous fixed-layer moment reconstructs
its bits.  This gives an unconditional exact `S2+M1` fallback at this size.
It is not needed for (1), but identifies the anticipated failure mechanism:
raw counts can first lose subsets when a layer has four positions, because
two distinct two-subsets may have the same position sum.

The next finite target is consequently `n=9`: test whether its new
`tau=n-1` count layer still removes all lower-node collisions before appending
moments.  Static injectivity still does not prove continuation completeness;
THM-796's non-lumpability obstruction remains.

## Tournament Analysis and preservation boundary

The seven refinement stages in the table are the Tournament Analysis
vertices.  The pairwise observable is the number of unordered literal-line
pairs separated, the switches are raw retention and retention per encoded
bit, and the displayed recursive order is the tie Hamiltonian path.

The useful mathematical vertices are marked complement lines, gap-contracted
faces, and mirror layers—not unmarked tournament vertices.  `Lambda` preserves
the literal line within the audited size and all declared lower projections.
It destroys the upper class, automorphism/path stabilizer, owner-labelled LRC
state, metric gaps, wall chronology, and carry.  Its exactness is therefore a
finite recursive-codec theorem, not an LRC(14) proof or an all-size Markov
theorem.
