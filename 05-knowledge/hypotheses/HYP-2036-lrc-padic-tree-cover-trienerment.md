---
id: HYP-2036
status: SUPPORTED
source: codex-2026-06-01-S546
related:
  - HYP-2017
  - HYP-2019
  - HYP-2029
  - HYP-2031
  - HYP-2032
  - HYP-2033
  - HYP-2035
  - THM-369
---

# HYP-2036: p-adic zero-branch covers form a tree trienerment for LRC sieve channels

**Claim.** The p-adic tree layer becomes computationally sharp when its
vertices are product zero-branch obligations

```text
q = 2,3,...,n,
z_q(V) = #{v in V : q divides v}.
```

The node `q` is the observer residue ball `0 mod q` in the product of the
p-adic trees for the primes dividing `q`.  If `z_q=0` and `q<n`, then `t=1/q`
is an open THM-369 sieve witness.  If `z_n=0`, it is the compactified AP-style
wall witness at `t=1/n`.

The associated tree trienerment has vertices `q=2..n`.  Lower zero-branch mass
wins; equal-mass comparable nodes orient toward the deeper divisibility branch;
equal-mass incomparable nodes tie.  Thus ties record equal cover debt in
incomparable p-adic branches.

**Evidence.** `lrc_padic_tree_trienerment_s546.py` scans primitive boxes
`n=6,8,10,12,14,16,18` with `max_speed=n+2`.  It measures open sieve
survivors, compact survivors, singleton carriers, exact minimum set covers of
the open obligations, and tree-trienerment fingerprints.

```text
n   sets  open survivors  compact survivors  min open cover  mixed classes
6    56        25                 16                3              0
8   120        44                 30                4              0
10  220        70                 50                5              0
12  364       104                 77                6              0
14  560       147                112                7              0
16  816       200                156                8              0
18 1140       264                210                9              0
```

In every scanned `n`, all open survivors have exact minimum open cover size
`n/2`.  The tree-trienerment fingerprint classes have zero mixed
open-survivor/non-survivor fibers.  This is stronger than a raw channel-rank
summary: the product-tree cover profile itself is a pure certificate for
surviving all open THM-369 probes in these boxes.

The AP row has the same pattern throughout:

```text
AP: open_survivor=True, compact_survivor=False,
    open_empty=(), compact_empty=(n,),
    singleton nodes=(n/2,n/2+1,...,n-1),
    min_cover=n/2.
```

So the product p-adic tree sees AP as the canonical open-sieve survivor whose
only remaining sieve target is the compactified wall node `q=n`.

**Interpretation.** HYP-2035 says the coarse channel rank is `omega(n/2)`.
HYP-2036 adds a cover-core layer below that rank:

```text
rank / channel skeleton = which p-adic towers exist;
cover core = which speeds carry every zero branch q=2..n-1;
singleton carriers = forced boundary leaves of the product tree.
```

The fact that strict directed 3-cycles are always zero is not a failure of
Tournament Analysis.  It is the fingerprint of an ultrametric/product-poset
object: this p-adic relation is hierarchical rather than cyclic.  The
Hamiltonian-path count and ternary `B1` measure cover degeneracy and tie
freedom inside that hierarchy.

**Assumption challenge.** The session considered runners, gaps, fixed sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier/Gabor
cells, p-adic residue balls, composite CRT product nodes, and proof
obligations.

Chosen quotient:

```text
vertices = product p-adic zero-branch obligations q=2..n.
```

Predicate preserved:

```text
whether a THM-369 open sieve witness exists, and whether the q=n wall witness
remains after compactification.
```

Information destroyed:

```text
non-sieve LRC witnesses, exact circular gap shape, wall order, and most runner
identity beyond singleton carrier roles.
```

Challenged assumption:

```text
the p-adic tree only records channel rank or denominator shape.
```

The scan refines that: the tree also records a cover core with pure
open-sieve survivor fibers and forced singleton carriers.

**Predictions.**

1. The zero-mixed-fiber result persists for larger bounded boxes when the
   fingerprint includes `open_empty`/`compact_empty` counts and cover-core
   statistics.
2. For even `n` in initial-segment style boxes, every open survivor has minimum
   open cover size `n/2`; proving this would isolate the boundary-leaf burden.
3. AP/regular-polygon rows are exactly open-cover survivors with the remaining
   compact wall `q=n` and singleton carriers `n/2,...,n-1`.
4. The next symbolic-dynamics alphabet should include the first empty
   zero-branch `q`, plus the singleton carrier when `z_q=1`.
5. Gabor zero-column events should be indexed by these product-tree nodes:
   `q<n` gives an open observer-zero window at `t=1/q`, and `q=n` gives the
   compactified wall analogue.

**Files.** `04-computation/lrc_padic_tree_trienerment_s546.py`;
`05-knowledge/results/lrc_padic_tree_trienerment_s546.out`;
`07-reflections/lrc-padic-tree-cover-trienerment-s546.md`.
