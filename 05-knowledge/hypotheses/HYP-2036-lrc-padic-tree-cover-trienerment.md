---
id: HYP-2036
status: SUPPORTED
source: codex-2026-06-01-S546
related:
  - HYP-1811
  - HYP-1812
  - HYP-1992
  - HYP-2017
  - HYP-2019
  - HYP-2025
  - HYP-2026
  - HYP-2029
  - HYP-2031
  - HYP-2032
  - HYP-2033
  - HYP-2035
  - HYP-2037
  - HYP-2038
  - HYP-2039
  - THM-369
  - THM-390
  - THM-391
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

**Prime-power endpoint-core supplement.** The companion audit
`lrc_prime_power_zero_branch_core_s546b.py` restricts to prime powers
`q=p^d <= n` and asks whether a covered zero branch can itself hold a local
endpoint-protection core.  For each q-divisible speed, it builds the danger
intervals centered at every killed unit point `u/q`.

```text
Audited local branch cores: 86, nonempty=0
Audited full cover cores: 9, nonempty=0
```

This refines the product-tree cover scan: a branch with `z_q>0` really does kill
the THM-369 unit witness, but the local intervals around that unit point are
nested stars and peel to empty.  The zero branch is a gate/debt carrier, not a
local counterexample core.  In n=18, this holds even while the `n*=9` modular
full-support zero-flow counts remain enormous (`~1.88e14` to `~2.50e14`).

**Formalized sublayer.** THM-390 proves the exact sieve semantics behind this
coordinate: an empty `q<n` node gives the open witness `t=1/q`, an empty
`q=n` node gives the compactified wall witness `t=1/n`, and the AP row
`{1,...,n-1}` has a unique minimum open cover core

```text
U_n = {u in {1,...,n-1} : 2u >= n}
```

of size `floor(n/2)`.  The even-denominator scans in this hypothesis are the
special case where the proven AP core has size `n/2` and singleton carriers
`n/2,...,n-1`.

This leaves the genuinely computational part of HYP-2036 intact: whether the
zero-mixed tree-trienerment fibers persist beyond the bounded scan, and whether
the observed `n/2` minimum-cover law extends from AP to every open survivor in
the tested even boxes for structural rather than accidental reasons.

**Entropy integration.** Incoming HYP-2037/HYP-2038 add a complementary
entropy reading: AP/regular rows are the high-entropy or critical boundary
cases, while THM-390 identifies the exact finite denominator leaves that force
the AP wall core.  Thus the p-adic entropy layer should distinguish total
channel spread from the local zero-branch occupancy proved here: entropy can
say the row is spread or critical, but the `z_q=0` obligation is the exact
sieve gate that creates an explicit witness.

**Defect-transport integration.** HYP-2039 reframes LRC as transporting the
guaranteed large hole to the observer.  In this denominator quotient, an empty
zero branch `z_q=0` is a pinned transport success: the guaranteed hole is
already at the observer at `t=1/q`.  Covered branches are the cases where the
hole must be moved by endpoint descendants, owners, or cross-branch coupling.

**Star-core update (THM-391).** S548 also proves the local statement in full
generality.  The proof does not use q being a prime power: for any
`2 <= q <= n`, any nonzero q-grid centers, and any speeds all divisible by q,
the centered q-grid star has empty strict endpoint-protection core.  Prime
powers matter because they are p-adic branch labels, not because their local
interval geometry is special.

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

HYP-2037/HYP-2038's entropy threads are complementary rather than
contradictory: p-adic/H entropy and lonely-set box-dimension measure global
spread, mixing, and phase boundaries on tree-like objects, while this HYP asks
whether the local zero branch keeps a nonpeeling proof core.  The S546b answer
is negative for bare prime-power branches, so entropy becomes proof-useful only
after endpoint descendants, event owners, critical-wall labels, or product-tree
coupling are retained.

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
6. Nontrivial cyclic branch trienerments require endpoint descendants, event
   owners, or cross-prime product-tree coupling; bare prime-power branch debt is
   transitive in the S546b audit.

**Files.** `04-computation/lrc_padic_tree_trienerment_s546.py`;
`05-knowledge/results/lrc_padic_tree_trienerment_s546.out`;
`04-computation/lrc_prime_power_zero_branch_core_s546b.py`;
`05-knowledge/results/lrc_prime_power_zero_branch_core_s546b.out`;
`01-canon/theorems/THM-391-lrc-zero-branch-star-core-peeling.md`;
`04-computation/lrc_zero_branch_star_theorem_s548.py`;
`05-knowledge/results/lrc_zero_branch_star_theorem_s548.out`;
`07-reflections/lrc-padic-tree-cover-trienerment-s546.md`;
`07-reflections/lrc-prime-power-zero-branch-endpoint-core-s546b.md`.
