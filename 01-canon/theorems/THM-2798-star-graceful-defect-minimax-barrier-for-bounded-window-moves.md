---
id: THM-2798
title: "Star graceful-defect minimax barrier for bounded-window moves"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For a labelled m-edge star whose centre occupies position c, the
  missing-difference defect is exactly min(c,m-c).  Between the two graceful
  extreme orders, moves supported on at most r+1 consecutive positions have
  exact minimax barrier max(0,ceil((m-r)/2)).  Thus adjacent C2 swaps cost
  floor(m/2), while adjoining consecutive C3 rotations costs
  floor((m-1)/2); every fixed-width local grammar has unbounded defect.
source: root/star-defect-barrier-2026-07-28
depends_on:
  - THM-2795-tree-star-circuit-dwell-time-v4-diamonds-and-ternary-move-boundary
related:
  - THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion
  - THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction
  - THM-2768-modular-c2-c3-a4-s4-bass-serre-quotient
script: 04-computation/tree_star_graceful_defect_barrier_thm2798.py
output: 05-knowledge/results/tree_star_graceful_defect_barrier_thm2798.out
script_sha256: b139709e4a820336e572caf995f01dda400867de883cbc6d53d2a5dadf47dc10
output_sha256: 9b2bbabfc376ee8f929d475af3c73d4845a5b38ac6273b6b2778f08542e7da6d
hash_basis: LF-normalized bytes
---

# THM-2798 -- star graceful-defect minimax barrier for bounded-window moves

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2795 proves that the binary swap and ternary rotation grammar is locally
an affine `V4`/`S3` geometry but does not connect the graceful locus even on
stars.  The present theorem measures the exact price of leaving that locus.
The answer is linear in the size of the star.  In particular, adjoining the
`C3` generator to the `C2` generator lowers the barrier by at most one.

## 1. Star defect

Let `K_(1,m)` have centre `o` and labelled leaves `v_1,...,v_m`.  An order
of its `m+1` vertices gives the labels `0,...,m`.  If the centre occupies
position `c`, the set of absolute edge differences is

```text
{1,...,c} union {1,...,m-c}
  ={1,...,max(c,m-c)}.                              (1)
```

Define the missing-difference defect

```text
delta_m(c)
 =m-#{distinct edge differences}
 =min(c,m-c).                                       (2)
```

This is a collision defect, not a new graph invariant on arbitrary trees.
For the star it vanishes exactly at

```text
c=0 or c=m,                                         (3)
```

so the two centre positions are precisely the two graceful extremes.

## 2. Bounded-window move graph

Fix `r>=1`.  A width-`r` local move arbitrarily permutes the entries in one
consecutive window of at most `r+1` positions.  Such a move changes the
centre position by at most `r`.

Fix the two orders

```text
pi_L=(o,v_1,...,v_m),
pi_R=(v_1,...,v_m,o),                               (4)
```

which preserve the relative order of all leaves.  For a move path `P`
between them, put

```text
cost(P)=max_(pi in P) delta_m(position_pi(o)).       (5)
```

Then the exact minimax cost is

```text
B_r(m)
 =min_P cost(P)
 =max(0,ceil((m-r)/2)).                             (6)
```

### Lower bound

If every order on a path has defect at most `B`, then `(2)` confines the
centre to the two banks

```text
[0,B] union [m-B,m].                                (7)
```

A path from the left bank to the right bank contains a move crossing their
gap.  Since one move displaces the centre by at most `r`,

```text
m-2B<=r.                                            (8)
```

For `r<m`, this gives `B>=ceil((m-r)/2)`; for `r>=m`, the tautological
lower bound is zero.  This is `(6)`.

### Matching path

Let `B=B_r(m)`.  Move the centre from `0` to `B` in increments of at most
`r`, cross once from `B` to `m-B`, and then move from `m-B` to `m` in
increments of at most `r`.  Inequality `(8)` holds by the definition of
`B`.  A rightward move of `s<=r` is implemented by left-rotating the
consecutive block

```text
(o,w_1,...,w_s)  ->  (w_1,...,w_s,o).               (9)
```

It preserves the relative leaf order.  Every centre position before the
crossing lies at most `B`, and every position after it lies at least
`m-B`; hence every defect is at most `B`.  This proves equality in `(6)`.

## 3. The binary and ternary generators

For adjacent transpositions alone, `r=1`, and `(6)` gives

```text
B_C2(m)=ceil((m-1)/2)=floor(m/2).                   (10)
```

For adjacent transpositions together with cyclic rotations of consecutive
triples, the centre can move by one or two positions and the construction
`(9)` uses only those generators.  Therefore

```text
B_(C2,C3)(m)=ceil((m-2)/2)=floor((m-1)/2).          (11)
```

Consequently,

```text
B_C2(m)-B_(C2,C3)(m)
 =1  if m is even,
 =0  if m is odd.                                   (12)
```

The special roles of `2` and `3` are exact here: they are the orders of the
local swap and rotation, and their supports allow centre displacements one
and two.  The result is not a free `PSL_2(Z)=C2*C3` action.  THM-2795
already identifies the extra local `S3` relation; `(10)--(12)` show that
even the combined local grammar remains macroscopically trapped.

More generally, for every fixed `r`,

```text
B_r(m)=m/2+O_r(1).                                  (13)
```

Thus no fixed-width, single-seed local-move argument can connect the two
graceful star extremes while keeping a uniformly bounded number of missing
differences.  A uniform defect allowance `B` requires a move radius at
least

```text
r>=m-2B.                                            (14)
```

This isolates the missing operation: a successful connected-move program
must use a window growing linearly with `m`, change the underlying tree,
use multiple seeds, or carry an additional invariant that permits a
different induction.

## 4. Exact controls

The companion script verifies:

1. `(1)--(2)` for every centre position through `m=80`;
2. the minimax formula on all `3,320` pairs
   `1<=m<=80`, `1<=r<=m+1`;
3. the explicit relative-leaf-order-preserving path `(9)` in every case;
4. the full permutation move graph, rather than only its centre projection,
   for all `15` pairs through `m=6` and `r<=3`; and
5. normal, optimized, and stored transcript identity.

Run:

```text
python 04-computation/tree_star_graceful_defect_barrier_thm2798.py
python -O 04-computation/tree_star_graceful_defect_barrier_thm2798.py
```

The finite controls are not used in the proof.

## 5. Scope and stopping boundary

The theorem is an exact obstruction to one particular proof architecture.
It does **not** disprove the Graceful Tree Conjecture, obstruct several
starting graceful labelings, or rule out a leaf-extension induction such as
THM-2787.  Nor does it say that every move path through arbitrary trees has
the star barrier.  It says precisely that the most literal fixed-tree
binary/ternary local grammar cannot move between the complementary graceful
star components without accumulating an unbounded collision defect.
