---
id: THM-3135
title: "Projected forest-boundary parity and augmentation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  complete graph K_n and 0<=r<=n-2, the projection of the oriented
  (r+1)-forest boundary onto rank-r component partitions is surjective over
  Q when r is even.  When r is odd, its image is exactly the augmentation-zero
  hyperplane.  Thus THM-3117's K8/K9 rank-four surjectivity is an instance of
  a sharp all-rank parity law rather than a sporadic modular computation.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  An independent hostile audit rederived the endpoint-swap cancellation,
  realization of every one-point transfer, fixed-rank-layer connectivity,
  augmentation-zero pair span, parity ranks over F_2, characteristic-zero
  rank inference, signed augmentation obstruction, and boundary-lift scope.
  It independently replayed normal and optimized modes against the stored
  output, matched the declared LF hashes, and confirmed the exact relation
  to THM-3117 without importing its separate positivity no-go.
depends_on: []
related:
  - THM-3117-projected-five-forest-boundary-surjectivity-and-signed-holotopy-lift
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
script: 04-computation/gmc_projected_forest_boundary_parity_thm3135.py
output: 05-knowledge/results/gmc_projected_forest_boundary_parity_thm3135.out
script_sha256: d4866d44c4d980e38aaeefb1de5d990e4b4096b646b20ef1744548937a13d8de
output_sha256: c12dc652a314dc1c31ce2260ede74b185e0f5b2685e0bae5cd9bdb6452fc2f9c
hash_basis: LF-normalized bytes
---

# THM-3135 -- projected forest-boundary parity and augmentation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3117 proves by exact modular rank computation that every rational
rank-four current on the actual `K_8` and `K_9` product-Gamma macro spaces is
the projection of a rational forest boundary.  The mechanism is not special
to rank four or to those two vertex counts.  A local attachment swap produces
every difference between adjacent set partitions, and parity decides whether
one remaining augmentation direction is present.

The resulting theorem is a complete classification.  It also marks a sharp
holotopy boundary: even-rank macro currents have no linear obstruction to a
signed forest-boundary lift, while an odd-rank current has exactly one, namely
its total mass.

## 1. The projected boundary

Fix `n>=2` and a total order on the edges of `K_n`.  For `s>=0`, let
`C_s(n;Z)` be the free abelian group on the canonically ordered `s`-edge
forests.  Use the ordinary deletion boundary

```text
partial_s[e_1,...,e_s]
 =sum_(j=1)^s (-1)^(j-1)[e_1,...,omit e_j,...,e_s].          (1)
```

Let `Pi_(n,r)` be the set of partitions of `[n]` of atomic rank `r`, meaning
that they have `n-r` blocks, and let

```text
V_(n,r;Z)=Z[Pi_(n,r)].                                      (2)
```

Every `r`-forest has a rank-`r` component partition.  Forget its internal
tree and orientation by

```text
P_r:C_r(n;Z)->V_(n,r;Z),       P_r[F]=[pi(F)].               (3)
```

The map of interest is

```text
T_(n,r)=P_r partial_(r+1):C_(r+1)(n;Z)->V_(n,r;Z).           (4)
```

We assume

```text
0<=r<=n-2,                                                   (5)
```

so every `r`-forest has at least two components and an additional forest
edge can exist.  The excluded top rank `r=n-1` is genuinely different:
there is no `n`-edge forest, so `(4)` is the zero map.

## 2. Endpoint-swap pair lemma over `F_2`

Reduce `(4)` modulo two and write the result as `Tbar_(n,r)`.  Let `A` be an
`r`-edge forest, let

```text
a_0={x,y} in A,                                              (6)
```

and let `d` lie in a component of `A` different from the component containing
`x,y`.  The two edge sets

```text
F_x=A union {{x,d}},          F_y=A union {{y,d}}             (7)
```

are `(r+1)`-forests.  Compare their projected boundaries modulo two.

- Deleting the added edge gives `pi(A)` in both columns, so those terms
  cancel.
- If `a!=a_0` is deleted from `A`, the vertices `x,y` remain connected by
  `a_0`.  Attaching `d` at `x` or at `y` therefore gives the same component
  partition, so those terms cancel as well.
- After deleting `a_0`, the two columns attach the `d`-component to opposite
  sides of the split made by `a_0`.

Consequently

```text
Tbar[F_x]+Tbar[F_y]=[pi_x]+[pi_y],                            (8)
```

where `pi_x,pi_y` are two distinct rank-`r` partitions.

Every elementary one-point transfer is realized by `(8)`.  Indeed, suppose
`pi_x` has blocks

```text
X union {d},       Y,       and the unchanged remaining blocks,            (9)
```

with `X,Y` nonempty, and `pi_y` is obtained by moving `d` from the first
block to `Y`.  Start with spanning forests on the finer blocks

```text
X, {d}, Y, and the unchanged blocks,                         (10)
```

join `X` to `Y` by `a_0={x,y}`, and then apply `(7)`.  Deleting `a_0`
produces exactly the two partitions in `(9)`.

## 3. Connectivity of the fixed-rank layer

For `r>=1`, make a graph on `Pi_(n,r)` by joining two partitions when one is
obtained from the other by moving one point out of a block of size at least
two and into another block.  This graph is connected.

To see this, choose one block as an anchor.  From every other nonsingleton
block, move all but one point into the anchor.  The result has one block of
size `r+1` and `n-r-1` singleton blocks.  Any two such partitions are
connected by swaps of singleton identities: if `x` is in the large block and
`y` is a singleton, first move `x` into `{y}`, then move `y` back into the
large block.  Repeating reaches a fixed canonical large block.  Reversing
these moves connects any two vertices.

By `(8)`, the image of `Tbar_(n,r)` therefore contains

```text
[pi]+[pi']        for every pi,pi' in Pi_(n,r).               (11)
```

These pair sums span the augmentation-zero hyperplane

```text
A_(n,r)={z in F_2[Pi_(n,r)]:sum_pi z_pi=0}.                   (12)
```

For `r=0`, the target has one basis vector and every one-edge forest maps to
it, so the theorem below is immediate.

## 4. The parity classification

Every column of `Tbar_(n,r)` contains `r+1` boundary faces, counted with
their exact parity even when two component partitions coincide.  Hence:

- if `r` is odd, every column has even augmentation and lies in `(12)`;
- if `r` is even, every column has odd augmentation and lies outside `(12)`.

Combining this with `(11)--(12)` gives the exact finite-field rank law

```text
rank_F2 Tbar_(n,r)
 = |Pi_(n,r)|       if r is even,
 = |Pi_(n,r)|-1     if r is odd.                              (13)
```

Return now to characteristic zero.  An integer matrix has rational rank at
least its rank modulo two.  If `r` is even, `(13)` therefore makes `(4)`
surjective over `Q`.  If `r` is odd, the signed column augmentation is

```text
sum_(j=1)^(r+1)(-1)^(j-1)=0.                                 (14)
```

Thus the rational image lies in the augmentation-zero hyperplane.  Its rank
is at least `|Pi_(n,r)|-1` by `(13)`, so equality holds.  We have proved

```text
Im_Q T_(n,r)=V_(n,r;Q)                         for r even,
Im_Q T_(n,r)=ker(sum:V_(n,r;Q)->Q)             for r odd.     (15)
```

This conclusion is independent of the chosen total edge order.  The map
itself changes its signs with that choice, but the mod-two pair argument and
the signed augmentation bound force the same rational image dimension and
the same image in both parity cases.

## 5. Forest-boundary lifts and the product-Gamma consequence

Let `W` be a rational rank-`r` partition current satisfying the condition in
`(15)`; for even `r` there is no condition, and for odd `r` it is exactly
`sum W=0`.  Choose `Y` with

```text
T_(n,r)Y=W
```

and put `C=partial_(r+1)Y`.  Then

```text
partial_r C=0,                  P_r C=W.                       (16)
```

Thus `(15)` is also the exact existence criterion for a projected
**forest-boundary** lift.  It does not classify arbitrary forest cycles not
known to be boundaries.

At `r=4`, `(15)` applies for every `n>=6`.  In particular it recovers the
`K_8` and `K_9` surjectivity used by THM-3117 and explains its modular ranks
without enumerating `13,683` and `44,524` five-forests.  THM-3117's separate
one-sign propagation theorem remains fully load-bearing: `(15)` supplies a
signed rational lift, not a nonnegative or sign-preserving one.  It gives no
denominator, sparsity, locality, symmetry, insertion compatibility, or
positive Young-refinement transport.

The theorem also identifies the next topological boundary sharply.  At odd
rank, total mass is the only obstruction to a boundary lift.  At even rank,
there is no linear projected-boundary obstruction at all.  Any remaining GMC
failure must therefore live in positivity, locality, the macro-fibre gauge,
or compatibility with the response operation rather than in bare forest
homology.

## 6. Exact companion and scope

Reproduce with

```text
python 04-computation/gmc_projected_forest_boundary_parity_thm3135.py
python -O 04-computation/gmc_projected_forest_boundary_parity_thm3135.py
```

The companion independently enumerates every target rank for `K_3` through
`K_7`.  It checks `40,224` forest columns rank by rank, verifies full rank in
every even rank and codimension one in every odd rank, proves connectivity of
every finite transfer graph, and checks all `7,134` endpoint-swap identities
directly.  The exact printed forest-column total is retained rank by rank in
the stored output; no floating-point arithmetic is used.  Normal and
optimized runs byte-match that output.

This theorem is a general combinatorial holotopy theorem.  It does not prove
a positive forest lift, product-Gamma width-three positivity, the Gaussian
Moment Conjecture, NC2, LRC(14), JC(2), or DC(2).

QED.
