---
id: THM-513
title: The FKN neighborhood of the transitive tiling is root-packeted: one-flip atoms are the free interval roots, and two-flip triangle counts are controlled exactly by packet incidence
status: PROVED (one-flip score/c3/H laws, one-flip iso distinctness, two-flip c3 packet law). VERIFIED computationally for n=4..8 in fkn_root_packet_shells_codex_s11.py.
source: codex-2026-06-15-S11
depends_on:
  - THM-284   # one-flip H = 1 + 2^(skip-1)
  - THM-301   # same-end / cross-end quadratic sign rule
  - THM-302   # same-end quadratic formula
  - THM-511   # FKN/converse-parity framing
related:
  - THM-220   # near-transitive longest-root case
  - THM-287   # quadratic OCF decomposition
  - HYP-2537  # low-energy shells are packeted by root incidence
  - reflection: the-fkn-neighborhood-of-the-transitive-tiling-is-root-packeted-codex
---

# THM-513 — the local FKN geometry is the interval-root geometry

Fix the tiling model with base Hamiltonian path

```text
n -> n-1 -> ... -> 1
```

and free tiles `(x,y)` with `x >= y+2`. Write

```text
g(x,y) = x-y-1
```

for the number of vertices strictly between `y` and `x`.

The free tiles are exactly the non-simple interval roots of type `A_{n-1}`:
the interval `[y,x]` corresponds to the positive root
`alpha_y + ... + alpha_{x-1}` of height `x-y = g+1`.

## Statement

### 1. One-flip shell

Let `T_(x,y)` be the tournament obtained from the transitive ground state by
reversing the single free arc `x->y` to `y->x`. Then:

1. **Score defect.**

   ```text
   score(T_(x,y)) = score(transitive) + e_y - e_x.
   ```

   So a one-flip atom is a unit score transfer from `x` down to `y`.

2. **Triangle count.**

   ```text
   c3(T_(x,y)) = g(x,y).
   ```

3. **Hamiltonian path count.**

   ```text
   H(T_(x,y)) = 2^(g(x,y)) + 1.
   ```

4. **Iso-class distinctness.** Distinct free tiles give distinct isomorphism
   classes. Hence the one-flip shell contains exactly

   ```text
   m = C(n-1,2)
   ```

   pairwise non-isomorphic tournaments.

So the radius-1 FKN atoms are not just anonymous coordinates: they are the
free interval roots themselves.

### 2. Two-flip shell: exact packet law for c3

Let `t1, t2` be two distinct free tiles with gaps `g1, g2`.
There are three incidence types:

- `same-end`: they share upper or lower endpoint;
- `cross-end`: one tile's upper endpoint is the other's lower endpoint;
- `disjoint`: they share no endpoint.

For the tournament `T_{t1,t2}` obtained by reversing both tiles,

```text
c3(T_{t1,t2}) =
  g1 + g2 - 1     if {t1,t2} is same-end,
  g1 + g2 + 1     if {t1,t2} is cross-end,
  g1 + g2         if {t1,t2} is disjoint.
```

Thus the first nonlinear cyclic statistic is controlled exactly by the packet
incidence law of the two interval roots.

## Proof

### One flip: score, c3, H

The transitive tournament has score vector `(0,1,...,n-1)` (up to ordering by
labels), because vertex `k` beats exactly the lower labels.

Reversing `x->y` changes only those two vertices:

- `x` loses one win;
- `y` gains one win.

So the score defect is exactly `e_y-e_x`.

For `c3`, every new directed triangle must use the reversed arc `y->x`, because
all other arcs still follow the transitive order. The only possibilities are

```text
x -> z -> y -> x    with y < z < x,
```

one for each intermediate vertex `z`. There are exactly `x-y-1 = g(x,y)` such
vertices, proving the triangle formula.

The Hamiltonian-path formula is exactly THM-284:

```text
H(T_(x,y)) = 1 + 2^(x-y-1) = 1 + 2^g.
```

### One flip: iso distinctness

The score multiset of `T_(x,y)` is

```text
{0,1,...,n-1} \ {y-1, x-1} U {y, x-2}
```

as a multiset.

- If `x = y+2`, the unique duplicated score is `y`, while the missing scores are
  `y-1` and `y+1`, so the pair `(x,y)=(y+2,y)` is recovered.
- If `x > y+2`, the duplicated scores are `y` and `x-2`, so
  `y = min(duplicated)` and `x = max(duplicated)+2`.

Hence the score multiset determines `(x,y)`. Since score multiset is an
isomorphism invariant, distinct one-flip tiles are pairwise non-isomorphic.

### Two flips: c3 packet law

For a tile `t=(x,y)`, let

```text
F(t) = { {x,z,y} : y < z < x }.
```

This is exactly the family of directed triangles created by flipping `t`, so
`|F(t)| = g(t)`.

Now compare `F(t1)` and `F(t2)`.

1. **Same-end.** Suppose the tiles share their upper endpoint `x`, with
   `y1 < y2 < x` (the lower-end case is symmetric). Flipping `(x,y1)` alone
   creates `g1` triangles, one of which is

   ```text
   x -> y2 -> y1 -> x.
   ```

   After also flipping `(x,y2)`, that triangle is destroyed because the arc
   `x->y2` has been reversed. All other `g1-1` triangles from `(x,y1)` remain,
   and the second tile contributes its own `g2` triangles. Therefore

   ```text
   c3 = (g1-1) + g2 = g1 + g2 - 1.
   ```

2. **Cross-end.** Suppose `t1=(r,y)` and `t2=(x,r)` with `y<r<x`. Then
   `F(t1)` and `F(t2)` are disjoint, and in addition there is one new relay
   triangle

   ```text
   x -> r -> y -> x
   ```

   using both reversed arcs. Therefore

   ```text
   c3 = |F(t1)| + |F(t2)| + 1 = g1 + g2 + 1.
   ```

3. **Disjoint.** If the tiles share no endpoint, no triangle can contain both
   reversed arcs. So the two triangle families are disjoint and simply add:

   ```text
   c3 = |F(t1)| + |F(t2)| = g1 + g2.
   ```

This proves the packet law.

## Consequences

1. **The free dictator atoms are interval-addressed.** T823's "single-arc
   dictator atoms" are not all equivalent: there are `C(n-1,2)` of them, one
   for each free interval root.

2. **Different invariants forget different amounts.** At radius 1, the full
   iso class remembers the interval address, the score remembers the transfer
   `e_y-e_x`, and `H,c3` remember only the interval height/gap.

3. **The first nonlinear layer is packet incidence.** The two-flip triangle law
   is not a function of Hamming weight alone; it depends on whether the roots
   are same-end, relay-linked, or disjoint.

4. **Second-order H already knows this packet structure.** By THM-301/302,
   same-end pairs give a negative quadratic defect in `H`, while cross-end pairs
   give a positive one; the stored computation verifies nonnegative disjoint
   defect through `n=8`.
