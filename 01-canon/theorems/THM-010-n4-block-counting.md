# THM-010: n=4 Block-Counting Theorem

**Type:** Theorem (proved)
**Certainty:** 5 -- PROVED
**Status:** CANON
**Added by:** opus-2026-03-05-S1
**Source:** file.txt (inbox contribution)
**Tags:** #block-counting #3-cycle #n4 #bijection

---

## Statement

For any 4-vertex tournament T containing a directed 3-cycle C = (u -> w -> x -> u), exactly 3 Hamiltonian paths of T contain {u, w, x} as a contiguous block.

---

## Proof

Let z be the fourth vertex. The 3 cyclic orderings of the block are:
- (u, w, x): uses arcs u->w, w->x
- (w, x, u): uses arcs w->x, x->u
- (x, u, w): uses arcs x->u, u->w

For each ordering, z can appear before the block or after it (2 positions), giving 6 candidate paths total.

Paths with z before the block: z->(first of block)->...
- z -> u, w, x: valid iff z -> u (arc exists)
- z -> w, x, u: valid iff z -> w
- z -> x, u, w: valid iff z -> x
Count = [z->u] + [z->w] + [z->x] = outdeg(z in {u,w,x}) = |p(z)|

Paths with z after the block: ...(last of block) -> z
- u, w, x -> z: valid iff x -> z
- w, x, u -> z: valid iff u -> z
- x, u, w -> z: valid iff w -> z
Count = [x->z] + [u->z] + [w->z] = indeg(z from {u,w,x}) = 3 - |p(z)|

Total = |p(z)| + (3 - |p(z)|) = 3. QED

---

## Significance

This is a clean combinatorial result showing that at n = 4, the "contiguous block" count is exactly 3 regardless of how z relates to the 3-cycle. For general n, the analogous count H_C^+(T) = sum_Q f_C(Q) depends on the structure of the complement tournament T[V \ V(C)]. See THM-011 for the general formula.
